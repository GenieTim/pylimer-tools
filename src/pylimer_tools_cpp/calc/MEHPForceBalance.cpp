#include "MEHPForceBalance.h"
#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <nlopt.hpp>
#include <string>
#include <tuple>
#include <vector>

namespace pylimer_tools {
namespace calc {
  namespace mehp {

    /**
     * FORCE RELAXATION
     */
    void MEHPForceBalance::runForceRelaxation(
      long int maxNrOfSteps, // default: 10000
      double xtol,
      long int innerMaxNrOfSteps,
      double innerXtol,
      double innerAlphaTol)
    {
      this->simulationHasRun = true;
      double stress[3][3];

      for (size_t j = 0; j < 3; j++) {
        for (size_t k = 0; k < 3; k++) {
          stress[j][k] = 0.;
        }
      }

      ForceBalanceNetwork net = this->initialConfig;
      const int M = this->universe.getMolecules(crosslinkerType).size();
      const int N = this->universe.getMeanStrandLength(crosslinkerType) + 1;
      const double bM = this->universe.computeMeanBondLength();
      const int f =
        this->universe.determineFunctionalityPerType()[crosslinkerType];
      bool is2D = this->is2D;

      /* array allocation */
      Eigen::VectorXd u = Eigen::VectorXd::Zero(3 * net.nrOfLinks);

      /* force relaxation */
      double maxDistanceMoved = 0.0;
      size_t iterationsDone = 0;
      size_t totalInnerIterationsDone = 0;
      do {
        maxDistanceMoved = 0.0;
        // first, place cross-links
        for (size_t link_idx = 0; link_idx < net.nrOfNodes; ++link_idx) {
          double distanceMoved =
            this->displaceToMeanPosition(&net, u, link_idx);
          maxDistanceMoved = std::max(maxDistanceMoved, distanceMoved);
        }
        // then, place slip-link
        for (size_t link_idx = net.nrOfNodes; link_idx < net.nrOfLinks;
             ++link_idx) {
          size_t innerIterationsDone = 0;
          double displacementDone = 0.0;
          double parametrisationChange = 0.0;
          do {
            displacementDone = this->displaceToMeanPosition(&net, u, link_idx);
            parametrisationChange =
              this->updateSpringPartition(&net, u, link_idx);
            innerIterationsDone += 1;
          } while (displacementDone > innerXtol &&
                   innerIterationsDone < innerMaxNrOfSteps &&
                   parametrisationChange > innerAlphaTol);
        }
        iterationsDone += 1;
        // std::cout << "Iteration " << iterationsDone << " " <<
        // maxDistanceMoved
        //           << std::endl;
      } while (maxDistanceMoved > xtol && iterationsDone < maxNrOfSteps);

      // query solution & exit reason
      assert(u.size() == 3 * net.nrOfLinks);
      this->currentDisplacements = u;
      this->currentSpringDistances =
        this->evaluateSpringDistances(&net, this->currentDisplacements, is2D);

      this->exitReason = iterationsDone == maxNrOfSteps
                           ? ExitReason::MAX_STEPS
                           : ExitReason::X_TOLERANCE;
      this->nrOfStepsDone += iterationsDone;
      // TODO: instead of changing parametrisation in-place,
      // introduce a separate store
      this->finalConfig = net;
    }

    Eigen::Vector3d MEHPForceBalance::evaluateDistanceBetween(
      const ForceBalanceNetwork* net,
      const Eigen::VectorXd& u,
      const size_t linkIndexA,
      const size_t linkIndexB,
      const bool is2D)
    {
      Eigen::Vector3d distances;
      distances << net->coordinates[3 * linkIndexA] + u[3 * linkIndexA] -
                     (net->coordinates[3 * linkIndexB] + u[3 * linkIndexB]),
        net->coordinates[3 * linkIndexA + 1] + u[3 * linkIndexA + 1] -
          (net->coordinates[3 * linkIndexB + 1] + u[3 * linkIndexB + 1]),
        net->coordinates[3 * linkIndexA + 2] + u[3 * linkIndexA + 2] -
          (net->coordinates[3 * linkIndexB + 2] + u[3 * linkIndexB + 2]);

      // Possibly improvable PBC
      for (size_t j = 0; j < 3; ++j) {
        int iterations = 0;
        while (distances[j] > net->boxHalfs[j % 3]) {
          distances[j] -= net->L[j % 3];
          iterations++;
          if (iterations > 10) {
            throw std::runtime_error(
              "Too many iterations in PBC from " + std::to_string(linkIndexA) +
              " to " + std::to_string(linkIndexB) + ", currently at " +
              std::to_string(distances[j]));
          }
        }
        iterations = 0;
        while (distances[j] < -net->boxHalfs[j % 3]) {
          distances[j] += net->L[j % 3];
          iterations++;
          if (iterations > 10) {
            throw std::runtime_error(
              "Too many iterations in PBC from " + std::to_string(linkIndexA) +
              " to " + std::to_string(linkIndexB) + ", currently at " +
              std::to_string(distances[j]));
          }
        }
      }

      if (is2D) {
        distances[2] = 0.0;
      }

      return distances;
    }

    Eigen::VectorXd MEHPForceBalance::evaluateSpringDistances(
      const ForceBalanceNetwork* net,
      const Eigen::VectorXd& u,
      const bool is2D)
    {
      // first, the distances
      assert(u.size() == net->coordinates.size());
      Eigen::VectorXd springDistances =
        Eigen::VectorXd::Zero(3 * net->nrOfSprings);
      Eigen::VectorXd actualCoordinates = net->coordinates + u;

      for (size_t i = 0; i < net->nrOfSprings; ++i) {
        std::vector<size_t> springsPartners = net->linkIndicesOfSprings[i];
        for (size_t partner_idx = 0; partner_idx < springsPartners.size() - 1;
             ++partner_idx) {
          // add partial distance to the total distance
          Eigen::Vector3d partialDistance =
            MEHPForceBalance::evaluateDistanceBetween(
              net,
              u,
              springsPartners[partner_idx],
              springsPartners[partner_idx + 1],
              is2D);
          for (size_t distance_idx = 0; distance_idx < 3; ++distance_idx) {
            springDistances[i * 3 + distance_idx] +=
              partialDistance[distance_idx];
            assert(!isnan(partialDistance[distance_idx]));
          }
        }
      }

      if (is2D) {
        // springDistances(Eigen::seq(2, Eigen::last, Eigen::fix<3>)) =
        //   Eigen::VectorXd::Zero(net->nrOfSprings / 3);
        for (size_t i = 2; i < 3 * net->nrOfSprings; i += 3) {
          springDistances[i] = 0.0;
        }
      }
      assert(springDistances.size() == net->nrOfSprings * 3);

      return springDistances;
    }

    /**
     * FORCE RELAXATION DATA ACCESS
     */
    pylimer_tools::entities::Universe MEHPForceBalance::getCrosslinkerVerse(
      int newCrosslinkerType) const
    {
      // convert nodes & springs back to a universe
      pylimer_tools::entities::Universe xlinkUniverse =
        pylimer_tools::entities::Universe(this->universe.getBox());
      std::vector<long int> ids;
      std::vector<int> types = pylimer_tools::utils::initializeWithValue(
        this->initialConfig.nrOfNodes, crosslinkerType);
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> z;
      std::vector<int> zeros = pylimer_tools::utils::initializeWithValue(
        this->initialConfig.nrOfNodes, 0);
      ids.reserve(this->initialConfig.nrOfNodes);
      x.reserve(this->initialConfig.nrOfNodes);
      y.reserve(this->initialConfig.nrOfNodes);
      z.reserve(this->initialConfig.nrOfNodes);
      for (int i = 0; i < this->initialConfig.nrOfNodes; ++i) {
        x.push_back(this->initialConfig.coordinates[3 * i + 0] +
                    this->currentDisplacements[3 * i + 0]);
        y.push_back(this->initialConfig.coordinates[3 * i + 1] +
                    this->currentDisplacements[3 * i + 1]);
        z.push_back(this->initialConfig.coordinates[3 * i + 2] +
                    this->currentDisplacements[3 * i + 2]);
        ids.push_back(this->initialConfig.oldAtomIds[i]);
      }
      xlinkUniverse.addAtoms(ids, types, x, y, z, zeros, zeros, zeros);
      std::vector<long int> bondFrom;
      std::vector<long int> bondTo;
      bondFrom.reserve(this->initialConfig.nrOfSprings);
      bondTo.reserve(this->initialConfig.nrOfSprings);
      for (int i = 0; i < this->initialConfig.nrOfSprings; ++i) {
        bondFrom.push_back(
          this->initialConfig.oldAtomIds[this->initialConfig.springIndexA[i]]);
        bondTo.push_back(
          this->initialConfig.oldAtomIds[this->initialConfig.springIndexB[i]]);
      }
      xlinkUniverse.addBonds(
        bondFrom.size(),
        bondFrom,
        bondTo,
        pylimer_tools::utils::initializeWithValue(bondFrom.size(), 1),
        false,
        false); // disable simplify to keep the self-loops etc.
      return xlinkUniverse;
    }

    void MEHPForceBalance::addSlipLinks(const std::vector<size_t> strandIdx1,
                                        const std::vector<size_t> strandIdx2,
                                        const std::vector<double> x,
                                        const std::vector<double> y,
                                        const std::vector<double> z,
                                        const std::vector<double> alpha)
    {
      size_t additionalLen = strandIdx1.size();
      size_t currentNrOfLinks = this->initialConfig.nrOfLinks;
      if (additionalLen != x.size() || additionalLen != y.size() ||
          additionalLen != z.size()) {
        throw std::invalid_argument("x, y and z must have the same dimensions");
      }
      if (additionalLen != strandIdx2.size() || additionalLen != alpha.size()) {
        throw std::invalid_argument(
          "Strand indices and alpha estimates must have the same length");
      }
      // actually start adding them
      this->initialConfig.nrOfLinks += additionalLen;
      // but first, indicate the resize
      this->initialConfig.coordinates.conservativeResize(
        3 * this->initialConfig.nrOfLinks);
      this->initialConfig.springIndicesOfLinks.reserve(
        this->initialConfig.nrOfLinks);
      this->initialConfig.linkIsSliplink.conservativeResize(
        this->initialConfig.nrOfLinks);
      for (size_t i = 0; i < additionalLen; ++i) {
        this->initialConfig.coordinates[3 * currentNrOfLinks + 3 * i] = x[i];
        this->initialConfig.coordinates[3 * currentNrOfLinks + 3 * i + 1] =
          y[i];
        this->initialConfig.coordinates[3 * currentNrOfLinks + 3 * i + 2] =
          z[i];
        this->initialConfig.linkIsSliplink[currentNrOfLinks + i] = true;
        std::vector<size_t> springIndices{ strandIdx1[i], strandIdx2[i] };
        this->initialConfig.springIndicesOfLinks[currentNrOfLinks + i] =
          springIndices;
        // add to the springs
        for (size_t springIndex : springIndices) {
          // detect the position in the spring
          std::vector<double> partitionsStrand =
            this->initialConfig.springPartitions[springIndex];
          bool wasAdded = false;
          size_t targetIndexInSpring = 0;
          for (size_t p_idx = 0; p_idx < partitionsStrand.size(); ++p_idx) {
            if (partitionsStrand[p_idx] > alpha[i]) {
              targetIndexInSpring = p_idx;
              wasAdded = true;
              break;
            }
          }
          if (!wasAdded) {
            targetIndexInSpring =
              this->initialConfig.springPartitions[springIndex].size();
          }
          this->initialConfig.springPartitions[springIndex].insert(
            this->initialConfig.springPartitions[springIndex].begin() +
              targetIndexInSpring,
            currentNrOfLinks + i);
          this->initialConfig.linkIndicesOfSprings[springIndex].insert(
            this->initialConfig.linkIndicesOfSprings[springIndex].begin() +
              targetIndexInSpring +
              1, // + 1 to compensate for the first cross-link
            currentNrOfLinks + i);
        }
      }
    };

    /**
     * @brief Get the Average Spring Length at the current step
     *
     * @return double
     */
    double MEHPForceBalance::getAverageSpringLength() const
    {
      double r2 = 0.0;
      for (int i = 0; i < this->initialConfig.nrOfSprings; i++) {
        double r2local = 0.0;
        for (int j = 0; j < 3; ++j) {
          r2local += this->currentSpringDistances[i * 3 + j] *
                     this->currentSpringDistances[i * 3 + j];
        }
        r2 += sqrt(r2local);
      }
      return r2 / this->initialConfig.nrOfSprings;
    }

    /**
     * @brief Compute the stress tensor
     *
     * @param net
     * @param u
     * @return std::array<std::array<double, 3>, 3>
     */
    std::array<std::array<double, 3>, 3> MEHPForceBalance::evaluateStressTensor(
      const Eigen::VectorXd& springDistances,
      const double volume) const
    {
      std::array<std::array<double, 3>, 3> stress;

      for (size_t i = 0; i < springDistances.size() / 3; ++i) {
        double s[3] = { springDistances[3 * i + 0],
                        springDistances[3 * i + 1],
                        springDistances[3 * i + 2] };
        /* spring contribution to the overall stress tensor */
        for (size_t j = 0; j < 3; j++) {
          for (size_t k = 0; k < 3; k++) {
            double contribution = this->kappa * s[j] * s[k];
            stress[j][k] += contribution;
          }
        }
      }

      for (size_t j = 0; j < 3; j++) {
        for (size_t k = 0; k < 3; k++) {
          stress[j][k] /= volume;
        }
      }

      return stress;
    }

    /**
     * @brief Compute the stress tensor
     *
     * @param net
     * @param u
     * @param loopTol
     * @return std::array<std::array<double, 3>, 3>
     */
    std::array<std::array<double, 3>, 3> MEHPForceBalance::evaluateStressTensor(
      ForceBalanceNetwork* net,
      const Eigen::VectorXd& u,
      const double loopTol) const
    {
      Eigen::VectorXd springDistances =
        this->evaluateSpringDistances(net, u, this->is2D);

      return this->evaluateStressTensor(springDistances, net->vol);
    }

    std::array<std::array<double, 3>, 3> MEHPForceBalance::getStressTensor()
      const
    {
      return this->evaluateStressTensor(this->currentSpringDistances,
                                        this->initialConfig.vol);
    }

    /**
     * @brief Get the Effective Functionality Of each node
     *
     * Returns the number of active springs connected to each atom, atomId
     * used as index
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return std::unordered_map<long int, int>
     */
    std::unordered_map<long int, int>
    MEHPForceBalance::getEffectiveFunctionalityOfAtoms(double tolerance) const
    {
      std::unordered_map<long int, int> results;
      results.reserve(this->initialConfig.nrOfNodes);

      Eigen::VectorXi nrOfActiveSpringsConnected =
        this->getNrOfActiveSpringsConnected(tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfNodes; i++) {
        results.emplace(this->initialConfig.oldAtomIds[i],
                        nrOfActiveSpringsConnected[i]);
      }
      return results;
    }

    /**
     * @brief Get the Ids Of active Nodes
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @param minimumNrOfActiveConnections the number of active springs
     * required for this node to qualify as active
     * @return std::vector<long int> the atom ids
     */
    std::vector<long int> MEHPForceBalance::getIdsOfActiveNodes(
      double tolerance,
      int minimumNrOfActiveConnections,
      int maximumNrOfActiveConnections) const
    {
      std::vector<long int> results;
      results.reserve(this->initialConfig.nrOfNodes);

      Eigen::VectorXi nrOfActiveSpringsConnected =
        this->getNrOfActiveSpringsConnected(tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfNodes; i++) {
        if (nrOfActiveSpringsConnected[i] >= minimumNrOfActiveConnections &&
            (maximumNrOfActiveConnections < 0 ||
             maximumNrOfActiveConnections >= nrOfActiveSpringsConnected[i])) {
          results.push_back(this->initialConfig.oldAtomIds[i]);
        }
      }

      return results;
    }

    /**
     * @brief Get the Nr Of Active Springs connected to each node
     *
     * @param tolerance the tolerance: springs under a certain length are
     * considered inactive
     * @return Eigen::VectorXi
     */
    Eigen::VectorXi MEHPForceBalance::getNrOfActiveSpringsConnected(
      double tolerance) const
    {
      Eigen::VectorXi nrOfActiveSpringsConnected =
        Eigen::VectorXi::Zero(this->initialConfig.nrOfNodes);
      ArrayXb springIsActive =
        this->findActiveSprings(this->currentSpringDistances, tolerance);
      for (size_t i = 0; i < this->initialConfig.nrOfSprings; i++) {
        if (springIsActive[i] == true) /* active spring */
        {
          int a = this->initialConfig.springIndexA[i];
          int b = this->initialConfig.springIndexB[i];
          ++(nrOfActiveSpringsConnected[a]);
          ++(nrOfActiveSpringsConnected[b]);
        }
      }
      return nrOfActiveSpringsConnected;
    }

    /**
     * @brief Get the Gamma Factor at the current step
     *
     * @param r02 the melt <R_0^2>, for phantom = Nb^2
     * @param nrOfChains the nr of chains to average over (can be different
     * from the nr of springs thanks to omitted free chains or primary loops)
     * @return double
     */
    double MEHPForceBalance::getGammaFactor(double r02, int nrOfChains) const
    {
      if (r02 < 0) {
        r02 = this->defaultR0Squared;
      }
      if (nrOfChains < 1) {
        nrOfChains = this->defaultNrOfChains;
      }

      return this->evaluateGammaFactor(
        this->currentSpringDistances, r02, nrOfChains);
    }

    bool MEHPForceBalance::validateNetwork(const ForceBalanceNetwork* net)
    {
      if (net == nullptr) {
        net = &this->initialConfig;
      }
      RUNTIME_EXP_IFN(net->coordinates.size() == net->nrOfLinks * 3,
                      "Invalid size of coordinates");
      RUNTIME_EXP_IFN(net->springsContourLength.size() == net->nrOfSprings,
                      "Invalid size of contour lengths");
      RUNTIME_EXP_IFN(net->springIndicesOfLinks.size() == net->nrOfLinks,
                      "Invalid size of spring indices of links");
      RUNTIME_EXP_IFN(net->linkIndicesOfSprings.size() == net->nrOfSprings,
                      "Invalid size of link indices of springs");
      RUNTIME_EXP_IFN(net->linkIsSliplink.size() == net->nrOfLinks,
                      "Invalid size of link is sliplink");
      RUNTIME_EXP_IFN(net->oldAtomIds.size() == net->nrOfNodes,
                      "Invalid size of old atom ids");
      RUNTIME_EXP_IFN(net->springCoordinateIndexA.size() ==
                        net->nrOfSprings * 3,
                      "Invalid size of springCoordinateIndexA");
      RUNTIME_EXP_IFN(net->springCoordinateIndexB.size() ==
                        net->nrOfSprings * 3,
                      "Invalid size of springCoordinateIndexB");
      RUNTIME_EXP_IFN(net->springIndexA.size() == net->nrOfSprings,
                      "Invalid size of springIndexA");
      RUNTIME_EXP_IFN(net->springIndexB.size() == net->nrOfSprings,
                      "Invalid size of springIndexB");
      RUNTIME_EXP_IFN(net->springIsActive.size() == net->nrOfSprings,
                      "Invalid size of springIsActive");
      for (size_t i = 0; i < net->nrOfSprings; ++i) {
        for (size_t j = 0; j < net->springPartitions[i].size(); ++j) {
          RUNTIME_EXP_IFN(net->springPartitions[i][j] <= 1. &&
                            net->springPartitions[i][j] >= 0.,
                          "Spring partitions must be between 0 & 1");
          if (net->springPartitions[i].size() >= 2 && j < net->springPartitions[i].size() - 1) {
            RUNTIME_EXP_IFN(net->springPartitions[i][j] <=
                              net->springPartitions[i][j + 1],
                            "Spring partitions must be sequential");
          }
        }

        RUNTIME_EXP_IFN(net->linkIndicesOfSprings[i].size() >= 2,
                        "Each spring requires at least two links");
        RUNTIME_EXP_IFN(
          net->linkIndicesOfSprings[i].size() ==
            net->springPartitions[i].size() + 2,
          "Spring partitions must coincide with nr of participating links");
        RUNTIME_EXP_IFN(
          net->linkIndicesOfSprings[i][0] <
            net->linkIndicesOfSprings[i]
                                     [net->linkIndicesOfSprings[i].size() - 1],
          "Springs must have increasing end-point indices");
      }
      return true;
    }
  }
}
}
