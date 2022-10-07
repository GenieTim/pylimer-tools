#ifndef MEHP_FORCE_BALANCE_H
#define MEHP_FORCE_BALANCE_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "MEHPForceEvaluator.h"
#include "MEHPUtilityStructures.h"
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

    // heavily inspired by Prof. Dr. Andrei Gusev's Code
    class MEHPForceBalance
    {

    public:
      MEHPForceBalance(const pylimer_tools::entities::Universe u,
                       int crosslinkerType = 2,
                       bool is2D = false,
                       double kappa = 1.0)
        : universe(u)
      {
        this->crosslinkerType = crosslinkerType;
        // interpret network already to be able to give early results
        ForceBalanceNetwork net;
        ConvertNetwork(&net, crosslinkerType);
        this->initialConfig = net;
        this->is2D = is2D;
        this->currentDisplacements =
          Eigen::VectorXd::Zero(net.coordinates.size());
        this->currentSpringDistances =
          this->evaluateSpringDistances(&net, this->currentDisplacements, is2D);
        this->defaultR0Squared =
          universe.computeMeanSquareEndToEndDistance(crosslinkerType);
        this->defaultNrOfChains =
          universe.getMolecules(this->crosslinkerType).size();
      };

      /**
       * @brief Actually do run the simulation
       *
       * @param algorithm
       * @param maxNrOfSteps
       * @param xtol
       * @param ftol
       */
      void runForceRelaxation(long int maxNrOfSteps = 50000, // default: 10000
                              double xtol = 1e-12,
                              long int innerMaxNrOfSteps = 100,
                              double innerXtol = 1e-12,
                              double innerAlphaTol = 1e-8);

      /**
       * @brief Get the universe consisting of cross-linkers only
       *
       * @param newCrosslinkerType the type to give the cross-linkers
       * @return pylimer_tools::entities::Universe
       */
      pylimer_tools::entities::Universe getCrosslinkerVerse(
        int newCrosslinkerType = 2) const;

      int getDefaultNrOfChains() const { return this->defaultNrOfChains; }

      double getDefaultR0Square() const { return this->defaultR0Squared; }

      double getVolume() const { return this->initialConfig.vol; }

      int getNrOfNodes() const { return this->initialConfig.nrOfNodes; }

      int getNrOfLinks() const { return this->initialConfig.nrOfLinks; }

      int getNrOfSprings() const { return this->initialConfig.nrOfSprings; }

      /**
       * @brief Get the Nr Of Active Nodes
       *
       * @param tolerance  the tolerance: springs under a certain length are
       * considered inactive
       * @return int
       */
      int getNrOfActiveNodes(double tolerance = 0.1,
                             int minimumNrOfActiveConnections = 2,
                             int maximumNrOfActiveConnections = -1) const
      {
        return this
          ->getIdsOfActiveNodes(tolerance,
                                minimumNrOfActiveConnections,
                                maximumNrOfActiveConnections)
          .size();
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
      std::unordered_map<long int, int> getEffectiveFunctionalityOfAtoms(
        double tolerance = 0.1) const;

      /**
       * @brief Get the Ids Of active Nodes
       *
       * @param tolerance the tolerance: springs under a certain length are
       * considered inactive
       * @param minimumNrOfActiveConnections the number of active springs
       * required for this node to qualify as active
       * @return std::vector<long int> the atom ids
       */
      std::vector<long int> getIdsOfActiveNodes(
        double tolerance = 0.1,
        int minimumNrOfActiveConnections = 2,
        int maximumNrOfActiveConnections = -1) const;

      /**
       * @brief Get the Nr Of Active Springs connected to each node
       *
       * @param tolerance the tolerance: springs under a certain length are
       * considered inactive
       * @return Eigen::VectorXi
       */
      Eigen::VectorXi getNrOfActiveSpringsConnected(
        double tolerance = 0.1) const;

      Eigen::VectorXd getCurrentSpringDistances() const
      {
        return this->currentSpringDistances;
      }

      /**
       * @brief Get the Nr Of Active Springs object
       *
       * @param tol the tolerance: springs under a certain length are considered
       * inactive
       * @return int
       */
      int getNrOfActiveSprings(double tol = 0.1) const
      {
        return this->countNrOfActiveSprings(this->currentSpringDistances, tol);
      }

      /**
       * @brief Get the Average Spring Length at the current step
       *
       * @return double
       */
      double getAverageSpringLength() const;

      std::array<std::array<double, 3>, 3> getStressTensor() const;

      /**
       * @brief Get the Pressure
       *
       * @return double
       */
      double getPressure() const
      {
        return this->evaluatePressure(this->currentSpringDistances);
      }

      /**
       * @brief Get the Gamma Factor at the current step
       *
       * @param r02 the melt <R_0^2>, for phantom = Nb^2
       * @param nrOfChains the nr of chains to average over (can be different
       * from the nr of springs thanks to omitted free chains or primary loops)
       * @return double
       */
      double getGammaFactor(double r02 = -1.0, int nrOfChains = -1) const;

      int getNrOfIterations() const { return this->nrOfStepsDone; }

      ExitReason getExitReason() const { return this->exitReason; }

      void addSlipLinks(const std::vector<size_t> strandIdx1,
                        const std::vector<size_t> strandIdx2,
                        const std::vector<double> x,
                        const std::vector<double> y,
                        const std::vector<double> z)
      {
        std::vector<double> alphas;
        alphas.reserve(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
          alphas.push_back(0.5);
        }
        return this->addSlipLinks(strandIdx1, strandIdx2, x, y, z, alphas);
      }

      void addSlipLinks(const std::vector<size_t> strandIdx1,
                        const std::vector<size_t> strandIdx2,
                        const std::vector<double> x,
                        const std::vector<double> y,
                        const std::vector<double> z,
                        const std::vector<double> alpha);

      /**
       * @brief Compute the spring lenghts
       *
       * @param net the network to do the computation for
       * @param u the displacements on top of the network
       * @return Eigen::VectorXd
       */
      static Eigen::VectorXd evaluateSpringDistances(
        const ForceBalanceNetwork* net,
        const Eigen::VectorXd& u,
        const bool is2D);

      /**
       * @brief Compute one spring length
       *
       * @param net
       * @param linkIndexA
       * @param linkIndexB
       * @param is2D
       * @return Eigen::Vector3d
       */
      static Eigen::Vector3d evaluateDistanceBetween(
        const ForceBalanceNetwork* net,
        const Eigen::VectorXd& u,
        const size_t linkIndexA,
        const size_t linkIndexB,
        const bool is2D);

      bool validateNetwork(const ForceBalanceNetwork* net = nullptr);

    protected:
      /**
       * @brief Convert the universe to a network
       *
       * @param net the target network
       * @param crosslinkerType the atom type of the crosslinker
       * @return true
       * @return false
       */
      bool ConvertNetwork(ForceBalanceNetwork* net,
                          const int crosslinkerType = 2)
      {
        std::vector<pylimer_tools::entities::Atom> xlinkers =
          this->universe.getAtomsOfType(crosslinkerType);

        size_t nrOfXlinks = xlinkers.size();

        std::vector<pylimer_tools::entities::Molecule> crosslinkerChains =
          this->universe.getChainsWithCrosslinker(crosslinkerType);

        // need to include all but dangling and free chains in order to
        // model entanglement
        size_t nrOfSprings = 0;
        for (size_t i = 0; i < crosslinkerChains.size(); ++i) {
          std::vector<pylimer_tools::entities::Atom> endAtoms =
            crosslinkerChains[i].getAtomsOfType(crosslinkerType);
          if (crosslinkerChains[i].getType() ==
                pylimer_tools::entities::MoleculeType::NETWORK_STRAND ||
              crosslinkerChains[i].getType() ==
                pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
            nrOfSprings += 1;
          }
        }

        // crosslinkerUniverse.simplify();
        pylimer_tools::entities::Box box = this->universe.getBox();
        net->L[0] = box.getLx();
        net->L[1] = box.getLy();
        net->L[2] = box.getLz();
        net->boxHalfs[0] = 0.5 * net->L[0];
        net->boxHalfs[1] = 0.5 * net->L[1];
        net->boxHalfs[2] = 0.5 * net->L[2];
        net->nrOfNodes = nrOfXlinks;
        net->nrOfLinks = nrOfXlinks;
        net->nrOfSprings = nrOfSprings;
        net->coordinates = Eigen::VectorXd::Zero(3 * net->nrOfLinks);
        net->oldAtomIds = Eigen::ArrayXi::Zero(net->nrOfLinks);
        net->linkIsSliplink = ArrayXb::Constant(net->nrOfLinks, false);
        net->springIndicesOfLinks.reserve(net->nrOfLinks);
        for (size_t i = 0; i < net->nrOfLinks; ++i) {
          net->springIndicesOfLinks.push_back(std::vector<size_t>());
        }
        net->linkIndicesOfSprings.reserve(net->nrOfSprings);
        net->springPartitions.reserve(net->nrOfSprings);
        for (size_t i = 0; i < net->nrOfSprings; ++i) {
          net->linkIndicesOfSprings.push_back(std::vector<size_t>());
          net->springPartitions.push_back(std::vector<double>());
        }
        net->springIndexA = Eigen::ArrayXi::Zero(net->nrOfSprings);
        net->springIndexB = Eigen::ArrayXi::Zero(net->nrOfSprings);
        net->springCoordinateIndexA =
          Eigen::ArrayXi::Zero(3 * net->nrOfSprings);
        net->springCoordinateIndexB =
          Eigen::ArrayXi::Zero(3 * net->nrOfSprings);
        net->springIsActive = ArrayXb::Constant(net->nrOfSprings, false);
        net->springsContourLength = Eigen::VectorXd::Zero(net->nrOfSprings);

        // convert beads
        std::map<int, int> atomIdToNode;
        for (size_t i = 0; i < xlinkers.size(); ++i) {
          pylimer_tools::entities::Atom atom = xlinkers[i];
          atomIdToNode[atom.getId()] = i;
          net->oldAtomIds[i] = atom.getId();
          net->coordinates[3 * i + 0] = atom.getX();
          net->coordinates[3 * i + 1] = atom.getY();
          net->coordinates[3 * i + 2] = atom.getZ();
        }

        // convert springs
        size_t spring_idx = 0;
        for (size_t i = 0; i < crosslinkerChains.size(); ++i) {
          std::vector<pylimer_tools::entities::Atom> xlinkersOfChain =
            crosslinkerChains[i].getAtomsOfType(crosslinkerType);
          long int nodeIdxFrom;
          long int nodeIdxTo;
          bool addChain = false;
          if (crosslinkerChains[i].getType() ==
              pylimer_tools::entities::MoleculeType::NETWORK_STRAND) {
            assert(xlinkersOfChain.size() == 2);
            nodeIdxFrom = atomIdToNode.at(xlinkersOfChain[0].getId());
            nodeIdxTo = atomIdToNode.at(xlinkersOfChain[1].getId());
            if (nodeIdxFrom > nodeIdxTo) {
              std::swap(nodeIdxFrom, nodeIdxTo);
            }
            addChain = true;
          } else if (crosslinkerChains[i].getType() ==
                     pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
            assert(xlinkersOfChain.size() == 1 ||
                   (xlinkersOfChain.size() == 2 &&
                    xlinkersOfChain[0].getId() == xlinkersOfChain[1].getId()));

            nodeIdxFrom = atomIdToNode.at(xlinkersOfChain[0].getId());
            nodeIdxTo = nodeIdxFrom;
            addChain = true;
          }

          if (addChain) {
            net->springIndicesOfLinks[nodeIdxFrom].push_back(spring_idx);
            net->springIndicesOfLinks[nodeIdxTo].push_back(spring_idx);

            net->linkIndicesOfSprings[spring_idx].push_back(nodeIdxFrom);
            net->linkIndicesOfSprings[spring_idx].push_back(nodeIdxTo);

            net->springIndexA[spring_idx] = nodeIdxFrom;
            net->springIndexB[spring_idx] = nodeIdxTo;
            for (size_t j = 0; j < 3; j++) {
              net->springCoordinateIndexA[3 * spring_idx + j] =
                nodeIdxFrom * 3 + j;
              net->springCoordinateIndexB[3 * spring_idx + j] =
                nodeIdxTo * 3 + j;
            }

            net->springsContourLength[spring_idx] =
              crosslinkerChains[i].getNrOfAtoms(); // TODO: -2?
            spring_idx += 1;
          }
        }

        // box volume
        net->vol = net->L[0] * net->L[1] * net->L[2];
        net->meanSpringContourLength = net->springsContourLength.mean();

        return spring_idx == net->nrOfSprings;
      };

      double updateSpringPartition(ForceBalanceNetwork* net,
                                   Eigen::VectorXd& u,
                                   const size_t linkIdx) const
      {
        std::vector<size_t> springIndices = net->springIndicesOfLinks[linkIdx];
        double maxDiff = 0.0;
        for (size_t springIndex : springIndices) {
          std::vector<size_t> springsPartners =
            net->linkIndicesOfSprings[springIndices[springIndex]];
          for (size_t partner_idx = 1; partner_idx < springsPartners.size() - 1;
               ++partner_idx) {
            if (springsPartners[partner_idx] == linkIdx) {
              // found position of this link in this spring
              // want to find the ideal value for
              // net->springPartitions[springIndex][partner_idx-1]
              double distanceBack = ((MEHPForceBalance::evaluateDistanceBetween(
                                        net,
                                        u,
                                        springsPartners[partner_idx - 1],
                                        springsPartners[partner_idx],
                                        this->is2D))
                                       .squaredNorm());
              double distanceForward =
                ((MEHPForceBalance::evaluateDistanceBetween(
                    net,
                    u,
                    springsPartners[partner_idx],
                    springsPartners[partner_idx + 1],
                    this->is2D))
                   .squaredNorm());
              double idealValue =
                1. / (1. + sqrt(distanceBack / distanceForward));
              maxDiff =
                std::max(net->springPartitions[springIndex][partner_idx - 1] -
                           idealValue,
                         maxDiff);
              net->springPartitions[springIndex][partner_idx - 1] =
                partner_idx == 1
                  ? idealValue
                  : idealValue +
                      net->springPartitions[springIndex][partner_idx - 2];
            }
          }
        }
        return maxDiff;
      }

      /**
       * @brief Displace one link to the mean of all connected neighbours
       *
       * @param net the force balance network
       * @param u the current displacements, wherein the resulting coordinates
       * shall be stored
       * @param linkIdx the idx of the link to displace
       * @return double, the distance (squared norm) displaced
       */
      double displaceToMeanPosition(const ForceBalanceNetwork* net,
                                    Eigen::VectorXd& u,
                                    const size_t linkIdx) const
      {
        std::vector<size_t> springIndices = net->springIndicesOfLinks[linkIdx];
        Eigen::Vector3d currentDisplacement;
        currentDisplacement << u[3 * linkIdx], u[3 * linkIdx + 1],
          u[3 * linkIdx + 2];
        Eigen::Vector3d objectiveDisplacement =
          Eigen::Vector3d::Zero(); // = remainingDisplacement.array();
        double objectiveDisplacementContributors = 0.0;
        for (size_t spring_index = 0; spring_index < springIndices.size();
             ++spring_index) {
          // compute partial distances & total distance of this spring
          std::vector<size_t> springsPartners =
            net->linkIndicesOfSprings[springIndices[spring_index]];
          for (size_t partner_idx = 0; partner_idx < springsPartners.size() - 1;
               ++partner_idx) {
            if (springsPartners[partner_idx] == linkIdx ||
                springsPartners[partner_idx + 1] == linkIdx) {
              // add partial distance to the total distance
              Eigen::Vector3d partialDistance =
                MEHPForceBalance::evaluateDistanceBetween(
                  net,
                  u,
                  springsPartners[partner_idx],
                  springsPartners[partner_idx + 1],
                  this->is2D);
              // add to displacement
              double prefix =
                springsPartners[partner_idx] == linkIdx ? -1. : 1.;
              double contourLengthFraction = 1.0;
              if (springsPartners.size() > 2) {
                // TODO: this must be done more beautiful
                if (partner_idx == 0) {
                  // only the case if linkIdx is the first cross-linker in the
                  // strand
                  contourLengthFraction =
                    net->springPartitions[springIndices[spring_index]][0];
                } else if (partner_idx == springsPartners.size() - 2) {
                  // only the case if linkIdx is the last cross-linker in the
                  // strand
                  assert(springsPartners[partner_idx + 1] == linkIdx);
                  contourLengthFraction =
                    (1.0 - net->springPartitions[springIndices[spring_index]]
                                                [partner_idx]);
                } else {
                  if (springsPartners[partner_idx] == linkIdx) {
                    // fraction towards the "lower" end
                    contourLengthFraction =
                      net->springPartitions[springIndices[spring_index]]
                                           [partner_idx] -
                      net->springPartitions[springIndices[spring_index]]
                                           [partner_idx - 1];

                  } else {
                    assert(springsPartners[partner_idx + 1] == linkIdx);
                    // fraction towards the "higher" end
                    contourLengthFraction =
                      net->springPartitions[springIndices[spring_index]]
                                           [partner_idx + 1] -
                      net->springPartitions[springIndices[spring_index]]
                                           [partner_idx];
                  }
                }
              }
              double oneOverContourLengthFraction =
                1.0 / (net->springsContourLength[springIndices[spring_index]] *
                       contourLengthFraction);
              objectiveDisplacement +=
                prefix *
                (partialDistance)*oneOverContourLengthFraction; // /
                                                                // totalDistance.array());
              objectiveDisplacementContributors += oneOverContourLengthFraction;
            }
          }
        }
        // take mean for displacement
        u[3 * linkIdx] +=
          objectiveDisplacement[0] / objectiveDisplacementContributors;
        u[3 * linkIdx + 1] +=
          objectiveDisplacement[1] / objectiveDisplacementContributors;
        u[3 * linkIdx + 2] +=
          objectiveDisplacement[2] / objectiveDisplacementContributors;

        double dist = objectiveDisplacement.squaredNorm();
        // if (dist > 500) {
        //   std::cout << "Moving " << linkIdx << " for " << dist
        //             << " with displacements " << currentDisplacement[0] << ",
        //             "
        //             << currentDisplacement[1] << ", " <<
        //             currentDisplacement[2]
        //             << std::endl;
        //   std::cout << "For objective displacements "
        //             << objectiveDisplacement[0] << ", "
        //             << objectiveDisplacement[1] << ", "
        //             << objectiveDisplacement[2] << ", for "
        //             << objectiveDisplacementContributors << "." << std::endl;
        // }
        return dist;
      }

      /**
       * @brief Compute the gamma factor from certain spring distances
       *
       * @param springDistances
       * @param r02 the melt <R_0^2>, for phantom = Nb^2
       * @param nrOfChains the nr of chains to average over (can be different
       * from the nr of springs thanks to omitted free chains or primary loops)
       * @return double
       */
      double evaluateGammaFactor(const Eigen::VectorXd& springDistances,
                                 double r02,
                                 int nrOfChains) const
      {
        return springDistances.squaredNorm() /
               (static_cast<double>(nrOfChains) * r02);
      }

      /**
       * @brief Evaluate the pressure of the network at specific spring
       * distances
       *
       * @param springDistances the spring distances
       * @return double
       */
      double evaluatePressure(const Eigen::VectorXd& springDistances) const
      {
        auto stressTensor =
          this->evaluateStressTensor(springDistances, this->initialConfig.vol);
        return this->evaluatePressure(stressTensor);
      }

      /**
       * @brief Evaluate the pressure of the network at specific displacements
       *
       * @param net the network to evaluate the pressure for
       * @param u the displacements
       * @return double
       */
      double evaluatePressure(ForceBalanceNetwork* net,
                              const Eigen::VectorXd& u) const
      {
        auto stressTensor = this->evaluateStressTensor(net, u, -1);
        return this->evaluatePressure(stressTensor);
      }

      /**
       * @brief Evaluate the pressure from the stress tensor
       *
       * @param stressTensor
       * @return double
       */
      double evaluatePressure(
        const std::array<std::array<double, 3>, 3>& stressTensor) const
      {
        return (stressTensor[0][0] + stressTensor[1][1] + stressTensor[2][2]) /
               3.0;
      }

      /**
       * @brief Compute the stress tensor
       *
       * @param net
       * @param u
       * @return std::array<std::array<double, 3>, 3>
       */
      std::array<std::array<double, 3>, 3> evaluateStressTensor(
        const Eigen::VectorXd& springDistances,
        const double volume) const;

      /**
       * @brief Compute the stress tensor
       *
       * @param net
       * @param u
       * @param loopTol
       * @return std::array<std::array<double, 3>, 3>
       */
      std::array<std::array<double, 3>, 3> evaluateStressTensor(
        ForceBalanceNetwork* net,
        const Eigen::VectorXd& u,
        const double loopTol) const;

      /**
       * @brief Count how many of the springs are active (length > tolerance)
       *
       * @param springDistances
       * @param tolerance
       * @return int
       */
      int countNrOfActiveSprings(const Eigen::VectorXd& springDistances,
                                 const double tolerance = 0.1) const
      {
        return (this->findActiveSprings(springDistances, tolerance) == true)
          .count();
      }

      /**
       * @brief Iterate all spring distances, mark active ones (length >
       * tolerance)
       *
       * @param springDistances
       * @param tolerance
       * @return ArrayXb
       */
      ArrayXb findActiveSprings(const Eigen::VectorXd& springDistances,
                                const double tolerance = 0.1) const
      {
        ArrayXb result = ArrayXb::Constant(springDistances.size() / 3, false);
        for (size_t i = 0; i < springDistances.size() / 3; ++i) {
          result[i] =
            sqrt(springDistances[3 * i + 0] * springDistances[3 * i + 0] +
                 springDistances[3 * i + 1] * springDistances[3 * i + 1] +
                 springDistances[3 * i + 2] * springDistances[3 * i + 2]) >
            tolerance;
        }
        return result;
      }

    private:
      pylimer_tools::entities::Universe universe;
      bool is2D = false;
      double kappa = 1.0;
      bool simulationHasRun = false;
      int stepOutputFrequency = 0;
      int defaultNrOfChains = 0;
      double defaultR0Squared = 0.0;
      std::string stepOutputFile;
      bool outputEndNodes = false;
      std::string endNodesFile;
      ForceBalanceNetwork initialConfig;
      ForceBalanceNetwork finalConfig;
      Eigen::VectorXd currentDisplacements;
      Eigen::VectorXd currentSpringDistances;
      int crosslinkerType;
      int nrOfStepsDone = 0;
      ExitReason exitReason = ExitReason::UNSET;
    };
  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
