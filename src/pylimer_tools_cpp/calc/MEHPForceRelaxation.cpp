#include "MEHPForceRelaxation.h"
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
    void MEHPForceRelaxation::runForceRelaxation(
      const char* algorithm,
      long int maxNrOfSteps, // default: 10000
      double xtol,
      double ftol)
    {
      this->simulationHasRun = true;
      this->forceEvaluator->setNetwork(this->initialConfig);
      this->forceEvaluator->setIs2D(this->is2D);
      this->forceEvaluator->prepareForEvaluations();
      double stress[3][3];

      for (size_t j = 0; j < 3; j++) {
        for (size_t k = 0; k < 3; k++) {
          stress[j][k] = 0.;
        }
      }

      Network net = this->initialConfig;
      const int M = this->universe.getMolecules(crosslinkerType).size();
      const int N = this->universe.getMeanStrandLength(crosslinkerType) + 1;
      const double bM = this->universe.computeMeanBondLength();
      const int f =
        this->universe.determineFunctionalityPerType()[crosslinkerType];
      bool is2D = this->is2D;

      /* array allocation */
      std::vector<double> u0 =
        pylimer_tools::utils::initializeWithValue(3 * net.nrOfNodes, 0.0);
      Eigen::VectorXd u = Eigen::VectorXd::Zero(3 * net.nrOfNodes);

      /* force relaxation */
      nlopt::opt opt(algorithm, 3 * net.nrOfNodes);

      nlopt::func objectiveF =
        [](unsigned n, const double* x, double* grad, void* f_data) -> double {
        MEHPForceEvaluator* fEvaluator =
          static_cast<MEHPForceEvaluator*>(f_data);
        return fEvaluator->evaluateForceSetGradient(n, x, grad, f_data);
      };
      opt.set_min_objective(objectiveF, this->forceEvaluator);
      // set constraints to support more algorithms
      std::vector<double> upperBounds;
      upperBounds.reserve(3 * net.nrOfNodes);
      std::vector<double> lowerBounds;
      lowerBounds.reserve(3 * net.nrOfNodes);
      for (size_t i = 0; i < net.nrOfNodes; ++i) {
        for (size_t dir = 0; dir < 3; ++dir) {
          upperBounds.push_back(net.L[dir] * 0.5);
          lowerBounds.push_back(-net.L[dir] * 0.5);
        }
      }
      opt.set_upper_bounds(upperBounds);
      opt.set_lower_bounds(lowerBounds);
      // set exit conditions
      opt.set_xtol_rel(xtol);
      opt.set_ftol_rel(ftol);
      opt.set_ftol_abs(0.0);
      opt.set_maxeval(maxNrOfSteps);
      // opt.set_param("verbosity", 1.0);
      // start/set/run minimization
      double minf;
      nlopt::result res;
      std::exception_ptr nloptException = nullptr;
      try {
        res = opt.optimize(u0, minf);
      } catch (...) {
        nloptException = std::current_exception();
      }

      // query solution & exit reason
      assert(u0.size() == 3 * net.nrOfNodes);
      u = Eigen::Map<Eigen::VectorXd>(u0.data(), u0.size());
      this->currentDisplacements = u;
      this->currentSpringDistances =
        this->evaluateSpringDistances(&net, this->currentDisplacements, is2D);

      this->exitReason = ExitReason::OTHER;
      if (nloptException != nullptr) {
        this->exitReason = ExitReason::FAILURE;
      } else if (res == nlopt::result::FTOL_REACHED) {
        this->exitReason = ExitReason::F_TOLERANCE;
      } else if (res == nlopt::result::XTOL_REACHED) {
        this->exitReason = ExitReason::X_TOLERANCE;
      } else if (res == nlopt::result::MAXEVAL_REACHED) {
        this->exitReason = ExitReason::MAX_STEPS;
      }
      this->nrOfStepsDone += opt.get_numevals();
    }

    Eigen::VectorXd MEHPForceRelaxation::evaluateSpringDistances(
      const Network* net,
      const Eigen::VectorXd& u,
      const bool is2D)
    {
      double boxHalfs[3];
      boxHalfs[0] = 0.5 * net->L[0];
      boxHalfs[1] = 0.5 * net->L[1];
      boxHalfs[2] = 0.5 * net->L[2];
      // first, the distances
      assert(u.size() == net->coordinates.size());
      Eigen::VectorXd actualCoordinates = net->coordinates + u;
      // It *could* be more efficient to index u instead of the coordinates
      Eigen::VectorXd coordinatesSpringEndA =
        actualCoordinates(net->springCoordinateIndexA);
      Eigen::VectorXd coordinatesSpringEndB =
        actualCoordinates(net->springCoordinateIndexB);
      Eigen::VectorXd springDistances =
        (coordinatesSpringEndA - coordinatesSpringEndB);

      if (is2D) {
        // springDistances(Eigen::seq(2, Eigen::last, Eigen::fix<3>)) =
        //   Eigen::VectorXd::Zero(net->nrOfSprings / 3);
        for (size_t i = 2; i < 3 * net->nrOfSprings; i += 3) {
          springDistances[i] = 0.0;
        }
      }
      assert(springDistances.size() == net->nrOfSprings * 3);

      // Possibly improvable PBC
      for (size_t j = 0; j < 3 * net->nrOfSprings; ++j) {
        int iterations = 0;
        while (springDistances[j] > boxHalfs[j % 3]) {
          springDistances[j] -= net->L[j % 3];
          iterations++;
          if (iterations > 10) {
            throw std::runtime_error(
              "Too many iterations in PBC from " +
              std::to_string(coordinatesSpringEndA[j]) + " to " +
              std::to_string(coordinatesSpringEndB[j]) + ", currently at " +
              std::to_string(springDistances[j]));
          }
        }
        iterations = 0;
        while (springDistances[j] < -boxHalfs[j % 3]) {
          springDistances[j] += net->L[j % 3];
          iterations++;
          if (iterations > 10) {
            throw std::runtime_error(
              "Too many iterations in PBC from " +
              std::to_string(coordinatesSpringEndA[j]) + " to " +
              std::to_string(coordinatesSpringEndB[j]) + ", currently at " +
              std::to_string(springDistances[j]));
          }
        }
      }

      return springDistances;
    }

    /**
     * FORCE RELAXATION DATA ACCESS
     */
    pylimer_tools::entities::Universe MEHPForceRelaxation::getCrosslinkerVerse(
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

    /**
     * @brief Get the Average Spring Length at the current step
     *
     * @return double
     */
    double MEHPForceRelaxation::getAverageSpringLength() const
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
    std::array<std::array<double, 3>, 3>
    MEHPForceRelaxation::evaluateStressTensor(
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
            double contribution = this->forceEvaluator->evaluateStressContribution(s, j, k);
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
    std::array<std::array<double, 3>, 3>
    MEHPForceRelaxation::evaluateStressTensor(Network* net,
                                              const Eigen::VectorXd& u,
                                              const double loopTol) const
    {
      Eigen::VectorXd springDistances =
        this->evaluateSpringDistances(net, u, this->is2D);

      return this->evaluateStressTensor(springDistances, net->vol);
    }

    std::array<std::array<double, 3>, 3> MEHPForceRelaxation::getStressTensor()
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
    MEHPForceRelaxation::getEffectiveFunctionalityOfAtoms(
      double tolerance) const
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
    std::vector<long int> MEHPForceRelaxation::getIdsOfActiveNodes(
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
    Eigen::VectorXi MEHPForceRelaxation::getNrOfActiveSpringsConnected(
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
     * @brief Get the residuals (gradient) at the current step
     *
     * @return Eigen::VectorXd
     */
    Eigen::VectorXd MEHPForceRelaxation::getResiduals() const
    {
      double* r = new double[3 * this->initialConfig.nrOfNodes];
      for (size_t i = 0; i < this->initialConfig.nrOfNodes * 3; ++i) {
        r[i] = 0.0;
      }
      try {
        this->forceEvaluator->evaluateForceSetGradient(
          3 * this->initialConfig.nrOfNodes,
          this->currentSpringDistances,
          this->currentDisplacements,
          r);
      } catch (const std::exception& e) {
        delete[](r);
        throw e;
      }

      Eigen::VectorXd results =
        Eigen::VectorXd::Zero(this->initialConfig.nrOfNodes * 3);
      for (size_t i = 0; i < this->initialConfig.nrOfNodes * 3; ++i) {
        results[i] = r[i];
      }
      delete[](r);
      return results;
    }

    /**
     * @brief Get the Residual Norm at the current step
     *
     * @return double
     */
    double MEHPForceRelaxation::getResidualNorm() const
    {
      return this->getResiduals().norm();
    }

    /**
     * @brief Get the Force at the current step
     *
     * @return double
     */
    double MEHPForceRelaxation::getForce() const
    {
      return this->forceEvaluator->evaluateForceSetGradient(
        3 * this->initialConfig.nrOfNodes,
        this->currentSpringDistances,
        this->currentDisplacements,
        nullptr);
    }

    /**
     * @brief Get the Gamma Factor at the current step
     *
     * @param r02 the melt <R_0^2>, for phantom = Nb^2
     * @param nrOfChains the nr of chains to average over (can be different
     * from the nr of springs thanks to omitted free chains or primary loops)
     * @return double
     */
    double MEHPForceRelaxation::getGammaFactor(double r02, int nrOfChains) const
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
  }
}
}
