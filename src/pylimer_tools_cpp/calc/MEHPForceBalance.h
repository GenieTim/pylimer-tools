#ifndef MEHP_FORCE_BALANCE_H
#define MEHP_FORCE_BALANCE_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/NeighbourList.h"
#include "../entities/Universe.h"
#include "MEHPForceEvaluator.h"
#include "MEHPUtilityStructures.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <nlopt.hpp>
#include <random>
#include <string>
#include <tuple>
#include <vector>

namespace pylimer_tools {
namespace calc {
  namespace mehp {

    enum BalanceRunMode
    {
      EIGEN_RANDOM,
      EIGEN_HEURISTIC,
      EIGEN_STRANDS,
      EIGEN_ALL,
      ITERATIVE
    };

    // heavily inspired by Prof. Dr. Andrei Gusev's Code
    class MEHPForceBalance
    {

    public:
      MEHPForceBalance(const pylimer_tools::entities::Universe u,
                       int crosslinkerType = 2,
                       bool is2D = false,
                       double kappa = 1.0,
                       bool remove2functionalCrosslinkers = false)
        : universe(u)
      {
        this->crosslinkerType = crosslinkerType;
        // interpret network already to be able to give early results
        ForceBalanceNetwork net;
        ConvertNetwork(net, crosslinkerType, remove2functionalCrosslinkers);
        this->initialConfig = net;
        this->is2D = is2D;
        this->currentDisplacements =
          Eigen::VectorXd::Zero(net.coordinates.size());
        this->currentSpringDistances =
          this->evaluateSpringDistances(net, this->currentDisplacements, is2D);
        this->currentPartialSpringDistances =
          this->evaluatePartialSpringDistances(
            net, this->currentDisplacements, is2D);
        this->defaultR0Squared =
          universe.computeMeanSquareEndToEndDistance(crosslinkerType);
        this->defaultNrOfChains =
          universe.getMolecules(this->crosslinkerType).size();
        this->validateNetwork();
      };

      /**
       * @brief Actually do run the simulation
       *
       * @param algorithm
       * @param maxNrOfSteps
       * @param xtol
       * @param ftol
       */
      void runForceRelaxation(
        BalanceRunMode mode = BalanceRunMode::ITERATIVE,
        double damping = 1.0,
        long int maxNrOfSteps = 50000, // default: 10000
        double xtol = 1e-9,
        long int innerMaxNrOfSteps = 100,
        double innerAlphaTol = 1e-15,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const int maxFlag = 7,
        const bool removeInactiveCrosslinks = false,
        const bool remove2functionalCrosslinkers = false);

      /**
       * @brief Compute the spring update residual
       *
       * @param link_idx
       * @param displacements
       * @param springPartitions
       * @param oneOverSpringPartitionUpperLimit
       * @return double
       */
      double computePartitionUpdateZeroResidual(
        const std::vector<size_t> involvedPartitions,
        const size_t link_idx,
        const Eigen::VectorXd& displacements,
        Eigen::VectorXd& springPartitions,
        double oneOverSpringPartitionUpperLimit = 1.0)
      {
        // TODO: revise, hard!
        assert(involvedPartitions.size() == 4);
        double firstMeanVal = 0.5 * (springPartitions[involvedPartitions[0]] +
                                     springPartitions[involvedPartitions[1]]);
        double secondMeanVal = 0.5 * (springPartitions[involvedPartitions[2]] +
                                      springPartitions[involvedPartitions[3]]);
        // Eigen::ArrayXi involvedCoordinateIndices = Eigen::ArrayXi(12);
        // for (size_t i = 0; i < 4; ++i) {
        //   involvedCoordinateIndices[3 * i] = 3 * involvedPartitions[i];
        //   involvedCoordinateIndices[3 * i + 1] = 3 * involvedPartitions[i] +
        //   1; involvedCoordinateIndices[3 * i + 2] = 3 * involvedPartitions[i]
        //   + 2;
        // }
        // Eigen::VectorXd displacementsBefore =
        //   displacements(involvedCoordinateIndices);
        Eigen::Vector4d partitionsBefore = springPartitions(involvedPartitions);
        springPartitions[involvedPartitions[0]] = firstMeanVal;
        springPartitions[involvedPartitions[1]] = firstMeanVal;
        springPartitions[involvedPartitions[2]] = secondMeanVal;
        springPartitions[involvedPartitions[3]] = secondMeanVal;
        // this->displaceToMeanPosition(
        //   this->initialConfig, displacements, springPartitions, link_idx);
        double retVal =
          this->updateSpringPartition(this->initialConfig,
                                      displacements,
                                      springPartitions,
                                      link_idx,
                                      oneOverSpringPartitionUpperLimit);
        // it seems to be faster to re-use memory rather than copying the whole
        // vectors
        springPartitions(involvedPartitions) = partitionsBefore;
        // displacements(involvedCoordinateIndices) = displacementsBefore;
        return retVal;
      }

      /**
       * @brief Remove cross-links which do not have any springs with a certain
       * minimum length
       *
       * @param net
       * @param displacements
       * @param springPartitions
       * @param tolerance
       */
      void removeInactiveCrosslinks(ForceBalanceNetwork& net,
                                    Eigen::VectorXd& displacements,
                                    Eigen::VectorXd& springPartitions,
                                    double tolerance) const;

      /**
       * @brief Replace the two springs traversinga a two-functional cross-links
       * with a single spring
       *
       * @param net
       * @param displacements
       * @param springPartitions
       */
      void removeTwofunctionalCrosslinks(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& displacements,
        Eigen::VectorXd& springPartitions) const;

      /**
       * @brief Add slip-links to this system
       *
       * @param sliplinkDensity
       * @param cutoff
       */
      size_t randomlyAddSliplinks(const size_t nrOfSliplinksToSample,
                                  const double cutoff = 2.0,
                                  const size_t minimumNrOfSliplinks = 0,
                                  const double sameStrandCutoff = 2.0,
                                  const bool excludeCrosslinks = false)
      {
        INVALIDARG_EXP_IFN(nrOfSliplinksToSample > minimumNrOfSliplinks,
                           "Maximum nr. should be larger than minimum, got " +
                             std::to_string(nrOfSliplinksToSample) + " and " +
                             std::to_string(minimumNrOfSliplinks) + ".");
        INVALIDARG_EXP_IFN(cutoff > 0.0,
                           "Expected a cutoff > 0.0, got " +
                             std::to_string(cutoff) + ".");
        // RUNTIME_EXP_IFN(this->initialConfig.nrOfLinks ==
        //                   this->initialConfig.nrOfNodes,
        //                 "Slip-links are only added randomly when no other "
        //                 "slip-links are in place yet.");
        INVALIDARG_EXP_IFN(minimumNrOfSliplinks <
                             this->universe.getNrOfAtoms() / 2,
                           "Minimum number of slip-links must be less than the "
                           "possible number of slip-links to place.");
        INVALIDARG_EXP_IFN(nrOfSliplinksToSample <
                             this->universe.getNrOfAtoms() / 2,
                           "Number of slip-links to place must be less than "
                           "the possible number of slip-links to place.");
        // query all the cross-linker chains we actually use to place
        // slip-links on
        std::vector<pylimer_tools::entities::Molecule> crosslinkerChains =
          this->universe.getChainsWithCrosslinker(crosslinkerType);
        std::vector<size_t> usableSpringIdxs;
        usableSpringIdxs.reserve(crosslinkerChains.size());
        // and also query all the corresponding atoms we use to place slip-links
        // on
        size_t nrOfEligibleAtoms = 0;
        std::vector<pylimer_tools::entities::Atom> eligibleAtoms;
        eligibleAtoms.reserve(this->universe.getNrOfAtoms());
        std::vector<bool> vertexIdxIsEligible =
          pylimer_tools::utils::initializeWithValue<bool>(
            this->universe.getNrOfAtoms(), false);
        std::unordered_map<size_t, size_t> atomToStrand;
        atomToStrand.reserve(this->universe.getNrOfAtoms());
        std::unordered_map<size_t, size_t> atomIdxInStrand;
        atomIdxInStrand.reserve(this->universe.getNrOfAtoms());
        size_t springId = 0;
        for (size_t i = 0; i < crosslinkerChains.size(); ++i) {
          pylimer_tools::entities::Molecule chain = crosslinkerChains[i];
          RUNTIME_EXP_IFN(chain.getType() !=
                            pylimer_tools::entities::MoleculeType::UNDEFINED,
                          "Couldn't determine molecule type.");
          if (chain.getType() ==
                pylimer_tools::entities::MoleculeType::PRIMARY_LOOP ||
              chain.getType() ==
                pylimer_tools::entities::MoleculeType::NETWORK_STRAND) {
            assert(i == this->initialConfig.springToMoleculeIds[springId]);
            // TODO: also check that this is not a higher order dangling strand
            usableSpringIdxs.push_back(i);
            nrOfEligibleAtoms +=
              crosslinkerChains[i].getNrOfAtoms() -
              crosslinkerChains[i].getAtomsOfType(this->crosslinkerType).size();
            std::vector<pylimer_tools::entities::Atom> atoms =
              crosslinkerChains[i].getAtomsLinedUp(this->crosslinkerType);
            for (size_t atomIdx = 0; atomIdx < atoms.size(); ++atomIdx) {
              pylimer_tools::entities::Atom atom = atoms[atomIdx];
              if (atom.getType() != this->crosslinkerType) {
                eligibleAtoms.push_back(atom);
                vertexIdxIsEligible[this->universe.getIdxByAtomId(
                  atom.getId())] = true;
                atomToStrand.emplace(atom.getId(), springId);
                atomIdxInStrand.emplace(atom.getId(), atomIdx);
              }
            }
            springId += 1;
          }
        }
        // build neighbourlist
        std::vector<pylimer_tools::entities::Atom> atomsForNeighbourList =
          this->universe.getAtoms();
        if (excludeCrosslinks) {
          // TODO: check whether it is faster to just only query the other ones
          atomsForNeighbourList.erase(
            std::remove_if(atomsForNeighbourList.begin(),
                           atomsForNeighbourList.end(),
                           [&](pylimer_tools::entities::Atom a) {
                             return a.getType() == this->crosslinkerType;
                           }));
        }
        pylimer_tools::entities::NeighbourList neighbourList =
          pylimer_tools::entities::NeighbourList(
            atomsForNeighbourList, this->universe.getBox(), cutoff);

        std::random_device rd{};
        std::mt19937 rng = std::mt19937(rd());
        // build list of random samples
        // this way is more performant than
        // sampling integers and checking whether they have been sampled already
        std::vector<size_t> toSampleFrom;
        std::vector<bool> isMasked;
        toSampleFrom.reserve(this->universe.getNrOfAtoms());
        isMasked.reserve(this->universe.getNrOfAtoms());
        for (size_t i = 0; i < this->universe.getNrOfAtoms(); ++i) {
          if (!excludeCrosslinks ||
              this->universe.getAtomByVertexIdx(i).getType() !=
                this->crosslinkerType) {
            toSampleFrom.push_back(i);
          }
          isMasked[i] = false;
        }
        std::shuffle(toSampleFrom.begin(), toSampleFrom.end(), rng);

        // the resulting vectors to fill
        std::vector<double> slipLinkXs;
        slipLinkXs.reserve(nrOfSliplinksToSample);
        std::vector<double> slipLinkYs;
        slipLinkYs.reserve(nrOfSliplinksToSample);
        std::vector<double> slipLinkZs;
        slipLinkZs.reserve(nrOfSliplinksToSample);
        std::vector<size_t> slipLinkStrandA;
        slipLinkStrandA.reserve(nrOfSliplinksToSample);
        std::vector<size_t> slipLinkStrandB;
        slipLinkStrandB.reserve(nrOfSliplinksToSample);
        std::vector<double> slipLinkStrandAlpha;
        slipLinkStrandAlpha.reserve(nrOfSliplinksToSample);
        std::vector<double> slipLinkStrandBeta;
        slipLinkStrandBeta.reserve(nrOfSliplinksToSample);

        pylimer_tools::entities::Box box = this->universe.getBox();

        size_t nrOfSlipLinksPlaced = 0;
        size_t nrOfAttempts = 0;
        size_t sampleIdx = 0;
        // the actual sampling loop
        while (nrOfSlipLinksPlaced < minimumNrOfSliplinks ||
               nrOfAttempts < nrOfSliplinksToSample) {
          // first, randomly sample an atom
          while (isMasked[toSampleFrom[sampleIdx]] &&
                 sampleIdx < toSampleFrom.size()) {
            sampleIdx += 1;
          }
          if (sampleIdx >= toSampleFrom.size()) {
            // this is a path that should barely ever be reached
            std::cerr << "Sample index exceeds number of samples." << std::endl;
            break;
          }
          size_t sampledVertexId = toSampleFrom[sampleIdx];
          nrOfAttempts += 1;
          sampleIdx += 1;
          pylimer_tools::entities::Atom a1 =
            this->universe.getAtomByVertexIdx(sampledVertexId);
          isMasked[sampledVertexId] = true;
          // then, find neighbouring atoms (but not from the same strand?!)
          std::vector<pylimer_tools::entities::Atom> neighbours =
            neighbourList.getAtomsCloseTo(a1);
          neighbourList.removeAtom(a1);
          // filter the neighbours to include only those from other strands
          // NOTE: this skews the whole thing a bit
          neighbours.erase(
            std::remove_if(
              neighbours.begin(),
              neighbours.end(),
              [&](pylimer_tools::entities::Atom a) -> bool {
                return (
                  atomToStrand[a.getId()] ==
                    atomToStrand[a1.getId()] // do not use "at", because not all
                                             // atoms in the neighbours have
                                             // been assigned a strand
                  &&
                  std::abs(static_cast<double>(atomIdxInStrand[a.getId()] -
                                               atomIdxInStrand[a1.getId()])) <
                    sameStrandCutoff);
              }),
            neighbours.end());
          if (neighbours.size() == 0) {
            std::cerr << "Not enough neighbours found." << std::endl;
            continue;
          }
          // then, randomly select one of them
          pylimer_tools::entities::Atom a2 = neighbours[0];
          if (neighbours.size() > 1) {
            size_t randomA2Idx =
              std::uniform_int_distribution<size_t>{ 0, neighbours.size() - 1 }(
                rng);
            pylimer_tools::entities::Atom a2 = neighbours[randomA2Idx];
          }
          // finally, remove them from the neighbour lists so that they are not
          // sampled more than once
          size_t sampledVertexId2 = this->universe.getIdxByAtomId(a2.getId());
          neighbourList.removeAtom(a2);
          isMasked[sampledVertexId2] = true;
          // it is actually quite a lot of expensive stuff done until we get to
          // this check but only this way we have the balance of removing atoms
          // to sample them only once
          if (!vertexIdxIsEligible[sampledVertexId] ||
              !vertexIdxIsEligible[sampledVertexId2]) {
            // std::cout << "Sampled vertices are not eligible" << std::endl;
            continue;
          }
          // take the mean and their index etc. to add as slip-link
          std::array<double, 3> meanPositions = a1.meanPositionWith(a2, &box);
          slipLinkXs.push_back(meanPositions[0]);
          slipLinkYs.push_back(meanPositions[1]);
          slipLinkZs.push_back(meanPositions[2]);
          slipLinkStrandA.push_back(atomToStrand.at(a1.getId()));
          slipLinkStrandB.push_back(atomToStrand.at(a2.getId()));
          slipLinkStrandAlpha.push_back(
            static_cast<double>(atomIdxInStrand.at(a1.getId())) /
            (static_cast<double>(
              this->initialConfig
                .springsContourLength[atomToStrand.at(a1.getId())])));
          slipLinkStrandBeta.push_back(
            static_cast<double>(atomIdxInStrand.at(a2.getId())) /
            (static_cast<double>(
              this->initialConfig
                .springsContourLength[atomToStrand.at(a2.getId())])));
          nrOfSlipLinksPlaced += 1;
        }

        RUNTIME_EXP_IFN(
          slipLinkStrandA.size() == slipLinkStrandB.size() &&
            slipLinkXs.size() == slipLinkYs.size() &&
            slipLinkZs.size() == slipLinkXs.size(),
          "Expect all slip-link relevant properties to have the same length");
        // with the data assembled, we can actually add them to the structure
        // and stuff
        this->addSlipLinks(slipLinkStrandA,
                           slipLinkStrandB,
                           slipLinkXs,
                           slipLinkYs,
                           slipLinkZs,
                           slipLinkStrandAlpha,
                           slipLinkStrandBeta,
                           false);
        return nrOfSlipLinksPlaced;
      }

      /**
       * @brief Investigate one parametrisation optimisation
       *
       * @param link_idx
       * @param displacements
       * @param springPartitions
       * @param innerMaxNrOfSteps
       * @param innerAlphaTol
       * @param distanceBackTolerance
       * @param residualNormSTolerance
       * @param innerMinNrOfSteps
       * @return std::tuple<Eigen::VectorXd, Eigen::VectorXd, size_t, double,
       * double, double>
       */
      std::
        tuple<Eigen::VectorXd, Eigen::VectorXd, size_t, double, double, double>
        inspectParametrisationOptimsationForLink(
          size_t link_idx,
          Eigen::VectorXd& displacements,
          Eigen::VectorXd& springPartitions,
          long int innerMaxNrOfSteps = 100,
          double innerAlphaTol = 1e-9,
          long int innerMinNrOfSteps = 1,
          const double oneOverSpringPartitionUpperLimit = 1.0)
      {
        size_t innerIterationsDone = 0;
        double displacementDone = 0.0;
        double rOverr0 = 0.0;
        double r2 = 0.0;
        double r02 = this->computePartitionUpdateZeroResidual(
          this->getSpringpartitionIndicesOfSliplink(this->initialConfig,
                                                    link_idx),
          link_idx,
          displacements,
          springPartitions,
          oneOverSpringPartitionUpperLimit);
        do {
          r2 = this->updateSpringPartition(this->initialConfig,
                                           displacements,
                                           springPartitions,
                                           link_idx,
                                           oneOverSpringPartitionUpperLimit);
          rOverr0 = r2 / r02;
          displacementDone =
            this->displaceToMeanPosition(this->initialConfig,
                                         displacements,
                                         springPartitions,
                                         link_idx,
                                         oneOverSpringPartitionUpperLimit);
          innerIterationsDone += 1;
        } while ((innerIterationsDone < innerMaxNrOfSteps &&
                  rOverr0 > innerAlphaTol && std::isfinite(rOverr0)) ||
                 innerIterationsDone < innerMinNrOfSteps);
        return std::make_tuple(displacements,
                               springPartitions,
                               innerIterationsDone,
                               displacementDone,
                               rOverr0,
                               r02);
      }

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

      Eigen::VectorXd getCurrentDisplacements() const
      {
        return this->currentDisplacements;
      }

      void setCurrentDisplacements(const Eigen::VectorXd displacements)
      {
        this->currentDisplacements = displacements;
      }

      void setSpringContourLengths(const Eigen::VectorXd springsContourLengths)
      {
        INVALIDARG_EXP_IFN(springsContourLengths.size() ==
                             this->initialConfig.springsContourLength.size(),
                           "Contour length must have the correct dimensions.");
        this->initialConfig.springsContourLength = springsContourLengths;
      }

      std::vector<Eigen::ArrayXi> getIndependentCoordinateSets(
        const ForceBalanceNetwork& net) const;

      std::pair<std::vector<Eigen::ArrayXi>, std::vector<Eigen::ArrayXi>>
      getHeuristicallyIndependentCoordinateSets(
        const ForceBalanceNetwork& net) const;

      std::vector<Eigen::ArrayXi> getRandomCoordinateSets(
        const ForceBalanceNetwork& net) const;

      /**
       * @brief Get the Nr Of Active Nodes
       *
       * @param tolerance  the tolerance: springs under a certain length are
       * considered inactive
       * @return int
       */
      int getNrOfActiveNodes(double tolerance = 0.1,
                             int minimumNrOfActiveConnections = 2,
                             int maximumNrOfActiveConnections = -1,
                             bool usePartial = false) const
      {
        return this
          ->getIdsOfActiveNodes(tolerance,
                                minimumNrOfActiveConnections,
                                maximumNrOfActiveConnections,
                                usePartial)
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
        int maximumNrOfActiveConnections = -1,
        bool usePartial = false) const;

      Eigen::VectorXd getCurrentSpringDistances() const
      {
        return this->currentSpringDistances;
      }

      Eigen::VectorXd getCurrentPartialSpringDistances() const
      {
        return this->evaluatePartialSpringDistances(
          this->initialConfig, this->currentDisplacements, this->is2D);
      }

      /**
       * @brief Get the Nr Of Active Springs connected to each node
       *
       * @param tolerance the tolerance: springs under a certain length are
       * considered inactive
       * @return Eigen::VectorXi
       */
      Eigen::VectorXi getNrOfActiveSpringsConnected(
        double tolerance = 0.1) const;

      /**
       * @brief Get the Nr Of Active Springs connected to each node
       *
       * @param tolerance the tolerance: springs under a certain length are
       * considered inactive
       * @return Eigen::VectorXi
       */
      Eigen::VectorXi getNrOfActivePartialSpringsConnected(
        double tolerance = 0.1) const;

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
       * @brief Get the Nr Of Active Springs object
       *
       * @param tol the tolerance: springs under a certain length are considered
       * inactive
       * @return int
       */
      int getNrOfActivePartialSprings(double tol = 0.1) const
      {
        return this->countNrOfActiveSprings(this->currentPartialSpringDistances,
                                            tol);
      }

      /**
       * @brief Get the Average Spring Length at the current step
       *
       * @return double
       */
      double getAverageSpringLength() const;

      std::array<std::array<double, 3>, 3> getStressTensor(
        const double oneOverSpringPartitionUpperLimit = -1.0) const;

      std::array<std::array<double, 3>, 3> getStressTensorLinkBased(
        const double oneOverSpringPartitionUpperLimit = -1.0,
        const bool xlinksOnly = false) const;

      /**
       * @brief Get the Pressure
       *
       * @return double
       */
      double getPressure() const
      {
        return this->evaluatePressure(this->initialConfig,
                                      this->currentDisplacements,
                                      this->currentSpringPartitionsVec);
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

      void addSlipLinks(const std::vector<size_t>& strandIdx1,
                        const std::vector<size_t>& strandIdx2,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& z)
      {
        std::vector<double> alphas;
        alphas.reserve(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
          alphas.push_back(0.5);
        }
        return this->addSlipLinks(
          strandIdx1, strandIdx2, x, y, z, alphas, alphas);
      }

      void addSlipLinks(const std::vector<size_t>& strandIdx1,
                        const std::vector<size_t>& strandIdx2,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& z,
                        const std::vector<double>& alpha1,
                        const std::vector<double>& alpha2,
                        bool clampAlpha = false);

      /**
       * @brief Compute the spring lenghts
       *
       * @param net the network to do the computation for
       * @param u the displacements on top of the network
       * @return Eigen::VectorXd
       */
      Eigen::VectorXd evaluateSpringDistances(const ForceBalanceNetwork& net,
                                              const Eigen::VectorXd& u,
                                              const bool is2D) const;

      Eigen::VectorXd evaluatePartialSpringDistances(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const bool is2D) const;

      /**
       * @brief Compute one spring length
       *
       * @param net
       * @param linkIndexA
       * @param linkIndexB
       * @param is2D
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d evaluateDistanceBetween(const ForceBalanceNetwork& net,
                                              const Eigen::VectorXd& u,
                                              const size_t linkIndexA,
                                              const size_t linkIndexB,
                                              const bool is2D) const;

      bool validateNetwork() const
      {
        return this->validateNetwork(this->initialConfig, this->currentDisplacements, this->currentSpringPartitionsVec);
      }

      bool validateNetwork(const ForceBalanceNetwork& net) const
      {
        return this->validateNetwork(net, this->currentDisplacements, this->currentSpringPartitionsVec);
      }

      bool validateNetwork(const ForceBalanceNetwork& net, const Eigen::VectorXd& u, const Eigen::VectorXd& springPartitions) const;

      ForceBalanceNetwork getNetwork() { return this->initialConfig; }

      Eigen::VectorXd getSpringPartitions()
      {
        return this->currentSpringPartitionsVec;
      }

      void setSpringPartitions(const Eigen::VectorXd newSpringPartitionsVec)
      {
        this->currentSpringPartitionsVec = newSpringPartitionsVec;
      }

      Eigen::Matrix3d getForceOn(
        const size_t index,
        double oneOverSpringPartitionUpperLimit = 1.0) const
      {
        Eigen::VectorXi debugNrSpringsVisited =
          Eigen::VectorXi::Zero(this->initialConfig.nrOfPartialSprings);
        return this->initialConfig.linkIsSliplink[index]
                 ? this->evaluateForceOnSlipLink(
                     index,
                     this->initialConfig,
                     this->currentDisplacements,
                     this->currentSpringPartitionsVec,
                     debugNrSpringsVisited,
                     1.0,
                     oneOverSpringPartitionUpperLimit)
                 : this->evaluateForceOnCrossLink(
                     index,
                     this->initialConfig,
                     this->currentDisplacements,
                     this->currentSpringPartitionsVec,
                     debugNrSpringsVisited,
                     1.0,
                     oneOverSpringPartitionUpperLimit);
      }

      /**
       * @brief Assemble all indices of partial springs for a particular
       * slip-link
       *
       * @param linkIdx
       * @return std::vector<size_t>
       */
      std::vector<size_t> getSpringpartitionIndicesOfSliplink(
        const ForceBalanceNetwork& net,
        const size_t linkIdx) const
      {
        std::vector<size_t> indices =
          pylimer_tools::utils::initializeWithValue<size_t>(4, 0);
        assert(indices.size() == 4);
        this->setSpringpartitionIndicesOfSliplink(indices, net, linkIdx);
        return indices;
      };

      /**
       * @brief Assemble all indices of partial springs for a particular
       * slip-link
       *
       * @param linkIdx
       * @return void
       */
      void setSpringpartitionIndicesOfSliplink(std::vector<size_t>& res_vec,
                                               const ForceBalanceNetwork& net,
                                               const size_t linkIdx) const;

      /**
       * @brief Updates the partition/parametrisation of a spring around one
       * link
       *
       */
      Eigen::VectorXd inspectSpringPartitionUpdate(const size_t linkIdx) const
      {
        Eigen::VectorXd springPartitions = this->currentSpringPartitionsVec;
        this->updateSpringPartition(this->initialConfig,
                                    this->currentDisplacements,
                                    springPartitions,
                                    linkIdx);
        return springPartitions;
      };

      /**
       * @brief Updates the partition/parametrisation of a spring around one
       * link
       *
       */
      double updateSpringPartition(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions, /* gives the parametrisation of N */
        const size_t linkIdx,
        double oneOverSpringPartitionUpperLimit = 1.0) const;

      /**
       * @brief Displace one link to the mean of all connected neighbours
       *
       * @param u the current displacements, wherein the resulting coordinates
       * shall be stored
       * @param linkIdx the idx of the link to displace
       * @return double, the distance (squared norm) displaced
       */
      Eigen::VectorXd inspectDisplacementToMeanPositionUpdate(
        const size_t linkIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0) const
      {
        Eigen::VectorXd displacements = this->currentDisplacements;
        this->displaceToMeanPosition(this->initialConfig,
                                     displacements,
                                     this->currentSpringPartitionsVec,
                                     linkIdx,
                                     oneOverSpringPartitionUpperLimit);
        return displacements;
      };

      /**
       * @brief Displace one link to the mean of all connected neighbours
       *
       * @param u the current displacements, wherein the resulting coordinates
       * shall be stored
       * @param linkIdx the idx of the link to displace
       * @return double, the distance (squared norm) displaced
       */
      Eigen::VectorXd inspectLinkDisplacementToMeanPositionUpdate(
        const size_t linkIdx,
        double damping = 1.0) const
      {
        Eigen::VectorXd displacements = this->currentDisplacements;
        Eigen::VectorXd oneOverSpringPartitions =
          this->assembleOneOverSpringPartition(
            this->initialConfig, this->currentSpringPartitionsVec);
        Eigen::Array3i mask;
        mask << 3 * linkIdx, 3 * linkIdx + 1, 3 * linkIdx + 2;
        this->displaceLinksToMeanPosition(this->initialConfig,
                                          displacements,
                                          oneOverSpringPartitions,
                                          mask,
                                          damping);
        return displacements;
      };

      /**
       * @brief Displace one link to the mean of all connected neighbours
       *
       * @param net the force balance network
       * @param u the current displacements, wherein the resulting coordinates
       * shall be stored
       * @param linkIdx the idx of the link to displace
       * @return double, the distance (squared norm) displaced
       */
      double displaceToMeanPosition(
        const ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        const Eigen::VectorXd& springPartitions,
        const size_t linkIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      /**
       * @brief Translate the spring partition vector to its 3*size
       *
       * @param net
       * @param springPartitions0
       * @return Eigen::VectorXd
       */
      Eigen::VectorXd assembleOneOverSpringPartition(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& springPartitions0,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      double displaceLinksToMeanPosition(const ForceBalanceNetwork& net,
                                         Eigen::VectorXd& u,
                                         Eigen::VectorXd& springPartitions0,
                                         double damping = 0.5) const;

      /**
       * @brief Displace one link to the mean of all connected neighbours
       *
       * @param net the force balance network
       * @param u the current displacements, wherein the resulting coordinates
       * shall be stored
       * @return double, the distance (squared norm) displaced
       */
      double displaceLinksToMeanPosition(
        const ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        const Eigen::VectorXd& oneOverSpringPartitions,
        const Eigen::ArrayXi& resultingCoordinateIndexMask,
        const double damping) const;

      /**
       * @brief Displace one link to the mean of all connected neighbours
       *
       * @param net the force balance network
       * @param u the current displacements, wherein the resulting coordinates
       * shall be stored
       * @return double, the distance (squared norm) displaced
       */
      double displaceLinksToMeanPosition(
        const ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        const Eigen::VectorXd& oneOverSpringPartitions,
        const Eigen::ArrayXi& involvedSpringPartCoordinateIndexMask,
        const Eigen::ArrayXi& resultingCoordinateIndexMask,
        const double damping = 0.5) const;

      double getDisplacementResidualNorm(double cutoff) const;

      double getDisplacementResidualNormFor(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& oneOverSpringPartitions) const;

      template<typename VectorType>
      void handlePBC(const ForceBalanceNetwork& net,
                     VectorType& distances) const
      {
        // possibly improveable PBC
        for (size_t j = 0; j < distances.size(); ++j) {
          int iterations = 0;
          assert(!std::isinf(distances[j]) && !std::isnan(distances[j]));
          while (distances[j] > net.boxHalfs[j % 3]) {
            distances[j] -= net.L[j % 3];
            iterations++;
            if (iterations > 50) {
              throw std::runtime_error(
                "Too many iterations in PBC at distance index " +
                std::to_string(j) + ", currently at " +
                std::to_string(distances[j]) + " of " +
                std::to_string(net.boxHalfs[j % 3]) + " after " +
                std::to_string(iterations) + " iterations");
            }
          }
          iterations = 0;
          while (distances[j] < -net.boxHalfs[j % 3]) {
            distances[j] += net.L[j % 3];
            iterations++;
            if (iterations > 50) {
              throw std::runtime_error(
                "Too many iterations in PBC at distance index " +
                std::to_string(j) + ", currently at " +
                std::to_string(distances[j]) + " of " +
                std::to_string(net.boxHalfs[j % 3]) + " after " +
                std::to_string(iterations) + " iterations");
            }
          }
        }
      }

    protected:
      /**
       * @brief Convert the universe to a network
       *
       * @param net the target network
       * @param crosslinkerType the atom type of the crosslinker
       * @return true
       * @return false
       */
      bool ConvertNetwork(ForceBalanceNetwork& net,
                          const int crosslinkerType,
                          bool remove2functionalCrosslinkers);

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
       * @brief Evaluate the pressure of the network at specific displacements
       *
       * @param net the network to evaluate the pressure for
       * @param u the displacements
       * @return double
       */
      double evaluatePressure(const ForceBalanceNetwork& net,
                              const Eigen::VectorXd& u,
                              const Eigen::VectorXd& springPartitions) const
      {
        auto stressTensor =
          this->evaluateStressTensor(net, u, springPartitions);
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
       * @param loopTol
       * @return std::array<std::array<double, 3>, 3>
       */
      std::array<std::array<double, 3>, 3> evaluateStressTensorLinkBased(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& springPartitions,
        const double kappa0 = 1.0,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const bool xlinksOnly = false) const;

      /**
       * @brief Compute the stress tensor
       *
       * @param net
       * @param u
       * @param loopTol
       * @return std::array<std::array<double, 3>, 3>
       */
      std::array<std::array<double, 3>, 3> evaluateStressTensor(
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& springPartitions,
        const double kappa0 = 1.0,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      /**
       * @brief Compute the force acting on a slip-link
       *
       * TODO: use "global" partial distances
       *
       * @param linkIdx
       * @param net
       * @param u
       * @param springPartitions
       * @param kappa0
       * @param minCutoff
       * @return Eigen::Vector3d
       */
      Eigen::Matrix3d evaluateForceOnSlipLink(
        const size_t linkIdx,
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& springPartitions,
        Eigen::VectorXi& debugNrSpringsVisited,
        const double kappa0 = 1.0,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      /**
       * @brief Compute the force acting on a cross-link
       *
       * TODO: use "global" partial distances
       *
       * @param linkIdx
       * @param net
       * @param u
       * @param springPartitions
       * @param kappa0
       * @param minCutoff
       * @return Eigen::Vector3d
       */
      Eigen::Matrix3d evaluateForceOnCrossLink(
        const size_t linkIdx,
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& u,
        const Eigen::VectorXd& springPartitions,
        Eigen::VectorXi& debugNrSpringsVisited,
        const double kappa0 = 1.0,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

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
            springDistances.segment(3 * i, 3).squaredNorm() > tolerance;
        }
        return result;
      }

      static std::pair<size_t, size_t> makeConnectivityKey(size_t i1, size_t i2)
      {
        return i1 > i2 ? std::make_pair(i1, i2) : std::make_pair(i2, i1);
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
      Eigen::VectorXd currentDisplacements;
      Eigen::VectorXd currentSpringDistances;
      Eigen::VectorXd currentPartialSpringDistances;
      Eigen::VectorXd
        currentSpringPartitionsVec; /* gives the parametrisation of N */
      int crosslinkerType;
      int nrOfStepsDone = 0;
      ExitReason exitReason = ExitReason::UNSET;
    };
  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
