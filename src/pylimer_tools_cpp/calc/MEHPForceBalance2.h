#ifndef MEHP_FORCE_BALANCE2_H
#define MEHP_FORCE_BALANCE2_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "MEHPForceEvaluator.h"
#include "MEHPUtilityStructures.h"
#include "OutputSupportingSimulation.h"
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
    class MEHPForceBalance2
      : public pylimer_tools::calc::OutputSupportingSimulation
    {

    private:
      // structure
      pylimer_tools::entities::Universe universe;
      ForceBalanceNetwork net;
      igraph_t graph;

      // config
      bool is2D = false;
      double kappa = 1.0;
      bool simulationHasRun = false;
      int stepOutputFrequency = 0;
      int defaultNrOfChains = 0;
      double defaultR0Squared = 0.0;
      std::string stepOutputFile;
      bool outputEndNodes = false;
      std::string endNodesFile;
      int crosslinkerType = 2;
      int splipLinkType = 3;
      int partialBondType = 2;
      int normalBondType = 1;
      bool assumeBoxLargeEnough = false;

      // cache
      Eigen::VectorXd currentSpringDistances;
      Eigen::VectorXd currentPartialSpringDistances;
      Eigen::VectorXd
        currentSpringPartitionsVec; /* gives the parametrisation of N */

      // state
      int nrOfStepsDone = 0;
      ExitReason exitReason = ExitReason::UNSET;

    public:
      MEHPForceBalance2(const pylimer_tools::entities::Universe u,
                        int crosslinkerType = 2,
                        bool is2D = false,
                        double kappa = 1.0)
        : universe(u)
      {
        this->crosslinkerType = crosslinkerType;
        // interpret network already to be able to give early results
        ForceBalanceNetwork net;
        ConvertNetwork(&net, crosslinkerType);
        this->net = net;
        this->is2D = is2D;
        this->currentSpringDistances = this->evaluateSpringDistances();
        this->currentPartialSpringDistances =
          this->evaluatePartialSpringDistances();
        this->defaultR0Squared =
          universe.computeMeanSquareEndToEndDistance(crosslinkerType);
        this->defaultNrOfChains =
          universe.getMolecules(this->crosslinkerType).size();
        this->validateNetwork();
      };

      // rule of three:
      // 1. destructor (to destroy the graph)
      ~MEHPForceBalance2();
      // 2. copy constructor
      MEHPForceBalance2(const MEHPForceBalance2& src);
      // 3. copy assignment operator
      MEHPForceBalance2& operator=(MEHPForceBalance2 src);

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
        const double initialResidualToUse = -1.0,
        const StructureSimplificationMode simplificationMode =
          StructureSimplificationMode::NO_SIMPLIFICATION,
        const double inactiveRemovalCutoff = -1.0,
        const int outputFrequency = 50,
        bool doInnerIterations = false,
        LinkSwappingMode allowSlipLinksToPassEachOther =
          LinkSwappingMode::NO_SWAPPING,
        const int swappingFrequency = 10,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const int nrOfCrosslinkSwapsAllowedPerSliplink = -1);

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
        const ForceBalanceNetwork& net,
        const std::vector<size_t>& involvedPartitions,
        const size_t link_idx,
        Eigen::VectorXd& springPartitions,
        double oneOverSpringPartitionUpperLimit = 1.0) const
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
        //   this->net, displacements, springPartitions, link_idx);
        double retVal = this->updateSpringPartition(
          net, springPartitions, link_idx, oneOverSpringPartitionUpperLimit);
        // it seems to be faster to re-use memory rather than copying the whole
        // vectors
        springPartitions(involvedPartitions) = partitionsBefore;
        // displacements(involvedCoordinateIndices) = displacementsBefore;
        return retVal;
      }

      /**
       * @brief Remove cross-linkers, springs and associated slip-links with the
       * scheme suggested by Andrei
       *
       * @param net
       * @param displacements
       * @param springPartitions
       * @param tolerance
       * @return size_t
       */
      size_t doRemovalAndreisWay(ForceBalanceNetwork& net,
                                 Eigen::VectorXd& displacements,
                                 Eigen::VectorXd& springPartitions,
                                 double tolerance) const;

      /**
       * @brief Remove cross-links which do not have any springs with a certain
       * minimum length
       *
       * @param net
       * @param displacements
       * @param springPartitions
       * @param tolerance
       */
      size_t removeInactiveCrosslinks(ForceBalanceNetwork& net,
                                      Eigen::VectorXd& displacements,
                                      Eigen::VectorXd& springPartitions,
                                      double tolerance) const;

      /**
       * @brief Remove double listed springs from cross-links
       *
       * @param net
       */
      void cleanupPrimaryLoopsInStructure(ForceBalanceNetwork& net);

      /**
       * @brief remove all vertices that don't have any connections
       */
      void removeOrphanedVertices()
      {
        if (!igraph_cattribute_GAB(&this->graph, "is_up_to_date")) {
          this->updateGraph();
        }

        igraph_vector_int_t degrees;
        igraph_vector_int_init(&degrees, igraph_vcount(&this->graph));
        igraph_degree(
          &this->graph, &degrees, igraph_vss_all(), IGRAPH_ALL, true);

        std::vector<size_t> vertexIds;

        for (size_t i = 0; i < igraph_vector_int_size(&degrees); ++i) {
          if (igraph_vector_int_get(&degrees, i) == 0) {
            vertexIds.push_back(i);
          }
        }

        if (vertexIds.size() == 0) {
          return;
        }

        igraph_vector_int_t vertexIdsVec;
        igraph_vector_int_init(&vertexIdsVec, vertexIds.size());
        pylimer_tools::utils::StdVectorToIgraphVectorT(vertexIds,
                                                       &vertexIdsVec);
        igraph_delete_vertices(&this->graph, igraph_vss_vector(&vertexIdsVec));

        this->net.isUpToDate = false;
      }

      /**
       * @brief Remove a spring (and all its parts, incl. slip-links) from the
       * structures
       *
       * @param net
       * @param springPartitions
       */
      void removeSpring(const size_t springIdx)
      {
        if (!igraph_cattribute_GAB(&this->graph, "is_up_to_date")) {
          this->updateGraph();
        }

        igraph_vector_t parentEdges;
        igraph_vector_init(&parentEdges, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "parent_edge",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &parentEdges);

        std::vector<size_t> edgeIdsToRemove;
        edgeIdsToRemove.reserve(4);
        for (size_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
          if (igraph_vector_get(&parentEdges, i) == springIdx) {
            edgeIdsToRemove.push_back(i);
          }
        }

        // actually remove all edges
        igraph_vector_int_t edgeIdsToRemoveVec;
        igraph_vector_int_init(&edgeIdsToRemoveVec, edgeIdsToRemove.size());
        pylimer_tools::utils::StdVectorToIgraphVectorT(edgeIdsToRemove,
                                                       &edgeIdsToRemoveVec);
        igraph_es_t iterator;
        igraph_delete_edges(&this->graph,
                            igraph_ess_vector(&edgeIdsToRemoveVec));

        // remove resulting edges to itself
        // NO, don't: problems if bondBoxOffset > 0.
        // igraph_attribute_combination_t comb;
        // igraph_attribute_combination_init(&comb);
        // igraph_simplify(&this->graph, false, true, &comb);
        // igraph_attribute_combination_add(
        //   &comb, NULL, IGRAPH_ATTRIBUTE_COMBINE_MEAN, NULL);
        // igraph_attribute_combination_destroy(&comb);

        igraph_vector_destroy(&parentEdges);
        igraph_vector_int_destroy(&edgeIdsToRemoveVec);

        // remove vertices that "got lost"
        this->removeOrphanedVertices();

        this->net.isUpToDate = false;
      };

      /**
       * @brief marks a certain "parent" spring as non-existing
       */
      void combineParentSprings(size_t springIdxBefore, size_t springIdxNow);

      /**
       * @brief Remove a certain, 2-functional link from the structures
       *
       * @param net
       * @param displacements
       * @param linkIdx
       */
      void remove2fLink(const size_t linkIdx);

      /**
       * @brief Given a vertex id and a rail edge, returns the other two edges
       * that are not part of the rail
       */
      std::vector<size_t> getOffRailConnectedEdgeIds(size_t vertexId,
                                                     size_t railEdge);

      /**
       * @brief Given a vertex and a connected edge, returns the edge in the
       * opposite direction
       *
       */
      size_t getOtherRailEdgeId(size_t vertexId, size_t railEdge);

      /**
       * @brief Merge two springs around a given cross-link
       *
       * @param net
       * @param springPartitions
       */
      void mergeSprings(ForceBalanceNetwork& net,
                        Eigen::VectorXd& springPartitions,
                        const size_t removedSpringIdx,
                        const size_t keptSpringIdx,
                        const size_t linkToReduce) const;

      /**
       * @brief Merge two springs around a given cross-link
       *
       * This does not require the resulting network to be valid.
       *
       * @param net
       * @param springPartitions
       */
      void mergePartialSprings(ForceBalanceNetwork& net,
                               Eigen::VectorXd& springPartitions,
                               const size_t removedSpringIdx,
                               const size_t keptSpringIdx,
                               const size_t linkToReduce,
                               bool skipEigenResize = false) const;

      /**
       * @brief Add a slip-link to a given partial spring
       *
       * This does not require the resulting network to be valid.
       *
       * @param net
       * @param springPartitions
       * @param partialSpringIdx
       * @param slipLinkIdx
       */
      size_t addSlipLinkToPartialSpring(ForceBalanceNetwork& net,
                                        Eigen::VectorXd& springPartitions,
                                        const size_t partialSpringIdx,
                                        const size_t slipLinkIdx,
                                        const double alpha) const;

      void relaxationLight(ForceBalanceNetwork& net,
                           Eigen::VectorXd& springPartitions,
                           Eigen::VectorXd& displacements,
                           const size_t linkIdx,
                           const double oneOverSpringPartitionUpperLimit = 1.0)
      {
        Eigen::VectorXd oneOverSpringPartitions = Eigen::VectorXd::Zero(0);
        this->relaxationLight(net,
                              springPartitions,
                              oneOverSpringPartitions,
                              displacements,
                              linkIdx,
                              oneOverSpringPartitionUpperLimit);
      };

      void relaxationLight(ForceBalanceNetwork& net,
                           Eigen::VectorXd& springPartitions,
                           Eigen::VectorXd& oneOverSpringPartitions,
                           Eigen::VectorXd& displacements,
                           const size_t linkIdx,
                           const double oneOverSpringPartitionUpperLimit = 1.0);

      /**
       * @brief Replace the two springs traversinga a two-functional cross-links
       * with a single spring
       *
       * @param net
       * @param displacements
       * @param springPartitions
       */
      size_t removeTwofunctionalCrosslinks()
      {
        assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));

        size_t numRemoved = 0;
        for (long int i = igraph_vcount(&this->graph); i >= 0; --i) {
          igraph_integer_t degree;
          igraph_degree_1(&this->graph, &degree, i, IGRAPH_ALL, false);
          if (degree == 2) {
            // remove this link
            numRemoved += 1;
            this->remove2fLink(i);
          }
        }

        return numRemoved;
      };

      /**
       * @brief Add slip-links to this system
       *
       * @param sliplinkDensity
       * @param cutoff
       * @return size_t the nr of actually added slip-links
       */
      size_t randomlyAddSliplinks(const size_t nrOfSliplinksToSample,
                                  const double cutoff = 2.0,
                                  const size_t minimumNrOfSliplinks = 0,
                                  const double sameStrandCutoff = 2.0,
                                  const bool excludeCrosslinks = false,
                                  const int seed = -1);

      /**
       * @brief Add slip-links to this system based on entangled loops
       *
       * @return size_t the nr of actually added slip-links
       */
      size_t addSliplinksBasedOnCycles(const int maxLoopLength = -1);

      /**
       * @brief Deform the system to match the specified box
       *
       * @param box
       */
      void deformTo(pylimer_tools::entities::Box& newBox)
      {
        this->universe.getBox().adjustCoordinatesTo(this->net.coordinates,
                                                    newBox);
        this->universe.setBox(newBox, true);
        this->net.L[0] = newBox.getLx();
        this->net.L[1] = newBox.getLy();
        this->net.L[2] = newBox.getLz();
        this->net.boxHalfs[0] = 0.5 * this->net.L[0];
        this->net.boxHalfs[1] = 0.5 * this->net.L[1];
        this->net.boxHalfs[2] = 0.5 * this->net.L[2];
        this->updateGraph();
      }

      /**
       * @brief Get the universe consisting of cross-linkers only
       *
       * @param newCrosslinkerType the type to give the cross-linkers
       * @return pylimer_tools::entities::Universe
       */
      pylimer_tools::entities::Universe getCrosslinkerVerse() const;

      int getDefaultNrOfChains() const { return this->defaultNrOfChains; }

      double getDefaultR0Square() const { return this->defaultR0Squared; }

      double getVolume() override { return this->net.vol; }

      int getNrOfNodes() const { return this->net.nrOfNodes; }

      int getNrOfLinks() const { return this->net.nrOfLinks; }

      size_t getNumBonds() override { return this->getNrOfSprings(); }

      size_t getNumExtraBonds() override { return 0; }

      long int getNumBondsToForm() override { return 0; }

      size_t getNumAtoms() override { return this->getNrOfNodes(); }

      size_t getNumExtraAtoms() override { return this->getNrOfLinks(); }

      int getNrOfSprings() const { return this->net.nrOfSprings; }

      void setSpringContourLengths(const Eigen::VectorXd springsContourLengths)
      {
        INVALIDARG_EXP_IFN(springsContourLengths.size() ==
                             this->net.springsContourLength.size(),
                           "Contour length must have the correct dimensions.");
        this->net.springsContourLength = springsContourLengths;
      }

      std::vector<Eigen::ArrayXi> getIndependentCoordinateSets(
        const ForceBalanceNetwork& net) const;

      std::pair<std::vector<Eigen::ArrayXi>, std::vector<Eigen::ArrayXi>>
      getHeuristicallyIndependentCoordinateSets(
        const ForceBalanceNetwork& net) const;

      std::vector<Eigen::ArrayXi> getRandomCoordinateSets(
        const ForceBalanceNetwork& net) const;

      void configAssumeBoxLargeEnough(bool assumption)
      {
        throw std::invalid_argument(
          "Assumption of a large enough box is not supported yet");
        this->assumeBoxLargeEnough = assumption;
      }

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
       * @brief Get the Soluble Weight Fraction
       *
       * @param tolerance
       * @return double
       */
      double getSolubleWeightFraction(double tolerance = 0.1)
      {
        return this->computeSolubleWeightFraction(
          &this->net, this->currentSpringDistances, tolerance);
      }

      /**
       * @brief Get the Dangling Weight Fraction
       *
       * @param tolerance
       * @return double
       */
      double getDanglingWeightFraction(double tolerance = 0.1)
      {
        return this->computeDanglingWeightFraction(
          &this->net, this->currentSpringDistances, tolerance);
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
       * @brief Compute the weight fraction of non-active springs
       *
       * @param net
       * @param springDistances
       * @param tolerance
       * @return double
       */
      double computeDanglingWeightFraction(
        ForceBalanceNetwork* net,
        const Eigen::VectorXd& springDistances,
        const double tolerance = 0.1) const
      {
        if (net->nrOfSprings * 3 != springDistances.size()) {
          throw std::invalid_argument(
            "Spring distances and network don't match");
        }
        if (net->nrOfSprings < 1) {
          return 1.;
        }
        // find all active springs
        ArrayXb activeSprings =
          this->findActiveSprings(springDistances, tolerance);
        if (activeSprings.count() == 0) {
          return 1.;
        }
        // as of now, the springsContourLength is equal to the number of bonds
        // from cross-link to cross-link. therefore, the number of atoms of each
        // of these springs is one less
        Eigen::ArrayXd allActiveAtomsPerChains =
          activeSprings.cast<double>() *
          (net->springsContourLength.array() -
           Eigen::ArrayXd::Ones(net->nrOfSprings));
        // finally, normalise by the number of atoms.
        // NOTE: currently, the weight of the atoms is ignored
        return 1. - ((allActiveAtomsPerChains).matrix().sum() +
                     this->getNrOfActiveNodes()) /
                      this->universe.getNrOfAtoms();
      }

      /**
       * @brief Compute the weight fraction of springs connected to active
       * springs (any depth)
       *
       * @param net
       * @param springDistances
       * @param tolerance
       * @return double
       */
      double computeSolubleWeightFraction(
        ForceBalanceNetwork* net,
        const Eigen::VectorXd& springDistances,
        const double tolerance = 0.1) const
      {
        if (net->nrOfSprings * 3 != springDistances.size()) {
          throw std::invalid_argument(
            "Spring distances and network don't match");
        }
        if (net->nrOfSprings < 1) {
          return 1.;
        }
        std::vector<long int> activeCrosslinkIds =
          this->getIdsOfActiveNodes(tolerance, 2, -1);
        double activeClusterWeightFraction =
          this->universe.computeWeightFractionOfClustersAssociatedWith(
            activeCrosslinkIds);

        return 1. - activeClusterWeightFraction;
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
        return this->evaluatePartialSpringDistances();
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

      std::array<std::array<double, 3>, 3> getStressTensorArray(
        const double oneOverSpringPartitionUpperLimit = -1.0) const;

      std::array<std::array<double, 3>, 3> getStressTensorArrayLinkBased(
        const double oneOverSpringPartitionUpperLimit = -1.0,
        const bool xlinksOnly = false) const;

      /**
       * @brief Get the Pressure
       *
       * @return double
       */
      double getPressure() const
      {
        return this->evaluatePressure(this->net,
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
                        bool clampAlpha = false)
      {
        std::vector<std::vector<size_t>> loops;
        std::vector<std::vector<size_t>> loopsOfSliplinks;
        return this->addSlipLinks(strandIdx1,
                                  strandIdx2,
                                  x,
                                  y,
                                  z,
                                  alpha1,
                                  alpha2,
                                  loops,
                                  loopsOfSliplinks,
                                  clampAlpha);
      }

      void addSlipLinks(const std::vector<size_t>& strandIdx1,
                        const std::vector<size_t>& strandIdx2,
                        const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& z,
                        const std::vector<double>& alpha1,
                        const std::vector<double>& alpha2,
                        std::vector<std::vector<size_t>> loops,
                        std::vector<std::vector<size_t>> loopsOfSliplinks,
                        bool clampAlpha = false);

      /**
       * @brief Compute the spring lenghts
       *
       * @param net the network to do the computation for
       * @param u the displacements on top of the network
       * @return Eigen::VectorXd
       */
      Eigen::VectorXd evaluateSpringDistances() const;

      Eigen::VectorXd evaluatePartialSpringDistances() const;

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
                                              const size_t linkIndexA,
                                              const size_t linkIndexB) const;

      bool validateNetwork() const
      {
        return this->validateNetwork(this->net,
                                     this->currentSpringPartitionsVec);
      }

      bool validateNetwork(const ForceBalanceNetwork& net) const
      {
        return this->validateNetwork(net, this->currentSpringPartitionsVec);
      }

      bool validateNetwork(const ForceBalanceNetwork& net,
                           const Eigen::VectorXd& springPartitions) const;

      ForceBalanceNetwork getNetwork() { return this->net; }

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
          Eigen::VectorXi::Zero(this->net.nrOfPartialSprings);
        return this->net.linkIsSliplink[index]
                 ? this->evaluateForceOnSlipLink(
                     index,
                     this->net,
                     this->currentSpringPartitionsVec,
                     debugNrSpringsVisited,
                     1.0,
                     oneOverSpringPartitionUpperLimit)
                 : this->evaluateForceOnCrossLink(
                     index,
                     this->net,
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
        this->updateSpringPartition(this->net, springPartitions, linkIdx);
        return springPartitions;
      };

      /**
       * @brief Updates the partition/parametrisation of a spring around one
       * link
       *
       */
      double updateSpringPartition(
        const ForceBalanceNetwork& net,
        Eigen::VectorXd& springPartitions, /* gives the parametrisation of N */
        const size_t linkIdx,
        double oneOverSpringPartitionUpperLimit = 1.0,
        bool allowSlipLinksToPassEachOther = false) const
      {
        Eigen::VectorXd oneOverSpringPartitions = Eigen::VectorXd::Zero(0);
        return this->updateSpringPartition(net,
                                           springPartitions,
                                           oneOverSpringPartitions,
                                           linkIdx,
                                           oneOverSpringPartitionUpperLimit,
                                           allowSlipLinksToPassEachOther);
      };

      /**
       * @brief Updates the partition/parametrisation of a spring around one
       * link
       *
       */
      double updateSpringPartition(
        const ForceBalanceNetwork& net,
        Eigen::VectorXd& springPartitions, /* gives the parametrisation of N */
        Eigen::VectorXd&
          oneOverSpringPartitions, /* gives the parametrisation of N */
        const size_t linkIdx,
        double oneOverSpringPartitionUpperLimit = 1.0,
        bool allowSlipLinksToPassEachOther = false) const;

      /**
       * @brief Loop all slip-links and move them if appropriate to other
       * springs
       *
       * @param net
       * @param u
       * @param springPartitions
       * @param oneOverSpringPartitionUpperLimit
       */
      void moveSlipLinksToTheirBestBranch(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions,
        const double oneOverSpringPartitionUpperLimit,
        const int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
        const bool respectLoops = true);

      /**
       * @brief Move a slip-link if appropriate to other springs
       *
       * @param net
       * @param u
       * @param springPartitions
       * @param oneOverSpringPartitionUpperLimit
       */
      void moveSlipLinkToItsBestBranch(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions,
        size_t slipLinkIdx,
        const double oneOverSpringPartitionUpperLimit,
        const int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
        const bool respectLoops = true);

      /**
       * @brief Loop all springs, swap slip-links on them if they are close
       * enough
       *
       * @param net
       */
      void swapSlipLinksInclXlinks(ForceBalanceNetwork& net,
                                   Eigen::VectorXd& springPartitions,
                                   double swappableCutoff,
                                   const bool respectLoops = true);

      /**
       * @brief Loop all springs, swap slip-links on them if they are close
       * enough
       *
       * @param net
       */
      void swapSlipLinks(ForceBalanceNetwork& net,
                         Eigen::VectorXd& springPartitions,
                         double swappableCutoff);

      /**
       *
       */
      void swapSlipLinks(ForceBalanceNetwork& net,
                         const size_t partialSpringIdx)
      {
        const size_t linkIdx1 = net.springPartIndexA[partialSpringIdx];
        const size_t linkIdx2 = net.springPartIndexB[partialSpringIdx];
        INVALIDARG_EXP_IFN(linkIdx1 != linkIdx2,
                           "Cannot swap link with itself.");

        INVALIDARG_EXP_IFN(
          net.linkIsSliplink[linkIdx1],
          "Only partial springs with only slip-links allow swapping.");
        INVALIDARG_EXP_IFN(
          net.linkIsSliplink[linkIdx2],
          "Only partial springs with only slip-links allow swapping.");

        const size_t springIdx = net.partialToFullSpringIndex[partialSpringIdx];

        igraph_vector_int_t edgesOfLink1;
        igraph_vector_int_init(&edgesOfLink1, 4);
        igraph_incident(&this->graph, &edgesOfLink1, linkIdx1, IGRAPH_ALL);

        igraph_vector_int_t edgesOfLink2;
        igraph_vector_int_init(&edgesOfLink2, 4);
        igraph_incident(&this->graph, &edgesOfLink2, linkIdx2, IGRAPH_ALL);

        // TODO: figure out which edge to actually cut.
        // difficult when same "parent" bond

        igraph_vector_int_destroy(&edgesOfLink1);
        igraph_vector_int_destroy(&edgesOfLink2);
        this->net.isUpToDate = false;
      };

      bool swapSlipLinkReversibly(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions,
        const size_t partialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
        const bool respectLoops = true);

      long int rotateSlipLinkAroundCrosslink(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& springPartitions,
        const size_t partialSpringIdx,
        double oneOverSpringPartitionUpperLimit = 1.0,
        const bool respectLoops = true);

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
        const Eigen::VectorXd& oneOverSpringPartitions) const;

      /**
       * @brief Get the Link Indices of all neighbours of a specified link
       *
       * @param net
       * @param linkIdx
       * @return std::vector<size_t>
       */
      std::vector<size_t> getNeighbourLinkIndices(const size_t linkIdx)
      {
        std::vector<size_t> results;
        results.reserve(4);
        assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));

        igraph_vector_int_list_t res;
        igraph_vector_int_list_init(&res, 4);
        igraph_neighborhood(
          &this->graph, &res, igraph_vss_1(linkIdx), 1, IGRAPH_ALL, 1);

        for (size_t i = 0; i < igraph_vector_int_list_size(&res); ++i) {
          igraph_vector_int_t* resI = igraph_vector_int_list_get_ptr(&res, i);
          for (size_t j = 0; j < igraph_vector_int_size(resI); ++j) {
            results.push_back(igraph_vector_int_get(resI, j));
          }
        }

        igraph_vector_int_list_destroy(&res);

        return results;
      }

      void writeRestartFile(std::string& file) override
      {
        throw std::runtime_error("Restart not supported yet");
      }

      double getTimestep() override
      {
        return 1.; // this->dt;
      }

      double getCurrentTime(double currentStep) override
      {
        return this->nrOfStepsDone;
      }

      Eigen::Matrix3d getStressTensor() override
      {
        std::array<std::array<double, 3>, 3> stressTensor =
          this->getStressTensorArray();
        Eigen::Matrix3d result = Eigen::Matrix3d::Zero();
        for (size_t i = 0; i < 3; ++i) {
          for (size_t j = 0; j < 3; ++j) {
            result(i, j) = stressTensor[i][j];
          }
        }
        return result;
      }
      int getNumShifts() override { return 0; }
      int getNumRelocations() override { return 0; }

      Eigen::VectorXd getBondLengths() override
      {
        Eigen::VectorXd lens =
          Eigen::VectorXd::Zero(this->currentSpringDistances.size() / 3);

        for (size_t i = 0; i < this->currentSpringDistances.size() / 3; ++i) {
          double b = lens.segment(3 * i, 3).norm();
          lens[i] = b;
        }

        return lens;
      }

      Eigen::VectorXd getCoordinates() override
      {
        return this->net.coordinates;
      }

      double getTemperature() override
      {
        return -1; // TODO: implement?
      }

      size_t getNumParticles() override { return this->net.nrOfNodes; }

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
                          bool remove2functionalCrosslinkers = false,
                          bool removeDanglingChains = false);

      /**
       * @brief Convert the internal graph representation to the internal
       * ForceBalanceNetwork representation
       *
       */
      void convertFromGraph()
      {
        if (!igraph_cattribute_GAB(&this->graph, "is_up_to_date")) {
          this->updateGraph();
        }

        this->net.nrOfPartialSprings = igraph_ecount(&this->graph);
        this->net.nrOfLinks = igraph_vcount(&this->graph);
        igraph_vector_t edgeTypes;
        igraph_vector_init(&edgeTypes, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "type",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &edgeTypes);
        size_t numPartialSprings = 0;
        for (size_t i = 0; i < igraph_vector_size(&edgeTypes); ++i) {
          numPartialSprings +=
            (static_cast<igraph_integer_t>(igraph_vector_get(&edgeTypes, i)) ==
             this->partialBondType);
        }
        this->net.nrOfSprings =
          this->net.nrOfPartialSprings - numPartialSprings;

        // reset & resize
        this->net.coordinates.resize(3 * this->net.nrOfLinks);
        this->net.springsContourLength.resize(this->net.nrOfSprings);
        this->net.springIndicesOfLinks.clear();
        this->net.linkIndicesOfSprings.clear();
        this->net.partialSpringIsPartial.setConstant(false);
        this->net.localToGlobalSpringIndex.clear();
        this->net.partialToFullSpringIndex.resize(this->net.nrOfPartialSprings);
        this->net.springPartCoordinateIndexA.resize(
          3 * this->net.nrOfPartialSprings);
        this->net.springPartCoordinateIndexB.resize(
          3 * this->net.nrOfPartialSprings);
        this->net.springPartIndexA.resize(this->net.nrOfPartialSprings);
        this->net.springPartIndexB.resize(this->net.nrOfPartialSprings);
        this->net.springPartBoxOffset.resize(this->net.nrOfPartialSprings);
        this->currentSpringPartitionsVec.resize(this->net.nrOfPartialSprings);

        // reset everything we cleared
        for (size_t i = 0; i < this->net.nrOfSprings; ++i) {
          std::vector<size_t> vec;
          this->net.linkIndicesOfSprings.push_back(vec);
        }
        for (size_t i = 0; i < this->net.nrOfLinks; ++i) {
          std::vector<size_t> vec;
          this->net.springIndicesOfLinks.push_back(vec);
        }

        // fetch other properties needed
        // springs / partial springs
        igraph_vector_int_t allEdges;
        igraph_vector_int_init(&allEdges, 2 * this->net.nrOfPartialSprings);
        if (igraph_edges(
              &this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &allEdges)) {
          throw std::runtime_error("Failed to get all edges");
        }

        igraph_vector_t parentEdges;
        igraph_vector_init(&parentEdges, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "parent_edge",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &parentEdges);

        igraph_vector_t partitionFraction;
        igraph_vector_init(&partitionFraction, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "partition_fraction",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &partitionFraction);

        igraph_vector_t bondBoxOffsetX;
        igraph_vector_init(&bondBoxOffsetX, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "bond_box_x",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &bondBoxOffsetX);
        igraph_vector_t bondBoxOffsetY;
        igraph_vector_init(&bondBoxOffsetY, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "bond_box_y",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &bondBoxOffsetY);
        igraph_vector_t bondBoxOffsetZ;
        igraph_vector_init(&bondBoxOffsetZ, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "bond_box_z",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &bondBoxOffsetZ);

        for (size_t i = 0; i < this->net.nrOfPartialSprings; ++i) {
          this->net.partialToFullSpringIndex(i) =
            static_cast<igraph_integer_t>(igraph_vector_get(&parentEdges, i));
          this->net.springPartIndexA(i) =
            igraph_vector_int_get(&allEdges, 2 * i);
          this->net.springPartIndexB(i) =
            igraph_vector_int_get(&allEdges, 2 * i + 1);

          for (size_t dir = 0; dir < 3; ++dir) {
            this->net.springPartCoordinateIndexA(3 * i + dir) =
              3 * this->net.springPartIndexA(i) + dir;
            this->net.springPartCoordinateIndexB(3 * i + dir) =
              3 * this->net.springPartIndexB(i) + dir;
          }

          this->net
            .springIndicesOfLinks[igraph_vector_int_get(&allEdges, 2 * i)]
            .push_back(i);
          this->net
            .springIndicesOfLinks[igraph_vector_int_get(&allEdges, 2 * i + 1)]
            .push_back(i);
          this->net.springPartBoxOffset(i * 3 + 0) =
            igraph_vector_get(&bondBoxOffsetX, i);
          this->net.springPartBoxOffset(i * 3 + 1) =
            igraph_vector_get(&bondBoxOffsetY, i);
          this->net.springPartBoxOffset(i * 3 + 2) =
            igraph_vector_get(&bondBoxOffsetZ, i);
          this->net.partialSpringIsPartial(i) =
            (igraph_vector_get(&edgeTypes, i) == this->partialBondType);
          this->currentSpringPartitionsVec(i) =
            igraph_vector_get(&partitionFraction, i);
        }

        // cleanup
        igraph_vector_destroy(&edgeTypes);
        igraph_vector_int_destroy(&allEdges);
        igraph_vector_destroy(&parentEdges);
        igraph_vector_destroy(&bondBoxOffsetX);
        igraph_vector_destroy(&bondBoxOffsetY);
        igraph_vector_destroy(&bondBoxOffsetZ);

        // same for per-link properties
        igraph_vector_t coordsX;
        igraph_vector_init(&coordsX, this->net.nrOfLinks);
        igraph_cattribute_VANV(&this->graph, "x", igraph_vss_all(), &coordsX);

        igraph_vector_t coordsY;
        igraph_vector_init(&coordsY, this->net.nrOfLinks);
        igraph_cattribute_VANV(&this->graph, "y", igraph_vss_all(), &coordsY);

        igraph_vector_t coordsZ;
        igraph_vector_init(&coordsZ, this->net.nrOfLinks);
        igraph_cattribute_VANV(&this->graph, "z", igraph_vss_all(), &coordsZ);

        igraph_vector_t linkType;
        igraph_vector_init(&linkType, this->net.nrOfLinks);
        igraph_cattribute_VANV(
          &this->graph, "type", igraph_vss_all(), &linkType);

        // actually write things
        for (size_t i = 0; i < this->net.nrOfLinks; ++i) {
          this->net.coordinates(3 * i + 0) = igraph_vector_get(&coordsX, i);
          this->net.coordinates(3 * i + 1) = igraph_vector_get(&coordsY, i);
          this->net.coordinates(3 * i + 2) = igraph_vector_get(&coordsZ, i);
          this->net.linkIsSliplink(i) =
            static_cast<igraph_integer_t>(igraph_vector_get(&linkType, i)) ==
            this->splipLinkType;
        }

        // cleanup
        igraph_vector_destroy(&linkType);
        igraph_vector_destroy(&coordsX);
        igraph_vector_destroy(&coordsY);
        igraph_vector_destroy(&coordsZ);

        // mark as done
        net.isUpToDate = true;
      };

      void updateGraph()
      {
        // write the current coordinates to the graph
        igraph_vector_t coordsX;
        igraph_vector_init(&coordsX, this->net.nrOfLinks);

        igraph_vector_t coordsY;
        igraph_vector_init(&coordsY, this->net.nrOfLinks);

        igraph_vector_t coordsZ;
        igraph_vector_init(&coordsZ, this->net.nrOfLinks);

        for (size_t i = 0; i < this->net.nrOfLinks; ++i) {
          igraph_vector_set(&coordsX, i, this->net.coordinates(3 * i + 0));
          igraph_vector_set(&coordsY, i, this->net.coordinates(3 * i + 1));
          igraph_vector_set(&coordsZ, i, this->net.coordinates(3 * i + 2));
        }

        igraph_cattribute_VAN_setv(&this->graph, "x", &coordsX);
        igraph_cattribute_VAN_setv(&this->graph, "y", &coordsY);
        igraph_cattribute_VAN_setv(&this->graph, "z", &coordsZ);

        igraph_vector_destroy(&coordsX);
        igraph_vector_destroy(&coordsY);
        igraph_vector_destroy(&coordsZ);

        // as well as the current spring partition
        igraph_vector_t partitionFraction;
        igraph_vector_init(&partitionFraction, this->net.nrOfPartialSprings);

        for (size_t i = 0; i < this->net.nrOfPartialSprings; ++i) {
          igraph_vector_set(
            &partitionFraction, i, this->currentSpringPartitionsVec(i));
        }

        igraph_cattribute_EAN_setv(
          &this->graph, "partition_fraction", &partitionFraction);

        igraph_vector_destroy(&partitionFraction);

        igraph_cattribute_GAB_set(&this->graph, "is_up_to_date", true);
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
       * @brief Evaluate the pressure of the network at specific displacements
       *
       * @param net the network to evaluate the pressure for
       * @param u the displacements
       * @return double
       */
      double evaluatePressure(const ForceBalanceNetwork& net,
                              const Eigen::VectorXd& springPartitions) const
      {
        auto stressTensor = this->evaluateStressTensor(net, springPartitions);
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

      Eigen::Matrix3d evaluateStressTensorForLinks(
        const std::vector<size_t> linkIndices,
        const ForceBalanceNetwork& net,
        const Eigen::VectorXd& springPartitions,
        const double kappa0 = 1.0,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

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

      bool swapSlipLinksReversibly(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions,
        const size_t partialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0);

      bool swapSlipLinkWithXlinkReversibly(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& u,
        Eigen::VectorXd& springPartitions,
        const size_t partialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const bool respectLoops = true);
    };
  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
