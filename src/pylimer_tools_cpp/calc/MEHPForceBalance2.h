#ifndef MEHP_FORCE_BALANCE2_H
#define MEHP_FORCE_BALANCE2_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/NeighbourList.h"
#include "../entities/Universe.h"
#include "EntanglementDetector.h"
#include "MEHPForceEvaluator.h"
#include "MEHPUtilityStructures.h"
#include "OutputSupportingSimulation.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <nlopt.hpp>
#include <random>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>
extern "C"
{
#include <igraph/igraph.h>
}

namespace pylimer_tools {
namespace calc {
  namespace mehp {
#ifndef CLAMP_ONE_OVER_SPRINGPARTITION
/**
 * @brief a macro for doing the clamping in the routines using kappa,
 * to prevent division by zero issues / multiplications by infinity
 */
#define CLAMP_ONE_OVER_SPRINGPARTITION(                                        \
  isPartialSpring, val, N, oneOverSpringPartitionUpperLimit)                   \
  ((!isPartialSpring)                                                          \
     ? val                                                                     \
     : std::clamp(val,                                                         \
                  (oneOverSpringPartitionUpperLimit > 0.)                      \
                    ? (1. / (N - 1. / oneOverSpringPartitionUpperLimit))       \
                    : (0.0),                                                   \
                  (oneOverSpringPartitionUpperLimit > 0.)                      \
                    ? (oneOverSpringPartitionUpperLimit)                       \
                    : (N)));
#endif

    static inline igraph_integer_t castToIgraphInt(igraph_real_t c)
    {
      return static_cast<igraph_integer_t>(std::lround(c));
    }

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
      int defaultNrOfChains = 0;
      double defaultR0Squared = 0.0;
      bool outputEndNodes = false;
      std::string endNodesFile;
      int crosslinkerType = 2;
      int slipLinkType = 3;
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

      MEHPForceBalance2(const pylimer_tools::entities::Universe& u,
                        int crosslinkerType = 2,
                        bool is2D = false,
                        double kappa = 1.0)
        : universe(u)
      {
        this->crosslinkerType = crosslinkerType;
        // interpret network already to be able to give early results
        ForceBalanceNetwork network;
        for (size_t dir = 0; dir < 3; ++dir) {
          network.L[dir] = universe.getBox().getL()[dir];
          network.boxHalfs[dir] = network.L[dir] * 0.5;
        }
        network.vol = universe.getBox().getVolume();
        network.isUpToDate = false;
        this->net = network;
        this->is2D = is2D;
        igraph_empty(&this->graph, 0, IGRAPH_UNDIRECTED);
        igraph_cattribute_GAB_set(&this->graph, "is_up_to_date", false);
      };

    public:
      //----------------------------------------------------------------
      // MARK: Constructors
      //----------------------------------------------------------------

      // rule of three:
      // 1. destructor (to destroy the graph)
      ~MEHPForceBalance2() { igraph_destroy(&this->graph); };
      // 2. copy constructor
      MEHPForceBalance2(const MEHPForceBalance2& src)
        : MEHPForceBalance2(src.universe,
                            src.crosslinkerType,
                            src.is2D,
                            src.kappa)
      {
        igraph_copy(&this->graph, &src.graph);

        // structure
        this->net = src.net;
        // config
        this->simulationHasRun = src.simulationHasRun;
        this->defaultNrOfChains = src.defaultNrOfChains;
        this->outputEndNodes = src.outputEndNodes;
        this->endNodesFile = src.endNodesFile;
        this->slipLinkType = src.slipLinkType;
        this->partialBondType = src.partialBondType;
        this->normalBondType = src.normalBondType;
        this->assumeBoxLargeEnough = src.assumeBoxLargeEnough;
        // cache
        this->currentSpringDistances = src.currentSpringDistances;
        this->currentPartialSpringDistances = src.currentPartialSpringDistances;
        this->currentSpringPartitionsVec = src.currentSpringPartitionsVec;
        this->nrOfStepsDone = src.nrOfStepsDone;
        this->exitReason = src.exitReason;
      };
      // 3. copy assignment operator
      MEHPForceBalance2& operator=(MEHPForceBalance2 src)
      {
        // sructure
        std::swap(this->graph, src.graph);
        std::swap(this->universe, src.universe);
        std::swap(this->net, src.net);
        // config
        std::swap(this->is2D, src.is2D);
        std::swap(this->kappa, src.kappa);
        std::swap(this->simulationHasRun, src.simulationHasRun);
        std::swap(this->defaultNrOfChains, src.defaultNrOfChains);
        std::swap(this->outputEndNodes, src.outputEndNodes);
        std::swap(this->endNodesFile, src.endNodesFile);
        std::swap(this->crosslinkerType, src.crosslinkerType);
        std::swap(this->slipLinkType, src.slipLinkType);
        std::swap(this->partialBondType, src.partialBondType);
        std::swap(this->normalBondType, src.normalBondType);
        std::swap(this->assumeBoxLargeEnough, src.assumeBoxLargeEnough);
        // cache
        std::swap(this->currentSpringDistances, src.currentSpringDistances);
        std::swap(this->currentPartialSpringDistances,
                  src.currentPartialSpringDistances);
        std::swap(this->currentSpringPartitionsVec,
                  src.currentSpringPartitionsVec);
        std::swap(this->nrOfStepsDone, src.nrOfStepsDone);
        std::swap(this->exitReason, src.exitReason);

        return *this;
      };

      //----------------------------------------------------------------
      // MARK: Constructors
      //----------------------------------------------------------------

      static MEHPForceBalance2 constructWithoutSlipLinks(
        const pylimer_tools::entities::Universe& universe,
        int crosslinkerType = 2,
        bool is2D = false,
        double kappa = 1.0)
      {
        return MEHPForceBalance2::constructWithRandomSlipLinks(
          universe, 0, 1.0, 0, 1, "", crosslinkerType, is2D, kappa);
      }

      static MEHPForceBalance2 constructWithSlipLinks(
        const pylimer_tools::entities::Universe& universe,
        const std::vector<size_t>& strandIdx1,
        const std::vector<size_t>& strandIdx2,
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& z,
        const std::vector<double>& alpha1,
        const std::vector<double>& alpha2,
        int crosslinkerType = 2,
        bool is2D = false,
        double kappa = 1.0,
        bool clampAlpha = true)
      {
        MEHPForceBalance2 fb = MEHPForceBalance2::constructWithoutSlipLinks(
          universe, crosslinkerType, is2D, kappa);
        fb.validateNetwork();
        fb.addSlipLinks(
          strandIdx1, strandIdx2, x, y, z, alpha1, alpha2, clampAlpha);
        // convert the graph to the network usable for simulations
        fb.finaliseInitialisation();
        return fb;
      }

      static MEHPForceBalance2 constructWithRandomSlipLinks(
        const pylimer_tools::entities::Universe& universe,
        const size_t nrOfSliplinksToSample,
        const double cutoff,
        const size_t minimumNrOfSliplinks,
        const double sameStrandCutoff,
        const std::string seed = "",
        int crosslinkerType = 2,
        bool is2D = false,
        double kappa = 1.0)
      {
        MEHPForceBalance2 fb =
          MEHPForceBalance2(universe, crosslinkerType, is2D, kappa);
        fb.net.isUpToDate = false;
        igraph_cattribute_GAB_set(&fb.graph, "is_up_to_date", true);

        // sample the "entanglements"
        pylimer_tools::calc::entanglement_detection::AtomPairEntanglements
          entanglements = pylimer_tools::calc::entanglement_detection::
            randomlyFindEntanglements(universe,
                                      nrOfSliplinksToSample,
                                      cutoff,
                                      minimumNrOfSliplinks,
                                      sameStrandCutoff,
                                      seed,
                                      crosslinkerType);

        std::vector<std::pair<size_t, size_t>> pairsOfAtoms =
          entanglements.pairsOfAtoms;
        std::vector<long int> pairOfAtom = entanglements.pairOfAtom;
        // std::cout << "Found " << pairsOfAtoms.size() << " random slip-links."
        //           << std::endl;

        std::vector<pylimer_tools::entities::Molecule> crosslinkerChains =
          universe.getChainsWithCrosslinker(crosslinkerType);

        // add ends of chains
        std::vector<igraph_integer_t> endAtomIdxToVertexId =
          pylimer_tools::utils::initializeWithValue<igraph_integer_t>(
            universe.getNrOfAtoms(), -1);
        igraph_integer_t currentVertexId = 0;
        for (size_t i = 0; i < crosslinkerChains.size(); ++i) {
          pylimer_tools::entities::Molecule chain = crosslinkerChains[i];
          if (chain.getLength() < 2) {
            continue;
          }
          std::vector<pylimer_tools::entities::Atom> linedUpAtoms =
            chain.getAtomsLinedUp(crosslinkerType, false, true);
          size_t firstAtomIdx =
            universe.getIdxByAtomId(linedUpAtoms[0].getId());
          if (endAtomIdxToVertexId[firstAtomIdx] == -1) {
            endAtomIdxToVertexId[firstAtomIdx] = currentVertexId;
            currentVertexId += 1;
          }
          size_t lastAtomIdx = universe.getIdxByAtomId(
            pylimer_tools::utils::last(linedUpAtoms).getId());
          if (endAtomIdxToVertexId[lastAtomIdx] == -1) {
            endAtomIdxToVertexId[lastAtomIdx] = currentVertexId;
            currentVertexId += 1;
          }
        }

        // create `currentVertexId` vertices for the chain-end atoms, and
        // `pairsOfAtoms.size()` vertices for the so many slip-links
        igraph_add_vertices(
          &fb.graph, currentVertexId + pairsOfAtoms.size(), nullptr);

        size_t parentEdgeId = 0;
        for (size_t chainIdx = 0; chainIdx < crosslinkerChains.size();
             ++chainIdx) {
          pylimer_tools::entities::Molecule chain = crosslinkerChains[chainIdx];
          if (chain.getLength() < 2) {
            continue;
          }

          std::vector<pylimer_tools::entities::Atom> linedUpAtoms =
            chain.getAtomsLinedUp(crosslinkerType, false, true);
          size_t previousIdx = 0;
          igraph_integer_t previousVertexId =
            endAtomIdxToVertexId[universe.getIdxByAtomId(
              linedUpAtoms[0].getId())];
          fb.setVertexPropertiesFromAtom(
            previousVertexId, linedUpAtoms[0], fb.crosslinkerType);
          pylimer_tools::entities::Atom lastAtom =
            pylimer_tools::utils::last(linedUpAtoms);
          fb.setVertexPropertiesFromAtom(
            endAtomIdxToVertexId[universe.getIdxByAtomId(lastAtom.getId())],
            lastAtom,
            fb.crosslinkerType);
          assert(linedUpAtoms.size() == chain.getLength() ||
                 linedUpAtoms.size() == chain.getLength() + 1);
          if (pairsOfAtoms.size() > 0) {
            for (size_t i = 1; i < linedUpAtoms.size() - 1; i++) {
              pylimer_tools::entities::Atom a = linedUpAtoms[i];
              if (pairOfAtom[universe.getIdxByAtomId(a.getId())] != -1) {
                igraph_integer_t thisVertexId =
                  currentVertexId +
                  pairOfAtom[universe.getIdxByAtomId(a.getId())];
                // set the mean x,y,z of the two involved atoms
                pylimer_tools::entities::Atom a1 = universe.getAtom(
                  pairsOfAtoms[pairOfAtom[universe.getIdxByAtomId(a.getId())]]
                    .first);
                pylimer_tools::entities::Atom a2 = universe.getAtom(
                  pairsOfAtoms[pairOfAtom[universe.getIdxByAtomId(a.getId())]]
                    .second);
                fb.setVertexPropertiesFromAtoms(thisVertexId, a1, a2);
                igraph_integer_t currentEdgeId = igraph_ecount(&fb.graph);
                igraph_add_edge(&fb.graph, previousVertexId, thisVertexId);
                fb.setBondPropertiesBasedOnChain(
                  chain, previousIdx, i, currentEdgeId, parentEdgeId);
                //
                previousIdx = i;
                previousVertexId = thisVertexId;
              }
            }
          }

          // close the chain
          igraph_integer_t currentEdgeId = igraph_ecount(&fb.graph);
          igraph_add_edge(
            &fb.graph,
            previousVertexId,
            endAtomIdxToVertexId[universe.getIdxByAtomId(lastAtom.getId())]);
          fb.setBondPropertiesBasedOnChain(chain,
                                           previousIdx,
                                           linedUpAtoms.size() - 1,
                                           currentEdgeId,
                                           parentEdgeId);

#ifndef NDEBUG
          fb.validateIgraphSpring(parentEdgeId);
#endif
          parentEdgeId += 1;
        }

// validate the creation
#ifndef NDEBUG
        for (size_t i = 0; i < endAtomIdxToVertexId.size(); ++i) {
          if (endAtomIdxToVertexId[i] < 0) {
            continue;
          }
          igraph_integer_t vertexId = endAtomIdxToVertexId[i];
          size_t savedAtomId = castToIgraphInt(
            igraph_cattribute_VAN(&fb.graph, "atom_id", vertexId));
          assert(savedAtomId == universe.getAtomIdByIdx(i));
          size_t degreeNow = fb.getVertexDegree(vertexId);
          size_t degreeBefore = universe.getVertexDegree(i);
          assert(degreeNow == degreeBefore);
        }
#endif

        // cleanup the graph
        // fb.removeSubfunctionalVertices();

        // convert the graph to the network usable for simulations
        fb.finaliseInitialisation();
        assert(fb.getNumExtraAtoms() == pairsOfAtoms.size());

        return fb;
      }

      //----------------------------------------------------------------
      // MARK: Simulation / Optimization Procedures
      //----------------------------------------------------------------

      /**
       * @brief Actually do run the simulation
       *
       * @param algorithm
       * @param maxNrOfSteps
       * @param xtol
       * @param ftol
       */
      void runForceRelaxation(
        long int maxNrOfSteps = 50000, // default: 10000
        double xtol = 1e-9,
        const double initialResidualToUse = -1.0,
        const StructureSimplificationMode simplificationMode =
          StructureSimplificationMode::NO_SIMPLIFICATION,
        const double inactiveRemovalCutoff = -1.0,
        bool doInnerIterations = false,
        LinkSwappingMode allowSlipLinksToPassEachOther =
          LinkSwappingMode::NO_SWAPPING,
        const int swappingFrequency = 10,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const int nrOfCrosslinkSwapsAllowedPerSliplink = -1);

      /**
       * @brief Publicly available method to convert graph to network
       *
       */
      void synchronise()
      {
        INVALIDARG_EXP_IFN(
          igraph_cattribute_GAB(&this->graph, "is_up_to_date"),
          "Graph must be up to date to be converted to network.");
        this->convertFromGraph();
      }

      /**
       * @brief Remove all "parent" springs that have no active "children"
       *
       * @param tolerance the acceptance tolerance, partial springs longer than
       * this are active
       */
      size_t removeInactiveParentEdges(double tolerance);

      /**
       * @brief Remove cross-linkers, springs and associated slip-links with the
       * scheme suggested by Andrei
       *
       * @param tolerance
       * @return size_t
       */
      size_t doRemovalAndreisWay(double tolerance);

      /**
       * @brief Remove chains that have two otherwise not connected ends
       * Mostly useful for phantom simulations, to compare methods
       *
       * The current algorithm to construct the network does keep free chains.
       *
       * @return size_t
       */
      size_t removeFreeChains();

      /**
       * @brief Remove chains with links with f = 1.
       *
       * NOTE: this function may result in even more chains with links with f
       * = 1.
       *
       * @return size_t
       */
      size_t removeDanglingChains();

      /**
       * @brief Remove cross-links which do not have any springs with a certain
       * minimum length
       *
       * @param net
       * @param displacements
       * @param springPartitions
       * @param tolerance
       */
      size_t removeInactiveCrosslinks(double tolerance = 1e-5);

      /**
       * @brief Remove all vertices (incl. edges!) with a functionality < 3 for
       * cross-links, < 4 for slip-links
       *
       * @param minCrosslinkFunctionalityToBeKept
       * @return size_t the number of removed vertices
       */
      size_t removeSubfunctionalVertices();

      /**
       * @brief Remove double listed springs from cross-links (if they have
       * length 0)
       *
       * @param net
       * @return size_t the nr of removed edges
       */
      size_t cleanupPrimaryLoopsInStructure();

      /**
       * @brief remove all vertices that don't have any connections
       */
      size_t removeOrphanedVertices();

      void relaxationLight(const size_t linkIdx,
                           const double oneOverSpringPartitionUpperLimit = 1.0)
      {
        assert(net.isUpToDate);
        if (net.linkIsSliplink[linkIdx]) {
          this->updateSpringPartition(linkIdx,
                                      oneOverSpringPartitionUpperLimit);
        }
        this->displaceToMeanPosition(linkIdx, oneOverSpringPartitionUpperLimit);
      };

      /**
       * @brief Replace the two springs traversinga a two-functional cross-links
       * with a single spring
       *
       *
       */
      size_t removeTwofunctionalLinks();

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
      void deformTo(const pylimer_tools::entities::Box& newBox)
      {
        if (!this->net.isUpToDate) {
          this->convertFromGraph();
        }
        this->universe.getBox().adjustCoordinatesTo(this->net.coordinates,
                                                    newBox);
        this->universe.getBox().adjustCoordinatesTo(
          this->net.springPartBoxOffset, newBox);
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
      pylimer_tools::entities::Universe getCrosslinkerVerse();

      int getDefaultNrOfChains() const { return this->defaultNrOfChains; }

      double getDefaultR0Square() const { return this->defaultR0Squared; }

      double getVolume() override
      {
        return this->universe.getBox().getVolume();
      }

      int getNrOfNodes() { return this->getNetwork().nrOfNodes; }

      int getNrOfLinks() { return igraph_vcount(&this->graph); }

      size_t getNumBonds() override { return this->getNrOfSprings(); }

      size_t getNumExtraBonds() override { return 0; }

      long int getNumBondsToForm() override { return 0; }

      size_t getNumAtoms() override { return this->getNrOfNodes(); }

      size_t getNumExtraAtoms() override
      {
        return this->getNrOfLinks() - this->getNrOfNodes();
      }

      int getNrOfSprings() { return this->getNetwork().nrOfSprings; }

      int getNrOfPartialSprings() { return igraph_ecount(&this->graph); }

      void setSpringContourLengths(const Eigen::VectorXd springsContourLengths)
      {
        INVALIDARG_EXP_IFN(springsContourLengths.size() ==
                             this->net.springsContourLength.size(),
                           "Contour length must have the correct dimensions.");
        if (!this->net.isUpToDate) {
          this->convertFromGraph();
        }
        this->net.springsContourLength = springsContourLengths;
        igraph_vector_t partialContour;
        igraph_vector_init(&partialContour, this->net.nrOfPartialSprings);
        igraph_vector_t fullContour;
        igraph_vector_init(&fullContour, this->net.nrOfPartialSprings);
        for (size_t i = 0; i < this->net.nrOfSprings; ++i) {
          for (size_t partialSpring : this->net.localToGlobalSpringIndex[i]) {
            igraph_vector_set(
              &partialContour,
              partialSpring,
              springsContourLengths[i] /
                static_cast<double>(
                  this->net.localToGlobalSpringIndex[i].size()));
            igraph_vector_set(
              &fullContour, partialSpring, springsContourLengths[i]);
          }
        }
        igraph_cattribute_EAN_setv(
          &this->graph, "local_contour_length", &partialContour);
        igraph_vector_destroy(&partialContour);
        igraph_cattribute_EAN_setv(
          &this->graph, "contour_length", &fullContour);
        igraph_vector_destroy(&fullContour);
      }

      std::vector<Eigen::ArrayXi> getIndependentCoordinateSets() const;

      std::pair<std::vector<Eigen::ArrayXi>, std::vector<Eigen::ArrayXi>>
      getHeuristicallyIndependentCoordinateSets() const;

      std::vector<Eigen::ArrayXi> getRandomCoordinateSets() const;

      void configAssumeBoxLargeEnough(bool assumption)
      {
        this->assumeBoxLargeEnough = assumption;
      }

      /**
       * @brief Get the Nr Of Active Nodes
       *
       * @param tolerance  the tolerance: springs under a certain length are
       * considered inactive
       * @return int
       */
      int getNrOfActiveNodes(double tolerance = 0.01,
                             int minimumNrOfActiveConnections = 2,
                             int maximumNrOfActiveConnections = -1,
                             bool usePartial = false)
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
      double getSolubleWeightFraction(double tolerance = 0.01)
      {
        return this->computeSolubleWeightFraction(this->currentSpringDistances,
                                                  tolerance);
      }

      /**
       * @brief Get the Dangling Weight Fraction
       *
       * @param tolerance
       * @return double
       */
      double getDanglingWeightFraction(double tolerance = 0.01)
      {
        return this->computeDanglingWeightFraction(this->currentSpringDistances,
                                                   tolerance);
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
        double tolerance = 0.01);

      /**
       * @brief Compute the weight fraction of non-active springs
       *
       * @param net
       * @param springDistances
       * @param tolerance
       * @return double
       */
      double computeDanglingWeightFraction(
        const Eigen::VectorXd& springDistances,
        const double tolerance = 0.01)
      {
        RUNTIME_EXP_IFN(
          this->getNetwork().isUpToDate,
          "Cannot compute dangling weight fraction with out-of-date network.");
        if (this->net.nrOfSprings * 3 != springDistances.size()) {
          throw std::invalid_argument(
            "Spring distances and network don't match");
        }
        if (this->net.nrOfSprings < 1) {
          return 1.;
        }
        // find all active springs
        Eigen::ArrayXb activeSprings =
          this->findActiveSprings(springDistances, tolerance);
        if (activeSprings.count() == 0) {
          return 1.;
        }
        // as of now, the springsContourLength is equal to the number of bonds
        // from cross-link to cross-link. therefore, the number of atoms of each
        // of these springs is one less
        Eigen::ArrayXd allActiveAtomsPerChains =
          activeSprings.cast<double>() *
          (this->net.springsContourLength.array() -
           Eigen::ArrayXd::Ones(this->net.nrOfSprings));
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
        const Eigen::VectorXd& springDistances,
        const double tolerance = 0.01)
      {
        if (this->net.nrOfSprings * 3 != springDistances.size()) {
          throw std::invalid_argument(
            "Spring distances and network don't match");
        }
        if (this->net.nrOfSprings < 1) {
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
        double tolerance = 0.01,
        int minimumNrOfActiveConnections = 2,
        int maximumNrOfActiveConnections = -1,
        bool usePartial = false);

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
      Eigen::VectorXi getNrOfActiveSpringsConnected(double tolerance = 0.01);

      /**
       * @brief Get the Nr Of Active Springs connected to each node
       *
       * @param tolerance the tolerance: springs under a certain length are
       * considered inactive
       * @return Eigen::VectorXi
       */
      Eigen::VectorXi getNrOfActivePartialSpringsConnected(
        double tolerance = 0.01);

      /**
       * @brief Get the Nr Of Active Springs object
       *
       * @param tolerance the tolerance: springs under a certain length are
       * considered inactive
       * @return int
       */
      int getNrOfActiveSprings(double tolerance = 0.01) const
      {
        assert(this->net.isUpToDate);
        Eigen::ArrayXb result = Eigen::ArrayXb::Constant(
          this->currentPartialSpringDistances.size() / 3, false);
        for (size_t i = 0; i < this->net.nrOfPartialSprings; ++i) {
          size_t springIdx = this->net.partialToFullSpringIndex[i];
          result[springIdx] =
            result[springIdx] ||
            this->currentPartialSpringDistances.segment(3 * i, 3)
                .squaredNorm() > tolerance;
        }

        return (result == true).count();
      }

      /**
       * @brief Get the Nr Of Active Springs object
       *
       * @param tolerance the tolerance: springs under a certain length are
       * considered inactive
       * @return int
       */
      int getNrOfActivePartialSprings(double tolerance = 0.01) const
      {
        return this->countNrOfActiveSprings(this->currentPartialSpringDistances,
                                            tolerance);
      }

      /**
       * @brief Get the Average Spring Length at the current step
       *
       * @return double
       */
      double getAverageSpringLength();

      /**
       * @brief Compute the stress tensor
       *
       * @param net
       * @param u
       * @param loopTol
       * @return std::array<std::array<double, 3>, 3>
       */
      std::array<std::array<double, 3>, 3> evaluateStressTensorLinkBased(
        const double kappa0 = 1.0,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const bool xlinksOnly = false) const;

      /**
       * @brief Get the Pressure
       *
       * @return double
       */
      double getPressure() const { return this->evaluatePressure(); }

      /**
       * @brief Get the Gamma Factor at the current step
       *
       * @param r02 the melt <R_0^2>, for phantom = Nb^2
       * @param nrOfChains the nr of chains to average over (can be different
       * from the nr of springs thanks to omitted free chains or primary loops)
       * @return double
       */
      double getGammaFactor(double r02 = -1.0, int nrOfChains = -1) const
      {
        if (r02 < 0) {
          r02 = this->defaultR0Squared;
        }
        if (nrOfChains < 1) {
          nrOfChains = this->defaultNrOfChains;
        }

        return this->evaluateGammaFactor(
          this->currentSpringDistances, r02, nrOfChains);
      };

      int getNrOfIterations() const { return this->nrOfStepsDone; }

      ExitReason getExitReason() const { return this->exitReason; }
      /**
       * @brief Compute the spring lenghts
       *
       * @param net the network to do the computation for
       * @return Eigen::VectorXd
       */
      Eigen::VectorXd evaluateSpringDistances();

      Eigen::VectorXd evaluatePartialSpringDistances() const;

      /**
       * @brief Compute the distance between two links
       *
       * @param net
       * @param linkIndexA
       * @param linkIndexB
       * @param is2D
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d evaluateDistanceBetween(const size_t linkIndexA,
                                              const size_t linkIndexB,
                                              bool is2D = false) const
      {
        Eigen::Vector3d distances = this->getVertexCoordinates(linkIndexB) -
                                    this->getVertexCoordinates(linkIndexA);

        // Possibly improvable PBC
        this->universe.getBox().handlePBC<Eigen::Vector3d>(distances);

        if (is2D) {
          distances[2] = 0.0;
        }

        return distances;
      };

      /**
       * @brief Compute one spring length
       *
       * @param net
       * @param springIdx
       * @param is2D
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d evaluatePartialSpringDistance(const size_t springIdx,
                                                    bool is2d = false) const
      {
        return this->evaluatePartialSpringDistanceTo(
          springIdx, this->net.springIndexB(springIdx), is2d);
      }

      /**
       * @brief Compute one spring length, in a specific direction
       *
       * @param net
       * @param springIdx
       * @param linkIdx the vector "target"
       * @param is2D
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d evaluatePartialSpringDistanceTo(const size_t springIdx,
                                                      const size_t linkIdx,
                                                      bool is2d = false) const
      {
        return this->computeEdgeDistanceTo(springIdx, linkIdx);
      }

      /**
       * @brief Compute one spring length, in a specific direction
       *
       * @param net
       * @param springIdx
       * @param linkIdx the vector "source"
       * @param is2D
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d evaluatePartialSpringDistanceFrom(const size_t springIdx,
                                                        const size_t linkIdx,
                                                        bool is2d = false) const
      {
        return -1. *
               this->evaluatePartialSpringDistanceTo(springIdx, linkIdx, is2d);
      }

      bool validateNetwork()
      {
        return this->validateNetwork(this->getNetwork(),
                                     this->currentSpringPartitionsVec);
      }

      bool validateNetwork(const ForceBalanceNetwork& net)
      {
        return this->validateNetwork(net, this->currentSpringPartitionsVec);
      }

      bool validateNetwork(const ForceBalanceNetwork& net,
                           const Eigen::VectorXd& springPartitions);

      void debugParentEdge(const size_t parentEdgeId)
      {
#ifndef VERBOSE_DEBUG
        return;
#endif
        igraph_vector_t parentEdges;
        igraph_vector_init(&parentEdges, 0);
        igraph_cattribute_EANV(&this->graph,
                               "parent_edge",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &parentEdges);
        igraph_vector_t partitionFraction;
        igraph_vector_init(&partitionFraction, 0);
        igraph_cattribute_EANV(&this->graph,
                               "partition_fraction",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &partitionFraction);
        igraph_vector_t contourLengths;
        igraph_vector_init(&contourLengths, 0);
        igraph_cattribute_EANV(&this->graph,
                               "contour_length",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &contourLengths);

        std::cout << "Debugging edge " << parentEdgeId << std::endl;
        for (igraph_integer_t i = 0; i < igraph_vector_size(&parentEdges);
             ++i) {
          if (castToIgraphInt(igraph_vector_get(&parentEdges, i)) ==
              parentEdgeId) {
            igraph_integer_t from, to;
            igraph_edge(&this->graph, i, &from, &to);
            std::cout << i << " (" << igraph_vector_get(&contourLengths, i)
                      << "):\t" << from << "\t" << to << "\t("
                      << igraph_vector_get(&partitionFraction, i) << ")"
                      << std::endl;
          }
        }
        std::cout << "That's all." << std::endl;

        igraph_vector_destroy(&parentEdges);
        igraph_vector_destroy(&partitionFraction);
        igraph_vector_destroy(&contourLengths);
      }

      bool validateIgraphSprings()
      {
        igraph_vector_t parentEdges;
        igraph_vector_init(&parentEdges, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "parent_edge",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &parentEdges);

        std::unordered_set<igraph_integer_t> visitedSprings;

        for (size_t i = 0; i < igraph_ecount(&this->graph); i++) {
          igraph_integer_t spring =
            castToIgraphInt(igraph_vector_get(&parentEdges, i));
          if (!visitedSprings.contains(spring)) {
            this->validateIgraphSpring(spring);
            visitedSprings.insert(spring);
          }
        }

        igraph_vector_destroy(&parentEdges);
        return true;
      }

      bool validateIgraphSpring(const size_t parentEdgeId)
      {
        this->debugParentEdge(parentEdgeId);
        igraph_vector_t parentEdges;
        igraph_vector_init(&parentEdges, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "parent_edge",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &parentEdges);

        double partitionSum = 0.;
        std::string partialSpringString = "";
        igraph_vector_int_t edgesOnPath;
        igraph_vector_int_init(&edgesOnPath, 0);
        for (igraph_integer_t i = 0; i < igraph_vector_size(&parentEdges);
             ++i) {
          if (castToIgraphInt(igraph_vector_get(&parentEdges, i)) ==
              parentEdgeId) {
            igraph_vector_int_push_back(&edgesOnPath, i);
            double thisPartition =
              igraph_cattribute_EAN(&this->graph, "partition_fraction", i);
            partitionSum += thisPartition;

            partialSpringString += std::to_string(thisPartition) + ", ";
            if (!APPROX_WITHIN(thisPartition, 0., 1., 1e-5)) {
              RUNTIME_EXP_IFN(APPROX_WITHIN(thisPartition, 0., 1., 1e-5),
                              "Expected partition to be between 1 and 0, got " +
                                std::to_string(thisPartition) + " for spring " +
                                std::to_string(i) + ".");
            }
          }
        }

        if (!(APPROX_EQUAL(partitionSum, 1.0, 1e-5))) {
          RUNTIME_EXP_IFN(
            APPROX_EQUAL(partitionSum, 1.0, 1e-5),
            "Expected partition sum to be closer to 1., got " +
              std::to_string(partitionSum) + " for parent spring " +
              std::to_string(parentEdgeId) +
              ". Partial Springs are: " + partialSpringString + ".");
        }

        // validate that we can make one chain
        igraph_t subgraph;
        igraph_empty(&subgraph, 0, IGRAPH_UNDIRECTED);
        igraph_subgraph_from_edges(
          &this->graph, &subgraph, igraph_ess_vector(&edgesOnPath), false);

        igraph_vector_int_t verticesOnPath;
        igraph_vector_int_init(&verticesOnPath, 0);
        if (igraph_eulerian_path(&subgraph, &edgesOnPath, &verticesOnPath)) {
          throw std::runtime_error(
            "Failed to find simple path for the spring " +
            std::to_string(parentEdgeId));
        }
        if (castToIgraphInt(igraph_cattribute_VAN(
              &subgraph, "type", igraph_vector_int_get(&verticesOnPath, 0))) !=
            this->crosslinkerType) {
          throw std::runtime_error("Parent edge " +
                                   std::to_string(parentEdgeId) +
                                   " does not start with a cross-link");
        }
        if (castToIgraphInt(igraph_cattribute_VAN(
              &subgraph,
              "type",
              igraph_vector_int_get(&verticesOnPath,
                                    igraph_vector_int_size(&verticesOnPath) -
                                      1))) != this->crosslinkerType) {
          throw std::runtime_error("Parent edge " +
                                   std::to_string(parentEdgeId) +
                                   " does not end with a cross-link");
        }
        this->printVerticesOnPath(&verticesOnPath);

        igraph_vector_int_destroy(&edgesOnPath);
        igraph_destroy(&subgraph);

        igraph_vector_destroy(&parentEdges);
        return true;
      }

      void printVerticesOnPath(igraph_vector_int_t* vertices)
      {
#ifndef VERBOSE_DEBUG
        return;
#endif
        std::cout << "\t";
        std::cout << igraph_vector_int_get(vertices, 0);
        for (igraph_integer_t i = 1; i < igraph_vector_int_size(vertices);
             ++i) {
          std::cout << " - " << igraph_vector_int_get(vertices, i);
        }
        std::cout << "\n";
      }

      ForceBalanceNetwork getNetwork()
      {
        if (!this->net.isUpToDate) {
          this->convertFromGraph();
        }
        return this->net;
      }

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
        return this->evaluateForceOnLink(
          index, debugNrSpringsVisited, 1.0, oneOverSpringPartitionUpperLimit);
      }

      /**
       * @brief Assemble all indices of partial springs for a particular
       * slip-link
       *
       * @param linkIdx
       * @return std::vector<size_t>
       */
      std::vector<size_t> getSpringpartitionIndicesOfSliplink(
        const size_t linkIdx)
      {
        INVALIDARG_EXP_IFN(
          this->linkIsSlipLink(linkIdx),
          "Only slip-links may be asked for their partition indices, link " +
            std::to_string(linkIdx) + " is not one.");
        std::vector<size_t> indices;
        indices.reserve(4);

        igraph_vector_int_t edgesOfLink;
        igraph_vector_int_init(&edgesOfLink, 2);
        igraph_incident(&this->graph, &edgesOfLink, linkIdx, IGRAPH_ALL);
        assert(igraph_vector_int_size(&edgesOfLink) == 4);

        for (size_t i = 0; i < igraph_vector_int_size(&edgesOfLink); ++i) {
          indices.push_back(igraph_vector_int_get(&edgesOfLink, i));
        }

        igraph_vector_int_destroy(&edgesOfLink);

        return indices;
      };

      /**
       * @brief Updates the partition/parametrisation of a spring around one
       * link
       *
       */
      double updateSpringPartition(
        const igraph_integer_t linkIdx,
        double oneOverSpringPartitionUpperLimit = 1.0,
        bool allowSlipLinksToPassEachOther = false)
      {
        INVALIDARG_EXP_IFN(linkIdx < igraph_vcount(&this->graph),
                           "Link to update needs to be in the list");
        INVALIDARG_EXP_IFN(this->linkIsSlipLink(linkIdx),
                           "Only slip-links may slip along a spring, link " +
                             std::to_string(linkIdx) + " is not one.");

        igraph_vector_int_t edgesOfVertex;
        igraph_vector_int_init(&edgesOfVertex, 1);
        igraph_incident(&this->graph, &edgesOfVertex, linkIdx, IGRAPH_ALL);

        double residualNorm = 0.0;
        int residualNormContributions = 0;
        std::unordered_set<igraph_integer_t> handledEdges;
        for (size_t i = 0; i < igraph_vector_int_size(&edgesOfVertex); ++i) {
          igraph_integer_t edgeId = igraph_vector_int_get(&edgesOfVertex, i);
          if (handledEdges.contains(edgeId)) {
            continue;
          }
          igraph_integer_t otherEdgeId =
            this->getOtherRailEdgeId(linkIdx, edgeId);
          if (handledEdges.contains(otherEdgeId)) {
            // primary loop edge /!
            igraph_integer_t from, to;
            igraph_edge(&this->graph, otherEdgeId, &from, &to);
            assert(from == to);
          }
          // remember so we don't do them twice
          handledEdges.insert(edgeId);
          handledEdges.insert(otherEdgeId);

          Eigen::Vector3d vecBack =
            this->evaluatePartialSpringDistanceFrom(edgeId, linkIdx);
          Eigen::Vector3d vecForward =
            this->evaluatePartialSpringDistanceFrom(otherEdgeId, linkIdx);

          double distanceBack = vecBack.squaredNorm();
          double distanceForward = vecForward.squaredNorm();
          double idealValue = 1. / (1. + sqrt(distanceForward / distanceBack));
          if (distanceBack <= 0.0) {
            idealValue = 0.0; // TODO: really?
          }

          double currentS =
            igraph_cattribute_EAN(&this->graph, "partition_fraction", edgeId);
          double nextS = igraph_cattribute_EAN(
            &this->graph, "partition_fraction", otherEdgeId);
          const double N =
            igraph_cattribute_EAN(&this->graph, "contour_length", edgeId);
          assert(N == igraph_cattribute_EAN(
                        &this->graph, "contour_length", otherEdgeId));

          const double l = (currentS + nextS);

          if (oneOverSpringPartitionUpperLimit > 0.) {
            // TODO: sketch theory why this should/not be necessary!!!
            const double limit =
              std::clamp(1. / (oneOverSpringPartitionUpperLimit *
                               (nextS + currentS) * (N)),
                         0.,
                         1.);
            idealValue = std::clamp(idealValue, limit, 1. - limit);
          }
          double newS = idealValue * l;
          double idealValueM1 = (1. - idealValue);
          double complementaryS = (1. - idealValue) * l;
          double localResidualNorm = 0.0;
          residualNormContributions += 2;
          if (idealValue > 0.0 && idealValue < 1.0) {
            double idealValue2 = idealValue * idealValue;
            double idealValueM12 = idealValueM1 * idealValueM1;
            // way too complicated expression to solve subtraction
            // truncation issues?
            localResidualNorm =
              (std::fma(distanceBack,
                        idealValueM12,
                        -1. * distanceForward * idealValue2)) /
              (idealValueM12 * idealValue2);
          }

          localResidualNorm /= (N * l);

          // lots of validation
          RUNTIME_EXP_IFN(
            APPROX_WITHIN(newS + complementaryS, 0., 1., 1e-9),
            "Require newS + complementaryS to be within 0, 1, got " +
              std::to_string(newS + complementaryS) + " from " +
              std::to_string(newS) + " and " + std::to_string(complementaryS) +
              " with ideal = " + std::to_string(idealValue) + " of " +
              std::to_string(nextS + currentS) + " for link " +
              std::to_string(linkIdx) +
              ". Diff: " + std::to_string(1. - (newS + complementaryS)) + ".");
          RUNTIME_EXP_IFN(
            APPROX_EQUAL(nextS + currentS, newS + complementaryS, 1e-9),
            "Require nextS + currentS == newS + complementaryS, got " +
              std::to_string(nextS + currentS) + " vs. " +
              std::to_string(newS + complementaryS) + " from " +
              std::to_string(nextS) + " and " + std::to_string(currentS) +
              ", " + std::to_string(newS) + " and " +
              std::to_string(complementaryS) + ". Diff: " +
              std::to_string((nextS + currentS) - (newS + complementaryS)) +
              ".");
          RUNTIME_EXP_IFN(
            APPROX_WITHIN(nextS + currentS, 0., 1., 1e-9),
            "Require nextS + currentS to be within 0, 1, got " +
              std::to_string(nextS + currentS) + " from " +
              std::to_string(nextS) + " and " + std::to_string(currentS) +
              ". Diff: " + std::to_string(1. - (nextS + currentS)) + ".");
          RUNTIME_EXP_IFN(nextS >= -0.00000000000001,
                          "nextS must be >= 0., got " + std::to_string(nextS) +
                            " from " + std::to_string(nextS) + " and " +
                            std::to_string(currentS) + ", " +
                            std::to_string(newS) + " and " +
                            std::to_string(complementaryS) + ".");
          RUNTIME_EXP_IFN(complementaryS >= -0.00000000000001,
                          "complementaryS must be >= 0., got " +
                            std::to_string(complementaryS) + " from " +
                            std::to_string(nextS) + " and " +
                            std::to_string(currentS) + ", " +
                            std::to_string(newS) + " and " +
                            std::to_string(complementaryS) + ".");
          residualNorm += localResidualNorm * localResidualNorm;

          // actually set the values
          igraph_cattribute_EAN_set(
            &this->graph, "partition_fraction", edgeId, newS);
          igraph_cattribute_EAN_set(
            &this->graph, "partition_fraction", otherEdgeId, complementaryS);
        }

        igraph_vector_int_destroy(&edgesOfVertex);
        // assert(residualNormContributions == 4 || cautionSecondaryLoop);
        return residualNorm;
      };

      /**
       * @brief Loop all slip-links and move them if appropriate to other
       * springs
       *
       * @param oneOverSpringPartitionUpperLimit
       */
      void moveSlipLinksToTheirBestBranch(
        const double oneOverSpringPartitionUpperLimit,
        const int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
        const bool respectLoops = true)
      {
        igraph_integer_t numVertices = igraph_vcount(&this->graph);
        igraph_vector_t linkType;
        igraph_vector_init(&linkType, this->net.nrOfLinks);
        igraph_cattribute_VANV(
          &this->graph, "type", igraph_vss_all(), &linkType);
        for (igraph_integer_t linkIdx = 0; linkIdx < numVertices; ++linkIdx) {
          // check this slip-link
          // std::cout << "Moving slip-link " << sliplinkIdx << " to its best
          // branch"
          //           << std::endl;
          if (castToIgraphInt(igraph_vector_get(&linkType, linkIdx)) ==
              this->slipLinkType) {
            this->moveSlipLinkToItsBestBranch(
              linkIdx,
              oneOverSpringPartitionUpperLimit,
              nrOfCrosslinkSwapsAllowedPerSliplink,
              respectLoops);
          }
          // this->validateNetwork(net, springPartitions);
        }
        this->validateNetwork();
      };

      /**
       * @brief Decide whether the spring fraction is short enough to have a
       * spring deserve swapping
       *
       * @param edgeId the edge to check for its fraction
       * @param oneOverSpringPartitionUpperLimit
       * @return true
       * @return false
       */
      bool partialSpringRequestsSwapping(
        const igraph_integer_t edgeId,
        const double oneOverSpringPartitionUpperLimit)
      {
        assert(edgeId < igraph_ecount(&this->graph));
        const double N =
          igraph_cattribute_EAN(&this->graph, "contour_length", edgeId);
        const double swappableCutoff =
          (oneOverSpringPartitionUpperLimit > 0.)
            ? 1. / (N - 1. / oneOverSpringPartitionUpperLimit)
            : 1e-12;

        igraph_integer_t from, to;
        igraph_edge(&this->graph, edgeId, &from, &to);

        return (igraph_cattribute_EAN(&this->graph,
                                      "partition_fraction",
                                      edgeId) < swappableCutoff) &&
               from != to;
      }

      /**
       * @brief Move a slip-link if appropriate to other springs
       *
       * @param oneOverSpringPartitionUpperLimit
       */
      bool moveSlipLinkToItsBestBranch(
        size_t slipLinkIdx,
        const double oneOverSpringPartitionUpperLimit,
        const int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
        const bool respectLoops = true)
      {
        if (castToIgraphInt(igraph_cattribute_VAN(
              &this->graph, "type", slipLinkIdx)) != this->slipLinkType) {
          INVALIDARG_EXP_IFN(castToIgraphInt(igraph_cattribute_VAN(
                               &this->graph, "type", slipLinkIdx)) ==
                               this->slipLinkType,
                             "Passed slip-link must be one.");
        }
        RUNTIME_EXP_IFN(igraph_cattribute_GAB(&this->graph, "is_up_to_date"),
                        "Should not move slip-link if graph is not up-to-date");

        igraph_vector_int_t edgesOfLink;
        igraph_vector_int_init(&edgesOfLink, 4);
        igraph_incident(&this->graph, &edgesOfLink, slipLinkIdx, IGRAPH_ALL);

        if (igraph_vector_int_size(&edgesOfLink) != 4) {
          igraph_vector_int_destroy(&edgesOfLink);
          return false;
        }

        bool didSwap = false;
        for (size_t i = 0; i < igraph_vector_int_size(&edgesOfLink); ++i) {
          const igraph_integer_t edgeId =
            igraph_vector_int_get(&edgesOfLink, i);
          if (this->partialSpringRequestsSwapping(
                edgeId, oneOverSpringPartitionUpperLimit)) {
            didSwap =
              this->swapSlipLinkReversibly(edgeId,
                                           oneOverSpringPartitionUpperLimit,
                                           nrOfCrosslinkSwapsAllowedPerSliplink,
                                           respectLoops);
          }
          if (didSwap) {
#ifndef NDEBUG
            this->validateIgraphSpring(castToIgraphInt(
              igraph_cattribute_EAN(&this->graph, "parent_edge", edgeId)));
#endif
            break;
          }
        }

        igraph_vector_int_destroy(&edgesOfLink);
        return didSwap;
      };

      /**
       * @brief Loop all springs, swap slip-links on them if they are close
       * enough
       *
       * @param net
       */
      void swapSlipLinksInclXlinks(double oneOverSpringPartitionUpperLimit,
                                   const bool respectLoops = true)
      {
        for (long int i = igraph_ecount(&this->graph) - 1; i >= 0; --i) {
          // problem: when moving slip-links, the edges ids change
          // "solution": edge ids increase -> if we remove one, we might
          // re-visit the previous, but otherwise, things should be fine
          if (this->partialSpringRequestsSwapping(
                i, oneOverSpringPartitionUpperLimit)) {
            if (this->springInvolvesCrossLink(i)) {
              this->rotateSlipLinkAroundCrosslink(
                i, oneOverSpringPartitionUpperLimit, respectLoops);
            } else {
              this->swapSlipLinksOfEdge(i);
            }
          }
        }
      };

      /**
       * @brief Loop all springs, swap slip-links on them if they are close
       * enough
       *
       * @param net
       */
      void swapSlipLinks(double oneOverSpringPartitionUpperLimit)
      {
        for (long int i = igraph_ecount(&this->graph) - 1; i >= 0; --i) {
          // problem: when moving slip-links, the edges ids change
          // "solution": edge ids increase -> if we remove one, we might
          // re-visit the previous, but otherwise, things should be fine
          if (this->partialSpringRequestsSwapping(
                i, oneOverSpringPartitionUpperLimit)) {
            if (!this->springInvolvesCrossLink(i)) {
              this->swapSlipLinksOfEdge(i);
            }
          }
        }
      };

      /**
       * @brief Swap the two links on one partial spring
       */
      void swapSlipLinksOfEdge(const size_t partialSpringIdx)
      {
        assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));

#ifndef NDEBUG
        this->validateIgraphSpring(castToIgraphInt(igraph_cattribute_EAN(
          &this->graph, "parent_edge", partialSpringIdx)));
#endif

        igraph_integer_t linkIdx1, linkIdx2;
        igraph_edge(&this->graph, partialSpringIdx, &linkIdx1, &linkIdx2);
        INVALIDARG_EXP_IFN(linkIdx1 != linkIdx2,
                           "Cannot swap link with itself: got " +
                             std::to_string(linkIdx1) + " and " +
                             std::to_string(linkIdx2) + ".");
        INVALIDARG_EXP_IFN(
          castToIgraphInt(igraph_cattribute_VAN(
            &this->graph, "type", linkIdx1)) == this->slipLinkType,
          "Only partial springs with only slip-links allow swapping.");
        INVALIDARG_EXP_IFN(
          castToIgraphInt(igraph_cattribute_VAN(
            &this->graph, "type", linkIdx2)) == this->slipLinkType,
          "Only partial springs with only slip-links allow swapping.");

        // figure out the path of the new edges
        igraph_integer_t otherEdge1 =
          this->getOtherRailEdgeId(linkIdx1, partialSpringIdx);
        igraph_integer_t otherEdge2 =
          this->getOtherRailEdgeId(linkIdx2, partialSpringIdx, otherEdge1);
        if (otherEdge1 == otherEdge2) {
          otherEdge1 =
            this->getOtherRailEdgeId(linkIdx1, partialSpringIdx, otherEdge2);
        }
        igraph_integer_t parentEdge = castToIgraphInt(
          igraph_cattribute_EAN(&this->graph, "parent_edge", partialSpringIdx));
        assert(parentEdge == castToIgraphInt(igraph_cattribute_EAN(
                               &this->graph, "parent_edge", otherEdge1)));
        assert(parentEdge == castToIgraphInt(igraph_cattribute_EAN(
                               &this->graph, "parent_edge", otherEdge2)));
        // TODO: this is an issue with loops that loop back in some way
        // there, the eulerian path cannot decide what is forward, what is
        // backward
        if (otherEdge1 == otherEdge2) {
          std::cerr << "WARNING: the connectivity could not be fully restored."
                    << std::endl;
          return;
        }
        // ...actually add the new edges
        igraph_integer_t newEdge1 = this->createEdge(
          linkIdx1,
          this->getOtherEdgePartner(otherEdge2, linkIdx2),
          parentEdge,
          igraph_cattribute_EAN(&this->graph, "partition_fraction", otherEdge2),
          igraph_cattribute_EAN(&this->graph, "contour_length", otherEdge2),
          this->getBondBoxOffsetForEdgeTo(otherEdge1, linkIdx1) +
            this->getBondBoxOffsetForEdgeFrom(partialSpringIdx, linkIdx1));

        igraph_integer_t newEdge2 = this->createEdge(
          linkIdx2,
          this->getOtherEdgePartner(otherEdge1, linkIdx1),
          parentEdge,
          igraph_cattribute_EAN(&this->graph, "partition_fraction", otherEdge1),
          igraph_cattribute_EAN(&this->graph, "contour_length", otherEdge1),
          this->getBondBoxOffsetForEdgeTo(otherEdge2, linkIdx2) +
            this->getBondBoxOffsetForEdgeFrom(partialSpringIdx, linkIdx2));

        // actually "cut" by removing the old edge
        std::vector<igraph_integer_t> eToRemove = { otherEdge1, otherEdge2 };
        this->removePartialSprings(eToRemove);

#ifndef NDEBUG
        this->validateIgraphSpring(parentEdge);
#endif

        this->net.isUpToDate = false;
      };

      /**
       * @brief Swap slip- or cross-links along a partial spring, iff the
       * slip-link may still swap, and iff the swap is MC favourable
       *
       * @param partialSpringIdx edge idx
       * @param oneOverSpringPartitionUpperLimit
       * @param nrOfCrosslinkSwapsAllowedPerSliplink
       * @param respectLoops
       * @return true if the swap was successful
       * @return false if the swap was not successful
       */
      bool swapSlipLinkReversibly(
        const size_t partialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const int nrOfCrosslinkSwapsAllowedPerSliplink = -1,
        const bool respectLoops = true)
      {
        // analyse spring
        if (this->springInvolvesCrossLink(partialSpringIdx)) {
          igraph_integer_t from, to;
          igraph_edge(&this->graph, partialSpringIdx, &from, &to);
          igraph_integer_t slipLinkIdx = from;
          if (castToIgraphInt(igraph_cattribute_VAN(
                &this->graph, "type", from)) == this->crosslinkerType) {
            slipLinkIdx = to;
          }
          assert(castToIgraphInt(igraph_cattribute_VAN(
                   &this->graph, "type", slipLinkIdx)) == this->slipLinkType);

          // first check if allowed.
          if ((nrOfCrosslinkSwapsAllowedPerSliplink < 0) ||
              (castToIgraphInt(igraph_cattribute_VAN(
                 &this->graph, "num_link_swaps", slipLinkIdx)) <
               nrOfCrosslinkSwapsAllowedPerSliplink)) {
            bool didSwap = this->swapSlipLinkWithXlinkReversibly(
              partialSpringIdx, oneOverSpringPartitionUpperLimit, respectLoops);
            return didSwap;
          }
          return false;
        } else {
          return this->swapSlipLinksReversibly(
            partialSpringIdx, oneOverSpringPartitionUpperLimit);
        }
      };

      /**
       * @brief Move a slip-link from one spring attached to a cross-link to
       * another spring attached to the same cross-link
       *
       * @param partialSpringIdx
       */
      void rotateSlipLinkAroundCrosslink(
        const size_t partialSpringIdx,
        double oneOverSpringPartitionUpperLimit = 1.0,
        const bool respectLoops = true)
      {
        igraph_integer_t from, to;
        igraph_edge(&this->graph, partialSpringIdx, &from, &to);

        igraph_integer_t xlinkIdx =
          castToIgraphInt(igraph_cattribute_VAN(&this->graph, "type", from)) ==
              this->crosslinkerType
            ? from
            : to;
        igraph_integer_t slipLinkIdx = xlinkIdx == from ? to : from;
        assert(castToIgraphInt(igraph_cattribute_VAN(
                 &this->graph, "type", xlinkIdx)) == this->crosslinkerType);

        this->unlinkSlipLinkFromRail(slipLinkIdx, partialSpringIdx);

        igraph_vector_int_t edgesOfCrossLink;
        igraph_vector_int_init(&edgesOfCrossLink, 4);
        igraph_incident(&this->graph, &edgesOfCrossLink, xlinkIdx, IGRAPH_ALL);
        RUNTIME_EXP_IFN(igraph_vector_int_size(&edgesOfCrossLink) > 0,
                        "Expected to find more than one edge on slip-link");
        int nextEdgeIndex = igraph_vector_int_get(&edgesOfCrossLink, 0);
        // use the smallest edge index as the new rail
        for (size_t i = 1; i < igraph_vector_int_size(&edgesOfCrossLink); ++i) {
          if (igraph_vector_int_get(&edgesOfCrossLink, i) < nextEdgeIndex) {
            nextEdgeIndex = igraph_vector_int_get(&edgesOfCrossLink, i);
          }
        }

        igraph_vector_int_destroy(&edgesOfCrossLink);

        this->insertSlipLinkIntoRail(slipLinkIdx, nextEdgeIndex);
      };

      /**
       * @brief Displace all links to their mean position
       *
       * @param oneOverSpringPartitionUpperLimit
       */
      void displaceLinksToMeanPosition(
        const double oneOverSpringPartitionUpperLimit = 1.0)
      {
        for (igraph_integer_t i = 0; i < igraph_vcount(&this->graph); ++i) {
          this->displaceToMeanPosition(i, oneOverSpringPartitionUpperLimit);
        }
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
      double displaceToMeanPosition(
        const size_t linkIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0);

      double getDisplacementResidualNorm(
        double oneOverSpringPartitionUpperLimit = 1.0) const;

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
        return this->getStressTensor(1.0);
      }

      Eigen::Matrix3d getStressTensor(double oneOverSpringPartitionUpperLimit)
      {
        std::array<std::array<double, 3>, 3> stressTensor =
          this->evaluateStressTensor(this->kappa,
                                     oneOverSpringPartitionUpperLimit);
        Eigen::Matrix3d result = Eigen::Matrix3d::Zero();
        for (size_t i = 0; i < 3; ++i) {
          for (size_t j = 0; j < 3; ++j) {
            result(i, j) = stressTensor[i][j];
          }
        }
        return result;
      }

      Eigen::Matrix3d getStressTensorLinkBased(
        double oneOverSpringPartitionUpperLimit = 1.0)
      {
        std::array<std::array<double, 3>, 3> stressTensor =
          this->evaluateStressTensorLinkBased(this->kappa,
                                              oneOverSpringPartitionUpperLimit);
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
       * @brief
       *
       * CAUTION: make sure the graph is up to date!
       *
       * @param vertexId
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d getVertexCoordinates(igraph_integer_t vertexId) const
      {
        assert(vertexId < igraph_vcount(&this->graph));
        Eigen::Vector3d coordinates;
        coordinates << igraph_cattribute_VAN(&this->graph, "x", vertexId),
          igraph_cattribute_VAN(&this->graph, "y", vertexId),
          igraph_cattribute_VAN(&this->graph, "z", vertexId);
        return coordinates;
      }

      /**
       * @brief Set the Vertex Coordinates in the graph
       *
       * @param vertexId
       * @param coordinates
       */
      void setVertexCoordinates(const igraph_integer_t vertexId,
                                const Eigen::Vector3d& coordinates)
      {
        assert(vertexId < igraph_vcount(&this->graph));
        igraph_cattribute_VAN_set(&this->graph, "x", vertexId, coordinates[0]);
        igraph_cattribute_VAN_set(&this->graph, "y", vertexId, coordinates[1]);
        igraph_cattribute_VAN_set(&this->graph, "z", vertexId, coordinates[2]);
      }

      /**
       * @brief Check whether a given vertex is a slip-link
       *
       * @param vertexId
       * @return true
       * @return false
       */
      bool linkIsSlipLink(const igraph_integer_t vertexId) const
      {
        assert(vertexId < igraph_vcount(&this->graph));
        return (castToIgraphInt(igraph_cattribute_VAN(
                  &this->graph, "type", vertexId)) == this->slipLinkType);
      }

      /**
       * @brief Check whether a given vertex is a slip-link
       *
       * @param vertexId
       * @return true
       * @return false
       */
      bool linkIsCrossLink(const igraph_integer_t vertexId) const
      {
        assert(vertexId < igraph_vcount(&this->graph));
        return (castToIgraphInt(igraph_cattribute_VAN(
                  &this->graph, "type", vertexId)) == this->crosslinkerType);
      }

      /**
       * @brief Check whether a given edge has one end that is a cross-link
       *
       *
       * @param edgeId
       * @return true
       * @return false
       */
      bool springInvolvesCrossLink(const igraph_integer_t edgeId) const
      {
        assert(edgeId < igraph_ecount(&this->graph));
        igraph_integer_t from, to;
        igraph_edge(&this->graph, edgeId, &from, &to);

        return this->linkIsCrossLink(from) || this->linkIsCrossLink(to);
      }

      /**
       * @brief Check whether a given edge is between two cross-links
       *
       * @param edgeId
       * @return true
       * @return false
       */
      bool springIsBetweenCrossLinks(const igraph_integer_t edgeId) const
      {
        assert(edgeId < igraph_ecount(&this->graph));
        igraph_integer_t from, to;
        igraph_edge(&this->graph, edgeId, &from, &to);

        return this->linkIsCrossLink(from) && this->linkIsCrossLink(to);
      }

      /**
       * @brief Shorthand to query the degree of a vertex
       *
       * @param vertexId
       * @param loops
       * @return igraph_integer_t
       */
      igraph_integer_t getVertexDegree(const igraph_integer_t vertexId,
                                       const bool loops = true) const
      {
        assert(vertexId >= 0);
        assert(vertexId < igraph_vcount(&this->graph));
        igraph_integer_t degree;
        igraph_degree_1(&this->graph, &degree, vertexId, IGRAPH_ALL, loops);
        return degree;
      }

      /**
       * @brief Set the Vertex Properties From Atom object (for a cross-link) in
       * the graph
       *
       * @param vertexId
       * @param atom
       */
      void setVertexPropertiesFromAtom(
        const igraph_integer_t vertexId,
        const pylimer_tools::entities::Atom& atom,
        const int atomType)
      {
        assert(vertexId >= 0);
        assert(vertexId < igraph_vcount(&this->graph));
        assert(atomType == this->crosslinkerType ||
               atomType == this->slipLinkType);
        // assert(atom.getType() == this->crosslinkerType);
        igraph_cattribute_VAN_set(
          &this->graph, "atom_id", vertexId, atom.getId());
        igraph_cattribute_VAN_set(&this->graph, "type", vertexId, atomType);
        igraph_cattribute_VAN_set(&this->graph, "num_link_swaps", vertexId, 0);
        Eigen::Vector3d coords = atom.getCoordinates();
        this->setVertexCoordinates(vertexId, coords);
      }

      /**
       * @brief Set the Vertex Properties of a slip-link based on the two
       * involved atoms
       *
       * @param vertexId
       * @param atom1
       * @param atom2
       */
      void setVertexPropertiesFromAtoms(const igraph_integer_t vertexId,
                                        pylimer_tools::entities::Atom& atom1,
                                        pylimer_tools::entities::Atom& atom2)
      {
        assert(vertexId >= 0);
        assert(vertexId < igraph_vcount(&this->graph));
        assert(atom1.getType() != this->crosslinkerType &&
               atom2.getType() != this->crosslinkerType);
        if (atom1.getId() > atom2.getId()) {
          // make sure a second call to this function would result in same
          // result
          std::swap(atom1, atom2);
        }
        igraph_cattribute_VAN_set(
          &this->graph, "type", vertexId, this->slipLinkType);
        igraph_cattribute_VAN_set(
          &this->graph, "atom_id", vertexId, atom1.getId());
        igraph_cattribute_VAN_set(
          &this->graph, "atom_1_id", vertexId, atom1.getId());
        igraph_cattribute_VAN_set(
          &this->graph, "atom_2_id", vertexId, atom2.getId());
        igraph_cattribute_VAN_set(&this->graph, "num_link_swaps", vertexId, 0);
        Eigen::Vector3d coords = atom1.getCoordinates();
        Eigen::Vector3d dist = atom2.getCoordinates() - coords;
        this->universe.getBox().handlePBC(dist);
        coords += 0.5 * dist;
        this->setVertexCoordinates(vertexId, coords);
      }

      /**
       * @brief Create a new edge in the graph with the specified properties
       *
       * @param from
       * @param to
       * @param parent
       * @param partitionFraction
       * @param boxOffset
       * @return igraph_integer_t
       */
      igraph_integer_t createEdge(igraph_integer_t from,
                                  igraph_integer_t to,
                                  long int parent,
                                  double partitionFraction,
                                  double contourLength,
                                  Eigen::Vector3d boxOffset)
      {
        // igraph, in undirected graphs, orders the edges from small to large
        if (from > to) {
          std::swap(from, to);
          boxOffset *= -1.;
        }
        assert(from < igraph_vcount(&this->graph));
        assert(to < igraph_vcount(&this->graph));
        igraph_integer_t newEdgeId = igraph_ecount(&this->graph);
        igraph_add_edge(&this->graph, from, to);
        igraph_integer_t newFrom, newTo;
        igraph_edge(&this->graph, newEdgeId, &newFrom, &newTo);
        assert(newFrom == from && newTo == to);
        igraph_cattribute_EAN_set(
          &this->graph, "partition_fraction", newEdgeId, partitionFraction);
        this->setBondBoxOffsetForEdge(newEdgeId, boxOffset);
        igraph_cattribute_EAN_set(
          &this->graph, "contour_length", newEdgeId, contourLength);
        igraph_cattribute_EAN_set(
          &this->graph, "parent_edge", newEdgeId, parent);
        return newEdgeId;
      }

      /**
       * @brief Set the Bond Box Offset For an edge
       *
       * @param edgeId
       * @param bondBoxOffset
       */
      void setBondBoxOffsetForEdge(const igraph_integer_t edgeId,
                                   const Eigen::Vector3d& bondBoxOffset)
      {
        assert(edgeId < igraph_ecount(&this->graph));
        assert(bondBoxOffset.array().isFinite().all());
        igraph_cattribute_EAN_set(
          &this->graph, "bond_box_x", edgeId, bondBoxOffset[0]);
        igraph_cattribute_EAN_set(
          &this->graph, "bond_box_y", edgeId, bondBoxOffset[1]);
        igraph_cattribute_EAN_set(
          &this->graph, "bond_box_z", edgeId, bondBoxOffset[2]);
      }

      /**
       * @brief Copy one specific bond property from one edge to another
       *
       * @param name
       * @param from
       * @param to
       */
      void copyBondProperty(const char* name,
                            const igraph_integer_t from,
                            const igraph_integer_t to)
      {
        igraph_cattribute_EAN_set(
          &this->graph,
          name,
          to,
          igraph_cattribute_EAN(&this->graph, name, from));
      }

      /**
       * @brief Copy all relevant properties from one edge to another
       *
       * @param from source edge id
       * @param to target edge id
       */
      void copyBondProperties(const igraph_integer_t from,
                              const igraph_integer_t to)
      {
        for (std::string property : { "bond_box_x",
                                      "bond_box_y",
                                      "bond_box_z",
                                      "parent_edge",
                                      "partition_fraction",
                                      "contour_length" }) {
          this->copyBondProperty(property.c_str(), from, to);
        }
      }

      /**
       * @brief Copy some relevant properties from one edge to another, but
       * scaled
       *
       * @param from
       * @param to
       * @param scaleFactor
       */
      void copyScaleBondProperties(const igraph_integer_t from,
                                   const igraph_integer_t to,
                                   double scaleFactor = 0.5)
      {
        for (std::string property :
             { "partition_fraction",
               "local_contour_length" }) { //, "contour_length"
          igraph_cattribute_EAN_set(
            &this->graph,
            property.c_str(),
            to,
            igraph_cattribute_EAN(&this->graph, property.c_str(), from) *
              scaleFactor);
        }
      }

      /**
       * @brief Use an existing chain to set the relevant edge properties
       *
       * @param chain
       * @param atom1Idx
       * @param atom2Idx
       * @param edgeId
       * @param chainIdx
       */
      void setBondPropertiesBasedOnChain(
        const pylimer_tools::entities::Molecule& chain,
        const size_t atom1Idx,
        const size_t atom2Idx,
        const igraph_integer_t edgeId,
        const size_t chainIdx)
      {
        assert(atom1Idx < atom2Idx);
        if (chain.getType() ==
            pylimer_tools::entities::MoleculeType::PRIMARY_LOOP) {
          assert(APPROX_WITHIN(atom1Idx, 0, chain.getLength() - 1, 1e-2));
          assert(APPROX_WITHIN(atom2Idx, 1, chain.getLength(), 1e-2));
        } else {
          assert(APPROX_WITHIN(atom1Idx, 0, chain.getLength() - 2, 1e-2));
          assert(APPROX_WITHIN(atom2Idx, 1, chain.getLength() - 1, 1e-2));
        }
        igraph_integer_t from, to;
        igraph_edge(&this->graph, edgeId, &from, &to);
        std::vector<pylimer_tools::entities::Atom> linedUpAtoms =
          chain.getAtomsLinedUp(crosslinkerType, false, true);
        igraph_cattribute_EAN_set(&this->graph,
                                  "partition_fraction",
                                  edgeId,
                                  static_cast<double>(atom2Idx - atom1Idx) /
                                    static_cast<double>(chain.getNrOfBonds()));
        igraph_cattribute_EAN_set(
          &this->graph, "parent_edge", edgeId, chainIdx);
        igraph_cattribute_EAN_set(
          &this->graph, "contour_length", edgeId, chain.getNrOfBonds());
        igraph_cattribute_EAN_set(
          &this->graph,
          "local_contour_length",
          edgeId,
          chain.getNrOfBondsFromTo(linedUpAtoms[atom1Idx].getId(),
                                   linedUpAtoms[atom2Idx].getId(),
                                   crosslinkerType));
        // use the actual position of the vertices!
        Eigen::Vector3d expectedDistance =
          chain.getOverallBondSumFromTo(linedUpAtoms[atom1Idx].getId(),
                                        linedUpAtoms[atom2Idx].getId(),
                                        crosslinkerType);
        Eigen::Vector3d additionalDistance1 =
          linedUpAtoms[atom1Idx].getCoordinates() -
          this->getVertexCoordinates(from);
        this->universe.getBox().handlePBC(additionalDistance1);
        Eigen::Vector3d additionalDistance2 =
          this->getVertexCoordinates(to) -
          linedUpAtoms[atom2Idx].getCoordinates();
        this->universe.getBox().handlePBC(additionalDistance2);
        expectedDistance += additionalDistance1 + additionalDistance2;
        if (this->is2D) {
          expectedDistance[2] = 0.0;
        }

        this->setBondBoxOffsetForEdge(edgeId, Eigen::Vector3d::Zero());
        Eigen::Vector3d actualDistance = this->computeEdgeDistance(edgeId);
        this->setBondBoxOffsetForEdge(edgeId,
                                      expectedDistance - actualDistance);
        assert(this->universe.getBox().isValidOffset(expectedDistance -
                                                     actualDistance));
#ifndef NDEBUG
        Eigen::Vector3d newActualDistance = this->computeEdgeDistance(edgeId);
        assert(newActualDistance.isApprox(expectedDistance));
#endif
      }

      /**
       * @brief Returns the box offset for a given edge
       *
       * CAUTION: make sure the graph is up to date!
       *
       * @param edgeId
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d getBondBoxOffsetForEdgeFrom(igraph_integer_t edgeId,
                                                  igraph_integer_t vertexIdx)
      {
        Eigen::Vector3d boxOffset = this->getBondBoxOffsetForEdge(edgeId);

        igraph_integer_t from, to;
        igraph_edge(&this->graph, edgeId, &from, &to);
        assert(from == vertexIdx || to == vertexIdx);

        if (from == vertexIdx) {
          return boxOffset;
        } else {
          return -1. * boxOffset;
        }
      }

      Eigen::Vector3d getBondBoxOffsetForEdgeTo(igraph_integer_t edgeId,
                                                igraph_integer_t vertexIdx)
      {
        return -1. * this->getBondBoxOffsetForEdgeFrom(edgeId, vertexIdx);
      }

      /**
       * @brief Returns the box offset for a given edge
       *
       * CAUTION: make sure the graph is up to date!
       *
       * @param edgeId
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d getBondBoxOffsetForEdge(igraph_integer_t edgeId) const
      {
        Eigen::Vector3d bondBoxOffset;
        bondBoxOffset << igraph_cattribute_EAN(
          &this->graph, "bond_box_x", edgeId),
          igraph_cattribute_EAN(&this->graph, "bond_box_y", edgeId),
          igraph_cattribute_EAN(&this->graph, "bond_box_z", edgeId);

        assert(bondBoxOffset.allFinite());
        return bondBoxOffset;
      }
      /**
       * @brief Compute the length of an edge based on the graph
       *
       * CAUTION: make sure the graph is up to date!
       *
       * @param edgeId
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d computeEdgeDistance(igraph_integer_t edgeId) const
      {
        assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
        igraph_integer_t from, to;
        igraph_edge(&this->graph, edgeId, &from, &to);

        Eigen::Vector3d dist = this->getVertexCoordinates(to) -
                               this->getVertexCoordinates(from) +
                               this->getBondBoxOffsetForEdge(edgeId);
        if (this->assumeBoxLargeEnough) {
          this->universe.getBox().handlePBC(dist);
        }
        if (this->is2D) {
          dist[2] = 0.;
        }
        assert(dist.allFinite());
        return dist;
      }

      /**
       * @brief Shortcut to query the edge's partition fraction
       *
       * @param edgeId
       * @return double
       */
      double getEdgeFraction(igraph_integer_t edgeId) const
      {
        assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
        return igraph_cattribute_EAN(
          &this->graph, "partition_fraction", edgeId);
      }

      /**
       * @brief Shortcut to compute 1/(N*fraction)
       *
       * @param edgeId
       * @param oneOverSpringPartitionUpperLimit
       * @param isPartialSpring
       * @return double
       */
      double getEdgeDenominator(igraph_integer_t edgeId,
                                const double oneOverSpringPartitionUpperLimit,
                                bool isPartialSpring) const
      {
        assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
        const double N =
          igraph_cattribute_EAN(&this->graph, "contour_length", edgeId);
        const double fraction =
          igraph_cattribute_EAN(&this->graph, "partition_fraction", edgeId);
        double denominator = 1. / (fraction * N);
        if (oneOverSpringPartitionUpperLimit > 0. ||
            !std::isfinite(denominator)) {
          denominator = CLAMP_ONE_OVER_SPRINGPARTITION(
            isPartialSpring, denominator, N, oneOverSpringPartitionUpperLimit);
        }

        assert(std::isfinite(denominator));
        return denominator;
      }

      double getEdgeDenominator(
        igraph_integer_t edgeId,
        const double oneOverSpringPartitionUpperLimit) const
      {
        return this->getEdgeDenominator(
          edgeId,
          oneOverSpringPartitionUpperLimit,
          !this->springIsBetweenCrossLinks(edgeId));
      }

      /**
       * @brief Compute the length of an edge based on the graph in a certain
       * direction
       *
       * CAUTION: make sure the graph is up to date!
       *
       * @param edgeId
       * @param vertexId
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d computeEdgeDistanceFrom(igraph_integer_t edgeId,
                                              igraph_integer_t vertexId) const
      {
        igraph_integer_t from, to;
        igraph_edge(&this->graph, edgeId, &from, &to);
        assert(from == vertexId || to == vertexId);

        Eigen::Vector3d dist = this->computeEdgeDistance(edgeId);
        return dist * (vertexId == from ? 1. : -1.);
      }

      /**
       * @brief Compute the length of an edge based on the graph in a certain
       * direction
       *
       * CAUTION: make sure the graph is up to date!
       *
       * @param edgeId
       * @param vertexId the target vertex id
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d computeEdgeDistanceTo(igraph_integer_t edgeId,
                                            igraph_integer_t vertexId) const
      {
        return -1. * this->computeEdgeDistanceFrom(edgeId, vertexId);
      }

      /**
       * @brief Compute the length of an edge based on the graph
       *
       * CAUTION: make sure the graph is up to date!
       * Also, this function is O(|E|), since all edges are iterated to check
       * their parent edge id.
       *
       * @param parentEdgeId
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d computeParentEdgeLength(size_t parentEdgeId) const
      {
        Eigen::Vector3d dist = Eigen::Vector3d::Zero();

        igraph_vector_t parentEdges;
        igraph_vector_init(&parentEdges, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "parent_edge",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &parentEdges);

        std::vector<size_t> edgeIdsToRemove;
        edgeIdsToRemove.reserve(4);
        for (igraph_integer_t i = 0; i < igraph_vector_size(&parentEdges);
             ++i) {
          if (igraph_vector_get(&parentEdges, i) == parentEdgeId) {
            dist += this->computeEdgeDistance(i);
          }
        }

        igraph_vector_destroy(&parentEdges);

        return dist;
      }

      /**
       * @brief Remove a spring (and all its parts, incl. slip-links) from the
       * structures
       *
       * NOTE: this may result in slip-links with f = 2.
       * They are not automatically removed in order to preserve vertex
       * iterations.
       *
       * @param net
       * @param springPartitions
       */
      void removeParentSpring(const size_t springIdx);

      /**
       * @brief Remove a set of edges from the graph
       *
       * @param edgeIdsToRemove
       */
      void removePartialSprings(std::vector<igraph_integer_t>& edgeIdsToRemove);

      /**
       * @brief marks a certain "parent" spring as non-existing
       */
      void combineParentSprings(size_t springIdxBefore,
                                size_t springIdxNow,
                                double contourLengthBefore = -1.);

      /**
       * @brief When springs have been removed, it is possible that the
       * numbering is not sequential anymore. This function fixes that.
       *
       */
      void renumberParentSprings();

      /**
       * @brief When spring has been removed, it is possible that the
       * sum of the partitions don't add up to one anymore. This function fixes
       * that.
       *
       */
      void renormalizePartitions();

      /**
       * @brief Combine two partial springs to be only one
       *
       * @param edge1Id
       * @param edge2Id
       */
      void combinePartialSprings(const igraph_integer_t edge1Id,
                                 const igraph_integer_t edge2Id,
                                 const igraph_integer_t centralLink);

      /**
       * @brief Remove a slip-link and combine the two edges corresponding to
       * the rail
       *
       * @param vertexId
       * @param railEdgeId
       */
      void unlinkSlipLinkFromRail(const igraph_integer_t vertexId,
                                  const igraph_integer_t railEdgeId);

      /**
       * @brief Inserts the given slip-link into a partial spring
       *
       * @param vertexId the slip-link to insert into the spring
       * @param railEdgeId the spring to be halfed
       */
      void insertSlipLinkIntoRail(const igraph_integer_t vertexId,
                                  const igraph_integer_t railEdgeId);

      /**
       * @brief Move a slip-link from one rail to another
       *
       * @param vertexId the slip-link to remove from one and insert into
       * another spring
       * @param sourceRailEdgeId the spring to be twiced
       * @param targetRailEdgeId the spring to be halfed
       */
      void moveSlipLinkFromRailToRail(const igraph_integer_t vertexId,
                                      const igraph_integer_t sourceRailEdgeId,
                                      const igraph_integer_t targetRailEdgeId);

      /**
       * @brief List the edges and vertices of one spring, in order
       *
       * @param springIdx
       * @param vertices
       * @param edges
       */
      void findEdgesAndVerticesOfSpring(size_t springIdx,
                                        igraph_vector_int_t* vertices,
                                        igraph_vector_int_t* edges);
      /**
       * @brief List the edges and vertices of one spring, in order
       *
       * @param unorderedEdges
       * @param vertices
       * @param edges
       */
      void findEdgesAndVerticesOfSpring(igraph_vector_int_t* unorderedEdges,
                                        igraph_vector_int_t* vertices,
                                        igraph_vector_int_t* edges);

      /**
       * @brief Find the index of a partial spring, given the fractn of the
       * total spring to traverse
       *
       * @param springIdx
       * @param alpha
       * @return igraph_integer_t
       */
      igraph_integer_t findPartialSpringByFraction(size_t springIdx,
                                                   double alpha,
                                                   double& fractionTillThere);

      /**
       * @brief Remove a certain, 2-functional link from the structures,
       * combining the two strands
       */
      void remove2fLink(const size_t linkIdx);

      /**
       * @brief Remove a certain, 3-functional link from the structures,
       * combining the two strands
       */
      void remove3fLink(const size_t linkIdx);

      /**
       * @brief Remove a certain link from the structures, removing all
       * connections
       */
      void removeLink(const size_t linkIdx);

      igraph_integer_t getOtherEdgePartner(
        igraph_integer_t edgeId,
        igraph_integer_t wrongPartnerVertexIdx)
      {
        igraph_integer_t from, to;
        igraph_edge(&this->graph, edgeId, &from, &to);
        return (from == wrongPartnerVertexIdx) ? to : from;
      }

      /**
       * @brief Given a vertex id and a rail edge, returns the other two edges
       * that are not part of the rail
       */
      std::vector<igraph_integer_t> getOffRailConnectedEdgeIds(
        igraph_integer_t vertexId,
        igraph_integer_t railEdgeId,
        igraph_integer_t avoidEdgeId = -1);

      /**
       * @brief Given a vertex and a connected edge, returns the edge in the
       * opposite direction
       *
       */
      igraph_integer_t getOtherRailEdgeId(igraph_integer_t vertexId,
                                          igraph_integer_t railEdgeId,
                                          igraph_integer_t avoidEdgeId = -1);

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

      /**
       * @brief Concretely add new slip-links
       *
       * @param strandIdx1
       * @param strandIdx2
       * @param x
       * @param y
       * @param z
       * @param alpha1
       * @param alpha2
       * @param clampAlpha
       */
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
       * @brief Compute a few properties of the simulator
       *
       */
      void finaliseInitialisation()
      {
        igraph_cattribute_GAB_set(&this->graph, "is_up_to_date", true);
        this->convertFromGraph();
        this->defaultR0Squared =
          this->universe.computeMeanSquareEndToEndDistance(
            this->crosslinkerType);
        assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
        this->defaultNrOfChains =
          this->universe.getMolecules(this->crosslinkerType).size();
        assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
        this->validateNetwork();
        assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
      }

      /**
       * @brief Convert the internal graph representation to the internal
       * ForceBalanceNetwork representation
       *
       */
      void convertFromGraph()
      {
        assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));

        this->net.isUpToDate = false;
        this->net.nrOfPartialSprings = igraph_ecount(&this->graph);
        this->net.nrOfLinks = igraph_vcount(&this->graph);

        igraph_vector_t parentEdges;
        igraph_vector_init(&parentEdges, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "parent_edge",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &parentEdges);

        std::unordered_set<igraph_integer_t> parentEdgeIds;
        igraph_integer_t maxParentEdgeId = 0;
        for (size_t i = 0; i < igraph_vector_size(&parentEdges); ++i) {
          igraph_integer_t parentEdgeId =
            castToIgraphInt(igraph_vector_get(&parentEdges, i));
          parentEdgeIds.insert(parentEdgeId);
          maxParentEdgeId = std::max(parentEdgeId, maxParentEdgeId);
        }
        this->net.nrOfSprings = parentEdgeIds.size();
        assert(maxParentEdgeId == parentEdgeIds.size() - 1 ||
               parentEdgeIds.size() == 0);
        size_t numPartialSprings =
          this->net.nrOfPartialSprings - this->net.nrOfSprings;

        // reset & resize
        this->net.coordinates.resize(3 * this->net.nrOfLinks);
        this->net.springsContourLength.resize(this->net.nrOfSprings);
        this->net.springIndicesOfLinks.clear();
        this->net.linkIndicesOfSprings.clear();
        this->net.partialSpringIsPartial.resize(this->net.nrOfPartialSprings);
        this->net.localToGlobalSpringIndex.clear();
        this->net.partialToFullSpringIndex.resize(this->net.nrOfPartialSprings);
        this->net.linkIsSliplink.resize(this->net.nrOfLinks);
        this->net.springPartCoordinateIndexA.resize(
          3 * this->net.nrOfPartialSprings);
        this->net.springPartCoordinateIndexB.resize(
          3 * this->net.nrOfPartialSprings);
        this->net.springPartIndexA.resize(this->net.nrOfPartialSprings);
        this->net.springPartIndexB.resize(this->net.nrOfPartialSprings);
        this->net.springPartBoxOffset.resize(3 * this->net.nrOfPartialSprings);
        this->currentSpringPartitionsVec.resize(this->net.nrOfPartialSprings);

        this->net.springIsActive.resize(this->net.nrOfSprings);
        this->net.springIndexA.resize(this->net.nrOfSprings);
        this->net.springIndexB.resize(this->net.nrOfSprings);
        this->net.springCoordinateIndexA.resize(3 * this->net.nrOfSprings);
        this->net.springCoordinateIndexB.resize(3 * this->net.nrOfSprings);

        // reset everything we cleared
        for (size_t i = 0; i < this->net.nrOfSprings; ++i) {
          std::vector<size_t> vec;
          this->net.linkIndicesOfSprings.push_back(vec);
          this->net.localToGlobalSpringIndex.push_back(vec);
        }
        for (size_t i = 0; i < this->net.nrOfLinks; ++i) {
          std::vector<size_t> vec;
          this->net.springIndicesOfLinks.push_back(vec);
        }

        // fetch other properties needed
        // springs / partial springs
        igraph_vector_int_t allEdges;
        igraph_vector_int_init(&allEdges, igraph_ecount(&this->graph) * 2);
        if (igraph_edges(
              &this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &allEdges)) {
          throw std::runtime_error("Failed to get all edges");
        }

        igraph_vector_t partitionFraction;
        igraph_vector_init(&partitionFraction, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "partition_fraction",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &partitionFraction);
        igraph_vector_t contourLength;
        igraph_vector_init(&contourLength, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "contour_length",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &contourLength);
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

        this->net.springsContourLength.setZero();
        for (size_t i = 0; i < this->net.nrOfPartialSprings; ++i) {
          this->net.partialToFullSpringIndex(i) =
            castToIgraphInt(igraph_vector_get(&parentEdges, i));
          assert(std::isfinite(igraph_vector_get(&contourLength, i)));
          this->net.springsContourLength(this->net.partialToFullSpringIndex(
            i)) = igraph_vector_get(&contourLength, i);
          this->net
            .localToGlobalSpringIndex[castToIgraphInt(
              igraph_vector_get(&parentEdges, i))]
            .push_back(i);
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

          this->net.springPartBoxOffset(i * 3 + 0) =
            igraph_vector_get(&bondBoxOffsetX, i);
          this->net.springPartBoxOffset(i * 3 + 1) =
            igraph_vector_get(&bondBoxOffsetY, i);
          this->net.springPartBoxOffset(i * 3 + 2) =
            igraph_vector_get(&bondBoxOffsetZ, i);
          this->currentSpringPartitionsVec(i) =
            igraph_vector_get(&partitionFraction, i);
        }
        if (this->net.nrOfSprings > 0) {
          this->net.meanSpringContourLength =
            this->net.springsContourLength.mean();
        } else {
          this->net.meanSpringContourLength = 0.;
        }

        // then, for the "parent" springs, we need to know the order of the
        // partial springs
        // NOTE: this order may change. But, since we convert back an forth,
        // this is not a problem
        for (size_t i = 0; i < this->net.nrOfSprings; ++i) {
          igraph_vector_int_t edgesOnPath;
          igraph_vector_int_init(&edgesOnPath, 0);
          igraph_vector_int_t verticesOnPath;
          igraph_vector_int_init(&verticesOnPath, 0);

          this->findEdgesAndVerticesOfSpring(i, &verticesOnPath, &edgesOnPath);
          assert(igraph_vector_int_size(&edgesOnPath) > 0);
          assert(igraph_vector_int_size(&verticesOnPath) ==
                 igraph_vector_int_size(&edgesOnPath) + 1);

          std::vector<size_t> linkIndicesOfThisSpring;
          pylimer_tools::utils::igraphVectorTToStdVector(
            &verticesOnPath, linkIndicesOfThisSpring);
          assert(linkIndicesOfThisSpring.size() ==
                 igraph_vector_int_size(&verticesOnPath));
          assert(linkIndicesOfThisSpring.size() > 1);
          this->net.linkIndicesOfSprings[i] = linkIndicesOfThisSpring;
          for (size_t linkIdx : linkIndicesOfThisSpring) {
            this->net.springIndicesOfLinks[linkIdx].push_back(i);
          }
          std::vector<size_t> edgeIdsOfThisSpring;
          pylimer_tools::utils::igraphVectorTToStdVector(&edgesOnPath,
                                                         edgeIdsOfThisSpring);
          this->net.localToGlobalSpringIndex[i] = edgeIdsOfThisSpring;
          this->net.springIndexA[i] = this->net.linkIndicesOfSprings[i][0];
          this->net.springIndexB[i] =
            pylimer_tools::utils::last(this->net.linkIndicesOfSprings[i]);
          for (size_t dir = 0; dir < 3; ++dir) {
            this->net.springCoordinateIndexA[3 * i + dir] =
              3 * this->net.springIndexA[i] + dir;
            this->net.springCoordinateIndexB[3 * i + dir] =
              3 * this->net.springIndexB[i] + dir;
          }

          igraph_vector_int_destroy(&edgesOnPath);
          igraph_vector_int_destroy(&verticesOnPath);
        }

        for (size_t i = 0; i < this->net.nrOfPartialSprings; ++i) {
          this->net.partialSpringIsPartial(i) =
            this->net
              .linkIndicesOfSprings[this->net.partialToFullSpringIndex(i)]
              .size() > 2;
        }

        // cleanup
        igraph_vector_int_destroy(&allEdges);
        igraph_vector_destroy(&parentEdges);
        igraph_vector_destroy(&contourLength);
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

        igraph_vector_t numLinkSwaps;
        igraph_vector_init(&numLinkSwaps, this->net.nrOfLinks);
        igraph_cattribute_VANV(
          &this->graph, "num_link_swaps", igraph_vss_all(), &numLinkSwaps);

        igraph_vector_t linkType;
        igraph_vector_init(&linkType, this->net.nrOfLinks);
        igraph_cattribute_VANV(
          &this->graph, "type", igraph_vss_all(), &linkType);

        igraph_vector_t linkAtomId;
        igraph_vector_init(&linkAtomId, this->net.nrOfLinks);
        igraph_cattribute_VANV(
          &this->graph, "atom_id", igraph_vss_all(), &linkAtomId);

        // actually write things
        size_t numCrosslinks = 0;
        for (size_t i = 0; i < this->net.nrOfLinks; ++i) {
          this->net.coordinates(3 * i + 0) = igraph_vector_get(&coordsX, i);
          this->net.coordinates(3 * i + 1) = igraph_vector_get(&coordsY, i);
          this->net.coordinates(3 * i + 2) = igraph_vector_get(&coordsZ, i);
          int alinkType = castToIgraphInt(igraph_vector_get(&linkType, i));
          assert(alinkType == this->slipLinkType ||
                 alinkType == this->crosslinkerType);
          this->net.linkIsSliplink(i) = alinkType == this->slipLinkType;
          if (alinkType == this->crosslinkerType) {
            numCrosslinks += 1;
          }
        }
        this->net.nrOfNodes = numCrosslinks;

        this->net.nrOfCrosslinkSwapsEndured.resize(this->net.nrOfLinks -
                                                   this->net.nrOfNodes);
        this->net.oldAtomIds.resize(numCrosslinks);
        size_t slipLinkIdx = 0;
        size_t crossLinkIdx = 0;
        for (size_t i = 0; i < this->net.nrOfLinks; ++i) {
          if (castToIgraphInt(igraph_vector_get(&linkType, i)) !=
              this->crosslinkerType) {
            this->net.nrOfCrosslinkSwapsEndured(slipLinkIdx) =
              castToIgraphInt(igraph_vector_get(&numLinkSwaps, i));
            slipLinkIdx += 1;
          } else {
            this->net.oldAtomIds[crossLinkIdx] =
              castToIgraphInt(igraph_vector_get(&linkAtomId, i));
            crossLinkIdx += 1;
          }
        }
        assert(slipLinkIdx == this->net.nrOfLinks - this->net.nrOfNodes);
        assert(crossLinkIdx == this->net.nrOfNodes);

        // cleanup
        igraph_vector_destroy(&linkType);
        igraph_vector_destroy(&coordsX);
        igraph_vector_destroy(&coordsY);
        igraph_vector_destroy(&coordsZ);

        // mark as done
        this->net.isUpToDate = true;

        // compute other local properties
        this->currentSpringDistances = this->evaluateSpringDistances();
        assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
        this->currentPartialSpringDistances =
          this->evaluatePartialSpringDistances();
        assert(igraph_cattribute_GAB(&this->graph, "is_up_to_date"));
      };

      /**
       * @brief Write the new coordinates to the graph
       */
      void updateGraph()
      {
        RUNTIME_EXP_IFN(this->net.isUpToDate,
                        "Should not update graph with outdated network");
        RUNTIME_EXP_IFN(igraph_vcount(&this->graph) == this->net.nrOfLinks,
                        "Should not update graph with outdated network");
        RUNTIME_EXP_IFN(igraph_ecount(&this->graph) ==
                          this->net.nrOfPartialSprings,
                        "Should not update graph with outdated network");

        // write the current coordinates to the graph
        igraph_vector_t coordsX;
        igraph_vector_init(&coordsX, this->net.nrOfLinks);

        igraph_vector_t coordsY;
        igraph_vector_init(&coordsY, this->net.nrOfLinks);

        igraph_vector_t coordsZ;
        igraph_vector_init(&coordsZ, this->net.nrOfLinks);

        for (igraph_integer_t i = 0; i < this->net.nrOfLinks; ++i) {
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

        for (igraph_integer_t i = 0; i < this->net.nrOfPartialSprings; ++i) {
          igraph_vector_set(
            &partitionFraction, i, this->currentSpringPartitionsVec(i));
          this->setBondBoxOffsetForEdge(
            i, this->net.springPartBoxOffset.segment(3 * i, 3));
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
      double evaluatePressure() const
      {
        auto stressTensor = this->evaluateStressTensor();
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

      Eigen::Matrix3d evaluateStressTensorForLinks(
        const std::vector<size_t> linkIndices,
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
        const double kappa0 = 1.0,
        const double oneOverSpringPartitionUpperLimit = 1.0) const;

      /**
       * @brief Compute the force acting on a slip-link
       *
       * @param linkIdx
       * @param net
       * @param u
       * @param springPartitions
       * @param kappa0
       * @param minCutoff
       * @return Eigen::Vector3d
       */
      Eigen::Matrix3d evaluateForceOnLink(
        const size_t linkIdx,
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
                                 const double tolerance = 0.01) const
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
       * @return Eigen::ArrayXb
       */
      Eigen::ArrayXb findActiveSprings(const Eigen::VectorXd& springDistances,
                                       const double tolerance = 0.01) const
      {
        Eigen::ArrayXb result =
          Eigen::ArrayXb::Constant(springDistances.size() / 3, false);
        for (size_t i = 0; i < springDistances.size() / 3; ++i) {
          result[i] =
            springDistances.segment(3 * i, 3).squaredNorm() > tolerance;
        }
        return result;
      }

      /**
       * @brief Swap two slip-links along a partial spring, iff the move leads
       * to a smaller stress diagonal squared norm
       *
       * @param partialSpringIdx
       * @param oneOverSpringPartitionUpperLimit
       * @return true
       * @return false
       */
      bool swapSlipLinksReversibly(
        const size_t partialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0)
      {
#ifndef NDEBUG
        this->validateIgraphSpring(
          igraph_cattribute_EAN(&this->graph, "parent_edge", partialSpringIdx));
#endif

        igraph_integer_t from, to;
        igraph_edge(&this->graph, partialSpringIdx, &from, &to);
        assert(castToIgraphInt(igraph_cattribute_VAN(
                 &this->graph, "type", from)) == this->slipLinkType);
        assert(castToIgraphInt(igraph_cattribute_VAN(
                 &this->graph, "type", to)) == this->slipLinkType);

        igraph_integer_t otherRailFrom =
          this->getOtherRailEdgeId(from, partialSpringIdx);
        igraph_integer_t otherRailTo =
          this->getOtherRailEdgeId(to, partialSpringIdx);

        double thisSpringDenominator = this->getEdgeDenominator(
          partialSpringIdx, oneOverSpringPartitionUpperLimit);
        Eigen::Vector3d thisSpringDistance =
          this->computeEdgeDistance(partialSpringIdx);
        double otherRailFromDenominator = this->getEdgeDenominator(
          otherRailFrom, oneOverSpringPartitionUpperLimit);
        Eigen::Vector3d otherRailFromSpringDistance =
          this->computeEdgeDistanceTo(otherRailFrom, from);
        double otherRailToDenominator = this->getEdgeDenominator(
          otherRailTo, oneOverSpringPartitionUpperLimit);
        Eigen::Vector3d otherRailToSpringDistance =
          this->computeEdgeDistanceFrom(otherRailTo, to);

        double forceEstimateBefore =
          (-otherRailFromSpringDistance * otherRailFromDenominator +
           thisSpringDistance * thisSpringDenominator)
            .squaredNorm() +
          (-thisSpringDistance * thisSpringDenominator +
           otherRailToSpringDistance * otherRailToDenominator)
            .squaredNorm();
        double forceEstimateAfter =
          (thisSpringDistance * thisSpringDenominator +
           (otherRailToSpringDistance + thisSpringDistance) *
             otherRailToDenominator)
            .squaredNorm() +
          (-thisSpringDistance * thisSpringDenominator -
           (thisSpringDistance + otherRailFromSpringDistance) *
             otherRailFromDenominator)
            .squaredNorm();

        if (forceEstimateAfter <= forceEstimateBefore) {
          this->swapSlipLinksOfEdge(partialSpringIdx);
          this->net.isUpToDate = false;
          return true;
        }
        return false;
      };

      /**
       * @brief Swap a slip-link and a cross-link along a partial spring, iff
       * the move leads to a smaller stress diagonal squared norm
       *
       * @param partialSpringIdx
       * @param oneOverSpringPartitionUpperLimit
       * @return true
       * @return false
       */
      bool swapSlipLinkWithXlinkReversibly(
        const size_t partialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const bool respectLoops = true)
      {
        igraph_integer_t from, to;
        igraph_edge(&this->graph, partialSpringIdx, &from, &to);

        igraph_integer_t xlinkIdx =
          castToIgraphInt(igraph_cattribute_VAN(&this->graph, "type", from)) ==
              this->crosslinkerType
            ? from
            : to;
        igraph_integer_t slipLinkIdx = xlinkIdx == from ? to : from;
        assert(castToIgraphInt(igraph_cattribute_VAN(
                 &this->graph, "type", slipLinkIdx)) == this->slipLinkType);

        igraph_vector_int_t edgesOfCrossLink;
        igraph_vector_int_init(&edgesOfCrossLink, 4);
        igraph_incident(&this->graph, &edgesOfCrossLink, xlinkIdx, IGRAPH_ALL);
        RUNTIME_EXP_IFN(igraph_vector_int_size(&edgesOfCrossLink) > 0,
                        "Expected to find more than one edge on cross-link");

        igraph_integer_t otherRailPart =
          this->getOtherRailEdgeId(slipLinkIdx, partialSpringIdx);

        Eigen::Vector3d otherRailDistance =
          this->computeEdgeDistanceFrom(otherRailPart, slipLinkIdx);
        Eigen::Vector3d thisRailDistance =
          this->computeEdgeDistanceTo(partialSpringIdx, xlinkIdx);

        for (size_t i = 0; i < igraph_vector_int_size(&edgesOfCrossLink); ++i) {
          igraph_integer_t attemptedEdge =
            igraph_vector_int_get(&edgesOfCrossLink, i);
          if (attemptedEdge == partialSpringIdx) {
            continue;
          }
          Eigen::Vector3d attemptSpringDistance =
            this->computeEdgeDistanceFrom(attemptedEdge, xlinkIdx);

          double forceEstimateBefore =
            (otherRailDistance + thisRailDistance).squaredNorm() +
            (attemptSpringDistance).squaredNorm();
          double forceEstimateAfter =
            (-thisRailDistance - otherRailDistance).squaredNorm() +
            (attemptSpringDistance + thisRailDistance).squaredNorm();

          if (forceEstimateAfter < forceEstimateBefore) {
            igraph_cattribute_VAN_set(&this->graph,
                                      "num_link_swaps",
                                      slipLinkIdx,
                                      igraph_cattribute_VAN(&this->graph,
                                                            "num_link_swaps",
                                                            slipLinkIdx) +
                                        1);
            this->moveSlipLinkFromRailToRail(
              slipLinkIdx, partialSpringIdx, attemptedEdge);
            this->net.isUpToDate = false;
            return true;
          }
        }

        return false;
      };
    };
  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
