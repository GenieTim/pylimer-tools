#ifndef MEHP_FORCE_BALANCE2_H
#define MEHP_FORCE_BALANCE2_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/NeighbourList.h"
#include "../entities/Universe.h"
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
#include <vector>

namespace pylimer_tools {
namespace calc {
  namespace mehp {
#ifndef CLAMP_ONE_OVER_SPRINGPARTITION
/**
 * @brief a macro for doing the clamping in the routines using kappa,
 * to prevent deivision by zero issues / multiplications by infinity
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

    igraph_integer_t castToIgraphInt(igraph_real_t c)
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

      MEHPForceBalance2(const pylimer_tools::entities::Universe& u,
                        int crosslinkerType = 2,
                        bool is2D = false,
                        double kappa = 1.0)
        : universe(u)
      {
        this->crosslinkerType = crosslinkerType;
        // interpret network already to be able to give early results
        ForceBalanceNetwork network;
        this->net = network;
        this->is2D = is2D;
        igraph_empty(&this->graph, 0, IGRAPH_UNDIRECTED);
      };

    public:
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
        this->splipLinkType = src.splipLinkType;
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
        std::swap(this->splipLinkType, src.splipLinkType);
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

      static MEHPForceBalance2 constructWithoutSlipLinks(
        const pylimer_tools::entities::Universe& universe,
        int crosslinkerType = 2,
        bool is2D = false,
        double kappa = 1.0)
      {
        return MEHPForceBalance2::constructWithRandomSlipLinks(
          universe, 0, 1.0, 0, 1, 0, crosslinkerType, is2D, kappa);
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
        bool clampAlpha = false)
      {
        MEHPForceBalance2 fb =
          MEHPForceBalance2(universe, crosslinkerType, is2D, kappa);
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
        const int seed,
        int crosslinkerType = 2,
        bool is2D = false,
        double kappa = 1.0)
      {
        MEHPForceBalance2 fb =
          MEHPForceBalance2(universe, crosslinkerType, is2D, kappa);

        INVALIDARG_EXP_IFN(minimumNrOfSliplinks <
                             fb.universe.getNrOfAtoms() / 2,
                           "Minimum number of slip-links must be less than the "
                           "possible number of slip-links to place.");
        INVALIDARG_EXP_IFN(nrOfSliplinksToSample <
                             fb.universe.getNrOfAtoms() / 2,
                           "Number of slip-links to place must be less than "
                           "the possible number of slip-links to place.");
        INVALIDARG_EXP_IFN(nrOfSliplinksToSample > minimumNrOfSliplinks,
                           "Maximum nr. should be larger than minimum, got " +
                             std::to_string(nrOfSliplinksToSample) + " and " +
                             std::to_string(minimumNrOfSliplinks) + ".");
        INVALIDARG_EXP_IFN(cutoff > 0.0,
                           "Expected a cutoff > 0.0, got " +
                             std::to_string(cutoff) + ".");

        std::vector<pylimer_tools::entities::Molecule> crosslinkerChains =
          fb.universe.getChainsWithCrosslinker(crosslinkerType);
        std::vector<std::pair<size_t, size_t>> pairsOfAtoms;
        pairsOfAtoms.reserve(nrOfSliplinksToSample);
        std::vector<long int> pairOfAtom =
          pylimer_tools::utils::initializeWithValue<long int>(
            fb.universe.getNrOfAtoms(), -1);

        std::unordered_map<size_t, size_t> atomToStrand;
        atomToStrand.reserve(universe.getNrOfAtoms());
        std::unordered_map<size_t, size_t> atomIdxInStrand;
        atomIdxInStrand.reserve(universe.getNrOfAtoms());
        for (size_t i = 0; i < crosslinkerChains.size(); ++i) {
          pylimer_tools::entities::Molecule chain = crosslinkerChains[i];
          RUNTIME_EXP_IFN(chain.getType() !=
                            pylimer_tools::entities::MoleculeType::UNDEFINED,
                          "Couldn't determine molecule type.");
          std::vector<pylimer_tools::entities::Atom> atoms =
            crosslinkerChains[i].getAtomsLinedUp(crosslinkerType);
          for (size_t atomIdx = 0; atomIdx < atoms.size(); ++atomIdx) {
            pylimer_tools::entities::Atom atom = atoms[atomIdx];
            if (atom.getType() != crosslinkerType) {
              atomToStrand.emplace(atom.getId(), i);
              atomIdxInStrand.emplace(atom.getId(), atomIdx);
            }
          }
        }

        // filter, we don't want cross-links etc.
        std::vector<pylimer_tools::entities::Atom> atomsForNeighbourList =
          fb.universe.getAtomsOfDegree(2);
        std::random_device rd{};
        std::mt19937 rng = std::mt19937(seed > 0 ? seed : rd());
        std::shuffle(
          atomsForNeighbourList.begin(), atomsForNeighbourList.end(), rng);
        pylimer_tools::entities::NeighbourList neighbourList =
          pylimer_tools::entities::NeighbourList(
            atomsForNeighbourList, fb.universe.getBox(), cutoff);
        size_t numLinksFoundInIteration = 0;
        while (pairsOfAtoms.size() < minimumNrOfSliplinks &&
               numLinksFoundInIteration > 0) {
          for (pylimer_tools::entities::Atom a1 : atomsForNeighbourList) {
            size_t atomVertexIdx1 = universe.getIdxByAtomId(a1.getId());
            // find neighbors
            if (pairOfAtom[atomVertexIdx1] == -1) {
              continue;
            }
            // then, find neighbouring atoms (but not from the same strand?!)
            std::vector<pylimer_tools::entities::Atom> neighbours =
              neighbourList.getAtomsCloseTo(a1);
            neighbourList.removeAtom(
              a1, "After querying neighbours. Impossible case.");
            // filter the neighbours to include only those from other strands
            // NOTE: this skews the whole thing a bit
            neighbours.erase(
              std::remove_if(
                neighbours.begin(),
                neighbours.end(),
                [&](const pylimer_tools::entities::Atom& a) -> bool {
                  return (
                    (atomToStrand[a.getId()] ==
                       atomToStrand[a1.getId()] // do not use "at", because not
                                                // all atoms in the neighbours
                                                // have been assigned a strand
                     && (std::abs(static_cast<double>(
                           atomIdxInStrand[a.getId()] -
                           atomIdxInStrand[a1.getId()])) < sameStrandCutoff)));
                }),
              neighbours.end());
            if (neighbours.size() == 0) {
              // std::cerr << "Not enough close neighbours found." << std::endl;
              continue;
            }
            // then, randomly select one of them
            pylimer_tools::entities::Atom a2 = neighbours[0];
            if (neighbours.size() > 1) {
              size_t randomA2Idx = std::uniform_int_distribution<size_t>{
                0, neighbours.size() - 1
              }(rng);
              a2 = neighbours[randomA2Idx];
            }

            size_t atomVertexIdx2 = universe.getIdxByAtomId(a2.getId());
            assert(pairOfAtom[atomVertexIdx2] == -1);
            pairOfAtom[atomVertexIdx2] = pairsOfAtoms.size();
            pairOfAtom[atomVertexIdx1] = pairsOfAtoms.size();
            pairsOfAtoms.push_back(std::make_pair(a1.getId(), a2.getId()));
            numLinksFoundInIteration += 1;
            if (pairsOfAtoms.size() >= nrOfSliplinksToSample) {
              break;
            }
          }
          if (pairsOfAtoms.size() >= nrOfSliplinksToSample) {
            break;
          }
        }

        // add ends of chains
        std::unordered_map<size_t, igraph_integer_t> endAtomIdToVertexId;
        igraph_integer_t currentVertexId = 0;
        for (size_t i = 0; i < crosslinkerChains.size(); ++i) {
          pylimer_tools::entities::Molecule chain = crosslinkerChains[i];
          if (chain.getLength() < 2) {
            continue;
          }
          std::vector<pylimer_tools::entities::Atom> linedUpAtoms =
            chain.getAtomsLinedUp(crosslinkerType);
          if (!pylimer_tools::utils::map_has_key(endAtomIdToVertexId,
                                                 linedUpAtoms[0].getId())) {
            endAtomIdToVertexId[linedUpAtoms[0].getId()] = currentVertexId;
            currentVertexId += 1;
          }
          if (!pylimer_tools::utils::map_has_key(
                endAtomIdToVertexId,
                pylimer_tools::utils::last(linedUpAtoms).getId())) {
            endAtomIdToVertexId[pylimer_tools::utils::last(linedUpAtoms)
                                  .getId()] = currentVertexId;
            currentVertexId += 1;
          }
        }

        // create `currentVertexId` vertices for the chain-end atoms, and
        // `pairsOfAtoms.size()` vertices for the so many slip-links
        igraph_add_vertices(
          &fb.graph, currentVertexId + pairsOfAtoms.size(), nullptr);

        for (size_t chainIdx = 0; chainIdx < crosslinkerChains.size();
             ++chainIdx) {
          pylimer_tools::entities::Molecule chain = crosslinkerChains[chainIdx];
          if (chain.getLength() < 2) {
            continue;
          }

          std::vector<pylimer_tools::entities::Atom> linedUpAtoms =
            chain.getAtomsLinedUp(crosslinkerType);
          size_t previousIdx = 0;
          igraph_integer_t previousVertexId =
            endAtomIdToVertexId.at(linedUpAtoms[0].getId());
          fb.setVertexPropertiesFromAtom(
            endAtomIdToVertexId.at(linedUpAtoms[0].getId()), linedUpAtoms[0]);
          pylimer_tools::entities::Atom lastAtom =
            pylimer_tools::utils::last(linedUpAtoms);
          fb.setVertexPropertiesFromAtom(
            endAtomIdToVertexId.at(lastAtom.getId()), lastAtom);
          for (size_t i = 1; i < linedUpAtoms.size(); i++) {
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
                chain, previousIdx, i, currentEdgeId, chainIdx);
              //
              previousIdx = i;
              previousVertexId = thisVertexId;
            }
          }

          // close the chain
          igraph_integer_t currentEdgeId = igraph_ecount(&fb.graph);
          igraph_add_edge(&fb.graph,
                          previousVertexId,
                          endAtomIdToVertexId.at(lastAtom.getId()));
          fb.setBondPropertiesBasedOnChain(
            chain, previousIdx, chain.getLength() - 1, currentEdgeId, chainIdx);
        }

        // cleanup the graph
        fb.removeSubfunctionalVertices();

        // convert the graph to the network usable for simulations
        fb.finaliseInitialisation();

        return fb;
      }

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
       * @brief
       *
       * CAUTION: make sure the graph is up to date!
       *
       * @param vertexId
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d getCoordinatesForVertex(igraph_integer_t vertexId)
      {
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
        igraph_cattribute_VAN_set(&this->graph, "x", vertexId, coordinates[0]);
        igraph_cattribute_VAN_set(&this->graph, "y", vertexId, coordinates[1]);
        igraph_cattribute_VAN_set(&this->graph, "z", vertexId, coordinates[2]);
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
        const pylimer_tools::entities::Atom& atom)
      {
        assert(atom.getType() == this->crosslinkerType);
        igraph_cattribute_VAN_set(
          &this->graph, "atom_id", vertexId, atom.getId());
        igraph_cattribute_VAN_set(
          &this->graph, "type", vertexId, atom.getType());
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
        assert(atom1.getType() != this->crosslinkerType &&
               atom2.getType() != this->crosslinkerType);
        if (atom1.getId() > atom2.getId()) {
          // make sure a second call to this function would result in same
          // result
          std::swap(atom1, atom2);
        }
        igraph_cattribute_VAN_set(
          &this->graph, "type", vertexId, this->splipLinkType);
        igraph_cattribute_VAN_set(
          &this->graph, "atom_1_id", vertexId, atom1.getId());
        igraph_cattribute_VAN_set(
          &this->graph, "atom_2_id", vertexId, atom2.getId());
        Eigen::Vector3d coords = atom1.getCoordinates();
        Eigen::Vector3d dist = atom2.getCoordinates() - coords;
        this->universe.getBox().handlePBC(dist);
        coords += 0.5 * dist;
        this->setVertexCoordinates(vertexId, coords);
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
        igraph_cattribute_EAN_set(
          &this->graph, "bond_box_x", edgeId, bondBoxOffset[0]);
        igraph_cattribute_EAN_set(
          &this->graph, "bond_box_y", edgeId, bondBoxOffset[1]);
        igraph_cattribute_EAN_set(
          &this->graph, "bond_box_z", edgeId, bondBoxOffset[2]);
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
        igraph_integer_t from, to;
        igraph_edge(&this->graph, edgeId, &from, &to);
        std::vector<pylimer_tools::entities::Atom> linedUpAtoms =
          chain.getAtomsLinedUp(crosslinkerType);
        igraph_cattribute_EAN_set(&this->graph,
                                  "partition_fraction",
                                  edgeId,
                                  (atom2Idx - atom1Idx + 1) /
                                    chain.getLength());
        igraph_cattribute_EAN_set(
          &this->graph, "parent_edge", edgeId, chainIdx);
        // use the actual position of the vertices!
        Eigen::Vector3d expectedDistance =
          chain.getOverallBondSumFromTo(linedUpAtoms[atom1Idx].getId(),
                                        linedUpAtoms[atom2Idx].getId(),
                                        crosslinkerType);
        Eigen::Vector3d additionalDistance1 =
          linedUpAtoms[atom1Idx].getCoordinates() -
          this->getCoordinatesForVertex(from);
        this->universe.getBox().handlePBC(additionalDistance1);
        Eigen::Vector3d additionalDistance2 =
          this->getCoordinatesForVertex(to) -
          linedUpAtoms[atom2Idx].getCoordinates();
        this->universe.getBox().handlePBC(additionalDistance2);
        expectedDistance += additionalDistance1 + additionalDistance2;

        Eigen::Vector3d actualDistance = this->computeEdgeLength(edgeId);
        this->setBondBoxOffsetForEdge(edgeId,
                                      actualDistance - expectedDistance);
        assert(this->universe.getBox().isValidOffset(actualDistance -
                                                     expectedDistance));
        assert(this->computeEdgeLength(edgeId).isApprox(expectedDistance));
      }

      /**
       * @brief Returns the box offset for a given edge
       *
       * CAUTION: make sure the graph is up to date!
       *
       * @param edgeId
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d getBondBoxOffsetForEdge(igraph_integer_t edgeId)
      {
        Eigen::Vector3d bondBoxOffset;
        bondBoxOffset << igraph_cattribute_EAN(
          &this->graph, "bond_box_x", edgeId),
          igraph_cattribute_EAN(&this->graph, "bond_box_y", edgeId),
          igraph_cattribute_EAN(&this->graph, "bond_box_z", edgeId);

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
      Eigen::Vector3d computeEdgeLength(igraph_integer_t edgeId)
      {
        igraph_integer_t from, to;
        igraph_edge(&this->graph, edgeId, &from, &to);

        Eigen::Vector3d dist = this->getCoordinatesForVertex(to) -
                               this->getCoordinatesForVertex(from) +
                               this->getBondBoxOffsetForEdge(edgeId);
        if (this->assumeBoxLargeEnough) {
          this->universe.getBox().handlePBC(dist);
        }
        return dist;
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
      Eigen::Vector3d computeParentEdgeLength(size_t parentEdgeId)
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
            dist += this->computeEdgeLength(i);
          }
        }

        igraph_vector_destroy(&parentEdges);

        return dist;
      }

      /**
       * @brief Remove all "parent" springs that have no active "children"
       *
       * @param tolerance the acceptance tolerance, partial springs longer than
       * this are active
       */
      size_t removeInactiveParentEdges(double tolerance)
      {
        size_t numRemoved = 0;

        igraph_vector_t parentEdges;
        igraph_vector_init(&parentEdges, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "parent_edge",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &parentEdges);

        std::unordered_map<igraph_integer_t, bool> isRemovalCandidate;
        for (igraph_integer_t i = 0; i < igraph_vector_size(&parentEdges);
             ++i) {
          igraph_integer_t parentEdgeid =
            castToIgraphInt(igraph_vector_get(&parentEdges, i));
          isRemovalCandidate[parentEdgeid] = true;
        }
        for (igraph_integer_t i = 0; i < igraph_vector_size(&parentEdges);
             ++i) {
          igraph_integer_t parentEdgeid =
            castToIgraphInt(igraph_vector_get(&parentEdges, i));
          if (isRemovalCandidate.at(parentEdgeid)) {
            isRemovalCandidate[parentEdgeid] =
              this->computeEdgeLength(i).squaredNorm() <= tolerance;
          }
        }

        igraph_vector_destroy(&parentEdges);

        std::vector<igraph_integer_t> edgeIdsToRemove;
        for (auto& [key, value] : isRemovalCandidate) {
          if (value) {
            edgeIdsToRemove.push_back(key);
          }
        }

        this->removePartialSprings(edgeIdsToRemove);

        return numRemoved;
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
      size_t doRemovalAndreisWay(Eigen::VectorXd& springPartitions,
                                 double tolerance);

      /**
       * @brief Remove cross-links which do not have any springs with a certain
       * minimum length
       *
       * @param net
       * @param displacements
       * @param springPartitions
       * @param tolerance
       */
      size_t removeInactiveCrosslinks(Eigen::VectorXd& springPartitions,
                                      double tolerance);

      /**
       * @brief Remove all vertices (incl. edges!) with a functionality < 3 for
       * cross-links, < 4 for slip-links
       *
       * @param minCrosslinkFunctionalityToBeKept
       * @return size_t the number of removed vertices
       */
      size_t removeSubfunctionalVertices()
      {
        size_t numRemovedTotal = 0;
        size_t numRemovedInIteration = 0;
        do {
          numRemovedInIteration = 0;

          igraph_vector_int_t degrees;
          igraph_vector_int_init(&degrees, igraph_vcount(&this->graph));
          igraph_vector_t types;
          igraph_vector_init(&types, igraph_vcount(&this->graph));
          igraph_cattribute_VANV(
            &this->graph, "type", igraph_vss_all(), &types);
          // we count self-loops here
          // only afterwards, we check, whether they actually have relevant bond
          // box offsets
          igraph_degree(
            &this->graph, &degrees, igraph_vss_all(), IGRAPH_ALL, true);

          igraph_vector_int_t indicesToRemove;
          igraph_vector_int_init(&indicesToRemove, 0);
          for (size_t i = 0; i < igraph_vcount(&this->graph); ++i) {
            if (igraph_vector_int_get(&degrees, i) < 2) {
              igraph_vector_int_push_back(&indicesToRemove, i);
            }
          }

          igraph_delete_vertices(&this->graph,
                                 igraph_vss_vector(&indicesToRemove));

          numRemovedInIteration = igraph_vector_int_size(&indicesToRemove);

          igraph_vector_int_destroy(&indicesToRemove);

          // the function may lead to slip-links being inconsistenly
          // linked (e.g., at the end of a chain)
          // -> here, we remove those dangling ones.
          // f = 1 and below have been removed
          // -> cleanup remaining f = 2 and f = 3
          if (numRemovedInIteration > 0) {
            // reload degree and type, since the vertex index changed
            igraph_degree(
              &this->graph, &degrees, igraph_vss_all(), IGRAPH_ALL, true);
            igraph_cattribute_VANV(
              &this->graph, "type", igraph_vss_all(), &types);
          }
          for (long int i = igraph_vcount(&this->graph); i >= 0; --i) {
            if (castToIgraphInt(igraph_vector_int_get(&degrees, i)) == 2) {
              this->remove2fLink(i);
              numRemovedInIteration += 1;
            } else if (castToIgraphInt(igraph_vector_int_get(&degrees, i)) ==
                         3 &&
                       castToIgraphInt(igraph_vector_get(&types, i)) ==
                         this->splipLinkType) {
              this->remove3fLink(i);
              numRemovedInIteration += 1;
            }
          }

          igraph_vector_int_destroy(&degrees);
          igraph_vector_destroy(&types);
          numRemovedTotal += numRemovedInIteration;
        } while (numRemovedInIteration > 0);

        if (numRemovedTotal > 0) {
          this->net.isUpToDate = false;
        }
        return numRemovedTotal;
      }

      /**
       * @brief Remove double listed springs from cross-links (if they have
       * length 0)
       *
       * @param net
       * @return size_t the nr of removed edges
       */
      size_t cleanupPrimaryLoopsInStructure()
      {
        igraph_vector_int_t allEdges;
        igraph_vector_int_init(&allEdges, 2 * this->net.nrOfPartialSprings);
        if (igraph_edges(
              &this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &allEdges)) {
          throw std::runtime_error("Failed to get all edges");
        }

        std::vector<igraph_integer_t> edgesToRemove;
        for (igraph_integer_t i = 0; i < igraph_ecount(&this->graph); ++i) {
          if (igraph_vector_int_get(&allEdges, 2 * i + 0) ==
              igraph_vector_int_get(&allEdges, 2 * i + 1)) {
            // check if this primary loop has length 0 -> remove
            if (this->getBondBoxOffsetForEdge(i).norm() <
                1e-3 * this->universe.getBox().getL().minCoeff()) {
              edgesToRemove.push_back(i);
            }
          }
        }

        this->removePartialSprings(edgesToRemove);
        return edgesToRemove.size();
      };

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
       * NOTE: this may result in slip-links with f = 2.
       * They are not automatically removed in order to preserve vertex
       * iterations.
       *
       * @param net
       * @param springPartitions
       */
      void removeParentSpring(const size_t springIdx)
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

        std::vector<igraph_integer_t> edgeIdsToRemove;
        edgeIdsToRemove.reserve(4);
        for (igraph_integer_t i = 0; i < igraph_vector_size(&parentEdges);
             ++i) {
          if (igraph_vector_get(&parentEdges, i) == springIdx) {
            edgeIdsToRemove.push_back(i);
          }
        }

        igraph_vector_destroy(&parentEdges);
        // actually remove all edges
        this->removePartialSprings(edgeIdsToRemove);
        // as the above may lead to slip-links with f = 2,
        // we want to remove those as well
        // this->removeTwofunctionalLinks();
      };

      /**
       * @brief Remove a set of edges from the graph
       *
       * @param edgeIdsToRemove
       */
      void removePartialSprings(std::vector<igraph_integer_t>& edgeIdsToRemove)
      {
        igraph_vector_int_t edgeIdsToRemoveVec;
        igraph_vector_int_init(&edgeIdsToRemoveVec, edgeIdsToRemove.size());
        pylimer_tools::utils::StdVectorToIgraphVectorT(edgeIdsToRemove,
                                                       &edgeIdsToRemoveVec);
        igraph_es_t iterator;
        igraph_delete_edges(&this->graph,
                            igraph_ess_vector(&edgeIdsToRemoveVec));
        igraph_vector_int_destroy(&edgeIdsToRemoveVec);
        // remove vertices that "got lost"
        this->removeOrphanedVertices();

        this->net.isUpToDate = false;
      }

      /**
       * @brief marks a certain "parent" spring as non-existing
       */
      void combineParentSprings(size_t springIdxBefore, size_t springIdxNow);

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
                                                   double alpha)
      {
        igraph_vector_int_t edges;
        igraph_vector_int_init(&edges, 0);

        this->findEdgesAndVerticesOfSpring(springIdx, nullptr, &edges);

        double currentAlpha = 0.0;
        igraph_integer_t previousEdgeId = igraph_vector_int_get(&edges, 0);
        for (size_t i = 1; i < igraph_vector_int_size(&edges); ++i) {
          igraph_integer_t edgeId = igraph_vector_int_get(&edges, i);
          currentAlpha +=
            igraph_cattribute_EAN(&this->graph, "partition_fraction", edgeId);
          if (currentAlpha > alpha) {
            break;
          }
          previousEdgeId = edgeId;
        }

        igraph_vector_int_destroy(&edges);

        return previousEdgeId;
      }

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

      void relaxationLight(ForceBalanceNetwork& net,
                           Eigen::VectorXd& springPartitions,
                           const size_t linkIdx,
                           const double oneOverSpringPartitionUpperLimit = 1.0)
      {
        Eigen::VectorXd oneOverSpringPartitions = Eigen::VectorXd::Zero(0);
        this->relaxationLight(net,
                              springPartitions,
                              oneOverSpringPartitions,
                              linkIdx,
                              oneOverSpringPartitionUpperLimit);
      };

      void relaxationLight(ForceBalanceNetwork& net,
                           Eigen::VectorXd& springPartitions,
                           Eigen::VectorXd& oneOverSpringPartitions,
                           const size_t linkIdx,
                           const double oneOverSpringPartitionUpperLimit = 1.0)
      {
        assert(net.isUpToDate);
        if (net.linkIsSliplink[linkIdx]) {
          this->updateSpringPartition(net,
                                      springPartitions,
                                      oneOverSpringPartitions,
                                      linkIdx,
                                      oneOverSpringPartitionUpperLimit);
        }
        this->displaceToMeanPosition(
          net, springPartitions, linkIdx, oneOverSpringPartitionUpperLimit);
      };

      /**
       * @brief Replace the two springs traversinga a two-functional cross-links
       * with a single spring
       *
       * @param net
       * @param displacements
       * @param springPartitions
       */
      size_t removeTwofunctionalLinks()
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
      void deformTo(const pylimer_tools::entities::Box& newBox)
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

      size_t getNumExtraAtoms() override
      {
        return this->getNrOfLinks() - this->getNrOfNodes();
      }

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
        return this->evaluatePartialSpringDistances(this->net, this->is2D);
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
      Eigen::VectorXd evaluateSpringDistances(const ForceBalanceNetwork& net,
                                              const bool is2D = false) const;

      Eigen::VectorXd evaluatePartialSpringDistances(
        const ForceBalanceNetwork& net,
        const bool is2D = false) const;

      /**
       * @brief Compute the distance between two links
       *
       * @param net
       * @param linkIndexA
       * @param linkIndexB
       * @param is2D
       * @return Eigen::Vector3d
       */
      Eigen::Vector3d evaluateDistanceBetween(const ForceBalanceNetwork& net,
                                              const size_t linkIndexA,
                                              const size_t linkIndexB,
                                              bool is2D = false) const
      {
        assert(net.isUpToDate);

        Eigen::Vector3d distances = net.coordinates.segment(3 * linkIndexB, 3) -
                                    net.coordinates.segment(3 * linkIndexA, 3);

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
      Eigen::Vector3d evaluatePartialSpringDistance(
        const ForceBalanceNetwork& net,
        const size_t springIdx,
        bool is2d = false) const
      {
        return this->evaluatePartialSpringDistanceTo(
          net, springIdx, net.springIndexB(springIdx), is2d);
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
      Eigen::Vector3d evaluatePartialSpringDistanceTo(
        const ForceBalanceNetwork& net,
        const size_t springIdx,
        const size_t linkIdx,
        bool is2d = false) const
      {
        assert(net.isUpToDate);
        assert(net.springIndexA(springIdx) == linkIdx ||
               net.springIndexB(springIdx) == linkIdx);

        Eigen::Vector3d dist =
          net.coordinates.segment(3 * net.springIndexB(springIdx), 3) -
          net.coordinates.segment(3 * net.springIndexA(springIdx), 3) +
          net.springPartBoxOffset.segment(3 * springIdx, 3);

        if (this->assumeBoxLargeEnough) {
          this->universe.getBox().handlePBC<Eigen::Vector3d>(dist);
        }

        if (is2D) {
          dist[2] = 0.0;
        }

        return dist * (net.springIndexA(springIdx) == linkIdx ? -1. : 1.);
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
      Eigen::Vector3d evaluatePartialSpringDistanceFrom(
        const ForceBalanceNetwork& net,
        const size_t springIdx,
        const size_t linkIdx,
        bool is2d = false) const
      {
        return -1. * this->evaluatePartialSpringDistanceTo(
                       net, springIdx, linkIdx, is2d);
      }

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
        bool allowSlipLinksToPassEachOther = false) const
      {

        INVALIDARG_EXP_IFN(linkIdx < net.springIndicesOfLinks.size(),
                           "Link to update needs to be in the list");
        INVALIDARG_EXP_IFN(net.linkIsSliplink[linkIdx],
                           "Only slip-links may slip along a spring, link " +
                             std::to_string(linkIdx) +
                             " is not one. Network has " +
                             std::to_string(net.nrOfNodes) + " cross- of " +
                             std::to_string(net.nrOfLinks) + " links.");
        std::vector<size_t> springIndices = net.springIndicesOfLinks[linkIdx];
        assert(springIndices.size() == 1 || springIndices.size() == 2);
        if (springIndices.size() == 2 && springIndices[0] == springIndices[1]) {
          springIndices.pop_back();
        }
        double residualNorm = 0.0;
        int residualNormContributions = 0;
        for (size_t springIndex : springIndices) {
          std::vector<size_t> springsPartners =
            net.linkIndicesOfSprings[springIndex];
          for (size_t partner_idx = 1; partner_idx < springsPartners.size() - 1;
               ++partner_idx) {
            if (springsPartners[partner_idx] == linkIdx) {
              size_t backSpringGlobalIdx =
                net.localToGlobalSpringIndex[springIndex][partner_idx - 1];
              size_t forwardSpringGlobalIdx =
                net.localToGlobalSpringIndex[springIndex][partner_idx];
              // found position of this link in this spring
              // want to find the ideal value for
              // net.springPartitions[springIndex][partner_idx-1]
              Eigen::Vector3d vecBack = this->evaluatePartialSpringDistanceTo(
                net,
                backSpringGlobalIdx,
                springsPartners[partner_idx - 1],
                this->is2D);
              double distanceBack = vecBack.squaredNorm();
              Eigen::Vector3d vecForward =
                this->evaluatePartialSpringDistanceTo(
                  net,
                  forwardSpringGlobalIdx,
                  springsPartners[partner_idx + 1],
                  this->is2D);
              double distanceForward = vecForward.squaredNorm();
              double idealValue =
                1. / (1. + sqrt(distanceForward / distanceBack));
              if (distanceBack <= 0.0) {
                idealValue = 0.0; // TODO: really?
              }
              size_t currentSpringGlobalIdx =
                net.localToGlobalSpringIndex[springIndex][partner_idx - 1];
              size_t neighbourSpringGlobalIdx =
                net.localToGlobalSpringIndex[springIndex][partner_idx];
              double currentS = springPartitions[currentSpringGlobalIdx];
              double nextS = springPartitions[neighbourSpringGlobalIdx];
              const double N = net.springsContourLength[springIndex];
              const double l = (currentS + nextS);
              if (oneOverSpringPartitionUpperLimit > 0.) {
                // TODO: sketch theory why this should/not be necessary!!!
                const double limit =
                  std::clamp(1. / (oneOverSpringPartitionUpperLimit *
                                   (nextS + currentS) * (N)),
                             0.,
                             1.);
                idealValue = std::clamp(idealValue, limit, 1. - limit);
                // double oneOverCurrent = 1. / (currentS * N);
                // double oneOverNext = 1. / (nextS * N);
                // double limitedOneOverCurrent =
                // CLAMP_ONE_OVER_SPRINGPARTITION(
                //   true, oneOverCurrent, N, oneOverSpringPartitionUpperLimit);
                // double limitedOneOverNext = CLAMP_ONE_OVER_SPRINGPARTITION(
                //   true, oneOverNext, N, oneOverSpringPartitionUpperLimit);
                // currentS = (1./N) * 1. / limitedOneOverCurrent;
                // nextS = (1./N) * 1. / limitedOneOverNext;
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
              // if ((1. - idealValue) != 0. && l > 0.) {
              //   double term1 = -
              //     (distanceForward / ((1. - idealValue) * (1. -
              //     idealValue))); localResidualNorm += term1;
              //     std::cout.precision(std::numeric_limits<double>::max_digits10);
              //     std::cout << "localResidualNorm term 1: " << term1 <<
              //     std::endl;

              // } else {
              //   std::cout << "localResidualNorm Case 1: " << l << " "
              //             << idealValue << std::endl;
              // }
              // if (idealValue != 0. && l > 0.) {
              //   double term2 =
              //     (distanceBack / (idealValue * idealValue));
              //     localResidualNorm += term2;
              //     std::cout << "localResidualNorm term 2: " << term2 <<
              //     std::endl;

              // } else {
              //   std::cout << "localResidualNorm Case 2: " << l << " "
              //             << idealValue << std::endl;
              // }
              localResidualNorm /= (N * l);
              // std::cout
              //   << "localResidualNorm val 2: "
              //   << localResidualNorm
              //                  << std::endl;

              RUNTIME_EXP_IFN(
                APPROX_WITHIN(newS + complementaryS, 0., 1., 1e-9),
                "Require newS + complementaryS to be within 0, 1, got " +
                  std::to_string(newS + complementaryS) + " from " +
                  std::to_string(newS) + " and " +
                  std::to_string(complementaryS) +
                  " with ideal = " + std::to_string(idealValue) + " of " +
                  std::to_string(nextS + currentS) + " for link " +
                  std::to_string(linkIdx) + ". Diff: " +
                  std::to_string(1. - (newS + complementaryS)) + ".");
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
              RUNTIME_EXP_IFN(
                nextS >= -0.00000000000001,
                "nextS must be >= 0., got " + std::to_string(nextS) + " from " +
                  std::to_string(nextS) + " and " + std::to_string(currentS) +
                  ", " + std::to_string(newS) + " and " +
                  std::to_string(complementaryS) + ".");
              RUNTIME_EXP_IFN(complementaryS >= -0.00000000000001,
                              "complementaryS must be >= 0., got " +
                                std::to_string(complementaryS) + " from " +
                                std::to_string(nextS) + " and " +
                                std::to_string(currentS) + ", " +
                                std::to_string(newS) + " and " +
                                std::to_string(complementaryS) + ".");

              // (complementaryS > residualNormSTolerance &&
              //  newS > residualNormSTolerance)
              //   ? ( -
              //      distanceBack / (newS * newS))
              //   : 0.0;
              // if (!(APPROX_EQUAL(newS, currentS, 0.2))) {
              //   std::cout << "Updating " << linkIdx << " to " << newS << "
              //   and
              //   "
              //             << complementaryS << " with global springs "
              //             << currentSpringGlobalIdx << " and "
              //             << neighbourSpringGlobalIdx << " from " << currentS
              //             << ", " << nextS << std::endl;
              // }
              // std::cout
              //   << "Contribution to "
              //   << linkIdx
              //           << " from global springs " <<
              //           currentSpringGlobalIdx
              //           << " (" << springsPartners[partner_idx - 1] <<
              //           ") "
              //           << vecBack[0] << ", " << vecBack[1] << ", " <<
              //           vecBack[2]
              //           << " and " << neighbourSpringGlobalIdx << " ("
              //           << springsPartners[partner_idx + 1] << ") "
              //           << vecForward[0] << ", " << vecForward[1] << ",
              //           "
              //           << vecForward[2] << "; "
              //           << " with " << currentS << ", " << nextS <<
              //           std::endl;
              //        std::cout
              // << "Distances are " << distanceForward
              // << ", "
              //           << distanceBack << " to get ideal value " <<
              //           idealValue
              //           << " for " << (nextS) << " , " << currentS <<
              //           std::endl;
              residualNorm += localResidualNorm * localResidualNorm;
              springPartitions[currentSpringGlobalIdx] = newS;
              springPartitions[neighbourSpringGlobalIdx] = complementaryS;
              if (oneOverSpringPartitions.size() > 0) {
                double primaryCorrectionMultiplierC = static_cast<double>(
                  net.springPartIndexA[currentSpringGlobalIdx] !=
                  net.springPartIndexB[currentSpringGlobalIdx]);
                double oneOverCurrent =
                  primaryCorrectionMultiplierC *
                  CLAMP_ONE_OVER_SPRINGPARTITION(
                    net.partialSpringIsPartial[currentSpringGlobalIdx],
                    (1.0 / (newS * N)),
                    N,
                    oneOverSpringPartitionUpperLimit);
                oneOverSpringPartitions.segment(3 * currentSpringGlobalIdx, 3) =
                  Eigen::Vector3d::Constant(oneOverCurrent);
                double primaryCorrectionMultiplierN = static_cast<double>(
                  net.springPartIndexA[neighbourSpringGlobalIdx] !=
                  net.springPartIndexB[neighbourSpringGlobalIdx]);
                double oneOverNeighbour =
                  primaryCorrectionMultiplierN *
                  CLAMP_ONE_OVER_SPRINGPARTITION(
                    net.partialSpringIsPartial[neighbourSpringGlobalIdx],
                    1.0 / (complementaryS * N),
                    N,
                    oneOverSpringPartitionUpperLimit);
                oneOverSpringPartitions.segment(3 * neighbourSpringGlobalIdx,
                                                3) =
                  Eigen::Vector3d::Constant(oneOverNeighbour);
              }
            }
          }
        }
        assert(residualNormContributions == 4);
        return residualNorm;
      };

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
        ForceBalanceNetwork& net,
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

      double displaceLinksToMeanPosition(ForceBalanceNetwork& net,
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
        ForceBalanceNetwork& net,
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
        ForceBalanceNetwork& net,
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
        return this->getStressTensor(1.0);
      }

      Eigen::Matrix3d getStressTensor(double oneOverSpringPartitionUpperLimit)
      {
        std::array<std::array<double, 3>, 3> stressTensor =
          this->evaluateStressTensor(this->net,
                                     this->currentSpringPartitionsVec,
                                     this->kappa,
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
          this->evaluateStressTensorLinkBased(this->net,
                                              this->currentSpringPartitionsVec,
                                              this->kappa,
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
       * @brief Compute a few properties of the simulator
       *
       */
      void finaliseInitialisation()
      {
        this->convertFromGraph();
        this->currentSpringDistances =
          this->evaluateSpringDistances(this->net, this->is2D);
        this->currentPartialSpringDistances =
          this->evaluatePartialSpringDistances(this->net, this->is2D);
        this->defaultR0Squared =
          this->universe.computeMeanSquareEndToEndDistance(
            this->crosslinkerType);
        this->defaultNrOfChains =
          this->universe.getMolecules(this->crosslinkerType).size();
        this->validateNetwork();
      }

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

        igraph_vector_t parentEdges;
        igraph_vector_init(&parentEdges, this->net.nrOfPartialSprings);
        igraph_cattribute_EANV(&this->graph,
                               "parent_edge",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &parentEdges);

        std::unordered_set<igraph_integer_t> parentEdgeIds;
        for (size_t i = 0; i < igraph_vector_size(&edgeTypes); ++i) {
          parentEdgeIds.insert(
            castToIgraphInt(igraph_vector_get(&parentEdges, i)));
        }
        this->net.nrOfSprings = parentEdgeIds.size();
        size_t numPartialSprings =
          this->net.nrOfPartialSprings - this->net.nrOfSprings;

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
            castToIgraphInt(igraph_vector_get(&parentEdges, i));
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

          std::vector<size_t> linkIndicesOfThisSpring;
          pylimer_tools::utils::igraphVectorTToStdVector(
            &verticesOnPath, linkIndicesOfThisSpring);
          this->net.linkIndicesOfSprings[i] = linkIndicesOfThisSpring;
          std::vector<size_t> edgeIdsOfThisSpring;
          pylimer_tools::utils::igraphVectorTToStdVector(&edgesOnPath,
                                                         edgeIdsOfThisSpring);
          this->net.localToGlobalSpringIndex[i] = edgeIdsOfThisSpring;

          igraph_vector_int_destroy(&edgesOnPath);
          igraph_vector_int_destroy(&verticesOnPath);
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
            castToIgraphInt(igraph_vector_get(&linkType, i)) ==
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

      /**
       * @brief Write the new coordinates to the graph
       */
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
        Eigen::VectorXd& springPartitions,
        const size_t partialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0);

      bool swapSlipLinkWithXlinkReversibly(
        ForceBalanceNetwork& net,
        Eigen::VectorXd& springPartitions,
        const size_t partialSpringIdx,
        const double oneOverSpringPartitionUpperLimit = 1.0,
        const bool respectLoops = true);
    };
  } // namespace mehp
} // namespace calc
} // namespace pylimer_tools
#endif
