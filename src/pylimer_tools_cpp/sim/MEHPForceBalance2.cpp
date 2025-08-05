#include "MEHPForceBalance2.h"
#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
// #include "../utils/MemoryUtil.h"
#include "../utils/ExtraEigenSolvers.h"
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <algorithm>
#include <array>
#include <cassert>
#include <future>
#include <iostream>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// #ifndef NDEBUG
// #define DEBUG_REMOVAL
// #endif

namespace pylimer_tools::sim::mehp {
MEHPForceBalance2::MEHPForceBalance2(
  const pylimer_tools::entities::Universe& u,
  const pylimer_tools::topo::entanglement_detection::AtomPairEntanglements&
    entanglements,
  const int crossLinkerType,
  const bool is2D,
  const bool entanglementsAsSprings)
{
  this->universe = u;
  this->is2D = is2D;
  this->crossLinkerType = crossLinkerType;
  this->box = u.getBox();

  // start to initialize the network
  this->initialConfig.L[0] = this->box.getL()[0];
  this->initialConfig.L[1] = this->box.getL()[1];
  this->initialConfig.L[2] = this->box.getL()[2];

  this->initialConfig.boxHalfs[0] = 0.5 * this->initialConfig.L[0];
  this->initialConfig.boxHalfs[1] = 0.5 * this->initialConfig.L[1];
  this->initialConfig.boxHalfs[2] = 0.5 * this->initialConfig.L[2];

  // move the entanglements to better accessible variables
  const std::vector<std::pair<size_t, size_t>> pairsOfAtoms =
    entanglements.pairsOfAtoms;
  std::vector<long int> pairOfAtom = entanglements.pairOfAtom;

  if (pairOfAtom.size() != u.getNrOfAtoms()) {
    assert(pairOfAtom.size() == 0);
    pairOfAtom =
      pylimer_tools::utils::initializeWithValue<long int>(u.getNrOfAtoms(), -1);
  }

  // use directly the graph for performance
  igraph_t graph = u.getCopyOfGraph();

  // what we do:
  // 1. find starting points for depth-first search: degree > 2
  // 2. do depth-first search from these points. Whenever we encounter a
  // vertex with f > 2, add it to the network.
  // Whenever we encounter a vertex that's an entanglement,
  // add it to the network,
  // Whenever we encounter a vertex that's been encountered,
  // make sure to adjust the offset of the bond.
  std::vector<int> vertexDegrees = u.getVertexDegrees();
  assert(vertexDegrees.size() == u.getNrOfAtoms());
  std::vector<igraph_integer_t> startingVertices;
  startingVertices.reserve(vertexDegrees.size());
  // we will want to start from all vertices, in a sense
  // but first, start from those with degree > 2
  for (size_t vIdx = 0; vIdx < vertexDegrees.size(); ++vIdx) {
    if (vertexDegrees[vIdx] > 2) {
      startingVertices.push_back(vIdx);
    }
  }
  // then, from the others, except entanglements.
  // this may be problematic in some ultra-high-entanglement-density cases
  // with free primary loops, but should be a reasonably improbable scenario.
  for (size_t vIdx = 0; vIdx < vertexDegrees.size(); ++vIdx) {
    if (vertexDegrees[vIdx] <= 2 && pairOfAtom[vIdx] < 0) {
      startingVertices.push_back(vIdx);
    }
  }
  assert(startingVertices.size() <= u.getNrOfAtoms());

  // some more validation
  for (size_t i = 0; i < pairOfAtom.size(); ++i) {
    if (pairOfAtom[i] >= 0) {
      INVALIDARG_EXP_IFN(vertexDegrees[i] == 2,
                         "Only 'normal' strand bead are currently allowed "
                         "entanglement targets.");
      INVALIDARG_EXP_IFN(pairsOfAtoms[pairOfAtom[i]].first !=
                           pairsOfAtoms[pairOfAtom[i]].second,
                         "");
    }
  }

  // the basic coordinates to use afterward
  Eigen::VectorXd vertexCoordinates =
    u.getUnwrappedVertexCoordinates(this->box);
  assert(vertexCoordinates.size() == u.getNrOfAtoms() * 3);
  Eigen::VectorXd effectiveCoordinates = vertexCoordinates;

  std::vector<long int> oldVertexIdToNewLinkId =
    pylimer_tools::utils::initializeWithValue<long int>(u.getNrOfAtoms(), -1);
  std::vector<bool> vertex_visited =
    pylimer_tools::utils::initializeWithValue(u.getNrOfAtoms(), false);
  std::vector<bool> edge_followed =
    pylimer_tools::utils::initializeWithValue(u.getNrOfBonds(), false);

  igraph_adjlist_t adjacencyList;
  igraph_adjlist_init(
    &graph, &adjacencyList, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);
  igraph_inclist_t incidenceList;
  igraph_inclist_init(&graph, &incidenceList, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

  size_t nVerticesWithRelevantDegree = (std::ranges::count_if(
    vertexDegrees, [](const int d) { return d > 0 && d != 2; }));
  size_t nLinksExpected = nVerticesWithRelevantDegree + pairsOfAtoms.size();
  size_t nStrandsExpected = nVerticesWithRelevantDegree / 2;
  size_t nSpringsExpected =
    nVerticesWithRelevantDegree / 2 + pairsOfAtoms.size();

  // data targets
  std::vector<size_t> springFromVertexIdx;
  springFromVertexIdx.reserve(nSpringsExpected);
  std::vector<size_t> springToVertexIdx;
  springToVertexIdx.reserve(nSpringsExpected);
  std::vector<size_t> springContourLength;
  springContourLength.reserve(nSpringsExpected);
  std::vector<Eigen::Vector3d> springBoxOffsets;
  springBoxOffsets.reserve(nSpringsExpected);
  std::vector<Eigen::Vector3d> newLinkPositions;
  std::vector<size_t> newLinkIdxToOldVertexIdx;
  newLinkIdxToOldVertexIdx.reserve(nLinksExpected);
  ArrayXArrayXi springIndicesOfStrand;
  springIndicesOfStrand.reserve(nStrandsExpected);
  ArrayXArrayXi linkIndicesOfStrand;

  long int nextLinkIdx = 0;
  long int currentStrandIdx = 0;
  long int currentSpringContourLength = 0;
  for (igraph_integer_t vIdx : startingVertices) {
    if (vertex_visited[vIdx]) {
      continue;
    }
    vertex_visited[vIdx] = true;
    // perform depth-first search from vIdx
    std::stack<igraph_integer_t> vertexStack;
    vertexStack.push(vIdx);
    std::stack<size_t> newLinkIdStack;
    std::stack<size_t> addedVertexIdStack;
    Eigen::Vector3d currentPosition = vertexCoordinates.segment<3>(3 * vIdx, 3);
    while (!vertexStack.empty()) {
      assert(addedVertexIdStack.size() == newLinkIdStack.size());
      igraph_integer_t currentVertex = vertexStack.top();

      Eigen::Vector3d distanceToCurrent =
        effectiveCoordinates.segment<3>(currentVertex * 3, 3) - currentPosition;
      this->box.handlePBC(distanceToCurrent);
      currentPosition += distanceToCurrent;

      igraph_vector_int_t* neighbors =
        igraph_adjlist_get(&adjacencyList, currentVertex);
      igraph_vector_int_t* incident_edges =
        igraph_inclist_get(&incidenceList, currentVertex);
      assert(igraph_vector_int_size(neighbors) ==
             igraph_vector_int_size(incident_edges));
      assert(igraph_vector_int_size(neighbors) == vertexDegrees[currentVertex]);

      // check whether we should remember this vertex
      // only if we don't remember it yet, ...
      if (oldVertexIdToNewLinkId[currentVertex] == -1) {
        bool addNewLink = false;
        bool addSpring = false;
        bool addNewStrand = false;
        // three cases on when to add
        if (newLinkIdStack.empty()) {
          // 1. first (end / central)
          addNewLink = true;
          assert(vertexStack.size() <= 1);
        } else if (pairOfAtom[currentVertex] >= 0) {
          // 2. entanglement
          auto [fst, snd] = pairsOfAtoms[pairOfAtom[currentVertex]];
          fst = this->universe.getIdxByAtomId(fst);
          snd = this->universe.getIdxByAtomId(snd);
          assert(fst == currentVertex || snd == currentVertex);
          addNewLink =
            entanglementsAsSprings || (oldVertexIdToNewLinkId[fst] < 0 &&
                                       oldVertexIdToNewLinkId[snd] < 0);
          if (!addNewLink) {
            assert(XOR(oldVertexIdToNewLinkId[fst] >= 0,
                       oldVertexIdToNewLinkId[snd] >= 0));
            // make this vertex identical to the other vertex in the
            // entanglement link
            oldVertexIdToNewLinkId[currentVertex] = std::max(
              oldVertexIdToNewLinkId[fst], oldVertexIdToNewLinkId[snd]);
          }
          // since the paired atom may have a different position, we need to
          // account for that here for PBC later to work.
          Eigen::Vector3d dist = effectiveCoordinates.segment(
                                   3 * (fst == currentVertex ? snd : fst), 3) -
                                 currentPosition;
          this->box.handlePBC(dist);
          currentPosition += dist * (addNewLink ? 0.5 : 1.);
          effectiveCoordinates.segment(3 * currentVertex, 3) = currentPosition;
          addSpring = true;
        } else if (vertexDegrees[currentVertex] > 2 ||
                   vertexDegrees[currentVertex] == 1) {
          // 3. crosslink/end
          addNewLink = true;
          addSpring = true;
          // close the previous strand, start a new strand
          addNewStrand = true;
        }

        if (addSpring) {
          assert(!newLinkIdStack.empty());
          // link this new link to the previous link
          if (springIndicesOfStrand.size() <= currentStrandIdx) {
            assert(currentStrandIdx == springIndicesOfStrand.size());
            springIndicesOfStrand.push_back(
              std::vector<size_t>{ springFromVertexIdx.size() });
            linkIndicesOfStrand.push_back(
              std::vector<size_t>{ newLinkIdStack.top() });
          } else {
            springIndicesOfStrand[currentStrandIdx].push_back(
              springFromVertexIdx.size());
          }
          springFromVertexIdx.push_back(addedVertexIdStack.top());

          if (oldVertexIdToNewLinkId[currentVertex] >= 0) {
            assert(!addNewLink);
            assert(!entanglementsAsSprings);
            assert(pairOfAtom[currentVertex] >= 0);
            effectiveCoordinates.segment(3 * currentVertex, 3) =
              currentPosition;
          } else {
            assert(addNewLink);
          }
          springToVertexIdx.push_back(currentVertex);
          linkIndicesOfStrand[currentStrandIdx].push_back(
            addNewLink ? nextLinkIdx : oldVertexIdToNewLinkId[currentVertex]);

          assert(linkIndicesOfStrand[currentStrandIdx].size() - 1 ==
                 springIndicesOfStrand[currentStrandIdx].size());
          assert(linkIndicesOfStrand[currentStrandIdx].size() >= 2);
          assert(springFromVertexIdx.size() == springToVertexIdx.size());

          // each entanglement link has two original bead positions,
          // which may differ by multiple boxes.
          // should be distance: The distance that is "correct" by following the
          // strand
          Eigen::Vector3d shouldBeDistance =
            currentPosition -
            effectiveCoordinates.segment<3>(3 * springFromVertexIdx.back(), 3);
          // is distance: The distance that would be calculated using the
          // coordinates of the beads
          Eigen::Vector3d isDistance =
            (addNewLink
               ? currentPosition
               : newLinkPositions[oldVertexIdToNewLinkId[currentVertex]]) -
            newLinkPositions
              [oldVertexIdToNewLinkId[springFromVertexIdx.back()]];

          Eigen::Vector3d offset = shouldBeDistance - isDistance;
          assert(this->box.isValidOffset(offset));
          springBoxOffsets.emplace_back(
            offset); // addNewLink ? Eigen::Vector3d::Zero()
          // : offset
          springContourLength.push_back(currentSpringContourLength);
          currentSpringContourLength = 0;

          // we have visited this, need to remember for next springs
          if (!addNewLink) {
            newLinkIdStack.push(
              oldVertexIdToNewLinkId[springToVertexIdx.back()]);
            addedVertexIdStack.push(currentVertex);
          }
        }
        if (addNewStrand) {
          assert(addSpring);
          assert(addNewLink);
          currentStrandIdx += 1;
        }
        if (addNewLink) {
          // add to network
          oldVertexIdToNewLinkId[currentVertex] = nextLinkIdx;
          newLinkIdxToOldVertexIdx.push_back(currentVertex);
          newLinkPositions.push_back(currentPosition);
          effectiveCoordinates.segment(3 * currentVertex, 3) = currentPosition;

          newLinkIdStack.push(nextLinkIdx);
          addedVertexIdStack.push(currentVertex);
          nextLinkIdx += 1;
          assert(newLinkIdxToOldVertexIdx.size() == nextLinkIdx);
          assert(newLinkIdxToOldVertexIdx.size() == newLinkPositions.size());
        }
      }

      assert(newLinkIdStack.top() ==
             oldVertexIdToNewLinkId[addedVertexIdStack.top()]);

      bool goDeeper = false;
      for (igraph_integer_t i = 0; i < igraph_vector_int_size(neighbors); ++i) {
        igraph_integer_t neighbor = igraph_vector_int_get(neighbors, i);
        igraph_integer_t edge = igraph_vector_int_get(incident_edges, i);

#ifndef NDEBUG
        // validate that the inclist & adjlist are actually compliant
        // by checking that such an edge exists in the graph
        igraph_integer_t a = 0;
        igraph_integer_t b = 0;
        igraph_edge(&graph, edge, &a, &b);
        assert(a == currentVertex || b == currentVertex);
        assert(a == neighbor || b == neighbor);
#endif

        assert(neighbor < u.getNrOfAtoms());

        // check if neighbor has been visited
        if (vertex_visited[neighbor] && edge_followed[edge]) {
          // this is the base case whenever we re-iterate this loop.
          // there are ways to avoid these additional iterations,
          // but they involve resizing the vectors,
          // which is generally for how many neighbours we have
          // not worth it, presumably.
          continue;
        } else if (vertex_visited[neighbor] == false) {
          // this is the case whenever we encounter a new vertex
          assert(edge_followed[edge] == false);
          // continue in this direction
          goDeeper = true;
          currentSpringContourLength += 1;
          vertex_visited[neighbor] = true;
          edge_followed[edge] = true;
          vertexStack.push(neighbor);
          // let's directly jump to this neighbor, search its neighbours
          break;
        } else {
          // this is the case for loops
          assert(vertex_visited[neighbor] == true);
          assert(edge_followed[edge] == false);
          assert(oldVertexIdToNewLinkId[neighbor] >= 0);
          // remember this bond (spring)
          if (springIndicesOfStrand.size() <= currentStrandIdx) {
            assert(currentStrandIdx == springIndicesOfStrand.size());
            springIndicesOfStrand.push_back(
              std::vector<size_t>{ springFromVertexIdx.size() });
            linkIndicesOfStrand.push_back(
              std::vector<size_t>{ newLinkIdStack.top() });
          } else {
            springIndicesOfStrand[currentStrandIdx].push_back(
              springFromVertexIdx.size());
          }
          linkIndicesOfStrand[currentStrandIdx].push_back(
            oldVertexIdToNewLinkId[neighbor]);
          assert(linkIndicesOfStrand[currentStrandIdx].size() - 1 ==
                 springIndicesOfStrand[currentStrandIdx].size());
          springFromVertexIdx.push_back(addedVertexIdStack.top());
          springToVertexIdx.push_back(neighbor);
          assert(linkIndicesOfStrand[currentStrandIdx].size() >= 2);
          assert(springFromVertexIdx.size() == springToVertexIdx.size());
          springContourLength.push_back(currentSpringContourLength + 1);
          edge_followed[edge] = true;
          //
          Eigen::Vector3d targetPos = effectiveCoordinates.segment(
            3 * neighbor,
            3); // newLinkPositions[oldVertexIdToNewLinkId[neighbor]]
          Eigen::Vector3d lastStepToTarget = currentPosition - targetPos;
          this->box.handlePBC(lastStepToTarget);
          Eigen::Vector3d actualTargetPos = currentPosition - lastStepToTarget;

          // should be distance: The distance that is "correct" by following the
          // strand
          Eigen::Vector3d shouldBeDistance =
            actualTargetPos -
            effectiveCoordinates.segment<3>(3 * springFromVertexIdx.back(), 3);
          // is distance: The distance that would be calculated using the
          // coordinates of the beads
          Eigen::Vector3d isDistance =
            newLinkPositions[oldVertexIdToNewLinkId[neighbor]] -
            newLinkPositions
              [oldVertexIdToNewLinkId[springFromVertexIdx.back()]];

          // compute box offset
          // Eigen::Vector3d distance = targetPos - previousPos;
          // springBoxOffsets.emplace_back(this->box.getOffset(-1. * distance));
          springBoxOffsets.emplace_back(shouldBeDistance - isDistance);
          assert(this->box.isValidOffset(springBoxOffsets.back()));

          // close the strand
          currentStrandIdx += 1;
          currentSpringContourLength = 0;
        }
      }

      if (!goDeeper) {
        // unroll to the previous junction
        if (newLinkIdStack.empty()) {
          break;
        }
        if (oldVertexIdToNewLinkId[currentVertex] == newLinkIdStack.top()) {
          newLinkIdStack.pop();
          addedVertexIdStack.pop();
          if (newLinkIdStack.empty()) {
            break;
          }
        }
        while (oldVertexIdToNewLinkId[currentVertex] != newLinkIdStack.top()) {
          vertexStack.pop();
          currentVertex = vertexStack.top();
        }
        currentPosition = newLinkPositions[newLinkIdStack.top()];
        // if we unroll, we must have reached the end of a strand
        assert(currentSpringContourLength == 0);
      }
    }
  }

  igraph_destroy(&graph);
  igraph_inclist_destroy(&incidenceList);
  igraph_adjlist_destroy(&adjacencyList);

  if (entanglementsAsSprings) {
    for (const auto& [fst, snd] : pairsOfAtoms) {
      size_t fstVertexIdx = this->universe.getIdxByAtomId(fst);
      size_t sndVertexIdx = this->universe.getIdxByAtomId(snd);
      const long int newLink1 = oldVertexIdToNewLinkId[fstVertexIdx];
      const long int newLink2 = oldVertexIdToNewLinkId[sndVertexIdx];
      assert(newLink1 >= 0);
      assert(newLink2 >= 0);
      springFromVertexIdx.push_back(fstVertexIdx);
      springToVertexIdx.push_back(sndVertexIdx);
      springContourLength.push_back(1.);
      springBoxOffsets.emplace_back(this->box.getOffset(
        newLinkPositions[newLink2] - newLinkPositions[newLink1]));
      linkIndicesOfStrand.push_back({ static_cast<unsigned long>(newLink1),
                                      static_cast<unsigned long>(newLink2) });
      springIndicesOfStrand.push_back({ springFromVertexIdx.size() - 1 });
      currentStrandIdx += 1;
      assert(springIndicesOfStrand.size() == currentStrandIdx);
    }
  }

  // the difficult part is done, now
  // we can fill the network data structure
  assert(springFromVertexIdx.size() == springToVertexIdx.size());
  assert(springContourLength.size() == springFromVertexIdx.size());
  assert(springFromVertexIdx.size() == springBoxOffsets.size());
  this->initialConfig.nrOfSprings = springFromVertexIdx.size();
  this->initialConfig.nrOfStrands = springIndicesOfStrand.size();
  this->initialConfig.springIndicesOfStrand = springIndicesOfStrand;
  this->initialConfig.linkIndicesOfStrand = linkIndicesOfStrand;

  // start with springs
  this->initialConfig.springIndexA =
    Eigen::ArrayXi(this->initialConfig.nrOfSprings);
  this->initialConfig.springIndexB =
    Eigen::ArrayXi(this->initialConfig.nrOfSprings);
  this->initialConfig.springCoordinateIndexA =
    Eigen::ArrayXi(this->initialConfig.nrOfSprings * 3);
  this->initialConfig.springCoordinateIndexB =
    Eigen::ArrayXi(this->initialConfig.nrOfSprings * 3);
  this->initialConfig.springContourLength =
    Eigen::VectorXd(this->initialConfig.nrOfSprings);
  this->initialConfig.springBoxOffset =
    Eigen::VectorXd::Zero(this->initialConfig.nrOfSprings * 3);
  this->initialConfig.springIsEntanglement =
    Eigen::ArrayXb::Constant(this->initialConfig.nrOfSprings, false);
  if (entanglementsAsSprings) {
    this->initialConfig.springIsEntanglement.tail(pairsOfAtoms.size()) = true;
  }

  // then allocate the data of the links
  assert(nextLinkIdx < u.getNrOfAtoms());
  this->initialConfig.nrOfLinks = nextLinkIdx;
  this->initialConfig.linkIsEntanglement =
    Eigen::ArrayXb::Constant(nextLinkIdx, false);
  this->initialConfig.oldAtomIds = Eigen::ArrayXi(nextLinkIdx);
  this->initialConfig.oldAtomTypes = Eigen::ArrayXi(nextLinkIdx);
  assert(newLinkIdxToOldVertexIdx.size() == nextLinkIdx);
  assert(newLinkPositions.size() == nextLinkIdx);

  // convert the springs
  for (size_t i = 0; i < this->initialConfig.nrOfSprings; ++i) {
    assert(APPROX_WITHIN(
      oldVertexIdToNewLinkId[springFromVertexIdx[i]], 0, nextLinkIdx - 1, 0.2));
    assert(APPROX_WITHIN(
      oldVertexIdToNewLinkId[springToVertexIdx[i]], 0, nextLinkIdx - 1, 0.2));
    this->initialConfig.springIndexA(i) =
      oldVertexIdToNewLinkId[springFromVertexIdx[i]];
    this->initialConfig.springIndexB(i) =
      oldVertexIdToNewLinkId[springToVertexIdx[i]];
    this->initialConfig.springContourLength(i) = springContourLength[i];
    for (size_t dir = 0; dir < 3; ++dir) {
      this->initialConfig.springCoordinateIndexA(3 * i + dir) =
        (this->initialConfig.springIndexA(i) * 3 + dir);
      this->initialConfig.springCoordinateIndexB(3 * i + dir) =
        (this->initialConfig.springIndexB(i) * 3 + dir);
    }

    this->initialConfig.springBoxOffset.segment(3 * i, 3) = springBoxOffsets[i];
  }

  // link between springs and strands
  this->initialConfig.strandIndexOfSpring =
    Eigen::ArrayXi::Constant(this->initialConfig.nrOfSprings, -1);
  for (size_t strandIdx = 0; strandIdx < springIndicesOfStrand.size();
       ++strandIdx) {
    std::vector<size_t>& springIndices = springIndicesOfStrand[strandIdx];
    for (size_t springIdx : springIndices) {
      this->initialConfig.strandIndexOfSpring(springIdx) = strandIdx;
    }
  }

  // link between strands and links
  this->initialConfig.strandIndicesOfLink =
    pylimer_tools::utils::initializeWithValue(this->initialConfig.nrOfLinks,
                                              std::vector<size_t>());
  for (size_t strandIdx = 0; strandIdx < linkIndicesOfStrand.size();
       ++strandIdx) {
    for (size_t linkIdx : linkIndicesOfStrand[strandIdx]) {
      pylimer_tools::utils::addIfNotContained(
        this->initialConfig.strandIndicesOfLink[linkIdx], strandIdx);
    }
  }

  // some more properties of the links
  this->initialConfig.coordinates = Eigen::VectorXd::Zero(nextLinkIdx * 3);

  for (size_t linkIdx = 0; linkIdx < nextLinkIdx; ++linkIdx) {
    this->initialConfig.oldAtomIds(linkIdx) =
      this->universe.getPropertyValue<long int>(
        "id", newLinkIdxToOldVertexIdx[linkIdx]);
    this->initialConfig.oldAtomTypes(linkIdx) =
      this->universe.getPropertyValue<int>("type",
                                           newLinkIdxToOldVertexIdx[linkIdx]);
    this->initialConfig.coordinates.segment(3 * linkIdx, 3) =
      newLinkPositions[linkIdx];
    this->initialConfig.linkIsEntanglement(linkIdx) =
      pairOfAtom[newLinkIdxToOldVertexIdx[linkIdx]] >= 0;
  }

  this->initialConfig.nrOfNodes =
    (this->initialConfig.linkIsEntanglement == false).count();

  if (!entanglementsAsSprings) {
    size_t contourLengthSum = this->initialConfig.springContourLength.sum();
    assert(contourLengthSum == this->universe.getNrOfBonds());
    assert(this->initialConfig.linkIsEntanglement.count() ==
           entanglements.pairsOfAtoms.size());
    assert(this->initialConfig.springIsEntanglement.count() == 0);
    assert(contourLengthSum - this->initialConfig.nrOfSprings +
             2 * this->initialConfig.nrOfLinks -
             this->initialConfig.nrOfNodes ==
           this->universe.getNrOfAtoms());
  } else {
    assert(this->initialConfig.springContourLength.sum() ==
           this->universe.getNrOfBonds() +
             this->initialConfig.springIsEntanglement.count());
  }

  this->completeInitialization();
}

MEHPForceBalance2::MEHPForceBalance2(
  const pylimer_tools::entities::Universe& universe,
  const Network& net1,
  const bool is2D)
{
  this->universe = universe;
  this->is2D = is2D;

  pylimer_tools::sim::mehp::ForceBalance2Network net2 = {};

  // Copy box properties
  for (size_t dir = 0; dir < 3; ++dir) {
    net2.L[dir] = net1.L[dir];
    net2.boxHalfs[dir] = 0.5 * net1.L[dir];
  }

  // Basic network properties
  net2.nrOfLinks = net1.nrOfNodes; // In Network, only nodes exist
  net2.nrOfNodes = net1.nrOfNodes;
  net2.nrOfSprings = net1.nrOfSprings;
  net2.nrOfStrands = net1.nrOfSprings; // Each spring becomes a strand

  // Copy coordinates and spring properties
  net2.coordinates = net1.coordinates;
  net2.springContourLength = net1.springsContourLength;

  // Convert spring connectivity
  net2.springIndexA = net1.springIndexA;
  net2.springIndexB = net1.springIndexB;

  // Set up coordinate indices
  net2.springCoordinateIndexA = Eigen::ArrayXi(net2.nrOfSprings * 3);
  net2.springCoordinateIndexB = Eigen::ArrayXi(net2.nrOfSprings * 3);
  for (size_t i = 0; i < net2.nrOfSprings; ++i) {
    for (size_t dir = 0; dir < 3; ++dir) {
      net2.springCoordinateIndexA[3 * i + dir] = 3 * net2.springIndexA[i] + dir;
      net2.springCoordinateIndexB[3 * i + dir] = 3 * net2.springIndexB[i] + dir;
    }
  }

  // Copy box offsets
  net2.springBoxOffset = net1.springBoxOffset;

  // Initialize strand-spring relationships
  net2.springIndicesOfStrand.resize(net2.nrOfStrands);
  net2.linkIndicesOfStrand.resize(net2.nrOfStrands);
  net2.strandIndexOfSpring = Eigen::ArrayXi(net2.nrOfSprings);

  for (size_t i = 0; i < net2.nrOfSprings; ++i) {
    // Each spring becomes its own strand
    net2.springIndicesOfStrand[i] = std::vector<size_t>{ i };
    net2.linkIndicesOfStrand[i] =
      std::vector<size_t>{ static_cast<size_t>(net2.springIndexA[i]),
                           static_cast<size_t>(net2.springIndexB[i]) };
    net2.strandIndexOfSpring[i] = static_cast<int>(i);
  }

  // Initialize strand-link relationships
  net2.strandIndicesOfLink.resize(net2.nrOfLinks);
  for (size_t i = 0; i < net2.nrOfSprings; ++i) {
    net2.strandIndicesOfLink[net2.springIndexA[i]].push_back(i);
    net2.strandIndicesOfLink[net2.springIndexB[i]].push_back(i);
  }

  // Initialize link and spring type arrays
  net2.linkIsEntanglement = Eigen::ArrayXb::Constant(net2.nrOfLinks, false);
  net2.springIsEntanglement = Eigen::ArrayXb::Constant(net2.nrOfSprings, false);

  // Copy atom information
  net2.oldAtomIds = net1.oldAtomIds;
  net2.oldAtomTypes = Eigen::ArrayXi::Zero(net2.nrOfLinks);

  // Set default atom types (since Network doesn't have this information)
  for (size_t i = 0; i < net2.nrOfLinks; ++i) {
    net2.oldAtomTypes[i] = this->crossLinkerType;
  }

  this->configAssumeNetworkIsComplete(net1.assumeComplete);
  this->initialConfig = net2;
  this->box = pylimer_tools::entities::Box(net1.L[0], net1.L[1], net1.L[2]);
  this->completeInitialization();
}

void
MEHPForceBalance2::completeInitialization()
{
  this->currentDisplacements =
    Eigen::VectorXd::Zero(3 * this->initialConfig.nrOfLinks);
  this->defaultBondLength = universe.computeMeanBondLength();
  RUNTIME_EXP_IFN(this->validateNetwork(), "Invalid internal state");
}

/**
 * FORCE RELAXATION
 */
/**
 * @brief Run force relaxation simulation with specified parameters
 *
 * @param simplificationMode the mode for structure simplification
 * @param inactiveRemovalCutoff the cutoff for removing inactive links
 * @param solverChoice the choice of solver for the linear equation system
 * @param residualReduction the target residual reduction for convergence
 * @param maxIterations the maximum number of iterations to perform
 * @param shouldInterrupt function to check if simulation should be interrupted
 * @param cleanupInterrupt function to call when simulation is interrupted
 */
void
MEHPForceBalance2::runForceRelaxation(
  const StructureSimplificationMode simplificationMode,
  const double inactiveRemovalCutoff,
  const SLESolver solverChoice,
  const double residualReduction,
  const int maxIterations,
  const std::function<bool()>& shouldInterrupt,
  const std::function<void()>& cleanupInterrupt)
{
  RUNTIME_EXP_IFN(this->validateNetwork(),
                  "Invalid internal state of the network.");
  // INVALIDARG_EXP_IFN(
  //   shouldRemoveInactiveCrosslinks == false &&
  //     remove2functionalCrosslinkers == true,
  //   "Removing 2-functional crosslinkers only makes sense when inactive "
  //   "crosslinkers may be removed too, during the procedure.");
  this->simulationHasRun = true;

  INVALIDARG_EXP_IFN(
    inactiveRemovalCutoff > 0.0 ||
      simplificationMode == StructureSimplificationMode::NO_SIMPLIFICATION,
    "Removal cut-off must be positive when simplification is enabled.");

  if (this->getNrOfStrands() == 0) {
    return;
  }

  if (simplificationMode != StructureSimplificationMode::NO_SIMPLIFICATION) {
    this->configAssumeNetworkIsComplete(false);
  }

  /* array allocation */
  Eigen::VectorXd oneOverSpringPartitions =
    this->assembleOneOverSpringPartition(this->initialConfig);
  const double initialResidual = this->getDisplacementResidualNormFor(
    this->initialConfig, this->currentDisplacements, oneOverSpringPartitions);
  std::cout << "Starting force balance procedure with a "
            << SLESolverNames[solverChoice] << " (" << solverChoice << ")"
            << " solver. Initial residual is " << initialResidual << std::endl;
  std::cout << "Simplification mode is " << simplificationMode << " ("
            << StructureSimplificationModeNames[simplificationMode] << ")"
            << std::endl;
  double currentResidual = initialResidual;
  size_t iterationsDone = 0;

  // best to start with zero as initial guess,
  // this leads to better results for free/melt chains and structures
  Eigen::VectorXd initialSolution =
    Eigen::VectorXd::Zero(3 * this->initialConfig.nrOfLinks);

  this->prepareAllOutputs();

  // actual loop
  bool wasInterrupted = false;
  size_t nRemovedInIteration;
  do {
    nRemovedInIteration = 0;
    std::vector<Eigen::Triplet<double>> triplets;
    // diagonal + the lower of the two components of each spring
    triplets.reserve(this->initialConfig.nrOfLinks * 3 +
                     this->initialConfig.nrOfSprings * 3 * 2);
    Eigen::VectorXd constants =
      Eigen::VectorXd::Zero(this->initialConfig.nrOfLinks * 3);
    // it's a bit more efficient to sum the diagonal ourselves
    Eigen::VectorXd diagonal =
      Eigen::VectorXd::Zero(this->initialConfig.nrOfLinks * 3);

    for (size_t springIdx = 0; springIdx < this->initialConfig.nrOfSprings;
         ++springIdx) {
      if (this->initialConfig.springIndexA[springIdx] ==
          this->initialConfig.springIndexB[springIdx]) {
        continue;
      }
      double oneOverContourLengthFraction =
        1.0 / (this->initialConfig.springContourLength[springIdx]);
      double multiplier =
        this->kappa *
        oneOverContourLengthFraction; // oneOverSpringPartitions(springIdx
      // * 3);
      assert(this->initialConfig.springIndexA[springIdx] <
             this->initialConfig.nrOfLinks);
      assert(this->initialConfig.springIndexB[springIdx] <
             this->initialConfig.nrOfLinks);
      // triplets will be summed up -> we can use the same indices multiple
      // times
      for (size_t dir = 0; dir < 3; ++dir) {
        // store only the lower part
        // if (this->initialConfig.springPartIndexA[springIdx] <
        //     this->initialConfig.springPartIndexB[springIdx]) {
        triplets.push_back(Eigen::Triplet<double>(
          this->initialConfig.springIndexA[springIdx] * 3 + dir,
          this->initialConfig.springIndexB[springIdx] * 3 + dir,
          1. * multiplier));
        // } else {
        triplets.push_back(Eigen::Triplet<double>(
          this->initialConfig.springIndexB[springIdx] * 3 + dir,
          this->initialConfig.springIndexA[springIdx] * 3 + dir,
          1. * multiplier));
        // }
      }
      diagonal.segment(3 * this->initialConfig.springIndexA[springIdx], 3) -=
        Eigen::Vector3d::Constant(multiplier);
      diagonal.segment(3 * this->initialConfig.springIndexB[springIdx], 3) -=
        Eigen::Vector3d::Constant(multiplier);

      // the constants, b in Ax = b
      constants.segment(3 * this->initialConfig.springIndexA[springIdx], 3) -=
        this->initialConfig.springBoxOffset.segment(3 * springIdx, 3) *
        multiplier;
      constants.segment(3 * this->initialConfig.springIndexB[springIdx], 3) +=
        this->initialConfig.springBoxOffset.segment(3 * springIdx, 3) *
        multiplier;
    }

    for (size_t linkIdx = 0; linkIdx < this->initialConfig.nrOfLinks;
         ++linkIdx) {
      for (size_t dir = 0; dir < 3; ++dir) {
        triplets.push_back(Eigen::Triplet<double>(
          linkIdx * 3 + dir, linkIdx * 3 + dir, diagonal(linkIdx * 3 + dir)));
      }
    }

    if (triplets.empty()) {
      break;
    }

    Eigen::SparseMatrix<double> sysMatrix(this->initialConfig.nrOfLinks * 3,
                                          this->initialConfig.nrOfLinks * 3);
    sysMatrix.setFromTriplets(triplets.begin(), triplets.end());
    sysMatrix.makeCompressed();

#ifndef NDEBUG
    // verify that the matrix fulfills the requirements of the gradient descent
    bool isSelfAdjoint = Eigen::isSelfAdjoint(sysMatrix);
    assert(isSelfAdjoint);
    // test whether the system matrix is positive definite (it's not
    // guaranteed!)
    // Eigen::VectorXd x =
    // Eigen::VectorXd::Random(sysMatrix.rows());
    // double result = x.transpose() * sysMatrix * x;
    // std::cout << "Quadratic test result: " << result << std::endl;
    // test whether the system matrix is full rank
    // std::cout << "Rank of matrix: " <<
    // Eigen::FullPivLU<Eigen::MatrixXd>(sysMatrix).rank()
    // << " for " << sysMatrix.rows() << "x" << sysMatrix.cols() << "matrix" <<
    // std::endl; Eigen::SparseLU<Eigen::SparseMatrix<double>> luSolver;
    // luSolver.compute(sysMatrix);
    // std::cout << "LU decomposition info: " << luSolver.info() << " had
    // determinant " << luSolver.determinant() << std::endl;
    // ...

    // Eigen::saveSparseMatrix(
    //   sysMatrix, "sysMatrix_" + std::to_string(constants.size()) + ".mtx");
    // Eigen::saveDenseVector(
    //   constants, "constants_" + std::to_string(constants.size()) + ".txt");
#endif

    Eigen::VectorXd finalCoordinates;

#define SOLVER_FAILED(action)                                                  \
  this->exitReason = ExitReason::FAILURE;                                      \
  RUNTIME_EXP(std::to_string(action) + " with solver " +                       \
              std::to_string(solverChoice) + " (" +                            \
              std::to_string(SLESolverNames[solverChoice]) +                   \
              "). Solver info: " + std::to_string(solver.info()));

#define COMPUTE_SOLVE(solver)                                                  \
  if (constants.isZero(1e-9)) {                                                \
    finalCoordinates = Eigen::VectorXd::Zero(constants.size());                \
  } else {                                                                     \
    solver.compute(sysMatrix);                                                 \
    if (solver.info() != Eigen::Success) {                                     \
      SOLVER_FAILED("System matrix computation failed");                       \
    }                                                                          \
    finalCoordinates = solver.solve(constants);                                \
    if (solver.info() != Eigen::Success) {                                     \
      SOLVER_FAILED("System could not be solved");                             \
    }                                                                          \
  }

#define ADD_ITERATIONS(solver) iterationsDone += solver.iterations();

#define SOLVE_ITERATIVE(solver)                                                \
  solver.setTolerance(residualReduction);                                      \
  solver.setMaxIterations(maxIterations);                                      \
  COMPUTE_SOLVE(solver);                                                       \
  ADD_ITERATIONS(solver);                                                      \
  this->exitReason = ExitReason::X_TOLERANCE;                                  \
  break;

#define SOLVE_DIRECT(solver)                                                   \
  COMPUTE_SOLVE(solver);                                                       \
  this->exitReason = ExitReason::X_TOLERANCE;                                  \
  break;

    // Common parameters for gradient descent solvers
    constexpr double gradientDescentLearningRate = 0.01;

    // Define common callback for gradient descent solvers to handle output and
    // interrupts
    auto createGradientDescentCallback = [&](size_t& iterationsDone) {
      return [&](int currentIteration,
                 const Eigen::VectorXd& currentSolution) -> bool {
        // Temporarily update the current displacements for handleOutput
        Eigen::VectorXd oldDisplacements = this->currentDisplacements;
        this->currentDisplacements =
          currentSolution - this->initialConfig.coordinates;

        this->handleOutput(iterationsDone + currentIteration);
        bool shouldStop = shouldInterrupt();

        // Restore the previous displacements if not stopping
        if (!shouldStop) {
          this->currentDisplacements = oldDisplacements;
        }

        return shouldStop;
      };
    };

    // Helper function to handle gradient descent solver results
    auto handleGradientDescentResult = [&](int solverIterations) {
      if (solverIterations >= maxIterations) {
        this->exitReason = ExitReason::MAX_STEPS;
      } else {
        this->exitReason = ExitReason::X_TOLERANCE;
      }
      iterationsDone += solverIterations;
    };

    switch (solverChoice) {
      // iterative solvers
      case SLESolver::CONJUGATE_GRADIENT:
      case SLESolver::CONJUGATE_GRADIENT_DIAGONALIZED: {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                 Eigen::Lower | Eigen::Upper,
                                 Eigen::DiagonalPreconditioner<double>>
          solver;
        SOLVE_ITERATIVE(solver);
      }
      case SLESolver::CONJUGATE_GRADIENT_IDENTITY: {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                 Eigen::Lower | Eigen::Upper,
                                 Eigen::IdentityPreconditioner>
          solver;
        SOLVE_ITERATIVE(solver);
      }
      case SLESolver::CONJUGATE_GRADIENT_INCOMPLETE_CHOLESKY: {
        Eigen::ConjugateGradient<
          Eigen::SparseMatrix<double>,
          Eigen::Lower | Eigen::Upper,
          Eigen::IncompleteCholesky<double, Eigen::Lower | Eigen::Upper>>
          solver;
        SOLVE_ITERATIVE(solver);
      }
      case SLESolver::LEAST_SQUARES_CONJUGATE_GRADIENT:
      case SLESolver::LEAST_SQUARES_CONJUGATE_GRADIENT_DIAGONALIZED: {
        Eigen::LeastSquaresConjugateGradient<
          Eigen::SparseMatrix<double>,
          Eigen::LeastSquareDiagonalPreconditioner<double>>
          solver;
        SOLVE_ITERATIVE(solver);
      }
      case SLESolver::LEAST_SQUARES_CONJUGATE_GRADIENT_IDENTITY: {
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>,
                                             Eigen::IdentityPreconditioner>
          solver;
        SOLVE_ITERATIVE(solver);
      }
      case SLESolver::BICGSTAB:
      case SLESolver::BICGSTAB_DIAGONALIZED: {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,
                        Eigen::DiagonalPreconditioner<double>>
          solver;
        SOLVE_ITERATIVE(solver);
      }
      case SLESolver::BICGSTAB_IDENTITY: {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,
                        Eigen::IdentityPreconditioner>
          solver;
        SOLVE_ITERATIVE(solver);
      }
      case SLESolver::BICGSTAB_INCOMPLETE_LU: {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>,
                        Eigen::IncompleteLUT<double>>
          solver;
        SOLVE_ITERATIVE(solver);
      }
      case SLESolver::GRADIENT_DESCENT: {
        int solverIterations = 0;
        auto callback = createGradientDescentCallback(iterationsDone);

        finalCoordinates = Eigen::gradientDescent(sysMatrix,
                                                  constants,
                                                  gradientDescentLearningRate,
                                                  residualReduction,
                                                  maxIterations,
                                                  solverIterations,
                                                  initialSolution,
                                                  // initial solution
                                                  initialResidual,
                                                  // initial residual
                                                  callback);

        handleGradientDescentResult(solverIterations);
        break;
      }
      case SLESolver::DEFAULT:
      case SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT:
      case SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG: {
        int solverIterations = 0;
        auto callback = createGradientDescentCallback(iterationsDone);

        finalCoordinates = Eigen::gradientDescentBarzilaiBorwein(
          sysMatrix,
          constants,
          gradientDescentLearningRate,
          residualReduction,
          maxIterations,
          solverIterations,
          solverChoice == SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT,
          initialSolution,
          // initial solution
          initialResidual,
          // initial residual
          callback);

        handleGradientDescentResult(solverIterations);
        break;
      }
      case SLESolver::GRADIENT_DESCENT_BARZILAI_BORWEIN_MOMENTUM: {
        int solverIterations = 0;
        auto callback = createGradientDescentCallback(iterationsDone);

        finalCoordinates = Eigen::gradientDescentHeavyBallBarzilaiBorwein(
          sysMatrix,
          constants,
          gradientDescentLearningRate,
          residualReduction,
          maxIterations,
          solverIterations,
          initialSolution,
          // initial solution
          initialResidual,
          // initial residual
          callback);

        handleGradientDescentResult(solverIterations);
        break;
      }
      // direct solvers
      case SLESolver::SIMPLICIAL_LLT: {
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        SOLVE_DIRECT(solver);
      }
      case SLESolver::SIMPLICIAL_LDLT: {
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        SOLVE_DIRECT(solver);
      }
      case SLESolver::SPARSE_LU: {
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        SOLVE_DIRECT(solver);
      }
      case SLESolver::SPARSE_QR: {
        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
          solver;
        SOLVE_DIRECT(solver);
      }
      default:
        throw std::runtime_error("This solver is not implemented.");
    }

#undef SOLVER_FAILED
#undef COMPUTE_SOLVE
#undef ADD_ITERATIONS
#undef SOLVE_ITERATIVE
#undef SOLVE_DIRECT

    this->currentDisplacements =
      finalCoordinates - this->initialConfig.coordinates;

    currentResidual = this->getDisplacementResidualNormFor(
      this->initialConfig, this->currentDisplacements);

    iterationsDone += 1;
    this->breakTooLongStrands(this->initialConfig, this->currentDisplacements);
    size_t nRemovedThisLoop = 0;

    do {
#ifndef NDEBUG
      RUNTIME_EXP_IFN(this->validateNetwork(), "Invalid internal state");
#endif
      nRemovedThisLoop = 0;
      if (simplificationMode == StructureSimplificationMode::INACTIVE_ONLY ||
          simplificationMode ==
            StructureSimplificationMode::INACTIVE_THEN_X2F) {
#ifdef DEBUG_REMOVAL
        std::cout << "Checking and possibly removing inactive crosslinks"
                  << std::endl;
#endif

        Eigen::ArrayXb springIsToDelete =
          this->markInactiveSpringsToDelete(this->initialConfig,
                                            this->currentDisplacements,
                                            inactiveRemovalCutoff);
        size_t toDeleteCount = springIsToDelete.count();
        if (toDeleteCount > 0) {
          std::vector<size_t> springIndicesToDelete;
          springIndicesToDelete.reserve(toDeleteCount);
          for (size_t i = springIsToDelete.size(); i > 0; --i) {
            if (springIsToDelete[i - 1]) {
              springIndicesToDelete.push_back(i - 1);
            }
          }
          // std::ranges::sort(springIndicesToDelete, std::greater<>());
          this->removeSprings(this->initialConfig,
                              this->currentDisplacements,
                              springIndicesToDelete);
          nRemovedThisLoop += toDeleteCount;
        }
      }
      if (simplificationMode == StructureSimplificationMode::X2F_ONLY ||
          simplificationMode ==
            StructureSimplificationMode::INACTIVE_THEN_X2F) {
#ifdef DEBUG_REMOVAL
        std::cout << "Checking and possibly removing crosslinks with f = 2"
                  << std::endl;
#endif
        nRemovedThisLoop += this->unlinkBifunctionalLinks(
          this->initialConfig, this->currentDisplacements);
      }
      if (simplificationMode ==
          StructureSimplificationMode::X1F_X2F_THEN_INACTIVE) {
        RUNTIME_EXP("This mode is not implemented.");
      }

      // cleanup some things
      if (simplificationMode !=
          StructureSimplificationMode::NO_SIMPLIFICATION) {
        this->validateNetwork(this->initialConfig, this->currentDisplacements);
        oneOverSpringPartitions =
          this->assembleOneOverSpringPartition(this->initialConfig);
      }

      nRemovedInIteration += nRemovedThisLoop;
    } while (nRemovedThisLoop > 0);

    // after running once, we can update the initial solution
    initialSolution =
      this->currentDisplacements + this->initialConfig.coordinates;

    // after removal, the residual changed, might even have increased
    // beyond initial
    // -> reset previous and current to prevent change to iterative
    // displacement
    if (nRemovedInIteration > 0) {
      oneOverSpringPartitions =
        this->assembleOneOverSpringPartition(this->initialConfig);
    }

    this->handleOutput(iterationsDone);

    if (shouldInterrupt()) {
      wasInterrupted = true;
      break;
    }
  } while (this->initialConfig.nrOfStrands > 0 && nRemovedInIteration > 0);

  // finish up
  this->closeAllOutputs();

  // some last output and additions
  this->nrOfStepsDone += iterationsDone;
  std::cout << iterationsDone << " steps done. "
            << "Current residual: " << currentResidual << ". "
            << "Initial residual: " << initialResidual << ". " << std::endl;

  assert(this->currentDisplacements.size() ==
         3 * this->initialConfig.nrOfLinks);
  RUNTIME_EXP_IFN(this->validateNetwork(), "Invalid internal state");
  if (wasInterrupted) {
    this->exitReason = ExitReason::INTERRUPT;
    cleanupInterrupt();
  }
}

/**
 * @brief Compute the displacement residual norm for the current
 * configuration
 *
 * @return double
 */
double
MEHPForceBalance2::getDisplacementResidualNorm() const
{
  const Eigen::VectorXd oneOverSpringPartitions =
    this->assembleOneOverSpringPartition(this->initialConfig);
  const Eigen::VectorXd displacements = this->currentDisplacements;
  return this->getDisplacementResidualNormFor(
    this->initialConfig, displacements, oneOverSpringPartitions);
}

/**
 * @brief Compute the displacement residual norm for a specific
 * configuration
 *
 * @param net the network to analyze
 * @param u the current displacements of the links
 * @return double the displacement residual norm
 */
double
MEHPForceBalance2::getDisplacementResidualNormFor(
  const ForceBalance2Network& net,
  const Eigen::VectorXd& u) const
{
  const Eigen::ArrayXd oneOverSpringPartitions =
    this->assembleOneOverSpringPartition(net).array();

  const Eigen::ArrayXd loopPartialSpringEliminator =
    (net.springCoordinateIndexA != net.springCoordinateIndexB).cast<double>();
  Eigen::ArrayXd forces = Eigen::ArrayXd::Zero(3 * net.nrOfLinks);
  const Eigen::ArrayXd distances =
    this->evaluateSpringVectors(net, u, this->is2D).array();
  forces(net.springCoordinateIndexA) +=
    (this->kappa * oneOverSpringPartitions * distances *
     loopPartialSpringEliminator);
  forces(net.springCoordinateIndexB) -=
    (this->kappa * oneOverSpringPartitions * distances *
     loopPartialSpringEliminator);

  // #ifndef NDEBUG
  //       Eigen::VectorXi debugNrSpringsVisited =
  //         Eigen::VectorXi::Zero(net.nrOfPartialSprings);
  //       for (size_t i = 0; i < net.nrOfLinks; ++i) {
  //         Eigen::Array3d forces2 =
  //           this
  //             ->evaluateForceOnLink(i,
  //                                   net,
  //                                   u,
  //
  //                                   debugNrSpringsVisited,
  //                                   )
  //             .array();
  //         double squareN1 = forces.segment(3 * i,
  //         3).matrix().squaredNorm(); double squareN2 =
  //         forces2.matrix().squaredNorm(); if (i == 129) {
  //           this->debugAtomVicinity(net.oldAtomIds[i]);
  //         }
  //         assert(APPROX_EQUAL(squareN1, squareN2, 1e-9));
  //         // assert(pylimer_tools::utils::vector_approx_equal(
  //         //   forces.segment(3 * i, 3), forces2, 1e-9, true));
  //       }
  //       assert((debugNrSpringsVisited.array() == 2).all());
  // #endif

  return forces.matrix().squaredNorm();
}

/**
 * @brief Compute the displacement residual norm for a specific
 * configuration with spring partitions
 *
 * @param net the network to analyze
 * @param u the current displacements of the links
 * @param oneOverSpringPartitions the reciprocal of spring partition values
 * @return double the displacement residual norm
 */
double
MEHPForceBalance2::getDisplacementResidualNormFor(
  const ForceBalance2Network& net,
  const Eigen::VectorXd& u,
  const Eigen::VectorXd& oneOverSpringPartitions) const
{
  Eigen::VectorXd displacedCoords = net.coordinates + u;

  Eigen::VectorXd relevantPartialDistances =
    (displacedCoords(net.springCoordinateIndexB) -
     displacedCoords(net.springCoordinateIndexA)) +
    net.springBoxOffset;

  if (this->is2D) {
    for (size_t i = 2; i < relevantPartialDistances.size(); i += 3) {
      relevantPartialDistances[i] = 0.;
    }
  }

#ifndef NDEBUG
  for (size_t i = 0; i < net.nrOfSprings; ++i) {
    Eigen::Vector3d dist = this->evaluateSpringVector(net, u, i);
    Eigen::Vector3d comparison = relevantPartialDistances.segment(3 * i, 3);
    assert(pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
      dist, comparison, 1e-9));
  }
#endif

  assert(relevantPartialDistances.size() == oneOverSpringPartitions.size());
  const Eigen::VectorXd partialDistancesOverSpringPartitions =
    (relevantPartialDistances.array() * oneOverSpringPartitions.array())
      .matrix();

  Eigen::VectorXd overallForces = Eigen::VectorXd::Zero(3 * net.nrOfLinks);
  overallForces(net.springCoordinateIndexB) -=
    partialDistancesOverSpringPartitions;
  overallForces(net.springCoordinateIndexA) +=
    partialDistancesOverSpringPartitions;

  return overallForces.squaredNorm();
}

/**
 * @brief Get the entanglement springs connected to a specific strand
 *
 * @param net the network to analyze
 * @param strandIdx the index of the strand to analyze
 * @return std::vector<size_t> indices of entanglement springs at the strand
 */
std::vector<size_t>
MEHPForceBalance2::getEntanglementSpringsAtStrand(
  const ForceBalance2Network& net,
  const size_t strandIdx) const
{
  std::vector<size_t> result;
  for (size_t linkI = 1; linkI < net.linkIndicesOfStrand[strandIdx].size() - 1;
       ++linkI) {
    const size_t linkIdx = net.linkIndicesOfStrand[strandIdx][linkI];
    const std::vector<size_t> springs =
      this->getPartialSpringIndicesOfLink(net, linkIdx);
    for (const size_t springIdx : springs) {
      if (net.springIsEntanglement[springIdx]) {
        result.push_back(springIdx);
      }
    }
  }
  return result;
}

/**
 * @brief Get the indices of neighboring links for a given link
 *
 * @param net the network to analyze
 * @param linkIdx the index of the link to find neighbors for
 * @return std::vector<size_t> indices of neighboring links
 */
std::vector<size_t>
MEHPForceBalance2::getNeighbourLinkIndices(const ForceBalance2Network& net,
                                           const size_t linkIdx)
{
  INVALIDINDEX_EXP_IFN(linkIdx < net.nrOfLinks,
                       "Link index " + std::to_string(linkIdx) +
                         " out of range, only " +
                         std::to_string(net.nrOfLinks) + " links are present.");
  std::vector<size_t> results;
  results.reserve(4);
  for (const size_t springIdx : net.strandIndicesOfLink[linkIdx]) {
    for (const size_t partialSpringIdx : net.springIndicesOfStrand[springIdx]) {
      if (net.springIndexA[partialSpringIdx] == linkIdx) {
        results.push_back(net.springIndexB[partialSpringIdx]);
      } //
      else if (net.springIndexB[partialSpringIdx] == linkIdx) {
        results.push_back(net.springIndexA[partialSpringIdx]);
      }
    }
  }
  return results;
}

/**
 * @brief Translate the spring partition vector to its 3*size
 * version with one over spring contour length values
 *
 * @param net the network to analyze
 * @return Eigen::VectorXd vector of one over spring contour length values
 */
Eigen::VectorXd
MEHPForceBalance2::assembleOneOverSpringPartition(
  const ForceBalance2Network& net)
{
  Eigen::VectorXd oneOverSpringPartitions =
    Eigen::VectorXd(3 * net.nrOfSprings);

  Eigen::ArrayXd primaryLoopCorrectionMultiplier =
    (net.springIndexA != net.springIndexB)
      .cast<double>(); // 0.0 for equal = primary loop, 1.0 otherwise

  for (size_t i = 0; i < net.nrOfSprings; ++i) {
    const double valueToSet = 1.0 / (net.springContourLength[i]);

    // if (springPartitions0[i] < 1e-9) {
    //   std::cout << "Got close call for partial spring " << i <<
    //   std::endl;
    // }
    oneOverSpringPartitions.segment(3 * i, 3) = Eigen::Vector3d::Constant(
      valueToSet * primaryLoopCorrectionMultiplier[i]);
  }

  return oneOverSpringPartitions;
}

/**
 * @brief Remove springs that exert a stress higher than
 * `this->springBreakingLength`
 *
 * @param net the network to modify
 * @param displacements the current displacements of the links
 * @return size_t the number of strands broken
 */
size_t
MEHPForceBalance2::breakTooLongStrands(ForceBalance2Network& net,
                                       Eigen::VectorXd& displacements) const
{
  if (this->springBreakingLength <= 0.) {
    return 0;
  }

  size_t numBroken = 0;

  std::vector<size_t> springIndicesToDelete;

  // iterate the springs, determine their distance, and determine if it
  // exceeds the breaking force
  for (size_t strandIdx = 0; strandIdx < net.nrOfStrands; ++strandIdx) {
    for (const size_t springIdx : net.springIndicesOfStrand[strandIdx]) {
      const double len =
        this->getWeightedSpringLength(net, displacements, springIdx);

      if (len > this->springBreakingLength) {
        // break this strand
        numBroken += 1;
        for (size_t i : net.springIndicesOfStrand[strandIdx]) {
          springIndicesToDelete.push_back(i);
        }
        pylimer_tools::utils::append(
          springIndicesToDelete,
          this->getEntanglementSpringsAtStrand(net, strandIdx));

        break;
      }
    }
  }

  std::ranges::sort(springIndicesToDelete, std::greater<>());
  this->removeSprings(net, displacements, springIndicesToDelete);

  return numBroken;
}

/**
 * @brief Decide for each spring if it should be removed,
 * remove them, then remove orphaned links
 *
 * @param net the network
 * @param displacements the current displacements
 * @param inactiveRemovalCutoff the cut-off for the activity tolerance criterion
 * @return the number of links removed
 */
size_t
MEHPForceBalance2::removeInactiveLinks(ForceBalance2Network& net,
                                       Eigen::VectorXd& displacements,
                                       const double inactiveRemovalCutoff) const
{
  Eigen::ArrayXb springIsToDelete = this->markInactiveSpringsToDelete(
    net, displacements, inactiveRemovalCutoff);
  const size_t toDeleteCount = springIsToDelete.count();
  std::vector<size_t> springIndicesToDelete;
  springIndicesToDelete.reserve(toDeleteCount);
  for (size_t i = springIsToDelete.size(); i > 0; --i) {
    if (springIsToDelete[i - 1]) {
      springIndicesToDelete.push_back(i - 1);
    }
  }
  assert(springIndicesToDelete.size() == toDeleteCount);
  // std::ranges::sort(springIndicesToDelete, std::greater<>());
  this->removeSprings(net, displacements, springIndicesToDelete);
  this->unlinkBifunctionalLinks(net, displacements);
  return this->removeZerofunctionalLinks(net, displacements);
};

/**
 * @brief Decide for each spring if it should be removed
 *
 * @param net the network
 * @param displacements the current displacements
 * @param inactiveRemovalCutoff the cut-off for the activity tolerance criterion
 * @return a vector of true/false for each spring, true if it should be
 * removed
 */
Eigen::ArrayXb
MEHPForceBalance2::markInactiveSpringsToDelete(
  const ForceBalance2Network& net,
  const Eigen::VectorXd& displacements,
  const double inactiveRemovalCutoff) const
{
  Eigen::ArrayXb result = Eigen::ArrayXb::Constant(net.nrOfSprings, false);
  Eigen::ArrayXb activeStrands =
    this->findActiveStrands(net, displacements, inactiveRemovalCutoff);

  for (size_t strandIdx = 0; strandIdx < net.nrOfStrands; ++strandIdx) {
    // if one end is inactive, we mark the whole strands as inactive.
    // an exception
    if (!activeStrands[strandIdx]) {
      for (const size_t springIdx : net.springIndicesOfStrand[strandIdx]) {
        result[springIdx] = true;
      }
      // also remove entanglement springs associated with any of the links on
      // the strand
      std::vector<size_t> entanglementSprings =
        this->getEntanglementSpringsAtStrand(net, strandIdx);
      for (const size_t springIdx : entanglementSprings) {
        result[springIdx] = true;
      }
    }
  }

  return result;
}

/**
 * @brief Replace the two springs traversing a two-functional crosslinkers
 * with a single spring
 *
 * Also handles entanglement beads
 *
 * @param net the network to modify
 * @param displacements the current displacements of the links
 * @return size_t the number of bifunctional links removed
 */
size_t
MEHPForceBalance2::unlinkBifunctionalLinks(ForceBalance2Network& net,
                                           Eigen::VectorXd& displacements) const
{
  size_t numRemoved = 0;
  std::vector<size_t> removedSprings = {};
  for (long int linkIdx = net.nrOfLinks - 1; linkIdx >= 0; --linkIdx) {
    std::vector<size_t> springIndices =
      this->getPartialSpringIndicesOfLink(net, linkIdx);
    if (springIndices.size() == 2 &&
        !this->isLoopingSpring(net, springIndices[0]) &&
        !this->isLoopingSpring(net, springIndices[1])) {
      assert(springIndices[0] != springIndices[1]);
      assert(!pylimer_tools::utils::contains<size_t>(removedSprings,
                                                     springIndices[0]));
      assert(!pylimer_tools::utils::contains<size_t>(removedSprings,
                                                     springIndices[1]));
      if (!net.linkIsEntanglement[linkIdx] &&
          net.strandIndexOfSpring[springIndices[0]] ==
            net.strandIndexOfSpring[springIndices[1]]) {
        // annoying primary loop we don't want to handle now.
        continue;
      }
      // merge these two springs,
      // unlink this link
      this->mergeSpringsWithoutRemoval(
        net, displacements, springIndices[0], springIndices[1], linkIdx);
      removedSprings.push_back(springIndices[0]);
      assert(net.strandIndexOfSpring[springIndices[0]] < 0);
      numRemoved += 1;
    }
  }

  std::ranges::sort(removedSprings, std::greater<>());
  this->removeSprings(net, displacements, removedSprings);

#ifndef NDEBUG
  assert(this->validateNetwork(net, displacements));
#endif
  return numRemoved;
}

/**
 * @brief Merge two springs without actually removing any spring from the
 * network data structure
 *
 * @param net the network to modify
 * @param u the current displacements of the links
 * @param removedSpringIdx the index of the spring to be logically removed
 * @param keptSpringIdx the index of the spring to be kept and updated
 * @param linkToReduce the index of the link to be reduced in functionality
 */
void
MEHPForceBalance2::mergeSpringsWithoutRemoval(ForceBalance2Network& net,
                                              const Eigen::VectorXd& u,
                                              const size_t removedSpringIdx,
                                              const size_t keptSpringIdx,
                                              const size_t linkToReduce) const
{
  INVALIDARG_EXP_IFN(removedSpringIdx < net.nrOfSprings,
                     "Removed spring index out of bounds.");
  INVALIDARG_EXP_IFN(keptSpringIdx < net.nrOfSprings,
                     "Kept spring index out of bounds.");
  INVALIDARG_EXP_IFN(keptSpringIdx != removedSpringIdx,
                     "Cannot merge one spring with the same one.");
  INVALIDARG_EXP_IFN(!net.springIsEntanglement[keptSpringIdx], "");
  INVALIDARG_EXP_IFN(!net.springIsEntanglement[removedSpringIdx], "");
  INVALIDARG_EXP_IFN(XOR(net.springIndexA[keptSpringIdx] == linkToReduce,
                         net.springIndexB[keptSpringIdx] == linkToReduce),
                     "");
  INVALIDARG_EXP_IFN(XOR(net.springIndexA[removedSpringIdx] == linkToReduce,
                         net.springIndexB[removedSpringIdx] == linkToReduce),
                     "");

#ifdef DEBUG_REMOVAL
  std::cout << "Merging springs " << removedPartialSpringIdx << " and "
            << keptPartialSpringIdx << " around " << linkToReduce << std::endl;
#endif

  // whether to only re-link the springs, or also update the strands relations
  const bool isStrandMerger = net.strandIndexOfSpring[keptSpringIdx] !=
                              net.strandIndexOfSpring[removedSpringIdx];

  // for validation, to keep the distances correct
  const Eigen::Vector3d totalDistance =
    this->evaluateSpringVectorTo(net, u, keptSpringIdx, linkToReduce) +
    this->evaluateSpringVectorFrom(net, u, removedSpringIdx, linkToReduce);

  // adjust obvious things
  // net.nrOfSprings -= 1;
  // net.nrOfStrands -= isStrandMerger;
  net.springContourLength[keptSpringIdx] +=
    net.springContourLength[removedSpringIdx];

  // more difficult stuff
  const bool keptIsA = net.springIndexA[keptSpringIdx] == linkToReduce;
  const size_t otherEndKept =
    keptIsA ? net.springIndexB[keptSpringIdx] : net.springIndexA[keptSpringIdx];
  const bool removedIsA = net.springIndexA[removedSpringIdx] == linkToReduce;
  const size_t otherEndRemoved = removedIsA
                                   ? net.springIndexB[removedSpringIdx]
                                   : net.springIndexA[removedSpringIdx];
  const size_t removedStrandIdx = net.strandIndexOfSpring[removedSpringIdx];
  assert(pylimer_tools::utils::contains(net.strandIndicesOfLink[linkToReduce],
                                        removedStrandIdx));
  assert(pylimer_tools::utils::contains(
    net.linkIndicesOfStrand[removedStrandIdx], linkToReduce));
  const size_t keptStrandIdx = net.strandIndexOfSpring[keptSpringIdx];
  assert(pylimer_tools::utils::contains(net.strandIndicesOfLink[linkToReduce],
                                        keptStrandIdx));
  assert(pylimer_tools::utils::contains(net.linkIndicesOfStrand[keptStrandIdx],
                                        linkToReduce));
  if (!isStrandMerger) {
    assert(XOR(keptIsA, removedIsA));
    assert(net.linkIndicesOfStrand[removedStrandIdx][0] != linkToReduce);
    assert(net.linkIndicesOfStrand[removedStrandIdx].back() != linkToReduce);
    assert(net.linkIndicesOfStrand[keptStrandIdx][0] != linkToReduce);
    assert(net.linkIndicesOfStrand[keptStrandIdx].back() != linkToReduce);
  }

  net.springIndexA[keptSpringIdx] = keptIsA ? otherEndRemoved : otherEndKept;
  net.springIndexB[keptSpringIdx] = keptIsA ? otherEndKept : otherEndRemoved;
  net.springBoxOffset.segment(3 * keptSpringIdx, 3) +=
    net.springBoxOffset.segment(3 * removedSpringIdx, 3) *
    (keptIsA == removedIsA ? -1. : 1.);
  for (size_t dir = 0; dir < 3; ++dir) {
    net.springCoordinateIndexA[3 * keptSpringIdx + dir] =
      3 * net.springIndexA[keptSpringIdx] + dir;
    net.springCoordinateIndexB[3 * keptSpringIdx + dir] =
      3 * net.springIndexB[keptSpringIdx] + dir;
  }

  assert(net.springIndicesOfStrand[keptStrandIdx].size() + 1 ==
         net.linkIndicesOfStrand[keptStrandIdx].size());
  // merge strands if needed
  if (isStrandMerger) {
    std::vector<size_t> removedStrandsLinks =
      net.linkIndicesOfStrand[removedStrandIdx];
    assert(removedStrandsLinks[0] == linkToReduce ||
           removedStrandsLinks.back() == linkToReduce);
    std::vector<size_t> removedStrandsSprings =
      net.springIndicesOfStrand[removedStrandIdx];
    assert(removedStrandsSprings[0] == removedSpringIdx ||
           removedStrandsSprings.back() == removedSpringIdx);
    // adjust link of springs to strands
    for (const size_t springIdx : removedStrandsSprings) {
      assert(net.strandIndexOfSpring[springIdx] == removedStrandIdx);
      net.strandIndexOfSpring[springIdx] = keptStrandIdx;
    }
    // adjust link of links to strands
    // need to remove duplicates since the strand may have more than one link to
    // the link (in case of a primary loop), while the link will still have only
    // one link to the strand
    std::vector<size_t> uniqueRemovedStrandsLinks = removedStrandsLinks;
    pylimer_tools::utils::sort_remove_duplicates(uniqueRemovedStrandsLinks);
    for (const size_t linkIdx : uniqueRemovedStrandsLinks) {
      assert(pylimer_tools::utils::contains(net.strandIndicesOfLink[linkIdx],
                                            removedStrandIdx));
      for (size_t& strandIdx : net.strandIndicesOfLink[linkIdx]) {
        if (strandIdx == removedStrandIdx) {
          strandIdx = keptStrandIdx;
          break;
        }
      }
      pylimer_tools::utils::sort_remove_duplicates(
        net.strandIndicesOfLink[linkIdx]);
      assert(!pylimer_tools::utils::contains(net.strandIndicesOfLink[linkIdx],
                                             removedStrandIdx));
    }

    // start with moving the links and springs to the new strand
    bool invertSpringDirections = false;
    // make sure to keep directionality (ordered listing of strands & links)
    // of the kept/extended strand
    if (!keptIsA && removedIsA) {
      // direction is already the same, "ascending"
      assert(net.linkIndicesOfStrand[keptStrandIdx].back() == linkToReduce);
      assert(removedStrandsLinks[0] == linkToReduce);
      net.linkIndicesOfStrand[keptStrandIdx].pop_back();
      pylimer_tools::utils::append(
        net.linkIndicesOfStrand[keptStrandIdx], removedStrandsLinks, 1);
      assert(net.springIndicesOfStrand[removedStrandIdx][0] ==
             removedSpringIdx);
      assert(net.springIndicesOfStrand[keptStrandIdx].back() == keptSpringIdx);
      pylimer_tools::utils::append(
        net.springIndicesOfStrand[keptStrandIdx], removedStrandsSprings, 1);
    } else if (keptIsA && !removedIsA) {
      // direction is the same, but in both cases "descending"
      assert(net.linkIndicesOfStrand[keptStrandIdx][0] == linkToReduce);
      assert(removedStrandsLinks.back() == linkToReduce);
      removedStrandsLinks.pop_back();
      net.linkIndicesOfStrand[keptStrandIdx][0] = removedStrandsLinks.back();
      pylimer_tools::utils::prepend(
        net.linkIndicesOfStrand[keptStrandIdx], removedStrandsLinks, 0, 1);
      assert(removedStrandsSprings.back() == removedSpringIdx);
      removedStrandsSprings.pop_back();
      pylimer_tools::utils::prepend(net.springIndicesOfStrand[keptStrandIdx],
                                    removedStrandsSprings);
    } else if (keptIsA && removedIsA) {
      // direction is different
      assert(net.linkIndicesOfStrand[keptStrandIdx][0] == linkToReduce);
      assert(removedStrandsLinks[0] == linkToReduce);

      net.linkIndicesOfStrand[keptStrandIdx][0] = removedStrandsLinks[1];
      pylimer_tools::utils::prepend_inverse(
        net.linkIndicesOfStrand[keptStrandIdx], removedStrandsLinks, 0, 2);
      pylimer_tools::utils::prepend_inverse(
        net.springIndicesOfStrand[keptStrandIdx], removedStrandsSprings, 0, 1);
      // need to swap some things in the strand indices
      invertSpringDirections = true;
    } else {
      assert(!keptIsA && !removedIsA);
      // direction is different
      assert(removedStrandsLinks.back() == linkToReduce);
      assert(net.linkIndicesOfStrand[keptStrandIdx].back() == linkToReduce);

      removedStrandsLinks.pop_back();
      net.linkIndicesOfStrand[keptStrandIdx].pop_back();
      removedStrandsSprings.pop_back();

      pylimer_tools::utils::append_inverse(
        net.linkIndicesOfStrand[keptStrandIdx], removedStrandsLinks);
      pylimer_tools::utils::append_inverse(
        net.springIndicesOfStrand[keptStrandIdx], removedStrandsSprings);
      // need to swap some things in the strand indices
      invertSpringDirections = true;
    }
    if (invertSpringDirections) {
      for (const size_t springIdx : removedStrandsSprings) {
        assert(springIdx != keptSpringIdx);
        if (springIdx == removedSpringIdx) {
          continue;
        }
        std::swap(net.springIndexA[springIdx], net.springIndexB[springIdx]);
        for (size_t dir = 0; dir < 3; ++dir) {
          std::swap(net.springCoordinateIndexA[3 * springIdx + dir],
                    net.springCoordinateIndexB[3 * springIdx + dir]);
        }
        net.springBoxOffset.segment(3 * springIdx, 3) *= -1.;
      }
    }

    // reset the removed strand
    net.linkIndicesOfStrand[removedStrandIdx].clear();
    net.springIndicesOfStrand[removedStrandIdx].clear();
    // and link
    assert(net.strandIndicesOfLink[linkToReduce].size() <= 2);
  }

  // reset the removed spring
  net.springIndexA[removedSpringIdx] = -1;
  net.springIndexB[removedSpringIdx] = -1;
  net.strandIndexOfSpring[removedSpringIdx] = -1;
  if (!isStrandMerger) {
    pylimer_tools::utils::removeIfContained(
      net.springIndicesOfStrand[keptStrandIdx], removedSpringIdx);
    // the removal of the link <> strand relationship is a bit more complex,
    // since it could have been involved in the same strand twice, as a primary
    // loop.
    assert(XOR(keptIsA, removedIsA));
    assert(keptStrandIdx == removedStrandIdx);
    assert(pylimer_tools::utils::contains(
      net.linkIndicesOfStrand[keptStrandIdx], linkToReduce));
    pylimer_tools::utils::removeIfContained(
      net.linkIndicesOfStrand[keptStrandIdx], linkToReduce);
    assert(net.strandIndicesOfLink[linkToReduce].size() == 1);
    assert(net.strandIndicesOfLink[linkToReduce][0] == keptStrandIdx);
  }
  net.strandIndicesOfLink[linkToReduce].clear();

  // some more validation checks
  assert(net.springIndicesOfStrand[keptStrandIdx].size() + 1 ==
         net.linkIndicesOfStrand[keptStrandIdx].size());
  assert(!pylimer_tools::utils::contains(
    net.springIndicesOfStrand[keptStrandIdx], removedSpringIdx));
  // end is relevant, as the unlinked link, in the case of a crosslink,
  // could have been involved in the same strand twice, as both ends.
  const bool linkIsEnd =
    net.linkIndicesOfStrand[keptStrandIdx].back() == linkToReduce ||
    net.linkIndicesOfStrand[keptStrandIdx][0] == linkToReduce;
  assert(!(linkIsEnd && !isStrandMerger));
  assert(!pylimer_tools::utils::contains(net.linkIndicesOfStrand[keptStrandIdx],
                                         linkToReduce) ||
         linkIsEnd);
  assert(!pylimer_tools::utils::contains(net.strandIndicesOfLink[linkToReduce],
                                         removedStrandIdx) ||
         linkIsEnd);
  assert(pylimer_tools::utils::contains(net.linkIndicesOfStrand[keptStrandIdx],
                                        otherEndKept));

  // finish with the validation of the distances
  assert(this->getOtherLinkOfSpring(net, keptSpringIdx, otherEndKept) ==
         otherEndRemoved);
  const Eigen::Vector3d distanceAfter =
    this->evaluateSpringVectorFrom(net, u, keptSpringIdx, otherEndKept);
  assert(pylimer_tools::utils::vector_approx_equal(
    totalDistance, distanceAfter, 1e-9));

#ifndef NDEBUG
  assert(this->validateNetwork(net, u));
#endif
}

/**
 * @brief Remove links without any springs
 *
 * @param net the network to be modified
 * @param displacements the current displacements, will be modified
 * @return the number of 0-functional links removed
 */
size_t
MEHPForceBalance2::removeZerofunctionalLinks(
  ForceBalance2Network& net,
  Eigen::VectorXd& displacements) const
{
  std::vector<size_t> linksToRemove;
  std::vector<size_t> linksWithDirsToRemove;
  for (long int i = net.nrOfLinks - 1; i >= 0; --i) {
    if (net.strandIndicesOfLink[i].empty()) {
      linksToRemove.push_back(i);
      for (int dir = 2; dir >= 0; --dir) {
        linksWithDirsToRemove.push_back(i * 3 + dir);
      }
    }
  }

  if (linksToRemove.empty()) {
    return 0;
  }

#define REMOVE_ROWS(vector)                                                    \
  pylimer_tools::utils::removeRows(vector, linksToRemove, true, true)

  // actually do remove these links
  REMOVE_ROWS(net.linkIsEntanglement);
  REMOVE_ROWS(net.strandIndicesOfLink);
  REMOVE_ROWS(net.oldAtomIds);
  REMOVE_ROWS(net.oldAtomTypes);
#undef REMOVE_ROWS

#define REMOVE_ROWS(vector)                                                    \
  pylimer_tools::utils::removeRows(vector, linksWithDirsToRemove, true, true)

  REMOVE_ROWS(displacements);
  REMOVE_ROWS(net.coordinates);

#undef REMOVE_ROWS

  // renumber
  const std::vector<long int> newIndices =
    pylimer_tools::utils::getMappingForRenumbering(linksToRemove,
                                                   net.nrOfLinks);
  assert(newIndices.size() == net.nrOfLinks);
  pylimer_tools::utils::renumberWithMapping(net.springIndexA, newIndices);
  pylimer_tools::utils::renumberWithMapping(net.springIndexB, newIndices);
  pylimer_tools::utils::renumberWithMapping(net.linkIndicesOfStrand,
                                            newIndices);

  std::vector<long int> newIndicesWithDir;
  newIndicesWithDir.reserve(newIndices.size() * 3);
  for (size_t i = 0; i < newIndices.size(); ++i) {
    for (size_t dir = 0; dir < 3; ++dir) {
      newIndicesWithDir.push_back(newIndices[i] * 3 + dir);
    }
  }

  pylimer_tools::utils::renumberWithMapping(net.springCoordinateIndexA,
                                            newIndicesWithDir);
  pylimer_tools::utils::renumberWithMapping(net.springCoordinateIndexB,
                                            newIndicesWithDir);

  // recount
  net.nrOfLinks -= linksToRemove.size();
  net.nrOfNodes = (net.linkIsEntanglement == false).count();

#ifndef NDEBUG
  assert(this->validateNetwork(net, displacements));
#endif

  return linksToRemove.size();
};

/**
 * @brief Delete the springs indicated by `toDelete` from the network
 *
 * @param net the network to be modified
 * @param displacements the current displacements, will be modified
 * @param springsToDeleteIndices a vector of indices of springs to be removed
 */
void
MEHPForceBalance2::removeSprings(
  ForceBalance2Network& net,
  Eigen::VectorXd& displacements,
  std::vector<size_t>& springsToDeleteIndices) const
{
  // make sure things are sorted
  for (long int i = springsToDeleteIndices.size() - 2; i >= 0; --i) {
    INVALIDARG_EXP_IFN(
      springsToDeleteIndices[i + 1] < springsToDeleteIndices[i],
      "Indices must be sorted descending and unique, got values " +
        std::to_string(springsToDeleteIndices[i]) + "@" + std::to_string(i) +
        " and " + std::to_string(springsToDeleteIndices[i + 1]) + "@" +
        std::to_string(i + 1) + ".");
  }

  // check which strands to remove
  std::vector<size_t> strandsToDelete;
  for (const size_t i : springsToDeleteIndices) {
    if (net.strandIndexOfSpring[i] >= 0) {
      strandsToDelete.push_back(net.strandIndexOfSpring[i]);
    }
  }
  for (size_t strandIdx = 0; strandIdx < net.nrOfStrands; ++strandIdx) {
    if (net.linkIndicesOfStrand[strandIdx].empty()) {
      assert(net.springIndicesOfStrand[strandIdx].empty());
      strandsToDelete.push_back(strandIdx);
    }
  }
  std::ranges::sort(strandsToDelete, std::greater<size_t>());
  // remove duplicate strands
  strandsToDelete.erase(std::ranges::unique(strandsToDelete).begin(),
                        strandsToDelete.end());
  // assert that all springs of those strands will be removed without lefovers
  for (const size_t strandIdx : strandsToDelete) {
    for (size_t springIdx : net.springIndicesOfStrand[strandIdx]) {
      INVALIDARG_EXP_IFN(
        std::ranges::binary_search(
          springsToDeleteIndices, springIdx, std::greater<size_t>()),
        "Spring " + std::to_string(springIdx) + " of strand " +
          std::to_string(strandIdx) +
          " should be removed as well, as all springs of this strand.");
    }
  }

  // unlink the strands to be deleted from the links
  for (size_t strandIdx : strandsToDelete) {
    for (const size_t linkIdx : net.linkIndicesOfStrand[strandIdx]) {
      pylimer_tools::utils::removeIfContained(net.strandIndicesOfLink[linkIdx],
                                              strandIdx);
    }
  }

#ifndef NDEBUG
  for (size_t springIdx = 0; springIdx < net.nrOfSprings; ++springIdx) {
    if (net.strandIndexOfSpring[springIdx] >= 0) {
      if (std::ranges::binary_search(
            springsToDeleteIndices, springIdx, std::greater<size_t>())) {
        RUNTIME_EXP_IFN(
          std::ranges::binary_search(strandsToDelete,
                                     net.strandIndexOfSpring[springIdx],
                                     std::greater<size_t>()),
          "Expected strand " +
            std::to_string(net.strandIndexOfSpring[springIdx]) + " of spring " +
            std::to_string(springIdx) + " to be removed, but it is not.");
      } else {
        RUNTIME_EXP_IFN(
          !std::ranges::binary_search(strandsToDelete,
                                      net.strandIndexOfSpring[springIdx],
                                      std::greater<size_t>()),
          "Expected not to remove strand " +
            std::to_string(net.strandIndexOfSpring[springIdx]) + " of spring " +
            std::to_string(springIdx) + ", but it is.");
      }
    } else {
      RUNTIME_EXP_IFN(std::ranges::binary_search(springsToDeleteIndices,
                                                 springIdx,
                                                 std::greater<size_t>()),
                      "Expected to delete spring " + std::to_string(springIdx) +
                        " since it's without strand");
    }
  }
#endif

  // actually remove some rows
#define REMOVE_ROWS(vector)                                                    \
  pylimer_tools::utils::removeRows(vector, springsToDeleteIndices, true, true)

  REMOVE_ROWS(net.springContourLength);
  REMOVE_ROWS(net.strandIndexOfSpring);
  REMOVE_ROWS(net.springIsEntanglement);
  REMOVE_ROWS(net.springIndexA);
  REMOVE_ROWS(net.springIndexB);
#undef REMOVE_ROWS

  // assemble the indices of the springs, but with directions
  std::vector<size_t> toDeleteWithDirs;
  toDeleteWithDirs.reserve(springsToDeleteIndices.size() * 3);
  for (const size_t i : springsToDeleteIndices) {
    toDeleteWithDirs.push_back(i * 3 + 2);
    toDeleteWithDirs.push_back(i * 3 + 1);
    toDeleteWithDirs.push_back(i * 3);
  }

  // remove those rows as well
#define REMOVE_DIRECTED_ROWS(vector)                                           \
  pylimer_tools::utils::removeRows(vector, toDeleteWithDirs, true, true)

  REMOVE_DIRECTED_ROWS(net.springBoxOffset);
  REMOVE_DIRECTED_ROWS(net.springCoordinateIndexA);
  REMOVE_DIRECTED_ROWS(net.springCoordinateIndexB);

#undef REMOVE_DIRECTED_ROWS

  // finally, remove the strands
  pylimer_tools::utils::removeRows(
    net.linkIndicesOfStrand, strandsToDelete, true, true);
  pylimer_tools::utils::removeRows(
    net.springIndicesOfStrand, strandsToDelete, true, true);

  // renumber relationships
  // first, change the indices of the strands
  const std::vector<long int> newStrandMapping =
    pylimer_tools::utils::getMappingForRenumbering(strandsToDelete,
                                                   net.nrOfStrands);

  pylimer_tools::utils::renumberWithMapping(net.strandIndicesOfLink,
                                            newStrandMapping);
  pylimer_tools::utils::renumberWithMapping(net.strandIndexOfSpring,
                                            newStrandMapping);

  net.nrOfStrands -= strandsToDelete.size();

  // then, change the indices of the springs
  const std::vector<long int> newSpringMapping =
    pylimer_tools::utils::getMappingForRenumbering(springsToDeleteIndices,
                                                   net.nrOfSprings);
  pylimer_tools::utils::renumberWithMapping(net.springIndicesOfStrand,
                                            newSpringMapping);

  net.nrOfSprings -= springsToDeleteIndices.size();

#ifndef NDEBUG
  this->validateNetwork(net, displacements);
#endif
};

/**
 * @brief Deform the system to match the specified box
 *
 * @param newBox the box to deform to
 */
void
MEHPForceBalance2::deformTo(const pylimer_tools::entities::Box& newBox)
{
  this->box.adjustCoordinatesTo(this->initialConfig.coordinates, newBox);
  this->box.adjustCoordinatesTo(this->currentDisplacements, newBox);
  this->box.adjustCoordinatesTo(this->initialConfig.springBoxOffset, newBox);
  this->box = newBox;
  this->universe.setBox(newBox, true);
  for (size_t i = 0; i < 3; ++i) {
    this->initialConfig.L[i] = this->box.getL(i);
    this->initialConfig.boxHalfs[i] = 0.5 * this->initialConfig.L[i];
  }
}

/**
 * @brief Displace all links to the mean of all connected neighbours
 *
 * @param net the force balance network
 * @param u the current displacements, wherein the resulting coordinates
 * shall be stored
 * @return double, the distance (squared norm) displaced
 */
double
MEHPForceBalance2::displaceToMeanPosition(const ForceBalance2Network& net,
                                          Eigen::VectorXd& u) const
{
  Eigen::ArrayXd objectiveDisplacement =
    Eigen::ArrayXd::Zero(3 * net.nrOfLinks);
  const Eigen::ArrayXd partialSpringDistances =
    this->evaluateSpringVectors(net, u, this->is2D).array();
  // assemble
  Eigen::ArrayXd oneOverSpringPartitions =
    Eigen::ArrayXd::Ones(3 * net.nrOfSprings);
  for (size_t i = 0; i < net.nrOfSprings; ++i) {
    for (size_t dir = 0; dir < 3; ++dir) {
      oneOverSpringPartitions(3 * i + dir) = 1. / net.springContourLength(i);
    }
  }
  // use
  objectiveDisplacement(net.springCoordinateIndexA) +=
    (oneOverSpringPartitions * partialSpringDistances);
  objectiveDisplacement(net.springCoordinateIndexB) -=
    (oneOverSpringPartitions * partialSpringDistances);

  Eigen::ArrayXd springPartWeightingFactor =
    Eigen::ArrayXd::Zero(net.nrOfLinks * 3);
  const Eigen::ArrayXd loopPartialSpringEliminator =
    (net.springCoordinateIndexA != net.springCoordinateIndexB).cast<double>();

  springPartWeightingFactor(net.springCoordinateIndexA) +=
    oneOverSpringPartitions * loopPartialSpringEliminator;
  springPartWeightingFactor(net.springCoordinateIndexB) +=
    oneOverSpringPartitions * loopPartialSpringEliminator;
  springPartWeightingFactor = springPartWeightingFactor.unaryExpr(
    [](const double v) { return v > 0. ? v : 1.0; });
  Eigen::ArrayXd remainingDisplacement =
    (objectiveDisplacement / springPartWeightingFactor);
#ifndef NDEBUG
  RUNTIME_EXP_IFN(
    pylimer_tools::utils::all_components_finite(remainingDisplacement),
    "Some displacements are not finite");
#endif
  // at this point, we have the ideal displacement if we were to do it
  // just one link at a time.
  // by doing all at once, as here, though, e.g. a pair of links would
  // oscillate back and forth to compensate for that:

  // NOTE: this stays mostly static, could be stored on the network
  Eigen::ArrayXd nSpringsPerLink = Eigen::ArrayXd::Zero(net.nrOfLinks * 3);
  // add a one for every partial spring that's not a primary loop
  nSpringsPerLink(net.springCoordinateIndexA) += loopPartialSpringEliminator;
  nSpringsPerLink(net.springCoordinateIndexB) += loopPartialSpringEliminator;
  nSpringsPerLink =
    nSpringsPerLink.unaryExpr([](const double v) { return v > 0. ? v : 1.0; });
  // make sure there are no infinite back-and-forth
  // and actually displace
  Eigen::ArrayXd backForthDisplacement =
    Eigen::ArrayXd::Zero(net.nrOfLinks * 3);
  backForthDisplacement(net.springCoordinateIndexA) +=
    loopPartialSpringEliminator *
    (remainingDisplacement(net.springCoordinateIndexB) /
     (nSpringsPerLink(net.springCoordinateIndexA) * 2.));
  backForthDisplacement(net.springCoordinateIndexB) +=
    loopPartialSpringEliminator *
    (remainingDisplacement(net.springCoordinateIndexA) /
     (nSpringsPerLink(net.springCoordinateIndexB) * 2.));
#ifndef NDEBUG
  RUNTIME_EXP_IFN(
    pylimer_tools::utils::all_components_finite(backForthDisplacement),
    "Some displacements are not finite");
#endif

  // actually displace
  const Eigen::VectorXd finalDisplacement =
    (remainingDisplacement + backForthDisplacement).matrix();
  RUNTIME_EXP_IFN(
    pylimer_tools::utils::all_components_finite(finalDisplacement),
    "Some displacements are not finite");
  // this->box.handlePBC(finalDisplacement);
  u += finalDisplacement;

  const double max_disp =
    pylimer_tools::utils::segmentwise_norm_max(finalDisplacement, 3);

  return max_disp;
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
double
MEHPForceBalance2::displaceToMeanPosition(const ForceBalance2Network& net,
                                          Eigen::VectorXd& u,
                                          const size_t linkIdx) const
{
#ifndef NDEBUG
  const Eigen::Vector3d forceBefore = this->getForceOn(net, u, linkIdx);
#endif
  // Eigen::Vector3d currentDisplacement = u.segment(3 * linkIdx, 3);
  Eigen::Vector3d objectiveDisplacement =
    Eigen::Vector3d::Zero(); // = remainingDisplacement.array();
  double objectiveDisplacementContributors = 0.0;

  const std::vector<size_t> partialSpringIndices =
    this->getPartialSpringIndicesOfLink(net, linkIdx);

  for (const size_t globalSpringIndex : partialSpringIndices) {
    assert(net.springIndexA[globalSpringIndex] == linkIdx ||
           net.springIndexB[globalSpringIndex] == linkIdx);
    if (net.springIndexA[globalSpringIndex] == linkIdx &&
        net.springIndexB[globalSpringIndex] == linkIdx) {
      // skip primary loops
      continue;
    }
    Eigen::Vector3d partialDistance =
      this->evaluateSpringVectorFrom(net, u, globalSpringIndex, linkIdx);
    double oneOverContourLengthFraction =
      1. / net.springContourLength[globalSpringIndex];

    if (std::isfinite(oneOverContourLengthFraction)) {
      objectiveDisplacement += (partialDistance)*oneOverContourLengthFraction;
      objectiveDisplacementContributors += oneOverContourLengthFraction;
    }
  }
  // take mean for displacement
  // prevent NaN from division by zero
  const double denominator = 1. / (objectiveDisplacementContributors == 0.0
                                     ? 1.0
                                     : objectiveDisplacementContributors);
  u.segment(3 * linkIdx, 3) += objectiveDisplacement * denominator;

#ifndef NDEBUG
  const Eigen::Vector3d forceAfter = this->getForceOn(net, u, linkIdx);

  // this is only true if we don't have "full" PBC
  assert((pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
    forceAfter, Eigen::Vector3d::Zero(), 0.01)));
  if (!pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
        forceBefore, Eigen::Vector3d::Zero(), 0.01)) {
    assert(forceBefore.squaredNorm() >= forceAfter.squaredNorm());
  }
#endif

  const double dist = (objectiveDisplacement * denominator).squaredNorm();
  // if (dist > 0.1) {
  //   std::cout << "Moving " << linkIdx << " for " << dist
  //             << " with displacements " << u.segment(3 * linkIdx, 3)[0]
  //             << ", " << u.segment(3 * linkIdx, 3)[1] << ", "
  //             << u.segment(3 * linkIdx, 3)[2] << std::endl;
  //   std::cout << "For objective displacements " <<
  //   objectiveDisplacement[0]
  //             << ", " << objectiveDisplacement[1] << ", "
  //             << objectiveDisplacement[2] << ", for "
  //             << objectiveDisplacementContributors << "." << std::endl;
  // }
  return dist;
}

/**
 * @brief Compute the stress tensor on one cross- or slip-link
 *
 * @param linkIdx the index of the link to compute stress for
 * @param net the force balance network
 * @param u the current displacements of the links
 * @param debugNrSpringsVisited debug counter for number of springs visited
 * @return Eigen::Matrix3d the stress tensor on the link
 */
Eigen::Matrix3d
MEHPForceBalance2::evaluateStressOnLink(
  const size_t linkIdx,
  const ForceBalance2Network& net,
  const Eigen::VectorXd& u,
  Eigen::VectorXi& debugNrSpringsVisited) const
{
  Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();

  const std::vector<size_t> partialSpringIndices =
    this->getPartialSpringIndicesOfLink(net, linkIdx);

  for (const size_t globalSpringIndex : partialSpringIndices) {
    Eigen::Vector3d partialDistance =
      this->evaluateSpringVectorFrom(net, u, globalSpringIndex, linkIdx);
    const double oneOverContourLengthFraction =
      1.0 / net.springContourLength[globalSpringIndex];

    double multiplier = this->kappa * oneOverContourLengthFraction;

    stress += multiplier * partialDistance * partialDistance.transpose();
    debugNrSpringsVisited[globalSpringIndex] += 1;

    // also account for primary loops.
    // they may have non-zero length thanks to assuming the box is not
    // large enough...
    if (net.springIndexA[globalSpringIndex] ==
        net.springIndexB[globalSpringIndex]) {
      stress +=
        multiplier * (-partialDistance) * (-partialDistance).transpose();

      debugNrSpringsVisited[globalSpringIndex] += 1;
    }
  }

  return stress;
}

/**
 * @brief Compute the force acting on one cross- or slip-link
 *
 * @param linkIdx the index of the link to compute force for
 * @param net the force balance network
 * @param u the current displacements of the links
 * @param debugNrSpringsVisited debug counter for number of springs visited
 * @return Eigen::Vector3d the force acting on the link
 */
Eigen::Vector3d
MEHPForceBalance2::evaluateForceOnLink(
  const size_t linkIdx,
  const ForceBalance2Network& net,
  const Eigen::VectorXd& u,

  Eigen::VectorXi& debugNrSpringsVisited) const
{
  Eigen::Vector3d force = Eigen::Vector3d::Zero();

  const std::vector<size_t> partialSpringIndices =
    this->getPartialSpringIndicesOfLink(net, linkIdx);

  for (const size_t globalSpringIndex : partialSpringIndices) {
    // partial spring's force goes both ways -> is zero anyway
    // but, as it would not be included twice in the list,
    // we have to skip them
    if (net.springIndexA[globalSpringIndex] == linkIdx &&
        net.springIndexB[globalSpringIndex] == linkIdx) {
      if (debugNrSpringsVisited.size() > 0) {
        debugNrSpringsVisited[globalSpringIndex] += 2;
      }
      continue;
    }
    Eigen::Vector3d partialDistance =
      this->evaluateSpringVectorFrom(net, u, globalSpringIndex, linkIdx);
    const double oneOverContourLengthFraction =
      1. / net.springContourLength[globalSpringIndex];

    force += this->kappa * oneOverContourLengthFraction * partialDistance;
    if (debugNrSpringsVisited.size() > 0) {
      debugNrSpringsVisited[globalSpringIndex] += 1;
    }
  }

  return force;
}

/**
 * @brief Count the number of intra-chain slip-links
 * i.e., slip-links that entangle a strand with itself
 *
 * @return int
 */
int
MEHPForceBalance2::getNumIntraChainSlipLinks() const
{
  int result = 0;
  for (size_t i = 0; i < this->initialConfig.nrOfLinks; ++i) {
    if (!this->initialConfig.linkIsEntanglement[i]) {
      continue;
    }
    if (this->initialConfig.strandIndicesOfLink[i].size() == 1) {
      result += 1;
    } else if (this->initialConfig.strandIndicesOfLink[i].size() == 2 &&
               this->initialConfig.strandIndicesOfLink[i][0] ==
                 this->initialConfig.strandIndicesOfLink[i][1]) {
      result += 1;
    }
  }

  return result;
};

Eigen::Vector3d
MEHPForceBalance2::evaluateSpringVector(const ForceBalance2Network& net,
                                        const Eigen::VectorXd& u,
                                        const size_t springIdx,
                                        const bool is2d)
{
  Eigen::Vector3d dist =
    ((net.coordinates.segment(3 * net.springIndexB(springIdx), 3) +
      u.segment(3 * net.springIndexB(springIdx), 3)) -
     (net.coordinates.segment(3 * net.springIndexA(springIdx), 3) +
      u.segment(3 * net.springIndexA(springIdx), 3))) +
    net.springBoxOffset.segment(3 * springIdx, 3);

  if (is2d) {
    dist[2] = 0.0;
  }

  return dist;
}

/**
 * @brief Evaluate the vectors between the two ends of all partial springs
 *
 * @param net
 * @param u
 * @param is2D
 * @param assumeLarge
 * @return Eigen::VectorXd
 */
Eigen::VectorXd
MEHPForceBalance2::evaluateSpringVectors(const ForceBalance2Network& net,
                                         const Eigen::VectorXd& u,
                                         const bool is2D)
{
  // first, the distances
  assert(u.size() == net.coordinates.size());

  Eigen::VectorXd displacedCoords = net.coordinates + u;
  Eigen::VectorXd partialDistances =
    (displacedCoords(net.springCoordinateIndexB) -
     displacedCoords(net.springCoordinateIndexA)) +
    net.springBoxOffset;

  // reset for 2D systems
  if (is2D) {
    // partialDistances(Eigen::seq(2, net.nrOfPartialSprings, 3)) =
    //   Eigen::VectorXd::Zero(net.nrOfPartialSprings);
    for (size_t i = 2; i < 3 * net.nrOfSprings; i += 3) {
      partialDistances[i] = 0.0;
    }
  }

  return partialDistances;
}

/**
 * FORCE BALANCE DATA ACCESS
 */
/**
 * @brief Convert the current network back into a universe, consisting
 * only of crosslinkers
 */
pylimer_tools::entities::Universe
MEHPForceBalance2::getCrosslinkerVerse() const
{
  // convert nodes & springs back to a universe
  pylimer_tools::entities::Universe xlinkUniverse =
    pylimer_tools::entities::Universe(this->box);
  std::vector<long int> ids;
  ids.reserve(this->initialConfig.nrOfNodes);
  std::vector<int> types;
  types.reserve(this->initialConfig.nrOfNodes);
  std::vector<double> x;
  x.reserve(this->initialConfig.nrOfNodes);
  std::vector<double> y;
  y.reserve(this->initialConfig.nrOfNodes);
  std::vector<double> z;
  z.reserve(this->initialConfig.nrOfNodes);
  const std::vector<int> zeros =
    pylimer_tools::utils::initializeWithValue(this->initialConfig.nrOfNodes, 0);
  for (int i = 0; i < this->initialConfig.nrOfLinks; ++i) {
    if (this->initialConfig.linkIsEntanglement[i]) {
      continue;
    }
    x.push_back(this->initialConfig.coordinates[3 * i + 0] +
                this->currentDisplacements[3 * i + 0]);
    y.push_back(this->initialConfig.coordinates[3 * i + 1] +
                this->currentDisplacements[3 * i + 1]);
    z.push_back(this->initialConfig.coordinates[3 * i + 2] +
                this->currentDisplacements[3 * i + 2]);
    ids.push_back(this->initialConfig.oldAtomIds[i]);
    // override type, since the types may be different from
    // crossLinkerType if converted with dangling chains
    types.push_back(this->initialConfig.oldAtomTypes[i]);
  }
  assert(ids.size() == this->initialConfig.nrOfNodes);
  xlinkUniverse.addAtoms(ids, types, x, y, z, zeros, zeros, zeros);
  std::vector<long int> bondFrom;
  std::vector<long int> bondTo;
  bondFrom.reserve(this->initialConfig.nrOfStrands);
  bondTo.reserve(this->initialConfig.nrOfStrands);
  for (int i = 0; i < this->initialConfig.nrOfStrands; ++i) {
    bondFrom.push_back(
      this->initialConfig
        .oldAtomIds[this->initialConfig.linkIndicesOfStrand[i][0]]);
    bondTo.push_back(this->initialConfig.oldAtomIds[pylimer_tools::utils::last(
      this->initialConfig.linkIndicesOfStrand[i])]);
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
double
MEHPForceBalance2::getAverageStrandLength() const
{
  Eigen::VectorXd partialSpringVectors = this->evaluateSpringVectors(
    this->initialConfig, this->currentDisplacements);

  Eigen::VectorXd springVectors =
    Eigen::VectorXd::Zero(3 * this->initialConfig.nrOfStrands);
  for (size_t i = 0; i < this->initialConfig.nrOfSprings; ++i) {
    springVectors.segment(3 * this->initialConfig.strandIndexOfSpring[i], 3) +=
      partialSpringVectors.segment(3 * i, 3);
  }

  const double result =
    pylimer_tools::utils::segmentwise_norm_mean(springVectors, 3);
  assert(result >= 0.0);
  return result;
}

Eigen::VectorXd
MEHPForceBalance2::evaluateStrandLengths(const ForceBalance2Network& net,
                                         Eigen::VectorXd u) const
{
  Eigen::VectorXd springVectors = this->evaluateStrandVectors(net, u);
  Eigen::VectorXd springLengths = Eigen::VectorXd::Zero(net.nrOfStrands);
  for (size_t i = 0; i < net.nrOfStrands; ++i) {
    springLengths[i] = springVectors.segment(3 * i, 3).norm();
  }

  return springLengths;
}

Eigen::Vector3d
MEHPForceBalance2::evaluateStrandVector(const ForceBalance2Network& net,
                                        const Eigen::VectorXd u,
                                        const size_t strandIdx) const
{
  INVALIDARG_EXP_IFN(strandIdx < net.nrOfStrands, "Invalid strand index.");
  Eigen::Vector3d result = Eigen::Vector3d::Zero();
  for (const size_t springIdx : net.springIndicesOfStrand[strandIdx]) {
    result += this->evaluateSpringVector(net, u, springIdx);
  }
  return result;
};

Eigen::VectorXd
MEHPForceBalance2::evaluateStrandVectors(const ForceBalance2Network& net,
                                         const Eigen::VectorXd u) const
{
  Eigen::VectorXd partialSpringVectors = this->evaluateSpringVectors(net, u);
  Eigen::VectorXd strandVectors = Eigen::VectorXd::Zero(3 * net.nrOfStrands);
  for (size_t i = 0; i < net.nrOfSprings; ++i) {
    strandVectors.segment(net.strandIndexOfSpring[i], 3) +=
      partialSpringVectors.segment(3 * i, 3);
  }

  return strandVectors;
}

/**
 * @brief Sum partial spring fractions to get total fraction up to a target link
 *
 * @param net the network to analyze
 * @param springPartition the spring partition values
 * @param springIdx the index of the spring (strand)
 * @param targetLink the index of the target link
 * @return double the total fraction summed for the target link
 */
double
MEHPForceBalance2::sumToTotalFraction(const ForceBalance2Network& net,
                                      Eigen::VectorXd springPartition,
                                      const size_t springIdx,
                                      const size_t targetLink)
{
  double alpha = 0.;
  for (size_t i = 0; i < net.springIndicesOfStrand[springIdx].size(); ++i) {
    const size_t currentPartialSpringIdx =
      net.springIndicesOfStrand[springIdx][i];
    if (net.springIndexA[currentPartialSpringIdx] == targetLink) {
      return alpha;
    }
    alpha += springPartition[currentPartialSpringIdx];
    if (net.springIndexB[currentPartialSpringIdx] == targetLink) {
      return alpha;
    }
  }
  throw std::runtime_error("Did not find target link in spring.");
};

/**
 * @brief Get the denominator for a specified partial spring
 *
 * @param net
 * @param partialSpringIdx
 * @return double
 */
double
MEHPForceBalance2::getDenominatorOfPartialSpring(
  const ForceBalance2Network& net,
  const size_t partialSpringIdx)
{
  const double denominator = 1. / net.springContourLength[partialSpringIdx];

  assert(std::isfinite(denominator));
  return denominator;
}

/**
 * @brief Compute the stress tensor
 *
 * @param linkIndices the indices of the links to respect
 * @param net
 * @param u the current displacements
 * @return std::array<std::array<double, 3>, 3>
 */
Eigen::Matrix3d
MEHPForceBalance2::evaluateStressTensorForLinks(
  const std::vector<size_t> linkIndices,
  const ForceBalance2Network& net,
  const Eigen::VectorXd& u) const
{
  Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();

  const double halfOverVolume = 0.5 / (net.L[0] * net.L[1] * net.L[2]);

  Eigen::VectorXi debugNrSpringsVisited =
    Eigen::VectorXi::Zero(net.nrOfSprings);

  for (const size_t linkIdx : linkIndices) {
    Eigen::Matrix3d force =
      this->evaluateStressOnLink(linkIdx, net, u, debugNrSpringsVisited);
    /* spring contribution to the overall stress tensor */
    RUNTIME_EXP_IFN(std::isfinite(force.squaredNorm()),
                    "Got non-finite force contribution to stress tensor: " +
                      std::to_string(force.squaredNorm()) + " at link " +
                      std::to_string(linkIdx) + "!");
    stress += force;
  }

  return halfOverVolume * stress;
};

/**
 * @brief Compute the stress tensor
 *
 * @param net
 * @param u
 * @return std::array<std::array<double, 3>, 3>
 */
std::array<std::array<double, 3>, 3>
MEHPForceBalance2::evaluateStressTensorLinkBased(
  const ForceBalance2Network& net,
  const Eigen::VectorXd& u,
  const bool crosslinksOnly) const
{
  Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();

  const double halfOverVolume = 0.5 / (net.L[0] * net.L[1] * net.L[2]);

  Eigen::VectorXi debugNrSpringsVisited =
    Eigen::VectorXi::Zero(net.nrOfSprings);

  for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
    if (crosslinksOnly && net.linkIsEntanglement[linkIdx]) {
      continue;
    }
    Eigen::Matrix3d stressOnLink =
      this->evaluateStressOnLink(linkIdx, net, u, debugNrSpringsVisited);
    /* spring contribution to the overall stress tensor */
    RUNTIME_EXP_IFN(std::isfinite(stressOnLink.squaredNorm()),
                    "Got non-finite force contribution to stress tensor: " +
                      std::to_string(stressOnLink.squaredNorm()) + " at link " +
                      std::to_string(linkIdx) + "!");
    stress += stressOnLink;
  }

  std::array<std::array<double, 3>, 3> stressA;
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      stressA[i][j] = halfOverVolume * stress(i, j);
    }
  }

  if (!crosslinksOnly) {
    RUNTIME_EXP_IFN(
      debugNrSpringsVisited.sum() == 2 * net.nrOfSprings,
      "Every spring must be visited twice, got min " +
        std::to_string(debugNrSpringsVisited.minCoeff()) + " and max " +
        std::to_string(debugNrSpringsVisited.maxCoeff()) + ". Sum is " +
        std::to_string(debugNrSpringsVisited.sum()) + " instead of " +
        std::to_string(2 * net.nrOfSprings) + ".");
    RUNTIME_EXP_IFN((debugNrSpringsVisited.array() == 2).all(),
                    "Every spring must be visited twice, got min " +
                      std::to_string(debugNrSpringsVisited.minCoeff()) +
                      " and max " +
                      std::to_string(debugNrSpringsVisited.maxCoeff()) + ".");
  }

  return stressA;
}

/**
 * @brief Compute the stress tensor
 *
 * @param net
 * @param u
 * @return std::array<std::array<double, 3>, 3>
 */
std::array<std::array<double, 3>, 3>
MEHPForceBalance2::evaluateStressTensor(const ForceBalance2Network& net,
                                        const Eigen::VectorXd& u) const
{
  std::array<std::array<double, 3>, 3> stress;
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      stress[i][j] = 0.0;
    }
  }

  const double oneOverVolume = 1. / (net.L[0] * net.L[1] * net.L[2]);

  Eigen::VectorXd displacedCoords = net.coordinates + u;
  Eigen::VectorXd relevantPartialDistancesA =
    (displacedCoords(net.springCoordinateIndexB) -
     displacedCoords(net.springCoordinateIndexA)) +
    net.springBoxOffset;

  if (this->is2D) {
    for (size_t i = 2; i < relevantPartialDistancesA.size(); i += 3) {
      relevantPartialDistancesA[i] = 0.;
    }
  }

  for (size_t partialSpringIdx = 0; partialSpringIdx < net.nrOfSprings;
       ++partialSpringIdx) {
    Eigen::Vector3d distance =
      relevantPartialDistancesA.segment(3 * partialSpringIdx, 3);
    const double oneOverContourLengthFraction =
      1. / net.springContourLength[partialSpringIdx];

    /* spring contribution to the overall stress tensor */
    for (size_t j = 0; j < 3; j++) {
      for (size_t k = 0; k < 3; k++) {
        const double contribution = distance[j] * distance[k] * this->kappa *
                                    oneOverContourLengthFraction;
        RUNTIME_EXP_IFN(
          std::isfinite(contribution),
          "Got non-finite contribution to stress tensor: " +
            std::to_string(contribution) + " at coordinates " +
            std::to_string(k) + ", " + std::to_string(j) +
            " for partial spring " + std::to_string(partialSpringIdx) +
            " from distances " + std::to_string(distance[j]) + ", " +
            std::to_string(distance[k]) + " and denominator " +
            std::to_string(oneOverContourLengthFraction) + ".");
        // if (std::isfinite(denominator) && std::isfinite(contribution))
        // {
        stress[j][k] += contribution;
        // }
      }
    }
  }

  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      stress[i][j] *= oneOverVolume;
      RUNTIME_EXP_IFN(std::isfinite(stress[i][j]),
                      "Got non-finite stress tensor component: " +
                        std::to_string(stress[i][j]) + " at coordinates " +
                        std::to_string(i) + ", " + std::to_string(j) +
                        " from denominator " + std::to_string(oneOverVolume) +
                        ".");
    }
  }

  return stress;
}

Eigen::Matrix3d
MEHPForceBalance2::getStressTensor()
{
  std::array<std::array<double, 3>, 3> res =
    this->evaluateStressTensor(this->initialConfig, this->currentDisplacements);

  // convert the array to an Eigen matrix
  Eigen::Matrix3d convertedRes = Eigen::Matrix3d::Zero();
  for (size_t i = 0; i < 3; ++i) {
    convertedRes.row(i) = Eigen::Vector3d::Map(res[i].data(), 3);
  }
  return convertedRes;
}

Eigen::Matrix3d
MEHPForceBalance2::getStressTensorLinkBased(const bool crosslinksOnly) const
{
  const std::array<std::array<double, 3>, 3> res =
    this->evaluateStressTensorLinkBased(
      this->initialConfig, this->currentDisplacements, crosslinksOnly);
  Eigen::Matrix3d convertedRes = Eigen::Matrix3d::Zero();
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      convertedRes(i, j) = res[i][j];
    }
  }
  return convertedRes;
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
MEHPForceBalance2::getEffectiveFunctionalityOfAtoms(
  const double tolerance) const
{
  std::unordered_map<long int, int> results;
  results.reserve(this->initialConfig.nrOfNodes);

  Eigen::VectorXi nrOfActiveSpringsConnected =
    this->getNrOfActiveStrandsConnected(tolerance);
  for (size_t i = 0; i < this->initialConfig.nrOfLinks; i++) {
    if (this->initialConfig.linkIsEntanglement[i]) {
      continue;
    }
    results.emplace(this->initialConfig.oldAtomIds[i],
                    nrOfActiveSpringsConnected[i]);
  }
  return results;
}

/**
 * @brief Compute the weight fraction of non-active springs
 *
 * We go the full route via active and soluble in order to compensate for
 * removed springs and atoms
 *
 * @param net the network to analyze
 * @param u the displacements
 * @param tolerance the tolerance for considering springs as active
 * @return double the dangling weight fraction
 */
double
MEHPForceBalance2::computeDanglingWeightFraction(ForceBalance2Network& net,
                                                 const Eigen::VectorXd& u,
                                                 const double tolerance) const
{
  const double activeWeightFraction =
    this->computeActiveWeightFraction(net, u, tolerance);
  RUNTIME_EXP_IFN(APPROX_WITHIN(activeWeightFraction, 0., 1., 1e-6),
                  "Expect active weight fraction to be between 0 and 1, got " +
                    std::to_string(activeWeightFraction) + ".");
  const double solubleWeightFraction =
    this->computeSolubleWeightFraction(net, u, tolerance);
  RUNTIME_EXP_IFN(APPROX_WITHIN(solubleWeightFraction, 0., 1., 1e-6),
                  "Expect soluble weight fraction to be between 0 and 1, got " +
                    std::to_string(solubleWeightFraction) + ".");
  RUNTIME_EXP_IFN(
    APPROX_WITHIN(activeWeightFraction + solubleWeightFraction, 0., 1., 1e-6),
    "Expect active and soluble weight fraction to add up to maximum 1, "
    "got " +
      std::to_string(activeWeightFraction + solubleWeightFraction) + ".");

  // TODO: currently, the weight of the atoms is ignored
  return 1. - activeWeightFraction - solubleWeightFraction;
}

/**
 * @brief Compute the weight fraction of active springs
 *
 * @param net the network to analyze
 * @param u the displacements
 * @param tolerance the tolerance for considering springs as active
 * @return double the active weight fraction
 */
double
MEHPForceBalance2::computeActiveWeightFraction(ForceBalance2Network& net,
                                               const Eigen::VectorXd& u,
                                               const double tolerance) const
{
  if (net.nrOfStrands < 1) {
    return 0.;
  }
  // find all active springs
  Eigen::ArrayXb activeStrands = this->findActiveStrands(net, u, tolerance);
  const double nActiveStrands = activeStrands.count();
  if (nActiveStrands == 0) {
    return 0.;
  }

  // we re-assign the active links and springs
  // this is due to some topological peculiarities:
  // for example, if we have an entanglement link a and a crosslink b,
  // it's possible to have a scenario …-a-b-a-…, in which case these two
  // springs would otherwise be considered inactive,
  // even though they should belong to the active fraction.
  Eigen::ArrayXb activeLinks = Eigen::ArrayXb::Constant(net.nrOfLinks, false);
  Eigen::ArrayXb activeSprings =
    Eigen::ArrayXb::Constant(net.nrOfSprings, false);

  for (size_t strandIdx = 0; strandIdx < net.nrOfStrands; ++strandIdx) {
    if (!activeStrands[strandIdx]) {
      continue;
    }

    for (const size_t linkIdx : net.linkIndicesOfStrand[strandIdx]) {
      activeLinks[linkIdx] = true;
    }
    for (const size_t springIdx : net.springIndicesOfStrand[strandIdx]) {
      activeSprings[springIdx] = true;
    }
  }

  const double nActiveSprings = activeSprings.count();
  // as of now, the springsContourLength is equal to the number of bonds
  // from link to link. therefore, the number of atoms of each
  // of these springs is one less
  const double nActiveAtomsFromSprings =
    ((net.springContourLength.array() - Eigen::ArrayXd::Ones(net.nrOfSprings)) *
     activeSprings.cast<double>())
      .sum();

  // TODO: currently, the weight of the atoms is ignored
  const double nActiveCrossLinks =
    (activeLinks * (net.linkIsEntanglement == false)).count();
  const double nActiveEntanglementLinks =
    (activeLinks * net.linkIsEntanglement).count();
  const double nEntanglementSprings = net.springIsEntanglement.count();
  // normalize by the number of atoms
  const double nAtoms = this->assumeNetworkIsComplete
                          ? (this->inferNrOfAtomsFromNetwork(net))
                          : this->universe.getNrOfAtoms();
  const double result =
    (nActiveAtomsFromSprings +
     (nEntanglementSprings > 0 ? 1. : 2.) * nActiveEntanglementLinks +
     nActiveCrossLinks) /
    nAtoms;
  assert(APPROX_WITHIN(result, 0., 1., 1e-6));
  return result;
}

/**
 * @brief Find whether springs and links are in any way connected to an
 * active spring
 *
 * @param net the network that includes the connectivity
 * @param u the current displacements of the links
 * @param tolerance the tolerance for considering springs as active
 * @return std::pair<Eigen::ArrayXb, Eigen::ArrayXb> indices of springs
 * (first) and links (second) connected in any way to active springs
 */
std::pair<Eigen::ArrayXb, Eigen::ArrayXb>
MEHPForceBalance2::findSpringsAndLinksClusteredToActive(
  const ForceBalance2Network& net,
  const Eigen::VectorXd& u,
  const double tolerance) const
{
  INVALIDARG_EXP_IFN(u.size() == net.nrOfLinks * 3, "Invalid sizes.");

  // find all active springs
  Eigen::ArrayXb strandIsActive = this->findActiveStrands(net, u, tolerance);
  Eigen::ArrayXb springIsActive =
    Eigen::ArrayXb::Constant(net.nrOfSprings, false);

  // then, iteratively walk along the springs to mark those as "active"
  // that are connected to active springs
  bool hadChanged = true;
  Eigen::ArrayXb linkIsActive = Eigen::ArrayXb::Constant(net.nrOfNodes, false);
  while (hadChanged) {
    hadChanged = false;
    for (Eigen::Index linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
      if (linkIsActive[linkIdx]) {
        continue;
      }
      for (const int strandIdx : net.strandIndicesOfLink[linkIdx]) {
        if (strandIsActive[strandIdx]) {
          hadChanged = true;
          linkIsActive[linkIdx] = true;
          // mark all strands of this cross-link as connected to active
          for (const int innerStrandIdx : net.strandIndicesOfLink[linkIdx]) {
            strandIsActive[innerStrandIdx] = true;
            // and all springs of this strand as active
            for (const int springIdx :
                 net.springIndicesOfStrand[innerStrandIdx]) {
              springIsActive[springIdx] = true;
            }
          }
          break;
        }
      }
    }
  }

  return std::make_pair(springIsActive, linkIsActive);
}

/**
 * @brief Find strands and nodes (crosslinks) that are clustered to active
 * springs
 *
 * @param net the network that includes the connectivity
 * @param u the current displacements of the links
 * @param tolerance the tolerance for considering springs as active
 * @return std::pair<Eigen::ArrayXb, Eigen::ArrayXb> indices of strands
 * (first) and nodes (second) connected in any way to active springs
 */
std::pair<Eigen::ArrayXb, Eigen::ArrayXb>
MEHPForceBalance2::findStrandsAndNodesClusteredToActive(
  const ForceBalance2Network& net,
  const Eigen::VectorXd& u,
  const double tolerance) const
{
  std::pair<Eigen::ArrayXb, Eigen::ArrayXb> clusteredToActive =
    this->findSpringsAndLinksClusteredToActive(net, u, tolerance);

  Eigen::ArrayXb strandIsActive =
    Eigen::ArrayXb::Constant(net.nrOfStrands, false);
  Eigen::ArrayXb nodeIsActive = Eigen::ArrayXb::Constant(net.nrOfNodes, false);

  for (size_t i = 0; i < net.nrOfSprings; ++i) {
    if (clusteredToActive.first[i]) {
      strandIsActive[net.strandIndexOfSpring[i]] = true;
    }
  }

  size_t xlinkIdx = 0;
  for (size_t i = 0; i < net.nrOfLinks; ++i) {
    if (net.linkIsEntanglement[i]) {
      // ignore entanglement links
      continue;
    }
    if (clusteredToActive.second[i]) {
      nodeIsActive[xlinkIdx] = true;
    }
    xlinkIdx += 1;
  }

  return std::make_pair(strandIsActive, nodeIsActive);
}

/**
 * @brief Count the number of atoms that can be considered part of an
 * active cluster, i.e., are somehow connected to an active spring
 *
 * @param net the network to analyze
 * @param u the current displacements of the links
 * @param tolerance the tolerance for considering springs as active
 * @return double the number of active clustered atoms
 */
double
MEHPForceBalance2::countActiveClusteredAtoms(ForceBalance2Network& net,
                                             const Eigen::VectorXd& u,
                                             const double tolerance) const
{
  INVALIDARG_EXP_IFN(net.nrOfLinks * 3 == u.size(),
                     "Link displacements and network don't match");
  if (net.nrOfStrands < 1) {
    return 0.;
  }

  // if no universe (an empty universe) is set,
  // we assume that no dangling chains had been omitted,
  // and the network is complete and can be used as a basis
  // for the clusters
  if (this->assumeNetworkIsComplete) {
    std::pair<Eigen::ArrayXb, Eigen::ArrayXb> clusteredToActive =
      this->findSpringsAndLinksClusteredToActive(net, u, tolerance);
    Eigen::ArrayXb clusteredToActiveSprings = clusteredToActive.first;
    Eigen::ArrayXb clusteredToActiveLinks = clusteredToActive.second;
    // actually count with what we have now
    double nAtomsInActiveStrands =
      (clusteredToActiveSprings.cast<double>() *
       // spring contour length is equal to the number of bonds
       (net.springContourLength.array() -
        Eigen::ArrayXd::Ones(static_cast<Eigen::Index>(net.nrOfSprings)))
         // we don't want to count entanglement springs,
         // but since those have a contour length of 1,
         // it does not matter whether we explicitly don't check as such:
         // (net.springIsEntanglement == false)
         .cast<double>())
        .sum();
    // then, account for the links, be it cross- or entanglement links
    double nActiveCrossLinks =
      ((clusteredToActiveLinks.cast<double>() *
        (net.linkIsEntanglement == false).cast<double>()))
        .sum();
    double nActiveEntanglementLinks = (clusteredToActiveLinks.cast<double>() *
                                       net.linkIsEntanglement.cast<double>())
                                        .sum();
    // since the entanglement links come from two atoms,
    // we have to count them twice
    // TODO: currently, the weight of the atoms is ignored
    return nAtomsInActiveStrands + nActiveCrossLinks +
           2 * nActiveEntanglementLinks;
  }

  RUNTIME_EXP_IFN(this->universe.getNrOfAtoms() >= net.nrOfLinks,
                  "Expected number of atoms in universe to be at least "
                  "equal to the number of links, got " +
                    std::to_string(this->universe.getNrOfAtoms()) + " < " +
                    std::to_string(net.nrOfLinks) + ".");

  const std::vector<pylimer_tools::entities::Universe> clusters =
    this->universe.getClusters();
  std::vector<long int> vertexIdxToClusterIdx(this->universe.getNrOfAtoms());
  for (size_t i = 0; i < clusters.size(); ++i) {
    for (const pylimer_tools::entities::Atom& atom : clusters[i].getAtoms()) {
      vertexIdxToClusterIdx[this->universe.getIdxByAtomId(atom.getId())] = i;
    }
  }

  std::vector<bool> clusterIsActive(clusters.size(), false);

  // find active atoms
  const std::vector<int> activeNodeIndices =
    this->getIndicesOfActiveNodes(net, u, tolerance);

  for (const long int& nodeIdx : activeNodeIndices) {
    const long int universeAtomIdx =
      this->universe.getIdxByAtomId(net.oldAtomIds[nodeIdx]);
    clusterIsActive[vertexIdxToClusterIdx[universeAtomIdx]] = true;
  }

  double nClusteredAtoms = 0.;
  for (size_t i = 0; i < clusters.size(); ++i) {
    if (clusterIsActive[i]) {
      nClusteredAtoms += clusters[i].getNrOfAtoms();
    }
  }

  return nClusteredAtoms;
}

/**
 * @brief Compute the weight fraction of springs connected to active
 * springs (any depth)
 *
 * @param net the network to analyze
 * @param u the current displacements of the links
 * @param tolerance the tolerance for considering springs as active
 * @return double the soluble weight fraction
 */
double
MEHPForceBalance2::computeSolubleWeightFraction(ForceBalance2Network& net,
                                                const Eigen::VectorXd& u,
                                                const double tolerance) const
{
  INVALIDARG_EXP_IFN(net.nrOfLinks * 3 == u.size(),
                     "Link displacements and network don't match");
  if (net.nrOfStrands < 1) {
    return 1.;
  }
  const double nActiveClusteredAtoms =
    this->countActiveClusteredAtoms(net, u, tolerance);
  const double nAtoms = this->assumeNetworkIsComplete
                          ? this->inferNrOfAtomsFromNetwork(net)
                          : static_cast<double>(this->universe.getNrOfAtoms());

#ifndef NDEBUG
  if (this->assumeNetworkIsComplete && this->universe.getNrOfAtoms() > 0) {
    assert(APPROX_EQUAL(this->inferNrOfAtomsFromNetwork(net),
                        this->universe.getNrOfAtoms(),
                        0.5));
  }
#endif

  RUNTIME_EXP_IFN(nAtoms > 0,
                  "Expected number of atoms to be positive, got " +
                    std::to_string(nAtoms) + ".");
  // finally, normalize by the number of atoms.
  // NOTE: currently, the weight of the atoms is ignored
  return 1. - (nActiveClusteredAtoms / (nAtoms));
}

/**
 * @brief Get the indices of active Nodes
 *
 * @param net the network to analyze
 * @param u the current displacements of the links
 * @param tolerance the tolerance: springs under a certain length are
 * considered inactive
 * @return std::vector<int> the indices of active nodes
 */
std::vector<int>
MEHPForceBalance2::getIndicesOfActiveNodes(const ForceBalance2Network& net,
                                           const Eigen::VectorXd& u,
                                           const double tolerance) const
{
  std::vector<int> results;
  results.reserve(net.nrOfNodes);

  // find all active springs
  Eigen::ArrayXb strandIsActive = this->findActiveStrands(net, u, tolerance);

  for (Eigen::Index i = 0; i < net.nrOfLinks; i++) {
    if (net.linkIsEntanglement[i]) {
      continue;
    }

    for (const int strandIdx : net.strandIndicesOfLink[i]) {
      if (strandIsActive[strandIdx]) {
        results.push_back(static_cast<int>(i));
        break;
      }
    }
  }

  return results;
};

/**
 * @brief Get the atom ids of the active crosslinks (not entanglement
 * beads/links)
 *
 * @param tolerance the tolerance: springs under a certain length are
 * considered inactive
 * @return std::vector<long int> the atom ids of active nodes
 */
std::vector<long int>
MEHPForceBalance2::getAtomIdsOfActiveNodes(const double tolerance) const
{
  std::vector<long int> results;
  // find all active springs
  const std::vector<int> activeNodes = this->getIndicesOfActiveNodes(
    this->initialConfig, this->currentDisplacements, tolerance);

  results.reserve(activeNodes.size());

  for (const int nodeIdx : activeNodes) {
    results.push_back(this->initialConfig.oldAtomIds[nodeIdx]);
  }

  return results;
}

/**
 * @brief Get the overall spring lengths
 *
 * @return std::vector<double> a vector with the sum of the norms of all the
 * springs per strand
 */
std::vector<double>
MEHPForceBalance2::getOverallSpringLengths() const
{
  const std::vector<double> partialSpringDistances =
    this->getCurrentSpringLengths();
  assert(partialSpringDistances.size() == this->initialConfig.nrOfSprings);
  std::vector<double> results =
    std::vector<double>(this->initialConfig.nrOfStrands, 0.);
  for (size_t i = 0; i < this->initialConfig.nrOfSprings; ++i) {
    results[this->initialConfig.strandIndexOfSpring[i]] +=
      partialSpringDistances[i];
  }

  return results;
}

/**
 * @brief Get the current spring distances
 *
 * @return Eigen::VectorXd the current spring distances
 */
Eigen::VectorXd
MEHPForceBalance2::getCurrentSpringDistances() const
{
  Eigen::VectorXd partialSpringVectors = this->evaluateSpringVectors(
    this->initialConfig, this->currentDisplacements);

  return partialSpringVectors;
}

/**
 * @brief Get the current spring lengths
 *
 * @return std::vector<double> the current spring lengths
 */
std::vector<double>
MEHPForceBalance2::getCurrentSpringLengths() const
{
  const Eigen::VectorXd vecs = this->evaluateSpringVectors(
    this->initialConfig, this->currentDisplacements);

  return pylimer_tools::utils::segmentwise_norm(vecs, 3);
}

/**
 * @brief Get the Nr Of Active Springs connected to each node
 *
 * @param tolerance the tolerance: springs under a certain length are
 * considered inactive
 * @return Eigen::VectorXi
 */
Eigen::VectorXi
MEHPForceBalance2::getNrOfActiveStrandsConnected(const double tolerance) const
{
  if (this->initialConfig.nrOfLinks != this->initialConfig.nrOfNodes) {
    RUNTIME_EXP("This method, `getNrOfActiveStrandsConnected`, is not yet "
                "implemented for entangled systems.");
  }
  Eigen::VectorXi nrOfActiveSpringsConnected =
    Eigen::VectorXi::Zero(this->initialConfig.nrOfNodes);
  Eigen::ArrayXb springIsActive = this->findActiveStrands(tolerance);
  for (size_t i = 0; i < this->initialConfig.nrOfStrands; i++) {
    if (springIsActive[i]) {
      /* active spring */
      const int a = this->initialConfig.linkIndicesOfStrand[i][0];
      const int b =
        pylimer_tools::utils::last(this->initialConfig.linkIndicesOfStrand[i]);
      ++(nrOfActiveSpringsConnected[a]);
      ++(nrOfActiveSpringsConnected[b]);
    }
  }
  return nrOfActiveSpringsConnected;
}

/**
 * @brief Get the Nr Of Active Springs connected to each node
 *
 * @param tolerance the tolerance: springs under a certain length are
 * considered inactive
 * @return Eigen::VectorXi
 */
Eigen::VectorXi
MEHPForceBalance2::getNrOfActiveSpringsConnected(const double tolerance) const
{
  if (this->initialConfig.nrOfLinks != this->initialConfig.nrOfNodes) {
    RUNTIME_EXP("This method, `getNrOfActiveStrandsConnected`, is not yet "
                "implemented for entangled systems.");
  }
  Eigen::VectorXi nrOfActivePartialSpringsConnected =
    Eigen::VectorXi::Zero(this->initialConfig.nrOfNodes);
  Eigen::ArrayXb partialSpringIsActive = this->findActiveSprings(tolerance);
  // translate this to the nodes
  for (Eigen::Index i = 0; i < this->initialConfig.nrOfSprings; ++i) {
    if (partialSpringIsActive[i]) {
      /* active spring */
      // size_t a =
      //   this->initialConfig
      //     .springIndexA[this->initialConfig.partialToFullSpringIndex[i]];
      // size_t b =
      //   this->initialConfig
      //     .springIndexB[this->initialConfig.partialToFullSpringIndex[i]];
      const size_t a = this->initialConfig.springIndexA[i];
      const size_t b = this->initialConfig.springIndexB[i];
      if (!this->initialConfig.linkIsEntanglement[a]) {
        ++(nrOfActivePartialSpringsConnected[a]);
      }

      if (!this->initialConfig.linkIsEntanglement[b]) {
        ++(nrOfActivePartialSpringsConnected[b]);
      }
    }
  }
  return nrOfActivePartialSpringsConnected;
}

/**
 * @brief Get the Gamma Factor at the current step
 *
 * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
 * computed as phantom = N<b^2>.
 * @param nrOfChains the nr of chains to average over (can be different
 * from the nr of springs thanks to omitted free chains or primary loops)
 * @return double
 */
double
MEHPForceBalance2::getGammaFactor(double b02, const int nrOfChains) const
{
  if (b02 < 0) {
    b02 = this->defaultBondLength * this->defaultBondLength;
  }

  if (this->getNrOfSprings() == 0) {
    return 0.;
  }

  const Eigen::VectorXd gammaFactors = this->getGammaFactors(b02);

  if (nrOfChains < 1) {
    return gammaFactors.mean();
  } else {
    return gammaFactors.sum() / static_cast<double>(nrOfChains);
  }
}

/**
 * @brief Get the per-(partial)-spring gamma factors
 *
 * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
 * computed as phantom = N<b^2>.
 * @return Eigen::VectorXd
 */
Eigen::VectorXd
MEHPForceBalance2::getGammaFactors(const double b02) const
{
  Eigen::VectorXd springVectors = this->evaluateSpringVectors(
    this->initialConfig, this->currentDisplacements);

  Eigen::VectorXd gammaFactors(springVectors.size() / 3);
  const double commonDenominator = 1. / b02;
  for (size_t i = 0; i < springVectors.size() / 3; ++i) {
    const double oneOverContourLengthFraction =
      1.0 / this->initialConfig.springContourLength[i];
    gammaFactors[i] = springVectors.segment(3 * i, 3).squaredNorm() *
                      commonDenominator * oneOverContourLengthFraction;
    RUNTIME_EXP_IFN(
      std::isfinite(gammaFactors[i]),
      "Non-finite gamma factor for partial spring " + std::to_string(i) +
        ", computed from 1/N = " +
        std::to_string(oneOverContourLengthFraction) +
        ", b02 = " + std::to_string(b02) + ", and squared distance = " +
        std::to_string(springVectors.segment(3 * i, 3).squaredNorm()) + ".");
  }
  return gammaFactors;
}

/**
 * @brief Get the per-(partial)-spring gamma factors
 *
 * @param b02 the melt <b^2>: mean bond length; vgl. the required <R_0^2>,
 * computed as phantom = N<b^2>.
 * @param dir the direction (0=x, 1=y, 2=z)
 * @return Eigen::VectorXd
 */
Eigen::VectorXd
MEHPForceBalance2::getGammaFactorsInDir(const double b02, const int dir) const
{
  INVALIDARG_EXP_IFN(dir >= 0 && dir <= 2, "Invalid direction.");

  if (this->initialConfig.nrOfSprings == 0) {
    return Eigen::VectorXd::Zero(0);
  }

  Eigen::VectorXd springVectors = this->evaluateSpringVectors(
    this->initialConfig, this->currentDisplacements);

  Eigen::VectorXd gammaFactors =
    Eigen::VectorXd::Zero(springVectors.size() / 3);
  const double commonDenominator = 1. / b02;
  for (size_t i = 0; i < springVectors.size() / 3; ++i) {
    const double oneOverContourLengthFraction =
      1.0 / this->initialConfig.springContourLength[i];
    gammaFactors[i] = SQUARE(springVectors[3 * i + dir]) * commonDenominator *
                      oneOverContourLengthFraction;
    RUNTIME_EXP_IFN(std::isfinite(gammaFactors[i]),
                    "Non-finite gamma factor for partial spring " +
                      std::to_string(i) + ", computed from 1/N = " +
                      std::to_string(oneOverContourLengthFraction) +
                      ", and squared distance = " +
                      std::to_string(SQUARE(springVectors[3 * i + dir])) +
                      " in dir " + std::to_string(dir) + ".");
  }
  return gammaFactors;
}

/**
 * @brief Get the Weighted Partial Spring Length
 *
 * @param net the network to analyze
 * @param u the current displacements of the links
 * @param partialSpringIdx the index of the partial spring
 * @return double the weighted spring length
 */
double
MEHPForceBalance2::getWeightedSpringLength(const ForceBalance2Network& net,
                                           const Eigen::VectorXd& u,
                                           const size_t partialSpringIdx) const
{
  const double oneOverContourLengthFraction =
    1. / net.springContourLength[partialSpringIdx];
  return this->evaluateSpringVector(net, u, partialSpringIdx).norm() *
         oneOverContourLengthFraction;
}

/**
 * @brief Get the indices of partial springs connected to a given link
 *
 * @param net the network to analyze
 * @param linkIdx the index of the link to find partial springs for
 * @return std::vector<size_t> indices of partial springs connected to the link
 */
std::vector<size_t>
MEHPForceBalance2::getPartialSpringIndicesOfLink(
  const ForceBalance2Network& net,
  const size_t linkIdx) const
{
  INVALIDARG_EXP_IFN(linkIdx < net.nrOfLinks,
                     "The requested link does not exist");
  std::vector<size_t> partialSpringIndices;

  const std::vector<size_t> strandIndices = net.strandIndicesOfLink[linkIdx];

  for (const size_t strandIdx : strandIndices) {
    std::vector<size_t> springs = net.springIndicesOfStrand[strandIdx];
    for (size_t springIdx : springs) {
      if (this->isPartOfSpring(net, linkIdx, springIdx)) {
        partialSpringIndices.push_back(springIdx);
      }
    }
  }
  return partialSpringIndices;
}

Eigen::VectorXd
MEHPForceBalance2::getBondLengths()
{
  return this->evaluateSpringVectors(this->initialConfig,
                                     this->currentDisplacements);
}

Eigen::VectorXd
MEHPForceBalance2::getInitialCoordinates()
{
  return this->initialConfig.coordinates;
}

Eigen::VectorXd
MEHPForceBalance2::getCoordinates()
{
  return this->initialConfig.coordinates + this->currentDisplacements;
}

double
MEHPForceBalance2::getTemperature()
{
  std::cerr << "Warning: Temperature is not a reasonable metric for this "
               "type of computation."
            << std::endl;
  return 0;
}

size_t
MEHPForceBalance2::getNumParticles()
{
  return this->initialConfig.nrOfNodes;
}

// LCOV_EXCL_START
/**
 * @brief Debug helper function to print connectivity information for a specific
 * atom
 *
 * @param atomId the ID of the atom to debug
 */
void
MEHPForceBalance2::debugAtomVicinity(const size_t atomId) const
{
  long int atomIdx = -1;
  for (size_t i = 0; i < this->initialConfig.oldAtomIds.size(); ++i) {
    if (this->initialConfig.oldAtomIds[i] == atomId) {
      atomIdx = i;
      break;
    }
  }
  RUNTIME_EXP_IFN(atomIdx >= 0, "Atom not found.");
  std::cout << "Atom " << atomIdx << " (" << atomId << ")"
            << " connectivity:" << std::endl;
  for (const long int parentSpringIdx :
       this->initialConfig.strandIndicesOfLink[atomIdx]) {
    std::vector<size_t> allSpringIndices =
      this->initialConfig.springIndicesOfStrand[parentSpringIdx];
    std::string prefix = "";
    for (const size_t springIdx : allSpringIndices) {
      prefix += "\t";
      std::cout << prefix << "Spring " << springIdx << " (";
      std::cout << this->initialConfig.springIndexA[springIdx] << " ⟷ "
                << this->initialConfig.springIndexB[springIdx];
      std::cout << prefix << "\t";

      for (const long int linkIdx :
           this->initialConfig.linkIndicesOfStrand[springIdx]) {
        std::cout << linkIdx << " ";
        if (linkIdx < this->initialConfig.nrOfNodes) {
          std::cout << "(" << this->initialConfig.oldAtomIds[linkIdx] << ") ";
        }
      }
      std::cout << std::endl;
    }
  }
}

// LCOV_EXCL_STOP
/**
 * @brief Validate the integrity and consistency of the network structure
 *
 * @param net the network to validate
 * @param u the current displacements of the links
 * @return bool true if the network is valid, false otherwise
 */
bool
MEHPForceBalance2::validateNetwork(const ForceBalance2Network& net,
                                   const Eigen::VectorXd& u) const
{
  // std::cout << "Validating network..." << std::endl;
  /**
   * First, test dimensions
   */
  RUNTIME_EXP_IFN(!std::isinf(net.L[0]) && !std::isnan(net.L[0]),
                  "Box direction x must be scalar");
  RUNTIME_EXP_IFN(!std::isinf(net.L[1]) && !std::isnan(net.L[1]),
                  "Box direction y must be scalar");
  RUNTIME_EXP_IFN(!std::isinf(net.L[2]) && !std::isnan(net.L[2]),
                  "Box direction z must be scalar");
  RUNTIME_EXP_IFN(net.coordinates.size() == net.nrOfLinks * 3,
                  "Invalid size of coordinates");
  RUNTIME_EXP_IFN(u.size() == net.nrOfLinks * 3,
                  "Invalid size of displacements");
  RUNTIME_EXP_IFN(u.size() == net.coordinates.size(),
                  "Invalid size of displacements or coordinates");
  RUNTIME_EXP_IFN(net.springIndicesOfStrand.size() == net.nrOfStrands,
                  "Invalid size of connectivity map, got " +
                    std::to_string(net.springIndicesOfStrand.size()) + " for " +
                    std::to_string(net.nrOfStrands) + " springs.");
  RUNTIME_EXP_IFN(net.springContourLength.size() == net.nrOfSprings,
                  "Invalid size of contour lengths, got " +
                    std::to_string(net.springContourLength.size()) + " for " +
                    std::to_string(net.nrOfSprings) + " springs.");
  RUNTIME_EXP_IFN(net.springIsEntanglement.size() == net.nrOfSprings,
                  "Invalid size of springs types, got " +
                    std::to_string(net.springIsEntanglement.size()) + " for " +
                    std::to_string(net.nrOfSprings) + " springs.");
  RUNTIME_EXP_IFN(net.strandIndicesOfLink.size() == net.nrOfLinks,
                  "Invalid size of spring indices of links, got " +
                    std::to_string(net.linkIndicesOfStrand.size()) + " for " +
                    std::to_string(net.nrOfStrands) + " springs.");
  RUNTIME_EXP_IFN(net.linkIndicesOfStrand.size() == net.nrOfStrands,
                  "Invalid size of link indices of springs, got " +
                    std::to_string(net.linkIndicesOfStrand.size()) + " for " +
                    std::to_string(net.nrOfStrands) + " springs.");
  RUNTIME_EXP_IFN(net.linkIsEntanglement.size() == net.nrOfLinks,
                  "Invalid size of link is sliplink");
  RUNTIME_EXP_IFN(
    net.linkIsEntanglement.count() == (net.nrOfLinks - net.nrOfNodes),
    "Nr of nodes plus nr of slp-links should give the total nr of links");
  RUNTIME_EXP_IFN(net.oldAtomIds.size() == net.nrOfLinks,
                  "Invalid size of old atom ids");
  RUNTIME_EXP_IFN(net.oldAtomTypes.size() == net.nrOfLinks,
                  "Invalid size of old atom types");
  RUNTIME_EXP_IFN(net.springCoordinateIndexA.size() == net.nrOfSprings * 3,
                  "Invalid size of springPartCoordinateIndexA");
  RUNTIME_EXP_IFN(net.springCoordinateIndexB.size() == net.nrOfSprings * 3,
                  "Invalid size of springPartCoordinateIndexB");
  RUNTIME_EXP_IFN(net.springIndexA.size() == net.nrOfSprings,
                  "Invalid size of springPartIndexA");
  RUNTIME_EXP_IFN(net.springBoxOffset.size() == net.nrOfSprings * 3,
                  "Invalid size of springPartBoxOffset: expected " +
                    std::to_string(net.nrOfSprings * 3) + ", got " +
                    std::to_string(net.springBoxOffset.size()) + " rows.");
  RUNTIME_EXP_IFN(net.springIndexB.size() == net.nrOfSprings,
                  "Invalid size of springPartIndexB");
  RUNTIME_EXP_IFN(
    net.strandIndexOfSpring.size() == net.nrOfSprings,
    "Every partial spring must be able to map to the full spring.");
  const bool hasEntanglementSprings = net.springIsEntanglement.any();

  /**
   * Test maximum values
   */
  if (net.nrOfStrands > 0) {
    RUNTIME_EXP_IFN(net.strandIndexOfSpring.maxCoeff() < net.nrOfStrands,
                    "Partial spring must map to full spring, which must have "
                    "a lower index.");
  }
  if (net.nrOfSprings > 0) {
    RUNTIME_EXP_IFN(net.springCoordinateIndexA.maxCoeff() < 3 * net.nrOfLinks,
                    "Part coordinates must map to coordinates.");
    RUNTIME_EXP_IFN(net.springCoordinateIndexB.maxCoeff() < 3 * net.nrOfLinks,
                    "Part coordinates must map to coordinates.");
    RUNTIME_EXP_IFN(net.springIndexA.maxCoeff() < net.nrOfLinks,
                    "Part indices must map to links.");
    RUNTIME_EXP_IFN(net.springIndexB.maxCoeff() < net.nrOfLinks,
                    "Part indices must map to links.");
  }

  /**
   * Test reversibility of link <-> spring mapping
   */
  for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
    std::vector<size_t> thisLinksStrands = net.strandIndicesOfLink[linkIdx];
    std::sort(thisLinksStrands.begin(), thisLinksStrands.end());
    auto last = std::unique(thisLinksStrands.begin(), thisLinksStrands.end());
    RUNTIME_EXP_IFN(last == thisLinksStrands.end(),
                    "Expect each link to only have one back-link to the "
                    "strand, found back-links " +
                      pylimer_tools::utils::join(thisLinksStrands.begin(),
                                                 thisLinksStrands.end(),
                                                 std::string("_")) +
                      " for link " + std::to_string(linkIdx) + ".");
    for (const size_t strandIdx : thisLinksStrands) {
      std::vector<size_t> thisSpringsLinks = net.linkIndicesOfStrand[strandIdx];
      RUNTIME_EXP_IFN(std::find(thisSpringsLinks.begin(),
                                thisSpringsLinks.end(),
                                linkIdx) != thisSpringsLinks.end(),
                      "Strand must have a connection to the link, too. Did "
                      "not find link " +
                        std::to_string(linkIdx) + " in strand " +
                        std::to_string(strandIdx) + ".");
    }
  }

  /**
   * Test the assumptions on slip-links
   */
  for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
    if (net.linkIsEntanglement[linkIdx]) {
      RUNTIME_EXP_IFN(
        net.strandIndicesOfLink[linkIdx].size() <= 2,
        "Expect each entanglement to be involved in at most two strands, "
        "got " +
          std::to_string(net.strandIndicesOfLink[linkIdx].size()) + ".");
      std::vector<size_t> springIndices =
        this->getPartialSpringIndicesOfLink(net, linkIdx);
      bool anyIsPrimaryLoop = false;
      int nEntanglementSprings = 0;
      for (const size_t springIndex : springIndices) {
        anyIsPrimaryLoop =
          anyIsPrimaryLoop || this->isLoopingSpring(net, springIndex);
        nEntanglementSprings += (net.springIsEntanglement[springIndex]);
      }
      if (springIndices.size() - nEntanglementSprings == 2) {
        RUNTIME_EXP_IFN(
          net.strandIndicesOfLink[linkIdx].size() - nEntanglementSprings == 1,
          "Bifunctional entanglement link must be involved in "
          "only one strand. Got " +
            std::to_string(net.strandIndicesOfLink[linkIdx].size()) +
            " for link " + std::to_string(linkIdx) + ".");
      }
      RUNTIME_EXP_IFN((springIndices.size() - nEntanglementSprings) % 2 == 0 ||
                        anyIsPrimaryLoop,
                      "Entanglements must be either zero-, bi- or "
                      "tetrafunctional links. Got f = " +
                        std::to_string(springIndices.size()) + " with " +
                        std::to_string(nEntanglementSprings) +
                        " entanglement springs " + " for link " +
                        std::to_string(linkIdx) + ".");
    }
  }

  /**
   * Test the validitiy of springs and their mapping
   */
  Eigen::ArrayXi nrOfLinkMentionsInStrands =
    Eigen::ArrayXi::Zero(net.nrOfLinks);
  for (size_t i = 0; i < net.nrOfStrands; ++i) {
    RUNTIME_EXP_IFN(net.linkIndicesOfStrand[i].size() >= 2 ||
                      net.linkIndicesOfStrand[i].size() == 0,
                    "Each strand requires at least two links, got " +
                      std::to_string(net.linkIndicesOfStrand[i].size()) +
                      " at i = " + std::to_string(i) + ".");
    RUNTIME_EXP_IFN(
      (net.springIndicesOfStrand[i].size() ==
       net.linkIndicesOfStrand[i].size() - 1) ||
        (net.linkIndicesOfStrand[i].empty() &&
         net.springIndicesOfStrand[i].empty()),
      "Require a global index for each local one, got " +
        std::to_string(net.springIndicesOfStrand[i].size()) +
        " != " + std::to_string(net.linkIndicesOfStrand[i].size() - 1) +
        " for spring " + std::to_string(i) + ".");
    const std::vector<size_t> springIndices = net.springIndicesOfStrand[i];
    for (size_t idxInSpringIndices = 0;
         idxInSpringIndices < springIndices.size();
         ++idxInSpringIndices) {
      const size_t springIdx = net.springIndicesOfStrand[i][idxInSpringIndices];
      const size_t partner0 = net.linkIndicesOfStrand[i][idxInSpringIndices];
      const size_t partner1 =
        net.linkIndicesOfStrand[i][idxInSpringIndices + 1];
      RUNTIME_EXP_IFN(
        ((net.springIndexA[springIdx] == partner0 &&
          net.springIndexB[springIdx] == partner1)),
        "Expect linkIndicesOfStrand and springIndicesOfStrand "
        "ordering to correspond. Got partner0 = " +
          std::to_string(partner0) +
          ", partner1 = " + std::to_string(partner1) +
          " vs. those of the spring " + std::to_string(springIdx) + ": " +
          std::to_string(net.springIndexA[springIdx]) + " and " +
          std::to_string(net.springIndexB[springIdx]) + " in strand " +
          std::to_string(i) + " with spring indices " +
          pylimer_tools::utils::join(net.springIndicesOfStrand[i].begin(),
                                     net.springIndicesOfStrand[i].end(),
                                     std::string(", ")) +
          " and links " +
          pylimer_tools::utils::join(net.linkIndicesOfStrand[i].begin(),
                                     net.linkIndicesOfStrand[i].end(),
                                     std::string(", ")) +
          ".");
    }
    // the following is not guaranteed anymore with the removal of links
    // while running RUNTIME_EXP_IFN(
    //   net.linkIndicesOfSprings[i][0] <=
    //     net.linkIndicesOfSprings[i][net.linkIndicesOfSprings[i].size()
    //     - 1],
    //   "Springs must have increasing end-point indices");
    std::vector<size_t> links = net.linkIndicesOfStrand[i];
    const bool isEntanglementStrand =
      (springIndices.size() == 1 && net.springIsEntanglement[springIndices[0]]);
    for (size_t j = 0; j < links.size(); ++j) {
      const size_t link_idx = links[j];
      nrOfLinkMentionsInStrands[link_idx] += 1;
      if (isEntanglementStrand) {
        RUNTIME_EXP_IFN(
          net.linkIsEntanglement[link_idx],
          "Expected ends of entanglement spring to be entanglement links.");
      } else {
        RUNTIME_EXP_IFN(net.linkIsEntanglement[link_idx] ==
                          ((j != 0) && (j != (links.size() - 1))),
                        "Crosslinks must be first and last in a spring, "
                        "entanglements in-between. Found discrepancy at " +
                          std::to_string(j) + "/" +
                          std::to_string(links.size()) + " in spring " +
                          std::to_string(i) + ".")
      }
      std::vector<size_t> thisLinksStrands = net.strandIndicesOfLink[link_idx];
      RUNTIME_EXP_IFN(
        std::ranges::find(thisLinksStrands, i) != thisLinksStrands.end(),
        "Link must have a connection to the strand, too. Did not find "
        "strand " +
          std::to_string(i) + " in link " + std::to_string(link_idx) + ".");
    }
  }
  Eigen::ArrayXi nrOfLinkMentionsInSprings =
    Eigen::ArrayXi::Zero(net.nrOfLinks);
  for (size_t i = 0; i < net.nrOfSprings; ++i) {
    if (net.springIndexA[i] >= 0 && net.springIndexB[i] >= 0) {
      nrOfLinkMentionsInSprings[net.springIndexA[i]] += 1;
      nrOfLinkMentionsInSprings[net.springIndexB[i]] += 1;
    }
  }
  for (size_t i = 0; i < net.nrOfLinks; ++i) {
    if (net.linkIsEntanglement[i]) {
      RUNTIME_EXP_IFN(
        nrOfLinkMentionsInStrands[i] <= 2 + hasEntanglementSprings,
        "Expect each entanglement link to be mentioned at most twice in the "
        "links-of-springs mapping, but " +
          std::to_string(i) + " was mentioned " +
          std::to_string(nrOfLinkMentionsInStrands[i]) + " times.");
      if (hasEntanglementSprings) {
        RUNTIME_EXP_IFN(nrOfLinkMentionsInSprings[i] == 3 ||
                          nrOfLinkMentionsInSprings[i] == 2 ||
                          nrOfLinkMentionsInSprings[i] == 0,
                        "Expect each entanglement link to be mentioned 0, 2 or "
                        "3 times in the "
                        "spring end indices, but " +
                          std::to_string(i) + " was mentioned " +
                          std::to_string(nrOfLinkMentionsInSprings[i]) +
                          " times.");
      } else {
        RUNTIME_EXP_IFN(nrOfLinkMentionsInSprings[i] == 4 ||
                          nrOfLinkMentionsInSprings[i] == 2 ||
                          nrOfLinkMentionsInSprings[i] == 0,
                        "Expect each entanglement link to be mentioned 0, 2 or "
                        "4 times in the "
                        "spring end indices, but " +
                          std::to_string(i) + " was mentioned " +
                          std::to_string(nrOfLinkMentionsInSprings[i]) +
                          " times.");
      }
    }
  }

  /**
   * Test the validity of partial springs and their mapping
   */
  for (size_t i = 0; i < net.nrOfSprings; i++) {
    const long int partialEndA = net.springIndexA[i];
    const long int partialEndB = net.springIndexB[i];
    RUNTIME_EXP_IFN(net.springContourLength[i] > 0,
                    "Unexpected spring contour length, got " +
                      std::to_string(net.springContourLength[i]) +
                      " for spring " + std::to_string(i) + ".");
    RUNTIME_EXP_IFN(partialEndA < static_cast<long int>(net.nrOfLinks),
                    "Cannot have a spring (" + std::to_string(i) +
                      ") part larger " + std::to_string(partialEndA) +
                      " than the nr of links (" +
                      std::to_string(net.nrOfLinks) + ").");
    RUNTIME_EXP_IFN(partialEndB < static_cast<long int>(net.nrOfLinks),
                    "Cannot have a spring (" + std::to_string(i) +
                      ") part larger " + std::to_string(partialEndB) +
                      " than the nr of links (" +
                      std::to_string(net.nrOfLinks) + ").");
    if (partialEndA < 0) {
      RUNTIME_EXP_IFN(partialEndA == -1 && partialEndB == -1,
                      "Negative spring part indices may only be -1.");
      RUNTIME_EXP_IFN(net.strandIndexOfSpring[i] == -1,
                      "Expect disabled spring to also have a disabled strand.");
    } else {
      RUNTIME_EXP_IFN(partialEndB >= 0,
                      "Unexpected spring part index, got " +
                        std::to_string(partialEndB) + " for spring " +
                        std::to_string(i) + ".");
      RUNTIME_EXP_IFN(net.springCoordinateIndexA[3 * i] % 3 == 0,
                      "Expected spring part coordinates to be sequentially "
                      "built from spring parts.");
      RUNTIME_EXP_IFN(net.springCoordinateIndexB[3 * i] % 3 == 0,
                      "Expected spring part coordinates to be sequentially "
                      "built from spring parts.");
      for (int dir = 0; dir < 3; ++dir) {
        RUNTIME_EXP_IFN(
          net.springCoordinateIndexA[3 * i + dir] == 3 * partialEndA + dir,
          "Spring part index and coordinate index must match. Got " +
            std::to_string(net.springCoordinateIndexA[3 * i + dir]) +
            " but expected " + std::to_string(3 * partialEndA + dir) +
            " with dir = " + std::to_string(dir) + ".");
        RUNTIME_EXP_IFN(
          net.springCoordinateIndexB[3 * i + dir] == 3 * partialEndB + dir,
          "Spring part index and coordinate index must match. Got " +
            std::to_string(net.springCoordinateIndexB[3 * i + dir]) +
            " but expected " + std::to_string(3 * partialEndB + dir) +
            " with dir = " + std::to_string(dir) + ".");
      }
      RUNTIME_EXP_IFN(net.strandIndexOfSpring[i] >= 0,
                      "A spring must have an associated strand."
                      "This seems not to be the case for spring " +
                        std::to_string(i) + ".");
      if (net.springIsEntanglement[i]) {
        RUNTIME_EXP_IFN(
          net.linkIndicesOfStrand[net.strandIndexOfSpring[i]].size() == 2,
          "Entangled strands must have exactly two links, strand " +
            std::to_string(net.strandIndexOfSpring[i]) + " has " +
            std::to_string(
              net.linkIndicesOfStrand[net.strandIndexOfSpring[i]].size()) +
            " links.");
        RUNTIME_EXP_IFN(
          net.springContourLength[i] == 1,
          "Entangled springs must have a contour length of 1, got" +
            std::to_string(net.springContourLength[i]) + "for spring" +
            std::to_string(i) + ".");
        RUNTIME_EXP_IFN(
          net.springIndicesOfStrand[net.strandIndexOfSpring[i]].size() == 1,
          "Entangled strands must have exactly one spring, strand " +
            std::to_string(net.strandIndexOfSpring[i]) + " has " +
            std::to_string(
              net.springIndicesOfStrand[net.strandIndexOfSpring[i]].size()) +
            " springs.");
      }
    }
  }

  /**
   * Check that we do not have any nan or inf values in our vectors
   */
  for (long int coordI = 0; coordI < net.coordinates.size(); coordI++) {
    RUNTIME_EXP_IFN(std::isfinite(net.coordinates[coordI]),
                    "Coordinate component " + std::to_string(coordI) +
                      " must be finite, got " +
                      std::to_string(net.coordinates[coordI]) + ".");
    RUNTIME_EXP_IFN(std::isfinite(u[coordI]),
                    "Displacement component " + std::to_string(coordI) +
                      " must be finite, got " + std::to_string(u[coordI]) +
                      ".");
  }
  for (int dir = 0; dir < 3; ++dir) {
    RUNTIME_EXP_IFN(std::isfinite(net.L[dir]),
                    "Expected box size to be finite, got " +
                      std::to_string(net.L[dir]) + " in dir " +
                      std::to_string(dir) + ".");
    RUNTIME_EXP_IFN(net.L[dir] > 0.0,
                    "Expected box size to be positive, got " +
                      std::to_string(net.L[dir]) + " in dir " +
                      std::to_string(dir) + ".");
    RUNTIME_EXP_IFN(APPROX_EQUAL(net.boxHalfs[dir], 0.5 * net.L[dir], 1e-12),
                    "Expected box half to be half of box length");
  }

  // std::cout << "Validation passed." << std::endl;
  return true;
}

}
