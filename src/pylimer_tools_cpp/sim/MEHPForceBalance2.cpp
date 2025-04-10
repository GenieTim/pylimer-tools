#include "MEHPForceBalance2.h"
#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
// #include "../utils/MemoryUtil.h"
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <set>
#include <string>
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
  const std::vector<std::pair<size_t, size_t> > pairsOfAtoms =
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
  // then, from the others
  for (size_t vIdx = 0; vIdx < vertexDegrees.size(); ++vIdx) {
    if (vertexDegrees[vIdx] <= 2) {
      startingVertices.push_back(vIdx);
    }
  }
  assert(startingVertices.size() == u.getNrOfAtoms());

  Eigen::VectorXd vertexCoordinates =
    u.getUnwrappedVertexCoordinates(this->box);
  assert(vertexCoordinates.size() == u.getNrOfAtoms() * 3);

  std::vector<long int> oldVertexIdToNewLinkId =
    pylimer_tools::utils::initializeWithValue<long int>(u.getNrOfAtoms(), -1);
  std::vector<bool> vertex_visited =
    pylimer_tools::utils::initializeWithValue(u.getNrOfAtoms(), false);
  std::vector<bool> edge_followed =
    pylimer_tools::utils::initializeWithValue(u.getNrOfBonds(), false);

  igraph_adjlist_t adjacencyList;
  igraph_adjlist_init(
    &graph,
    &adjacencyList,
    IGRAPH_ALL,
    IGRAPH_LOOPS_TWICE,
    IGRAPH_MULTIPLE);
  igraph_inclist_t incidenceList;
  igraph_inclist_init(&graph, &incidenceList, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

  size_t nVerticesWithRelevantDegree = (std::ranges::count_if(
    vertexDegrees,
    [](int d) { return d > 0 && d != 2; }));
  size_t nLinksExpected = nVerticesWithRelevantDegree + pairsOfAtoms.size();
  size_t nStrandsExpected = nVerticesWithRelevantDegree / 2;
  size_t nSpringsExpected =
    nVerticesWithRelevantDegree / 2 + pairsOfAtoms.size();

  // data targets
  std::vector<size_t> springFrom;
  springFrom.reserve(nSpringsExpected);
  std::vector<size_t> springTo;
  springTo.reserve(nSpringsExpected);
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
    Eigen::Vector3d currentPosition = vertexCoordinates.segment<3>(3 * vIdx, 3);
    while (!vertexStack.empty()) {
      igraph_integer_t currentVertex = vertexStack.top();

      Eigen::Vector3d distanceToCurrent =
        vertexCoordinates.segment<3>(currentVertex * 3, 3) - currentPosition;
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
              oldVertexIdToNewLinkId[fst],
              oldVertexIdToNewLinkId[snd]);
          }
          addSpring = true;
        } else if (vertexDegrees[currentVertex] > 2 ||
                   vertexDegrees[currentVertex] == 1) {
          // 3. cross-link/end
          addNewLink = true;
          addSpring = true;
          // close the previous strand, start a new strand
          addNewStrand = true;
        }

        if (addSpring) {
          assert(newLinkIdStack.size() > 0);
          // link this new link to the previous link
          if (springIndicesOfStrand.size() <= currentStrandIdx) {
            assert(currentStrandIdx == springIndicesOfStrand.size());
            springIndicesOfStrand.push_back(
              std::vector<size_t>{ { springFrom.size() } });
            linkIndicesOfStrand.push_back(
              std::vector<size_t>{ newLinkIdStack.top() });
          } else {
            springIndicesOfStrand[currentStrandIdx].push_back(
              springFrom.size());
          }
          linkIndicesOfStrand[currentStrandIdx].push_back(
            addNewLink ? nextLinkIdx : oldVertexIdToNewLinkId[currentVertex]);
          springFrom.push_back(newLinkIdStack.top());

          if (oldVertexIdToNewLinkId[currentVertex] >= 0) {
            assert(!addNewLink);
            springTo.push_back(oldVertexIdToNewLinkId[currentVertex]);
          } else {
            assert(addNewLink);
            springTo.push_back(nextLinkIdx);
          }
          springBoxOffsets.emplace_back(Eigen::Vector3d::Zero());
          springContourLength.push_back(currentSpringContourLength);
          currentSpringContourLength = 0;
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

          newLinkIdStack.push(nextLinkIdx);
          nextLinkIdx += 1;
          assert(newLinkIdxToOldVertexIdx.size() == nextLinkIdx);
          assert(newLinkIdxToOldVertexIdx.size() == newLinkPositions.size());
        }
      }

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
              std::vector<size_t>{ { springFrom.size() } });
            linkIndicesOfStrand.push_back(
              std::vector<size_t>{ newLinkIdStack.top() });
          } else {
            springIndicesOfStrand[currentStrandIdx].push_back(
              springFrom.size());
          }
          linkIndicesOfStrand[currentStrandIdx].push_back(
            oldVertexIdToNewLinkId[neighbor]);
          springFrom.push_back(newLinkIdStack.top());
          springTo.push_back(oldVertexIdToNewLinkId[neighbor]);
          springContourLength.push_back(currentSpringContourLength + 1);
          edge_followed[edge] = true;
          // compute box offset
          Eigen::Vector3d distance =
            newLinkPositions[oldVertexIdToNewLinkId[neighbor]] -
            currentPosition;
          springBoxOffsets.emplace_back(this->box.getOffset(distance));

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
      assert(oldVertexIdToNewLinkId[this->universe.getIdxByAtomId(fst)] >= 0);
      springFrom.push_back(
        oldVertexIdToNewLinkId[this->universe.getIdxByAtomId(fst)]);
      assert(oldVertexIdToNewLinkId[this->universe.getIdxByAtomId(snd)] >= 0);
      springTo.push_back(
        oldVertexIdToNewLinkId[this->universe.getIdxByAtomId(snd)]);
      springContourLength.push_back(1.);
      springBoxOffsets.emplace_back(Eigen::Vector3d::Zero());
    }
  }

  // the difficult part is done, now
  // we can fill the network data structure
  assert(springFrom.size() == springTo.size());
  assert(springContourLength.size() == springFrom.size());
  assert(springFrom.size() == springBoxOffsets.size());
  this->initialConfig.nrOfSprings = springFrom.size();
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
    assert(springFrom[i] < nextLinkIdx);
    assert(springTo[i] < nextLinkIdx);
    this->initialConfig.springIndexA(i) = springFrom[i];
    this->initialConfig.springIndexB(i) = springTo[i];
    this->initialConfig.springContourLength(i) = springContourLength[i];
    for (size_t dir = 0; dir < 3; ++dir) {
      this->initialConfig.springCoordinateIndexA(3 * i + dir) =
        (springFrom[i] * 3 + dir);
      this->initialConfig.springCoordinateIndexB(3 * i + dir) =
        (springTo[i] * 3 + dir);
    }

    this->initialConfig.springBoxOffset.segment(3 * i, 3) = springBoxOffsets[i];
  }

  // link between springs and strands
  this->initialConfig.strandIdxOfSpring =
    Eigen::ArrayXi(this->initialConfig.nrOfSprings);
  for (size_t strandIdx = 0; strandIdx < springIndicesOfStrand.size();
       ++strandIdx) {
    std::vector<size_t>& springIndices = springIndicesOfStrand[strandIdx];
    for (size_t springIdx : springIndices) {
      this->initialConfig.strandIdxOfSpring(springIdx) = strandIdx;
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
        this->initialConfig.strandIndicesOfLink[linkIdx],
        strandIdx);
    }
  }

  // some more properties of the links
  this->initialConfig.coordinates = Eigen::VectorXd::Zero(nextLinkIdx * 3);

  for (size_t linkIdx = 0; linkIdx < nextLinkIdx; ++linkIdx) {
    this->initialConfig.oldAtomIds(linkIdx) =
      this->universe.getPropertyValue<long int>(
        "id",
        newLinkIdxToOldVertexIdx[linkIdx]);
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
    assert(this->initialConfig.springContourLength.sum() ==
      this->universe.getNrOfBonds());
  } else {
    assert(this->initialConfig.springContourLength.sum() ==
      this->universe.getNrOfBonds() +
      this->initialConfig.springIsEntanglement.count());
  }

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
void
MEHPForceBalance2::runForceRelaxation(
  const StructureSimplificationMode simplificationMode,
  const double inactiveRemovalCutoff,
  const SLESolver solverChoice,
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

  /* array allocation */
  Eigen::VectorXd oneOverSpringPartitions =
    this->assembleOneOverSpringPartition(this->initialConfig);
  const double initialResidual = this->getDisplacementResidualNormFor(
    this->initialConfig,
    this->currentDisplacements,
    oneOverSpringPartitions);
  std::cout << "Starting force balance procedure "
    << "with " << initialResidual << " as initial residual."
    << std::endl;
  std::cout << "Simplification mode is " << simplificationMode << std::endl;
  double currentResidual = initialResidual;
  size_t iterationsDone = 0;

  this->prepareAllOutputs();

  // actual loop
  bool wasInterrupted = false;
  size_t nRemovedInIteration;
  do {
    nRemovedInIteration = 0;
    std::vector<Eigen::Triplet<double> > triplets;
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
          linkIdx * 3 + dir,
          linkIdx * 3 + dir,
          diagonal(linkIdx * 3 + dir)));
      }
    }

    Eigen::SparseMatrix<double> sysMatrix(this->initialConfig.nrOfLinks * 3,
                                          this->initialConfig.nrOfLinks * 3);
    sysMatrix.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd finalCoordinates;

#define COMPUTE_SOLVE(solver)                                                  \
  solver.compute(sysMatrix);                                                   \
  RUNTIME_EXP_IFN(solver.info() == Eigen::Success,                             \
                  "System matrix computation failed. Solver info: " +          \
                    std::to_string(solver.info()));                            \
  finalCoordinates = solver.solve(constants);                                  \
  RUNTIME_EXP_IFN(solver.info() == Eigen::Success,                             \
                  "System could not be solved. Solver info: " +                \
                    std::to_string(solver.info()));

#define ADD_ITERATIONS(solver) iterationsDone += solver.iterations();

#define SOLVE_ITERATIVE(solver)                                                \
  COMPUTE_SOLVE(solver);                                                       \
  ADD_ITERATIONS(solver);                                                      \
  break;

#define SOLVE_DIRECT(solver)                                                   \
  COMPUTE_SOLVE(solver);                                                       \
  break;

    switch (solverChoice) {
      // iterative solvers
      case SLESolver::DEFAULT:
      case SLESolver::CONJUGATE_GRADIENT:
      case SLESolver::CONJUGATE_GRADIENT_DIAGONALIZED: {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                 Eigen::Lower,
                                 Eigen::DiagonalPreconditioner<double> >
          solver;
        SOLVE_ITERATIVE(solver);
      }
      case SLESolver::CONJUGATE_GRADIENT_IDENTITY: {
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                 Eigen::Lower,
                                 Eigen::IdentityPreconditioner>
          solver;
        SOLVE_ITERATIVE(solver);
      }
      case SLESolver::LEAST_SQUARES_CONJUGATE_GRADIENT: {
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> >
          solver;
        SOLVE_ITERATIVE(solver);
      }
      case SLESolver::BICGSTAB: {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
        SOLVE_ITERATIVE(solver);
      }
      // direct solvers
      case SLESolver::SIMPLICIAL_LLT: {
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > solver;
        SOLVE_DIRECT(solver);
      }
      case SLESolver::SIMPLICIAL_DLT: {
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
        SOLVE_DIRECT(solver);
      }
      case SLESolver::SPARSE_LU: {
        Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
        SOLVE_DIRECT(solver);
      }
      case SLESolver::SPARSE_QR: {
        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<
                          int> >
          solver;
        SOLVE_DIRECT(solver);
      }
      default:
        throw std::runtime_error("This solver is not implemented.");
    }

    this->currentDisplacements =
      finalCoordinates - this->initialConfig.coordinates;

    currentResidual = this->getDisplacementResidualNormFor(
      this->initialConfig,
      this->currentDisplacements);

    iterationsDone += 1;
    this->breakTooLongStrands(this->initialConfig, this->currentDisplacements);
    size_t nRemovedThisLoop = 0;

    do {
#ifndef NDEBUG
      RUNTIME_EXP_IFN(this->validateNetwork(), "Invalid internal state");
#endif
      nRemovedThisLoop = 0;
      if (simplificationMode == StructureSimplificationMode::INACTIVE_ONLY ||
          simplificationMode == StructureSimplificationMode::ALL_TIM) {
#ifdef DEBUG_REMOVAL
        std::cout << "Checking and possibly removing inactive cross-links"
                  << std::endl;
#endif
        nRemovedThisLoop +=
          this->removeInactiveCrosslinks(this->initialConfig,
                                         this->currentDisplacements,
                                         inactiveRemovalCutoff);
      }
      if (simplificationMode == StructureSimplificationMode::X2F_ONLY ||
          simplificationMode == StructureSimplificationMode::ALL_TIM) {
#ifdef DEBUG_REMOVAL
        std::cout << "Checking and possibly removing cross-links with f = 2"
                  << std::endl;
#endif
        nRemovedThisLoop += this->removeBifunctionalCrosslinks(
          this->initialConfig,
          this->currentDisplacements);
      }
      if (simplificationMode == StructureSimplificationMode::ALL_ANDREI) {
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

  // query solution & exit reason
  this->exitReason = ExitReason::X_TOLERANCE;
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
  Eigen::VectorXd oneOverSpringPartitions =
    this->assembleOneOverSpringPartition(this->initialConfig);
  Eigen::VectorXd displacements = this->currentDisplacements;
  return this->getDisplacementResidualNormFor(
    this->initialConfig,
    displacements,
    oneOverSpringPartitions);
}

/**
 * @brief Compute the displacement residual norm for a specific
 * configuration
 *
 * @param net
 * @param u
 * @return double
 */
double
MEHPForceBalance2::getDisplacementResidualNormFor(
  const ForceBalance2Network& net,
  const Eigen::VectorXd& u) const
{
  Eigen::ArrayXd oneOverSpringPartitions =
    this->assembleOneOverSpringPartition(net).array();

  Eigen::ArrayXd loopPartialSpringEliminator =
    (net.springCoordinateIndexA != net.springCoordinateIndexB).cast<double>();
  Eigen::ArrayXd forces = Eigen::ArrayXd::Zero(3 * net.nrOfLinks);
  Eigen::ArrayXd distances =
    this->evaluateSpringVectors(net, u, this->is2D, this->assumeBoxLargeEnough)
    .array();
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
 * configuration
 *
 * @param net
 * @param u
 * @param oneOverSpringPartitions
 * @return double
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

  if (this->assumeBoxLargeEnough) {
    this->box.handlePBC(relevantPartialDistances);
  }

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
  Eigen::VectorXd partialDistancesOverSpringPartitions =
    (relevantPartialDistances.array() * oneOverSpringPartitions.array())
    .matrix();

  Eigen::VectorXd overallForces = Eigen::VectorXd::Zero(3 * net.nrOfLinks);
  overallForces(net.springCoordinateIndexB) -=
    partialDistancesOverSpringPartitions;
  overallForces(net.springCoordinateIndexA) +=
    partialDistancesOverSpringPartitions;

  return overallForces.squaredNorm();
}

std::vector<size_t>
MEHPForceBalance2::getNeighbourLinkIndices(const ForceBalance2Network& net,
                                           const size_t linkIdx) const
{
  std::vector<size_t> results;
  results.reserve(4);
  for (size_t springIdx : net.strandIndicesOfLink[linkIdx]) {
    for (size_t partialSpringIdx : net.springIndicesOfStrand[springIdx]) {
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
 *
 * @param net
 * @return Eigen::VectorXd
 */
Eigen::VectorXd
MEHPForceBalance2::assembleOneOverSpringPartition(
  const ForceBalance2Network& net) const
{
  Eigen::VectorXd oneOverSpringPartitions =
    Eigen::VectorXd(3 * net.nrOfSprings);

  Eigen::ArrayXd primaryLoopCorrectionMultiplier =
    (net.springIndexA != net.springIndexB)
    .cast<double>(); // 0.0 for equal = primary loop, 1.0 otherwise

  for (size_t i = 0; i < net.nrOfSprings; ++i) {
    double valueToSet = 1.0 / (net.springContourLength[i]);

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
 * @brief Remove double listed springs from cross-linkers
 *
 * @param net
 */
void
MEHPForceBalance2::removeDuplicateListedStrandsFromLinks(
  ForceBalance2Network& net) const
{
  for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
    this->removeDuplicateListedStrandsFromLink(net, linkIdx);
  }

#ifndef NDEBUG
  assert(this->validateNetwork());
#endif
}

void
MEHPForceBalance2::removeDuplicateListedStrandsFromLink(
  ForceBalance2Network& net,
  size_t linkIdx,
  bool allowOnEntanglement) const
{
  INVALIDARG_EXP_IFN(linkIdx < net.nrOfLinks,
                     "Cannot remove duplicate spring indices of index "
                     "higher than nr. of links.");
  // remove duplicate mentions of the same spring index
  std::sort(net.strandIndicesOfLink[linkIdx].begin(),
            net.strandIndicesOfLink[linkIdx].end());
  auto last = std::unique(net.strandIndicesOfLink[linkIdx].begin(),
                          net.strandIndicesOfLink[linkIdx].end());
  if (last != net.strandIndicesOfLink[linkIdx].end()) {
#ifdef DEBUG_REMOVAL
    std::cout << "Removed duplicate spring indices from link " << linkIdx
              << std::endl;
#endif
    if (!allowOnEntanglement) {
      RUNTIME_EXP_IFN(
        !net.linkIsEntanglement[linkIdx],
        "Require entanglement beads to not form primary loops. Link " +
        std::to_string(linkIdx) + " is slip-link " +
        std::to_string(net.linkIsEntanglement[linkIdx]) + ".");
    }
    net.strandIndicesOfLink[linkIdx].erase(
      last,
      net.strandIndicesOfLink[linkIdx].end());
  }
}

size_t
MEHPForceBalance2::removePrimaryLoops(ForceBalance2Network& net,
                                      Eigen::VectorXd& displacements) const
{
  size_t numRemoved = 0;
  for (long int springIdx = net.nrOfStrands - 1; springIdx >= 0; --springIdx) {
    if (net.linkIndicesOfStrand[springIdx][0] ==
        pylimer_tools::utils::last(net.linkIndicesOfStrand[springIdx]) &&
        net.springIndicesOfStrand[springIdx].size() == 1) {
      this->removeStrandFollowingEntanglementLinks(net,
                                                   displacements,

                                                   springIdx);
      springIdx = std::min<long int>(springIdx, net.nrOfStrands - 1);
      numRemoved += 1;
    }
  }

#ifndef NDEBUG
  assert(this->validateNetwork());
#endif
  return numRemoved;
}

/**
 * @brief Remove crosslinkers which do not have any springs with a certain
 * minimum length
 *
 * @param net
 * @param displacements
 * @param tolerance
 */
size_t
MEHPForceBalance2::removeInactiveCrosslinks(ForceBalance2Network& net,
                                            Eigen::VectorXd& displacements,
                                            const double tolerance) const
{
  size_t numRemoved = 0;
  //        this->removePrimaryLoops(net, displacements);
  // RUNTIME_EXP_IFN(this->validateNetwork(net, displacements), "Invalid
  // internal network representation"); first, we remove all inactive springs
  for (long int springIdx = net.nrOfStrands - 1; springIdx >= 0; --springIdx) {
    if (springIdx >= net.nrOfStrands) {
      springIdx = net.nrOfStrands - 1;
      continue;
    }
    std::vector<size_t> involvedPartialSprings =
      net.springIndicesOfStrand[springIdx];
    assert(involvedPartialSprings.size() > 0);
    bool isActive = false;
    for (size_t partialSpringIdx : involvedPartialSprings) {
      RUNTIME_EXP_IFN(
        net.coordinates.size() == displacements.size(),
        "Expected coordinates and displacements to have same size, got " +
        std::to_string(net.coordinates.size()) + " and " +
        std::to_string(displacements.size()) + ".");
      // assert(net.coordinates.size() == displacements.size());
      Eigen::Vector3d distance =
        this->evaluateSpringVector(net, displacements, partialSpringIdx);
      const double contourLength =
        net.springContourLength[net.strandIdxOfSpring[partialSpringIdx]];
      if (!this->distanceIsWithinTolerance(
        distance,
        tolerance,
        contourLength)) {
        isActive = true;
        break;
      }
    }
    if (!isActive) {
      // remove this spring
#ifdef DEBUG_REMOVAL
      std::cout << "Removing inactive spring " << springIdx
                << " with all dependencies" << std::endl;
#endif
      this->removeStrandFollowingEntanglementLinks(net,
                                                   displacements,

                                                   springIdx);

#ifndef NDEBUG
      assert(this->validateNetwork(net, displacements));
#endif
      numRemoved += 1;
    }
  }

  // then, we remove all crosslinkers that are 0- or 1-functional
  for (long int crosslinkIdx = net.nrOfNodes - 1; crosslinkIdx >= 0;
       --crosslinkIdx) {
    assert(net.strandIndicesOfLink.size() > crosslinkIdx);
    if (net.strandIndicesOfLink[crosslinkIdx].size() == 0 // f = 0
    ) {
#ifdef DEBUG_REMOVAL
      std::cout << "Removing f = 0 x-link " << crosslinkIdx << std::endl;
#endif

      this->removeLink(net, displacements, crosslinkIdx);
      numRemoved += 1;
#ifndef NDEBUG
      assert(this->validateNetwork(net, displacements));
#endif
    } else if ( // or f = 1, NOT primary loop
      (net.strandIndicesOfLink[crosslinkIdx].size() == 1) &&
      (XOR(
        net.linkIndicesOfStrand[net.strandIndicesOfLink[crosslinkIdx][0]][0] ==
        crosslinkIdx,
        pylimer_tools::utils::last(
          net.linkIndicesOfStrand[net.strandIndicesOfLink[crosslinkIdx][0]]) ==
        crosslinkIdx))) {
#ifdef DEBUG_REMOVAL
      std::cout << "Removing f = 1 x-link " << crosslinkIdx << std::endl;
#endif
      // need to first remove the spring
      this->removeStrandFollowingEntanglementLinks(
        net,
        displacements,
        net.strandIndicesOfLink[crosslinkIdx][0]);
      numRemoved += 1;
      // to then remove the cross-link
      this->removeLink(net, displacements, crosslinkIdx);
#ifdef DEBUG_REMOVAL
      std::cout << "Effectively removed f = 1 x-link " << newCrosslinkIdx
                << std::endl;
#endif
      crosslinkIdx = std::min<long int>(crosslinkIdx, net.nrOfNodes - 1);
      // => we should ever only have 2-functional entanglement links that
      // could be merged after this.

#ifndef NDEBUG
      assert(this->validateNetwork(net, displacements));
#endif
    }
  }

#ifndef NDEBUG
  assert(this->validateNetwork(net, displacements));
#endif

  return numRemoved;
}

/**
 * @brief Remove springs that exert a stress higher than
 * `this->springBreakingLength`
 *
 * @param net
 * @param displacements
 * @return size_t the number of springs broken
 */
size_t
MEHPForceBalance2::breakTooLongStrands(ForceBalance2Network& net,
                                       Eigen::VectorXd& displacements) const
{
  if (this->springBreakingLength <= 0.) {
    return 0;
  }

  size_t numBroken = 0;

  // iterate the springs, determine their distance, and determine if it
  // exceeds the breaking force
  for (long int partialSpringIdx = net.nrOfSprings; partialSpringIdx >= 0;
       --partialSpringIdx) {
    if (partialSpringIdx >= net.nrOfSprings) {
      partialSpringIdx = net.nrOfSprings - 1;
    }
    double len = this->getWeightedSpringLength(net,
                                               displacements,

                                               partialSpringIdx);
    if (len > this->springBreakingLength) {
      // break this spring
      numBroken += 1;
      this->breakSpring(net,
                        displacements,

                        partialSpringIdx);
    }
  }

  return numBroken;
}

/**
 * @brief Remove a spring (and all its parts, incl. slip-links) from the
 * structures
 *
 * @param net
 * @param displacements
 * @param springIdx
 */
void
MEHPForceBalance2::removeStrand(ForceBalance2Network& net,
                                Eigen::VectorXd& displacements,
                                const size_t springIdx) const
{
#ifdef DEBUG_REMOVAL
  std::cout << "Starting to remove spring " << springIdx << std::endl;
#endif
  INVALIDARG_EXP_IFN(springIdx < net.nrOfStrands,
                     "Can only remove springs, not partial springs.");
#ifndef NDEBUG
  Eigen::VectorXd allTotalSpringDistancesBefore =
    this->evaluateStrandLengths(net, displacements);
#endif

  std::vector<size_t> affectedLinks = net.linkIndicesOfStrand[springIdx];
  std::vector<size_t> uniqueAffectedLinks = net.linkIndicesOfStrand[springIdx];
  std::sort(uniqueAffectedLinks.begin(), uniqueAffectedLinks.end());
  uniqueAffectedLinks.erase(
    std::unique(uniqueAffectedLinks.begin(), uniqueAffectedLinks.end()),
    uniqueAffectedLinks.end());

  // remove the link to the link, höhö
  for (size_t affectedLinkIdx : uniqueAffectedLinks) {
    RUNTIME_EXP_IFN(
      std::find(net.strandIndicesOfLink[affectedLinkIdx].begin(),
        net.strandIndicesOfLink[affectedLinkIdx].end(),
        springIdx) != net.strandIndicesOfLink[affectedLinkIdx].end(),
      "Link must have a connection to the spring, too. Did not find "
      "spring " +
      std::to_string(springIdx) + " in link " +
      std::to_string(affectedLinkIdx) + ", got " +
      pylimer_tools::utils::join(
        net.strandIndicesOfLink[affectedLinkIdx].begin(),
        net.strandIndicesOfLink[affectedLinkIdx].end(),
        std::string(", ")) +
      ".");
    if (net.linkIsEntanglement[affectedLinkIdx]) {
      RUNTIME_EXP_IFN(
        net.strandIndicesOfLink[affectedLinkIdx].size() <= 2,
        "Expect slip-link to be associated with 2 springs only, got " +
        pylimer_tools::utils::join(
          net.strandIndicesOfLink[affectedLinkIdx].begin(),
          net.strandIndicesOfLink[affectedLinkIdx].end(),
          std::string(", ")) +
        ".");
    }

    size_t found =
      std::erase(net.strandIndicesOfLink[affectedLinkIdx], springIdx);

    RUNTIME_EXP_IFN(found > 0,
                    "Expected to find spring " + std::to_string(springIdx) +
                    " in link " + std::to_string(affectedLinkIdx) +
                    " but did not, got " +
                    pylimer_tools::utils::join(
                      net.strandIndicesOfLink[affectedLinkIdx].begin(),
                      net.strandIndicesOfLink[affectedLinkIdx].end(),
                      std::string(", ")) +
                    ".");
    if (net.linkIsEntanglement[affectedLinkIdx]) {
      RUNTIME_EXP_IFN(net.strandIndicesOfLink[affectedLinkIdx].size() <= 1,
                      "Expect slip-link to be associated with 1 springs "
                      "only after removing one, got " +
                      pylimer_tools::utils::join(
                        net.strandIndicesOfLink[affectedLinkIdx].begin(),
                        net.strandIndicesOfLink[affectedLinkIdx].end(),
                        std::string(", ")) +
                      ".");
    }
  }

  std::vector<size_t> affectedPartialSprings =
    net.springIndicesOfStrand[springIdx];
  assert(affectedPartialSprings.size() > 0);
  net.nrOfStrands -= 1;
  net.nrOfSprings -= affectedPartialSprings.size();

  // actually spring remove stuff
  net.springIndicesOfStrand.erase(net.springIndicesOfStrand.begin() +
                                  springIdx);
  net.linkIndicesOfStrand.erase(net.linkIndicesOfStrand.begin() + springIdx);
  pylimer_tools::utils::removeRow(net.springContourLength, springIdx);
  pylimer_tools::utils::removeRow(net.springIsEntanglement, springIdx);

  // need to remove descending
  std::sort(affectedPartialSprings.begin(),
            affectedPartialSprings.end(),
            std::greater<size_t>());
  for (size_t partialSpringIdx : affectedPartialSprings) {
    pylimer_tools::utils::removeRow(net.springIndexA, partialSpringIdx);
    pylimer_tools::utils::removeRows(
      net.springCoordinateIndexA,
      3 * partialSpringIdx,
      3);
    pylimer_tools::utils::removeRow(net.springIndexB, partialSpringIdx);
    pylimer_tools::utils::removeRows(
      net.springCoordinateIndexB,
      3 * partialSpringIdx,
      3);
    pylimer_tools::utils::removeRows(
      net.springBoxOffset,
      3 * partialSpringIdx,
      3);
    pylimer_tools::utils::removeRow(net.strandIdxOfSpring, partialSpringIdx);
  }

  // renumber the remaining stuff
  // first, renumber the springs
  for (size_t linkIdx = 0; linkIdx < net.nrOfLinks; ++linkIdx) {
    for (size_t i = 0; i < net.strandIndicesOfLink[linkIdx].size(); ++i) {
      assert(net.strandIndicesOfLink[linkIdx][i] != springIdx);
      if (net.strandIndicesOfLink[linkIdx][i] > springIdx) {
        net.strandIndicesOfLink[linkIdx][i] -= 1;
      }
    }
  }

  // then, update the partial springs
  // decrease by one if larger then springIdx
  assert((net.strandIdxOfSpring != springIdx).all());
  net.strandIdxOfSpring -= (net.strandIdxOfSpring > springIdx).cast<int>();

  assert(net.springIndicesOfStrand.size() == net.nrOfStrands);
  for (size_t loopingSpringIdx = 0; loopingSpringIdx < net.nrOfStrands;
       ++loopingSpringIdx) {
    for (size_t i = 0; i < net.springIndicesOfStrand[loopingSpringIdx].size();
         ++i) {
      for (size_t partialSpringIdx : affectedPartialSprings) {
        if (net.springIndicesOfStrand[loopingSpringIdx][i] > partialSpringIdx) {
          net.springIndicesOfStrand[loopingSpringIdx][i] -= 1;
        }
      }
    }
  }

  //
  // remove the affected slip-links that are now only on one spring
  //
  std::vector<size_t> linksToRemove;
  linksToRemove.reserve(affectedLinks.size() -
                        2); // keep the first and last links
  for (size_t i = 1; i < affectedLinks.size() - 1; ++i) {
    linksToRemove.push_back(affectedLinks[i]);
  }

  // need to remove descending to remove need to renumber these as well
  std::sort(linksToRemove.begin(), linksToRemove.end(), std::greater<size_t>());
  linksToRemove.erase(std::unique(linksToRemove.begin(), linksToRemove.end()),
                      linksToRemove.end());
  if (linksToRemove.size() > 2) {
    assert(linksToRemove[0] > linksToRemove[1]);
  }

  Eigen::ArrayXb springIsAffected =
    Eigen::ArrayXb::Constant(net.nrOfStrands, false);
  for (size_t outermostI = 0; outermostI < linksToRemove.size(); ++outermostI) {
    size_t slipLinkIdx = linksToRemove[outermostI];
    assert(net.linkIsEntanglement[slipLinkIdx]);
    // first, merge the two other partial springs
    std::vector<size_t> springsOfLink = net.strandIndicesOfLink[slipLinkIdx];
    RUNTIME_EXP_IFN(springsOfLink.size() <= 1,
                    "Expected slip-link " + std::to_string(slipLinkIdx) +
                    " to have only 1 remaining spring, got " +
                    std::to_string(springsOfLink.size()) + " due to " +
                    pylimer_tools::utils::join(springsOfLink.begin(),
                      springsOfLink.end(),
                      std::string(", ")) +
                    ".");
    std::vector<size_t> involvedPartialSprings;
    involvedPartialSprings.reserve(2);
    for (int springInLinkIdx = springsOfLink.size() - 1; springInLinkIdx >= 0;
         --springInLinkIdx) {
      for (size_t partialSpringIdx :
           net.springIndicesOfStrand[springsOfLink[springInLinkIdx]]) {
        if (net.springIndexA[partialSpringIdx] == slipLinkIdx ||
            net.springIndexB[partialSpringIdx] == slipLinkIdx) {
          involvedPartialSprings.push_back(partialSpringIdx);
        }
      }
    }
    RUNTIME_EXP_IFN(
      involvedPartialSprings.size() >= springsOfLink.size(),
      "Expected more or equal number of partial springs (" +
      std::to_string(involvedPartialSprings.size()) + "; " +
      pylimer_tools::utils::join(involvedPartialSprings.begin(),
        involvedPartialSprings.end(),
        std::string(", ")) +
      ") than springs (" + std::to_string(springsOfLink.size()) + "; " +
      pylimer_tools::utils::join(
        springsOfLink.begin(), springsOfLink.end(), std::string(", ")) +
      ").");
    // RUNTIME_EXP_IFN(springsOfLink.size() % 2 == 0, "Expected link to
    // have an even number of components, got " +
    // std::to_string(springsOfLink.size()) + ".");
    if (involvedPartialSprings.size() > 0) {
      RUNTIME_EXP_IFN(
        involvedPartialSprings.size() == 2,
        "Expected only 2 involved partial springs, got: " +
        pylimer_tools::utils::join(involvedPartialSprings.begin(),
          involvedPartialSprings.end(),
          std::string(", ")) +
        " for springs " +
        pylimer_tools::utils::join(
          springsOfLink.begin(), springsOfLink.end(), std::string(", ")) +
        " when removing spring " + std::to_string(springIdx) + ".");
      assert(involvedPartialSprings.size() == 2);
      size_t partialSpringToKeep =
        std::min(involvedPartialSprings[0], involvedPartialSprings[1]);
      size_t partialSpringToRemove =
        std::max(involvedPartialSprings[0], involvedPartialSprings[1]);
      assert(partialSpringToKeep != partialSpringToRemove);
      assert(net.strandIdxOfSpring[partialSpringToKeep] ==
        net.strandIdxOfSpring[partialSpringToRemove]);
      // actually do the merge
      this->mergeSprings(net,
                         displacements,

                         partialSpringToRemove,
                         partialSpringToKeep,
                         slipLinkIdx);
      // total distance changes -> cannot use for checking the total
      // distance
      springIsAffected[net.strandIdxOfSpring[partialSpringToKeep]] = true;
    }

    assert(net.strandIndicesOfLink[slipLinkIdx].empty());

    // then, actually remove the slip-link
#ifdef DEBUG_REMOVAL
    std::cout << "Removing slip-link " << slipLinkIdx << std::endl;
#endif
    this->removeLink(net, displacements, slipLinkIdx);
  }
#ifdef DEBUG_REMOVAL
  std::cout << "Removed spring " << springIdx << std::endl;
#endif
#ifndef NDEBUG
  Eigen::VectorXd allTotalSpringDistancesAfter =
    this->evaluateStrandLengths(net, displacements);
  assert(allTotalSpringDistancesAfter.size() == net.nrOfStrands);
  for (size_t i = 0; i < net.nrOfStrands; ++i) {
    size_t correspondingOldIdx = i >= springIdx ? i + 1 : i;
    if (springIsAffected[i]) {
      // when slip-links are removed, the overall distance must reduce.
      RUNTIME_EXP_IFN(allTotalSpringDistancesBefore[correspondingOldIdx] +
                      1e-9 >=
                      allTotalSpringDistancesAfter[i],
                      "Expected that the total distances stay constant for "
                      "non-changed springs.");
    } else {
      RUNTIME_EXP_IFN(
        APPROX_EQUAL(allTotalSpringDistancesBefore[correspondingOldIdx],
          allTotalSpringDistancesAfter[i],
          1e-9),
        "Expected that the total distances stay constant for non-changed "
        "springs.");
    }
  }
#endif
}

/**
 * @brief Remove a spring, but also all springs that are connected to it
 * and are connected via entanglement links.
 *
 * @param net
 * @param displacements
 * @param springIdx
 */
void
MEHPForceBalance2::removeStrandFollowingEntanglementLinks(
  ForceBalance2Network& net,
  Eigen::VectorXd& displacements,
  const size_t springIdx) const
{
  std::vector<size_t> springsToRemove = { springIdx };
  std::vector<size_t> linksToRemove = net.linkIndicesOfStrand[springIdx];
  assert(springsToRemove.size() > 0);
  std::ranges::sort(springsToRemove, std::greater<>());
  for (const size_t springIdxToDelete : springsToRemove) {
    this->removeStrand(net,
                       displacements,

                       springIdxToDelete);
  }
  std::ranges::sort(linksToRemove, std::greater<>());
  for (size_t linkIdxToDelete : linksToRemove) {
    assert(net.strandIndicesOfLink[linkIdxToDelete].size() <= 1);
    if (net.strandIndicesOfLink[linkIdxToDelete].size() == 1) {
#ifdef DEBUG_REMOVAL
      std::cout << "Removing additional spring between entanglement links "
                << net.springIndicesOfLinks[linkIdxToDelete][0] << std::endl;
#endif

      this->removeStrand(
        net,
        displacements,
        net.strandIndicesOfLink[linkIdxToDelete][0]);
    }
    this->removeLink(net, displacements, linkIdxToDelete);
  }
};

/**
 * @brief break a spring, given its partial spring index
 *
 * @param net
 * @param displacements
 * @param partialSpringIdx
 */
void
MEHPForceBalance2::breakSpring(ForceBalance2Network& net,
                               Eigen::VectorXd& displacements,
                               const size_t partialSpringIdx) const
{
  this->removeStrandFollowingEntanglementLinks(
    net,
    displacements,
    net.strandIdxOfSpring[partialSpringIdx]);
};

/**
 * @brief remove a link from the network
 *
 * @param net
 * @param displacements
 * @param linkIdx
 */
void
MEHPForceBalance2::removeLink(ForceBalance2Network& net,
                              Eigen::VectorXd& displacements,
                              const size_t linkIdx) const
{
  INVALIDARG_EXP_IFN(net.strandIndicesOfLink[linkIdx].size() == 0,
                     "The springs have to be removed or re-linked before "
                     "removing the link.");
#ifdef DEBUG_REMOVAL
  std::cout << "Removing link " << linkIdx << std::endl;
#endif

  pylimer_tools::utils::removeRows(net.coordinates, linkIdx * 3, 3);
  pylimer_tools::utils::removeRows(displacements, linkIdx * 3, 3);

  if (!net.linkIsEntanglement[linkIdx]) {
    net.nrOfNodes -= 1;
    pylimer_tools::utils::removeRow(net.oldAtomIds, linkIdx);
    pylimer_tools::utils::removeRow(net.oldAtomTypes, linkIdx);
  }
  net.nrOfLinks -= 1;
  pylimer_tools::utils::removeRow(net.linkIsEntanglement, linkIdx);
  net.strandIndicesOfLink.erase(net.strandIndicesOfLink.begin() + linkIdx);

  // renumber the remaining links
  for (size_t i = 0; i < net.linkIndicesOfStrand.size(); ++i) {
    for (size_t j = 0; j < net.linkIndicesOfStrand[i].size(); ++j) {
#ifndef NDEBUG
      RUNTIME_EXP_IFN(
        net.linkIndicesOfStrand[i][j] != linkIdx,
        "Expected not to find link to remove " + std::to_string(linkIdx) +
        " in any spring, found in spring " + std::to_string(i) + ", " +
        pylimer_tools::utils::join(net.linkIndicesOfStrand[i].begin(),
          net.linkIndicesOfStrand[i].end(),
          std::string(", ")) +
        ".");
#endif
      if (net.linkIndicesOfStrand[i][j] > linkIdx) {
        net.linkIndicesOfStrand[i][j] -= 1;
      }
    }
  }
  //
  assert(net.springIndexA.size() == net.springIndexB.size());
  for (size_t i = 0; i < net.springIndexA.size(); ++i) {
#ifndef NDEBUG
    RUNTIME_EXP_IFN(
      net.springIndexA[i] != linkIdx,
      "Exected link " + std::to_string(linkIdx) +
      " to not be linked anywhere anymore, found in partial spring " +
      std::to_string(i) + ".");
#endif
    if (net.springIndexA[i] > linkIdx) {
      net.springIndexA[i] -= 1;
      net.springCoordinateIndexA[3 * i] -= 3;
      net.springCoordinateIndexA[3 * i + 1] -= 3;
      net.springCoordinateIndexA[3 * i + 2] -= 3;
    }
#ifndef NDEBUG
    RUNTIME_EXP_IFN(
      net.springIndexB[i] != linkIdx,
      "Exected link " + std::to_string(linkIdx) +
      " to not be linked anywhere anymore, found in partial spring " +
      std::to_string(i) + ".");
#endif
    if (net.springIndexB[i] > linkIdx) {
      net.springIndexB[i] -= 1;
      net.springCoordinateIndexB[3 * i] -= 3;
      net.springCoordinateIndexB[3 * i + 1] -= 3;
      net.springCoordinateIndexB[3 * i + 2] -= 3;
    }
  }
}

/**
 * @brief Merge two springs around a given cross-link
 *
 * This does not require the resulting network to be valid.
 *
 * @param net
 * @param u
 * @param removedPartialSpringIdx
 * @param keptPartialSpringIdx
 * @param linkToReduce
 * @param skipEigenResize
 */
void
MEHPForceBalance2::mergeSprings(ForceBalance2Network& net,
                                const Eigen::VectorXd& u,
                                const size_t removedPartialSpringIdx,
                                const size_t keptPartialSpringIdx,
                                const size_t linkToReduce,
                                bool skipEigenResize) const
{
  INVALIDARG_EXP_IFN(net.linkIsEntanglement[linkToReduce],
                     "The link to reduce must be a slip-link");
  INVALIDARG_EXP_IFN(keptPartialSpringIdx != removedPartialSpringIdx,
                     "Cannot merge one spring with the same one.");
  INVALIDARG_EXP_IFN(
    net.strandIdxOfSpring[keptPartialSpringIdx] ==
    net.strandIdxOfSpring[removedPartialSpringIdx],
    "The partial springs must be part of the same spring to merge them.");
  INVALIDARG_EXP_IFN(
    (net.springIndexA[keptPartialSpringIdx] == linkToReduce &&
      net.springIndexB[removedPartialSpringIdx] == linkToReduce) ||
    (net.springIndexB[keptPartialSpringIdx] == linkToReduce &&
      net.springIndexA[removedPartialSpringIdx] == linkToReduce),
    "Link to reduce must be part of the partial springs "
    "that are to be removed and kept, got " +
    std::to_string(net.springIndexA[keptPartialSpringIdx]) + " and " +
    std::to_string(net.springIndexB[keptPartialSpringIdx]) +
    " in kept spring, and got " +
    std::to_string(net.springIndexA[removedPartialSpringIdx]) + " and " +
    std::to_string(net.springIndexB[removedPartialSpringIdx]) +
    "in removed spring, instead of " + std::to_string(linkToReduce) + ".");
#ifdef DEBUG_REMOVAL
  std::cout << "Merging partial springs " << removedPartialSpringIdx << " and "
            << keptPartialSpringIdx << " around " << linkToReduce << std::endl;
#endif

  Eigen::Vector3d distanceBefore =
    this->evaluateSpringVector(
      net,
      u,
      removedPartialSpringIdx,
      this->is2D,
      false) +
    this->evaluateSpringVector(net, u, keptPartialSpringIdx, this->is2D, false);
  size_t fullSpringIdx = net.strandIdxOfSpring[keptPartialSpringIdx];
  // start with removal
  net.nrOfSprings -= 1;
  // tell the kept one their new end
  // NOTE: if is possible, if the removedPartialSpring is a primary loop,
  // that this procedure is ambiguous.
  bool removedIsA = net.springIndexA[removedPartialSpringIdx] == linkToReduce;
  size_t newEnd = removedIsA
                    ? net.springIndexB[removedPartialSpringIdx]
                    : net.springIndexA[removedPartialSpringIdx];
  if (removedIsA) {
    assert(net.springIndexB[keptPartialSpringIdx] == linkToReduce);
    net.springIndexB[keptPartialSpringIdx] = newEnd;
    for (size_t dir = 0; dir < 3; ++dir) {
      net.springCoordinateIndexB[3 * keptPartialSpringIdx + dir] =
        3 * newEnd + dir;
    }
  } else {
    assert(net.springIndexA[keptPartialSpringIdx] == linkToReduce);
    net.springIndexA[keptPartialSpringIdx] = newEnd;
    for (size_t dir = 0; dir < 3; ++dir) {
      net.springCoordinateIndexA[3 * keptPartialSpringIdx + dir] =
        3 * newEnd + dir;
    }
  }
  net.springBoxOffset.segment(3 * keptPartialSpringIdx, 3) +=
    net.springBoxOffset.segment(3 * removedPartialSpringIdx, 3);
  // remove the spring from the link
  // NOTE: currently, we allow it not to be present,
  // as it might be removed earlier already
  // It is anyway the case, that this function does not necessarily
  // keep the network valid
  int found = std::erase(net.strandIndicesOfLink[linkToReduce], fullSpringIdx);
  // TODO: check the origin of this assertion
  // assert(found == 1 || found == 0);
  RUNTIME_EXP_IFN(
    net.springIndicesOfStrand[fullSpringIdx].size() ==
    net.linkIndicesOfStrand[fullSpringIdx].size() - 1,
    "Require a global index for each local one, got " +
    std::to_string(net.springIndicesOfStrand[fullSpringIdx].size()) + " != " +
    std::to_string(net.linkIndicesOfStrand[fullSpringIdx].size() - 1) +
    " for spring " + std::to_string(fullSpringIdx) + ".");
  found = 0;
  // tell the spring of the removed link
  int removed = 0;
  for (int j = net.linkIndicesOfStrand[fullSpringIdx].size() - 1; j >= 0; --j) {
    if (net.linkIndicesOfStrand[fullSpringIdx][j] == linkToReduce) {
      if (removed == 0) {
        if ((j > 0 && net.springIndicesOfStrand[fullSpringIdx][j - 1] ==
             removedPartialSpringIdx) ||
            (j < net.springIndicesOfStrand[fullSpringIdx].size() &&
             net.springIndicesOfStrand[fullSpringIdx][j] ==
             removedPartialSpringIdx)) {
          net.linkIndicesOfStrand[fullSpringIdx].erase(
            net.linkIndicesOfStrand[fullSpringIdx].begin() + j);
          removed += 1;
        }
      }
      if (found == 1) {
        // we are dealing with a double -> re-add
        net.strandIndicesOfLink[linkToReduce].push_back(fullSpringIdx);
      } else if (found > 1 && removed > 0) {
        break; // required for certain cases... dangerous, somewhat.
      }
      found += 1;
    }
  }
  this->removeDuplicateListedStrandsFromLink(net, linkToReduce);
  assert(found >= 1 && removed == 1);
  found = 0;
  for (int j = net.springIndicesOfStrand[fullSpringIdx].size() - 1; j >= 0;
       --j) {
    if (net.springIndicesOfStrand[fullSpringIdx][j] ==
        removedPartialSpringIdx) {
      net.springIndicesOfStrand[fullSpringIdx].erase(
        net.springIndicesOfStrand[fullSpringIdx].begin() + j);
      found += 1;
    }
  }
  assert(found == 1);
  RUNTIME_EXP_IFN(
    net.springIndicesOfStrand[fullSpringIdx].size() ==
    net.linkIndicesOfStrand[fullSpringIdx].size() - 1,
    "Require a global index for each local one, got " +
    std::to_string(net.springIndicesOfStrand[fullSpringIdx].size()) + " != " +
    std::to_string(net.linkIndicesOfStrand[fullSpringIdx].size() - 1) +
    " for spring " + std::to_string(fullSpringIdx) + ".");
  // actually remove the rows
  pylimer_tools::utils::removeRow(
    net.strandIdxOfSpring,
    removedPartialSpringIdx,
    skipEigenResize);
  pylimer_tools::utils::removeRow(
    net.springIndexA,
    removedPartialSpringIdx,
    skipEigenResize);
  pylimer_tools::utils::removeRow(
    net.springIndexB,
    removedPartialSpringIdx,
    skipEigenResize);
  pylimer_tools::utils::removeRows(net.springCoordinateIndexA,
                                   3 * removedPartialSpringIdx,
                                   3,
                                   skipEigenResize);
  pylimer_tools::utils::removeRows(net.springCoordinateIndexB,
                                   3 * removedPartialSpringIdx,
                                   3,
                                   skipEigenResize);
  pylimer_tools::utils::removeRows(
    net.springBoxOffset,
    3 * removedPartialSpringIdx,
    3,
    skipEigenResize);
  // renumber stuff
  for (size_t loopSpringIdx = 0;
       loopSpringIdx < net.springIndicesOfStrand.size();
       ++loopSpringIdx) {
    for (size_t i = 0; i < net.springIndicesOfStrand[loopSpringIdx].size();
         ++i) {
      if (net.springIndicesOfStrand[loopSpringIdx][i] >
          removedPartialSpringIdx) {
        net.springIndicesOfStrand[loopSpringIdx][i] -= 1;
      }
    }
  }

  // validation
  size_t newSpringIdx =
    keptPartialSpringIdx +
    (keptPartialSpringIdx > removedPartialSpringIdx ? -1 : 0);
  Eigen::Vector3d newDistance =
    this->evaluateSpringVector(net, u, newSpringIdx, this->is2D, false);
  RUNTIME_EXP_IFN(pylimer_tools::utils::vector_approx_equal(
                    newDistance, distanceBefore, 1e-5),
                  "After merging two partial springs, the overall distance "
                  "is not consistent. Expected distance " +
                  std::to_string(distanceBefore) + ", but got " +
                  std::to_string(newDistance) + " for spring " +
                  std::to_string(newSpringIdx) + ".");
}

/**
 * @brief Merge two springs around a given cross-link
 *
 * @param net
 * @param u the current displacements
 * @param removedStrandIdx
 * @param keptStrandIdx
 * @param linkToReduce the index of the link to remove (combine the springs
 * around)
 */
void
MEHPForceBalance2::mergeStrands(ForceBalance2Network& net,
                                const Eigen::VectorXd& u,
                                const size_t removedStrandIdx,
                                const size_t keptStrandIdx,
                                const size_t linkToReduce) const
{
  INVALIDARG_EXP_IFN(removedStrandIdx < net.nrOfStrands &&
                     keptStrandIdx < net.nrOfStrands,
                     "Only full springs can be merged.");
  INVALIDARG_EXP_IFN(!net.linkIsEntanglement[linkToReduce],
                     "The link to reduce must be a cross-link");
  INVALIDARG_EXP_IFN(keptStrandIdx != removedStrandIdx,
                     "Cannot replace one spring with the same one.");

  // handle links
  const std::vector<size_t> removedSpringsLinks =
    net.linkIndicesOfStrand[removedStrandIdx];
  const std::vector<size_t> keptSpringsLinks =
    net.linkIndicesOfStrand[keptStrandIdx];

  const size_t removedPartialSpringIdx =
    (removedSpringsLinks[removedSpringsLinks.size() - 1] == linkToReduce)
      ? pylimer_tools::utils::last(net.springIndicesOfStrand[removedStrandIdx])
      : net.springIndicesOfStrand[removedStrandIdx][0];
  const size_t remainingPartialSpringIdx =
    (keptSpringsLinks[keptSpringsLinks.size() - 1] == linkToReduce)
      ? pylimer_tools::utils::last(net.springIndicesOfStrand[keptStrandIdx])
      : net.springIndicesOfStrand[keptStrandIdx][0];

  RUNTIME_EXP_IFN(net.springIndexA[removedPartialSpringIdx] == linkToReduce ||
                  net.springIndexB[removedPartialSpringIdx] == linkToReduce,
                  "Did not detect correct partial springs");
  RUNTIME_EXP_IFN(net.springIndexA[remainingPartialSpringIdx] == linkToReduce ||
                  net.springIndexB[remainingPartialSpringIdx] == linkToReduce,
                  "Did not detect correct partial springs");

  Eigen::Vector3d distanceBefore = this->evaluateSpringVector(
    net,
    u,
    removedPartialSpringIdx,
    this->is2D,
    false);
  Eigen::Vector3d distanceBeforeRemainingSpring = this->evaluateSpringVector(
    net,
    u,
    remainingPartialSpringIdx,
    this->is2D,
    false);

  net.nrOfStrands -= 1;
  net.nrOfSprings -= 1;

  net.linkIndicesOfStrand[keptStrandIdx].reserve(
    keptSpringsLinks.size() + removedSpringsLinks.size() - 2);
  net.springIndicesOfStrand[keptStrandIdx].reserve(
    keptSpringsLinks.size() + removedSpringsLinks.size() - 2);
  RUNTIME_EXP_IFN(net.springIndicesOfStrand[keptStrandIdx].size() ==
                  net.linkIndicesOfStrand[keptStrandIdx].size() - 1,
                  "Invalid sizes when merging springs");
  // tell the partial springs their new full spring
  for (size_t partialSpringIndex :
       net.springIndicesOfStrand[removedStrandIdx]) {
    net.strandIdxOfSpring[partialSpringIndex] = keptStrandIdx;
  }
  // std::cout << "Kept spring is "
  //           << pylimer_tools::utils::join(keptSpringsLinks.begin(),
  //                                         keptSpringsLinks.end(),
  //                                         std::string(", "))
  //           << std::endl;
  // actually merge the springs
  if (keptSpringsLinks[keptSpringsLinks.size() - 1] == linkToReduce) {
    // add to end...
    if (removedSpringsLinks[removedSpringsLinks.size() - 1] == linkToReduce) {
      // std::cout << "End end" << std::endl;
      // ...from end
      net.linkIndicesOfStrand[keptStrandIdx][keptSpringsLinks.size() - 1] =
        removedSpringsLinks[removedSpringsLinks.size() - 2];
      for (size_t i = 3; i <= removedSpringsLinks.size(); ++i) {
        net.linkIndicesOfStrand[keptStrandIdx].push_back(
          removedSpringsLinks[removedSpringsLinks.size() - i]);
      }
      for (int i = net.springIndicesOfStrand[removedStrandIdx].size() - 2;
           i >= 0;
           --i) {
        net.springIndicesOfStrand[keptStrandIdx].push_back(
          net.springIndicesOfStrand[removedStrandIdx][i]);
      }
      // invert the direction of these transferred partial springs
      for (size_t partialSpringIdxToInvert :
           net.springIndicesOfStrand[removedStrandIdx]) {
        if (partialSpringIdxToInvert == removedPartialSpringIdx) {
          continue;
        }
        std::swap(net.springIndexA[partialSpringIdxToInvert],
                  net.springIndexB[partialSpringIdxToInvert]);
        for (int dir = 0; dir < 3; ++dir) {
          std::swap(
            net.springCoordinateIndexA[3 * partialSpringIdxToInvert + dir],
            net.springCoordinateIndexB[3 * partialSpringIdxToInvert + dir]);
        }
        net.springBoxOffset.segment(3 * partialSpringIdxToInvert, 3) *= -1.;
      }
      distanceBefore -= distanceBeforeRemainingSpring;
      distanceBefore *= -1.;
    } else {
      // ...from start
      // std::cout << "End start" << std::endl;
      RUNTIME_EXP_IFN(removedSpringsLinks[0] == linkToReduce,
                      "Things don't make sense anymore.");
      net.linkIndicesOfStrand[keptStrandIdx][keptSpringsLinks.size() - 1] =
        removedSpringsLinks[1];
      for (size_t i = 2; i < removedSpringsLinks.size(); ++i) {
        net.linkIndicesOfStrand[keptStrandIdx].push_back(
          removedSpringsLinks[i]);
      }
      for (size_t i = 1; i < net.springIndicesOfStrand[removedStrandIdx].size();
           ++i) {
        net.springIndicesOfStrand[keptStrandIdx].push_back(
          net.springIndicesOfStrand[removedStrandIdx][i]);
      }
      distanceBefore += distanceBeforeRemainingSpring;
    }
  } else {
    RUNTIME_EXP_IFN(keptSpringsLinks[0] == linkToReduce, "How could this be?");
    // add to start...
    if (removedSpringsLinks[removedSpringsLinks.size() - 1] == linkToReduce) {
      // std::cout << "Start end" << std::endl;
      // from end
      net.linkIndicesOfStrand[keptStrandIdx][0] =
        removedSpringsLinks[removedSpringsLinks.size() - 2];
      for (size_t i = 3; i <= removedSpringsLinks.size(); ++i) {
        net.linkIndicesOfStrand[keptStrandIdx].insert(
          net.linkIndicesOfStrand[keptStrandIdx].begin(),
          removedSpringsLinks[removedSpringsLinks.size() - i]);
      }
      for (int i = net.springIndicesOfStrand[removedStrandIdx].size() - 2;
           i >= 0;
           --i) {
        net.springIndicesOfStrand[keptStrandIdx].insert(
          net.springIndicesOfStrand[keptStrandIdx].begin(),
          net.springIndicesOfStrand[removedStrandIdx][i]);
      }
      distanceBefore += distanceBeforeRemainingSpring;
    } else {
      // std::cout << "Start start" << std::endl;
      // from start
      RUNTIME_EXP_IFN(removedSpringsLinks[0] == linkToReduce,
                      "No way this exception is ever shown, right?");
      net.linkIndicesOfStrand[keptStrandIdx][0] = removedSpringsLinks[1];
      // have to insert it reverse order
      // happens automatically if we always insert the next the start
      for (size_t i = 2; i < removedSpringsLinks.size(); ++i) {
        net.linkIndicesOfStrand[keptStrandIdx].insert(
          net.linkIndicesOfStrand[keptStrandIdx].begin(),
          removedSpringsLinks[i]);
      }
      // skip the first (removed) partial spring
      for (size_t i = 1; i < net.springIndicesOfStrand[removedStrandIdx].size();
           ++i) {
        net.springIndicesOfStrand[keptStrandIdx].insert(
          net.springIndicesOfStrand[keptStrandIdx].begin(),
          net.springIndicesOfStrand[removedStrandIdx][i]);
      }

      // invert the direction of these transferred partial springs
      for (size_t partialSpringIdxToInvert :
           net.springIndicesOfStrand[removedStrandIdx]) {
        if (partialSpringIdxToInvert == removedPartialSpringIdx) {
          continue;
        }
        std::swap(net.springIndexA[partialSpringIdxToInvert],
                  net.springIndexB[partialSpringIdxToInvert]);
        for (int dir = 0; dir < 3; ++dir) {
          std::swap(
            net.springCoordinateIndexA[3 * partialSpringIdxToInvert + dir],
            net.springCoordinateIndexB[3 * partialSpringIdxToInvert + dir]);
        }
        net.springBoxOffset.segment(3 * partialSpringIdxToInvert, 3) *= -1.;
      }
      distanceBefore -= distanceBeforeRemainingSpring;
      distanceBefore *= -1.;
    }
  }
  RUNTIME_EXP_IFN(
    std::find(net.linkIndicesOfStrand[keptStrandIdx].begin(),
      net.linkIndicesOfStrand[keptStrandIdx].end(),
      linkToReduce) == net.linkIndicesOfStrand[keptStrandIdx].end(),
    "Link " + std::to_string(linkToReduce) +
    " to reduce should not be in the kept links anymore, found " +
    pylimer_tools::utils::join(net.linkIndicesOfStrand[keptStrandIdx].begin(),
      net.linkIndicesOfStrand[keptStrandIdx].end(),
      std::string(", ")) +
    ".");
  assert(net.springIndicesOfStrand[keptStrandIdx].size() ==
    net.linkIndicesOfStrand[keptStrandIdx].size() - 1);
  assert(net.linkIndicesOfStrand[keptStrandIdx].size() ==
    keptSpringsLinks.size() + removedSpringsLinks.size() - 2);

  // tell the links of their new spring index
  for (size_t linkOfRemovedSpring : removedSpringsLinks) {
    for (size_t i = 0; i < net.strandIndicesOfLink[linkOfRemovedSpring].size();
         ++i) {
      if (net.strandIndicesOfLink[linkOfRemovedSpring][i] == removedStrandIdx) {
        net.strandIndicesOfLink[linkOfRemovedSpring][i] = keptStrandIdx;
      }
    }
    this->removeDuplicateListedStrandsFromLink(net, linkOfRemovedSpring, true);
  }

  for (int i = net.strandIndicesOfLink[linkToReduce].size() - 1; i >= 0; --i) {
    if (net.strandIndicesOfLink[linkToReduce][i] == removedStrandIdx ||
        net.strandIndicesOfLink[linkToReduce][i] == keptStrandIdx) {
      net.strandIndicesOfLink[linkToReduce].erase(
        net.strandIndicesOfLink[linkToReduce].begin() + i);
    }
  }
  this->removeDuplicateListedStrandsFromLink(net, linkToReduce);

  net.linkIndicesOfStrand.erase(net.linkIndicesOfStrand.begin() +
                                removedStrandIdx);
  // partial springs

  bool removedIsA = net.springIndexA[removedPartialSpringIdx] == linkToReduce;
  size_t otherEndOfRemovedSpring =
    removedIsA
      ? net.springIndexB[removedPartialSpringIdx]
      : net.springIndexA[removedPartialSpringIdx];
  double offsetMultiplier = removedIsA ? -1. : 1.;
  if (net.springIndexA[remainingPartialSpringIdx] == linkToReduce) {
    net.springIndexA[remainingPartialSpringIdx] = otherEndOfRemovedSpring;
    for (size_t dir = 0; dir < 3; ++dir) {
      net.springCoordinateIndexA[3 * remainingPartialSpringIdx + dir] =
        3 * otherEndOfRemovedSpring + dir;
    };
  } else {
    RUNTIME_EXP_IFN(net.springIndexB[remainingPartialSpringIdx] == linkToReduce,
                    "");
    net.springIndexB[remainingPartialSpringIdx] = otherEndOfRemovedSpring;
    for (size_t dir = 0; dir < 3; ++dir) {
      net.springCoordinateIndexB[3 * remainingPartialSpringIdx + dir] =
        3 * otherEndOfRemovedSpring + dir;
    }
    offsetMultiplier *= -1.;
  }
  net.springBoxOffset.segment(3 * remainingPartialSpringIdx, 3) +=
    offsetMultiplier *
    net.springBoxOffset.segment(3 * removedPartialSpringIdx, 3);
  pylimer_tools::utils::removeRow(net.springIndexA, removedPartialSpringIdx);
  pylimer_tools::utils::removeRow(net.springIndexB, removedPartialSpringIdx);
  pylimer_tools::utils::removeRows(
    net.springCoordinateIndexA,
    3 * removedPartialSpringIdx,
    3);
  pylimer_tools::utils::removeRows(
    net.springCoordinateIndexB,
    3 * removedPartialSpringIdx,
    3);
  pylimer_tools::utils::removeRows(
    net.springBoxOffset,
    3 * removedPartialSpringIdx,
    3);

  pylimer_tools::utils::removeRow(net.strandIdxOfSpring,
                                  removedPartialSpringIdx);
  net.springIndicesOfStrand.erase(net.springIndicesOfStrand.begin() +
                                  removedStrandIdx);
  // renumber the remaining springs
  for (size_t i = 0; i < net.strandIndicesOfLink.size(); ++i) {
    for (size_t j = 0; j < net.strandIndicesOfLink[i].size(); ++j) {
      RUNTIME_EXP_IFN(net.strandIndicesOfLink[i][j] != removedStrandIdx,
                      "Removed spring found in spring indices of link " +
                      std::to_string(i) + ". Must not happen.");
      if (net.strandIndicesOfLink[i][j] > removedStrandIdx) {
        net.strandIndicesOfLink[i][j] -= 1;
      }
    }
  }

  // and the partial springs
  for (size_t i = 0; i < net.strandIdxOfSpring.size(); ++i) {
    RUNTIME_EXP_IFN(net.strandIdxOfSpring[i] != removedStrandIdx, "");
    if (net.strandIdxOfSpring[i] > removedStrandIdx) {
      net.strandIdxOfSpring[i] -= 1;
    }
  }

  for (size_t i = 0; i < net.springIndicesOfStrand.size(); ++i) {
    for (size_t j = 0; j < net.springIndicesOfStrand[i].size(); ++j) {
      RUNTIME_EXP_IFN(
        net.springIndicesOfStrand[i][j] != removedPartialSpringIdx,
        "");
      if (net.springIndicesOfStrand[i][j] > removedPartialSpringIdx) {
        net.springIndicesOfStrand[i][j] -= 1;
      }
    }
  }

  // handle contour lengths
  net.springContourLength[keptStrandIdx] +=
    net.springContourLength[removedStrandIdx];
  pylimer_tools::utils::removeRow(net.springContourLength, removedStrandIdx);
  pylimer_tools::utils::removeRow(net.springIsEntanglement, removedStrandIdx);
  RUNTIME_EXP_IFN(net.springContourLength.size() == net.nrOfStrands, "");

  size_t newKeptSpringIdx =
    (keptStrandIdx < removedStrandIdx) ? keptStrandIdx : (keptStrandIdx - 1);

  // std::cout << "Removed springs around " << linkToReduce << " with
  // spring
  // "
  //           << removedSpringIdx << " and partial "
  //           << removedPartialSpringIdx << ", keeping " << keptSpringIdx
  //           << " and " << remainingPartialSpringIdx << std::endl;
  // std::cout << "Spring partitions sum to " << springPartitions.sum()
  //           << " for " << net.nrOfSprings
  //           << " springs, contour length before was " <<
  //           contourLengthBefore
  //           << " and is now " <<
  //           net.springsContourLength[newKeptSpringIdx]
  //           << std::endl;

  // validation
  size_t newPartialSpringIdx =
    remainingPartialSpringIdx +
    (remainingPartialSpringIdx > removedPartialSpringIdx ? -1 : 0);
  Eigen::Vector3d newDistance =
    this->evaluateSpringVector(net, u, newPartialSpringIdx, this->is2D, false);
  RUNTIME_EXP_IFN(pylimer_tools::utils::vector_approx_equal(
                    newDistance, distanceBefore, 1e-5),
                  "After merging two springs, the overall distance "
                  "is not consistent. Expected distance " +
                  std::to_string(distanceBefore) + ", but got " +
                  std::to_string(newDistance) + " for spring " +
                  std::to_string(newPartialSpringIdx) + ".");

#ifndef NDEBUG
  // check that the ordering is correct
  for (size_t i = 0; i < net.springIndicesOfStrand[newKeptSpringIdx].size();
       ++i) {
    size_t partialSpringIdx = net.springIndicesOfStrand[newKeptSpringIdx][i];
    size_t endA = net.springIndexA[partialSpringIdx];
    size_t endB = net.springIndexB[partialSpringIdx];
    std::vector<size_t> linkIndices = net.linkIndicesOfStrand[newKeptSpringIdx];
    assert(endA == linkIndices[i]);
    assert(endB == linkIndices[i + 1]);
  }
#endif
}

/**
 * @brief Replace the two springs traversing a two-functional crosslinkers
 * with a single spring
 *
 * Also handles entanglement beads
 *
 * @param net
 * @param displacements
 */
size_t
MEHPForceBalance2::removeBifunctionalCrosslinks(
  ForceBalance2Network& net,
  Eigen::VectorXd& displacements) const
{
  size_t numRemoved = 0;
  for (long int crosslinkIdx = net.nrOfNodes - 1; crosslinkIdx >= 0;
       --crosslinkIdx) {
    if (net.strandIndicesOfLink[crosslinkIdx].size() == 2) {
      std::vector<size_t> springsToMerge =
        net.strandIndicesOfLink[crosslinkIdx];

      assert(springsToMerge.size() == 2);

      // second primary loop check for slip-link entanglements
      // check that it's not a primary loop in any way:
      if (springsToMerge[0] != springsToMerge[1] &&
          (XOR(net.linkIndicesOfStrand[springsToMerge[0]][0] == crosslinkIdx,
               pylimer_tools::utils::last(
                 net.linkIndicesOfStrand[springsToMerge[0]]) ==
               crosslinkIdx)) &&
          (XOR(net.linkIndicesOfStrand[springsToMerge[1]][0] == crosslinkIdx,
               pylimer_tools::utils::last(
                 net.linkIndicesOfStrand[springsToMerge[1]]) ==
               crosslinkIdx))) {
#ifdef DEBUG_REMOVAL
        std::cout << "Merging springs " << springsToMerge[0] << " and "
                  << springsToMerge[1] << " around " << crosslinkIdx
                  << std::endl;
#endif

        // let's remove this
        // TODO: this is inefficient shit, so much data being moved
        this->mergeStrands(net,
                           displacements,

                           springsToMerge[0],
                           springsToMerge[1],
                           crosslinkIdx);

        // RUNTIME_EXP_IFN(this->validateNetwork(net, displacements), "Invalid
        // internal network representation"); std::cout << "Removing link " <<
        // crosslinkIdx << std::endl;
        this->removeLink(net, displacements, crosslinkIdx);

        // std::cout << "Removed cross-link " << crosslinkIdx << std::endl;

#ifndef NDEBUG
        assert(this->validateNetwork(net, displacements));
#endif
        numRemoved += 1;
      }
      // else: TODO: decide
    }
  }
#ifndef NDEBUG
  assert(this->validateNetwork(net, displacements));
#endif
  return numRemoved;
}

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
 * @param oneOverSpringPartitions the 1/contour length weights for the springs
 * @return double, the distance (squared norm) displaced
 */
double
MEHPForceBalance2::displaceToMeanPosition(
  const ForceBalance2Network& net,
  Eigen::VectorXd& u,
  const Eigen::ArrayXd& oneOverSpringPartitions) const
{
  assert(oneOverSpringPartitions.size() == net.nrOfSprings * 3);
  Eigen::ArrayXd objectiveDisplacement =
    Eigen::ArrayXd::Zero(3 * net.nrOfLinks);
  Eigen::ArrayXd partialSpringDistances =
    this->evaluateSpringVectors(net, u, this->is2D, this->assumeBoxLargeEnough)
    .array();
  objectiveDisplacement(net.springCoordinateIndexA) +=
    (oneOverSpringPartitions * partialSpringDistances);
  objectiveDisplacement(net.springCoordinateIndexB) -=
    (oneOverSpringPartitions * partialSpringDistances);

  Eigen::ArrayXd springPartWeightingFactor =
    Eigen::ArrayXd::Zero(net.nrOfLinks * 3);
  Eigen::ArrayXd loopPartialSpringEliminator =
    (net.springCoordinateIndexA != net.springCoordinateIndexB).cast<double>();

  springPartWeightingFactor(net.springCoordinateIndexA) +=
    oneOverSpringPartitions * loopPartialSpringEliminator;
  springPartWeightingFactor(net.springCoordinateIndexB) +=
    oneOverSpringPartitions * loopPartialSpringEliminator;
  springPartWeightingFactor = springPartWeightingFactor.unaryExpr(
    [](double v) { return v > 0. ? v : 1.0; });
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
    nSpringsPerLink.unaryExpr([](double v) { return v > 0. ? v : 1.0; });
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
  Eigen::VectorXd finalDisplacement =
    (remainingDisplacement + backForthDisplacement).matrix();
  RUNTIME_EXP_IFN(
    pylimer_tools::utils::all_components_finite(finalDisplacement),
    "Some displacements are not finite");
  // this->box.handlePBC(finalDisplacement);
  u += finalDisplacement;

  double max_disp =
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

  const std::unordered_set<size_t> partialSpringIndices =
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
      objectiveDisplacement += (partialDistance) * oneOverContourLengthFraction;
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
  if (!this->assumeBoxLargeEnough) {
    const Eigen::Vector3d forceAfter = this->getForceOn(net, u, linkIdx);

    // this is only true if we don't have "full" PBC
    assert((pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
      forceAfter, Eigen::Vector3d::Zero(), 0.01)));
    if (!pylimer_tools::utils::vector_approx_equal<Eigen::Vector3d>(
      forceBefore,
      Eigen::Vector3d::Zero(),
      0.01)) {
      assert(forceBefore.squaredNorm() >= forceAfter.squaredNorm());
    }
  }
#endif

  double dist = (objectiveDisplacement * denominator).squaredNorm();
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
 * @param linkIdx
 * @param net
 * @param u
 * @param debugNrSpringsVisited
 * @return Eigen::Matrix3d
 */
Eigen::Matrix3d
MEHPForceBalance2::evaluateStressOnLink(
  const size_t linkIdx,
  const ForceBalance2Network& net,
  const Eigen::VectorXd& u,
  Eigen::VectorXi& debugNrSpringsVisited) const
{
  std::vector<size_t> springIndices = net.strandIndicesOfLink[linkIdx];
  Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();

  std::unordered_set<size_t> partialSpringIndices =
    this->getPartialSpringIndicesOfLink(net, linkIdx);

  for (const size_t globalSpringIndex : partialSpringIndices) {
    Eigen::Vector3d partialDistance =
      this->evaluateSpringVectorFrom(net, u, globalSpringIndex, linkIdx);
    double oneOverContourLengthFraction =
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
 * @param linkIdx
 * @param net
 * @param u
 * @param debugNrSpringsVisited
 * @return Eigen::Vector3d
 */
Eigen::Vector3d
MEHPForceBalance2::evaluateForceOnLink(
  const size_t linkIdx,
  const ForceBalance2Network& net,
  const Eigen::VectorXd& u,

  Eigen::VectorXi& debugNrSpringsVisited) const
{
  Eigen::Vector3d force = Eigen::Vector3d::Zero();

  std::unordered_set<size_t> partialSpringIndices =
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
    double oneOverContourLengthFraction =
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
  for (size_t i = this->initialConfig.nrOfNodes;
       i < this->initialConfig.nrOfLinks;
       ++i) {
    if (this->initialConfig.strandIndicesOfLink[i].size() < 2) {
      result += 1;
    }
    if (this->initialConfig.strandIndicesOfLink[i].size() == 2 &&
        this->initialConfig.strandIndicesOfLink[i][0] ==
        this->initialConfig.strandIndicesOfLink[i][1]) {
      result += 1;
    }
  }

  return result;
};

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
                                         const bool is2D,
                                         const bool assumeLarge) const
{
  // first, the distances
  assert(u.size() == net.coordinates.size());

  Eigen::VectorXd displacedCoords = net.coordinates + u;
  Eigen::VectorXd partialDistances =
    (displacedCoords(net.springCoordinateIndexB) -
     displacedCoords(net.springCoordinateIndexA)) +
    net.springBoxOffset;

  if (assumeLarge) {
    this->box.handlePBC(partialDistances);
  }

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
  std::vector<int> types = pylimer_tools::utils::initializeWithValue(
    this->initialConfig.nrOfNodes,
    crossLinkerType);
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<int> zeros =
    pylimer_tools::utils::initializeWithValue(this->initialConfig.nrOfNodes, 0);
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
    // override type, since the types may be different from
    // crossLinkerType if converted with dangling chains
    types[i] = this->initialConfig.oldAtomTypes[i];
  }
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
MEHPForceBalance2::getAverageSpringLength() const
{
  Eigen::VectorXd partialSpringVectors = this->evaluateSpringVectors(
    this->initialConfig,
    this->currentDisplacements);

  Eigen::VectorXd springVectors =
    Eigen::VectorXd::Zero(3 * this->initialConfig.nrOfStrands);
  for (size_t i = 0; i < this->initialConfig.nrOfSprings; ++i) {
    springVectors.segment(this->initialConfig.strandIdxOfSpring[i], 3) +=
      partialSpringVectors.segment(3 * i, 3);
  }

  return pylimer_tools::utils::segmentwise_norm_mean(springVectors, 3);
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
};

Eigen::VectorXd
MEHPForceBalance2::evaluateStrandVectors(const ForceBalance2Network& net,
                                         Eigen::VectorXd u) const
{
  Eigen::VectorXd partialSpringVectors = this->evaluateSpringVectors(net, u);
  Eigen::VectorXd springVectors = Eigen::VectorXd::Zero(3 * net.nrOfStrands);
  for (size_t i = 0; i < net.nrOfSprings; ++i) {
    springVectors.segment(net.strandIdxOfSpring[i], 3) +=
      partialSpringVectors.segment(3 * i, 3);
  }

  return springVectors;
}

double
MEHPForceBalance2::sumToTotalFraction(const ForceBalance2Network& net,
                                      Eigen::VectorXd springPartition,
                                      size_t springIdx,
                                      size_t targetLink) const
{
  double alpha = 0.;
  for (size_t i = 0; i < net.springIndicesOfStrand[springIdx].size(); ++i) {
    size_t currentPartialSpringIdx = net.springIndicesOfStrand[springIdx][i];
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
  const size_t partialSpringIdx) const
{
  double denominator = 1. / net.springContourLength[partialSpringIdx];

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

  double halfOverVolume = 0.5 / (net.L[0] * net.L[1] * net.L[2]);

  Eigen::VectorXi debugNrSpringsVisited =
    Eigen::VectorXi::Zero(net.nrOfSprings);

  for (size_t linkIdx : linkIndices) {
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

  double halfOverVolume = 0.5 / (net.L[0] * net.L[1] * net.L[2]);

  Eigen::VectorXi debugNrSpringsVisited =
    Eigen::VectorXi::Zero(net.nrOfSprings);

  size_t nrOfLinksToInspect = crosslinksOnly ? net.nrOfNodes : net.nrOfLinks;
  for (size_t linkIdx = 0; linkIdx < nrOfLinksToInspect; ++linkIdx) {
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

  double oneOverVolume = 1. / (net.L[0] * net.L[1] * net.L[2]);

  Eigen::VectorXd displacedCoords = net.coordinates + u;
  Eigen::VectorXd relevantPartialDistancesA =
    (displacedCoords(net.springCoordinateIndexB) -
     displacedCoords(net.springCoordinateIndexA)) +
    net.springBoxOffset;

  if (this->assumeBoxLargeEnough) {
    this->box.handlePBC(relevantPartialDistancesA);
  }

  if (this->is2D) {
    for (size_t i = 2; i < relevantPartialDistancesA.size(); i += 3) {
      relevantPartialDistancesA[i] = 0.;
    }
  }

  for (size_t partialSpringIdx = 0; partialSpringIdx < net.nrOfSprings;
       ++partialSpringIdx) {
    Eigen::Vector3d distance =
      relevantPartialDistancesA.segment(3 * partialSpringIdx, 3);
    double oneOverContourLengthFraction =
      1. / net.springContourLength[partialSpringIdx];

    /* spring contribution to the overall stress tensor */
    for (size_t j = 0; j < 3; j++) {
      for (size_t k = 0; k < 3; k++) {
        double contribution = distance[j] * distance[k] * this->kappa *
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
  std::array<std::array<double, 3>, 3> res =
    this->evaluateStressTensorLinkBased(
      this->initialConfig,
      this->currentDisplacements,
      crosslinksOnly);
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
MEHPForceBalance2::getEffectiveFunctionalityOfAtoms(double tolerance) const
{
  std::unordered_map<long int, int> results;
  results.reserve(this->initialConfig.nrOfNodes);

  Eigen::VectorXi nrOfActiveSpringsConnected =
    this->getNrOfActiveStrandsConnected(tolerance);
  for (size_t i = 0; i < this->initialConfig.nrOfNodes; i++) {
    results.emplace(this->initialConfig.oldAtomIds[i],
                    nrOfActiveSpringsConnected[i]);
  }
  return results;
}

double
MEHPForceBalance2::computeDanglingWeightFraction(ForceBalance2Network* net,
                                                 const Eigen::VectorXd& u,
                                                 const double tolerance) const
{
  double activeWeightFraction =
    this->computeActiveWeightFraction(net, u, tolerance);
  RUNTIME_EXP_IFN(APPROX_WITHIN(activeWeightFraction, 0., 1., 1e-6),
                  "Expect active weight fraction to be between 0 and 1, got " +
                  std::to_string(activeWeightFraction) + ".");
  double solubleWeightFraction =
    this->computeSolubleWeightFraction(net, u, tolerance);
  RUNTIME_EXP_IFN(APPROX_WITHIN(solubleWeightFraction, 0., 1., 1e-6),
                  "Expect soluble weight fraction to be between 0 and 1, got " +
                  std::to_string(solubleWeightFraction) + ".");
  RUNTIME_EXP_IFN(
    APPROX_WITHIN(activeWeightFraction + solubleWeightFraction, 0., 1., 1e-6),
    "Expect active and soluble weight fraction to add up to maximum 1, "
    "got " +
    std::to_string(activeWeightFraction + solubleWeightFraction) + ".");

  // finally, normalize by the number of atoms.
  // TODO: currently, the weight of the atoms is ignored
  return 1. - activeWeightFraction - solubleWeightFraction;
}

double
MEHPForceBalance2::computeActiveWeightFraction(ForceBalance2Network* net,
                                               const Eigen::VectorXd& u,
                                               const double tolerance) const
{
  INVALIDARG_EXP_IFN(net->nrOfLinks * 3 == u.size(),
                     "Link displacements and network don't match");
  if (net->nrOfStrands < 1) {
    return 0.;
  }
  // find all active springs
  Eigen::ArrayXb activeSprings = this->findActiveStrands(net, u, tolerance);
  if (activeSprings.count() == 0) {
    return 0.;
  }
  // as of now, the springsContourLength is equal to the number of bonds
  // from cross-link to cross-link. therefore, the number of atoms of each
  // of these springs is one less
  Eigen::ArrayXd allActiveAtomsPerChains =
    activeSprings.cast<double>() *
    (net->springContourLength.array() - Eigen::ArrayXd::Ones(net->nrOfStrands));

  // TODO: currently, the weight of the atoms is ignored
  // normalize by the number of atoms
  return (allActiveAtomsPerChains.matrix().sum() +
          this->getNrOfActiveNodes(tolerance)) /
         (static_cast<double>(this->universe.getNrOfAtoms()));
}

std::pair<Eigen::ArrayXb, Eigen::ArrayXb>
MEHPForceBalance2::findClusteredToActive(const ForceBalance2Network* net,
                                         const Eigen::VectorXd& u,
                                         const double tolerance) const
{
  INVALIDARG_EXP_IFN(u.size() == net->nrOfLinks * 3, "Invalid sizes.");

  // find all active springs
  Eigen::ArrayXb springIsActive = this->findActiveStrands(net, u, tolerance);

  // then, iteratively walk along the springs to mark those as "active"
  // that are connected to active springs
  bool hadChanged = true;
  Eigen::ArrayXb nodeIsActive = Eigen::ArrayXb::Constant(net->nrOfNodes, false);
  while (hadChanged) {
    hadChanged = false;
    for (size_t i = 0; i < net->nrOfNodes; ++i) {
      if (nodeIsActive(i)) {
        continue;
      }
      for (size_t springIdx : net->strandIndicesOfLink[i]) {
        if (springIsActive[springIdx]) {
          hadChanged = true;
          nodeIsActive(i) = true;
          for (size_t innerSpringIdx : net->strandIndicesOfLink[i]) {
            springIsActive[innerSpringIdx] = true;
          }
          break;
        }
      }
    }
  }

  return std::make_pair(springIsActive, nodeIsActive);
}

double
MEHPForceBalance2::countActiveClusteredAtoms(ForceBalance2Network* net,
                                             const Eigen::VectorXd& u,
                                             const double tolerance) const
{
  INVALIDARG_EXP_IFN(net->nrOfLinks * 3 == u.size(),
                     "Link displacements and network don't match");
  if (net->nrOfStrands < 1) {
    return 0.;
  }

  std::vector<pylimer_tools::entities::Universe> clusters =
    this->universe.getClusters();
  std::vector<long int> atomIdxToClusterIdx(this->universe.getNrOfAtoms());
  for (size_t i = 0; i < clusters.size(); ++i) {
    for (const pylimer_tools::entities::Atom& atom : clusters[i].getAtoms()) {
      atomIdxToClusterIdx[this->universe.getIdxByAtomId(atom.getId())] = i;
    }
  }

  std::vector<bool> clusterIsActive(clusters.size(), false);

  // find active atoms
  std::vector<long int> activeNodeIndices =
    this->getIndicesOfActiveNodes(net, u, tolerance);

  for (const long int& nodeIdx : activeNodeIndices) {
    long int universeAtomIdx =
      this->universe.getIdxByAtomId(net->oldAtomIds[nodeIdx]);
    clusterIsActive[atomIdxToClusterIdx[universeAtomIdx]] = true;
  }

  double nClusteredAtoms = 0.;
  for (size_t i = 0; i < clusters.size(); ++i) {
    if (clusterIsActive[i]) {
      nClusteredAtoms += clusters[i].getNrOfAtoms();
    }
  }

  return nClusteredAtoms;
}

double
MEHPForceBalance2::computeSolubleWeightFraction(ForceBalance2Network* net,
                                                const Eigen::VectorXd& u,
                                                const double tolerance) const
{
  INVALIDARG_EXP_IFN(net->nrOfLinks * 3 == u.size(),
                     "Link displacements and network don't match");
  if (net->nrOfStrands < 1) {
    return 1.;
  }
  double nActiveClusteredAtoms =
    this->countActiveClusteredAtoms(net, u, tolerance);
  // finally, normalize by the number of atoms.
  // NOTE: currently, the weight of the atoms is ignored
  return 1. - (nActiveClusteredAtoms /
               (static_cast<double>(this->universe.getNrOfAtoms())));
}

/**
 * @brief Get the indices of active Nodes
 *
 * @param tolerance the tolerance: springs under a certain length are
 * considered inactive
 * @return std::vector<long int> the atom ids
 */
std::vector<long int>
MEHPForceBalance2::getIndicesOfActiveNodes(const ForceBalance2Network* net,
                                           const Eigen::VectorXd& u,

                                           double tolerance) const
{
  std::vector<long int> results;
  results.reserve(net->nrOfNodes);

  // find all active springs
  Eigen::ArrayXb springIsActive = this->findActiveStrands(net, u, tolerance);

  for (size_t i = 0; i < net->nrOfNodes; i++) {
    std::vector<size_t> springIndices = net->strandIndicesOfLink[i];
    for (const size_t springIndex : springIndices) {
      if (springIsActive[springIndex]) {
        results.push_back(i);
        break;
      }
    }
  }

  return results;
};

/**
 * @brief Get the atom ids of the active cross-links (not entanglement
 * beads/links)
 *
 * @param tolerance the tolerance: springs under a certain length are
 * considered inactive
 * @return std::vector<long int> the atom ids
 */
std::vector<long int>
MEHPForceBalance2::getIdsOfActiveNodes(double tolerance) const
{
  std::vector<long int> results;
  // find all active springs
  std::vector<long int> activeNodes = this->getIndicesOfActiveNodes(
    &this->initialConfig,
    this->currentDisplacements,
    tolerance);

  results.reserve(activeNodes.size());

  for (long int nodeIdx : activeNodes) {
    results.push_back(this->initialConfig.oldAtomIds[nodeIdx]);
  }

  return results;
}

std::vector<double>
MEHPForceBalance2::getOverallSpringLengths() const
{
  std::vector<double> partialSpringDistances = this->getCurrentSpringLengths();
  assert(partialSpringDistances.size() == this->initialConfig.nrOfSprings);
  std::vector<double> results =
    std::vector<double>(this->initialConfig.nrOfStrands, 0.);
  for (size_t i = 0; i < this->initialConfig.nrOfSprings; ++i) {
    results[this->initialConfig.strandIdxOfSpring[i]] +=
      partialSpringDistances[i];
  }

  return results;
}

Eigen::VectorXd
MEHPForceBalance2::getCurrentSpringDistances() const
{
  Eigen::VectorXd partialSpringVectors = this->evaluateSpringVectors(
    this->initialConfig,
    this->currentDisplacements);

  return partialSpringVectors;
}

std::vector<double>
MEHPForceBalance2::getCurrentSpringLengths() const
{
  Eigen::VectorXd vecs = this->evaluateSpringVectors(
    this->initialConfig,
    this->currentDisplacements);

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
  Eigen::VectorXi nrOfActiveSpringsConnected =
    Eigen::VectorXi::Zero(this->initialConfig.nrOfNodes);
  Eigen::ArrayXb springIsActive = this->findActiveStrands(tolerance);
  for (size_t i = 0; i < this->initialConfig.nrOfStrands; i++) {
    if (springIsActive[i]) {
      /* active spring */
      int a = this->initialConfig.linkIndicesOfStrand[i][0];
      int b =
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
  Eigen::VectorXi nrOfActivePartialSpringsConnected =
    Eigen::VectorXi::Zero(this->initialConfig.nrOfNodes);
  Eigen::ArrayXb partialSpringIsActive = this->findActiveSprings(tolerance);
  // translate this to the nodes
  for (size_t i = 0; i < this->initialConfig.nrOfSprings; ++i) {
    if (partialSpringIsActive[i]) {
      /* active spring */
      // size_t a =
      //   this->initialConfig
      //     .springIndexA[this->initialConfig.partialToFullSpringIndex[i]];
      // size_t b =
      //   this->initialConfig
      //     .springIndexB[this->initialConfig.partialToFullSpringIndex[i]];
      size_t a = this->initialConfig.springIndexA[i];
      size_t b = this->initialConfig.springIndexB[i];
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
MEHPForceBalance2::getGammaFactor(double b02, int nrOfChains) const
{
  if (b02 < 0) {
    b02 = this->defaultBondLength * this->defaultBondLength;
  }

  Eigen::VectorXd gammaFactors = this->getGammaFactors(b02);

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
MEHPForceBalance2::getGammaFactors(double b02) const
{
  Eigen::VectorXd springVectors = this->evaluateSpringVectors(
    this->initialConfig,
    this->currentDisplacements);

  Eigen::VectorXd gammaFactors(springVectors.size() / 3);
  const double commonDenominator = 1. / b02;
  for (size_t i = 0; i < springVectors.size() / 3; ++i) {
    double oneOverContourLengthFraction =
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
MEHPForceBalance2::getGammaFactorsInDir(double b02, int dir) const
{
  INVALIDARG_EXP_IFN(dir >= 0 && dir <= 2, "Invalid direction.");
  Eigen::VectorXd springVectors = this->evaluateSpringVectors(
    this->initialConfig,
    this->currentDisplacements);

  Eigen::VectorXd gammaFactors(springVectors.size() / 3);
  const double commonDenominator = 1. / b02;
  for (size_t i = 0; i < springVectors.size() / 3; ++i) {
    double oneOverContourLengthFraction =
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
 * @return double
 */
double
MEHPForceBalance2::getWeightedSpringLength(const ForceBalance2Network& net,
                                           const Eigen::VectorXd& u,
                                           size_t partialSpringIdx) const
{
  double oneOverContourLengthFraction =
    1. / net.springContourLength[partialSpringIdx];
  return this->evaluateSpringVector(net, u, partialSpringIdx).norm() *
         oneOverContourLengthFraction;
}

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
                  "Invalid size of springPartBoxOffset");
  RUNTIME_EXP_IFN(net.springIndexB.size() == net.nrOfSprings,
                  "Invalid size of springPartIndexB");
  RUNTIME_EXP_IFN(
    net.strandIdxOfSpring.size() == net.nrOfSprings,
    "Every partial spring must be able to map to the full spring.");

  /**
   * Test maximum values
   */
  if (net.nrOfStrands > 0) {
    RUNTIME_EXP_IFN(net.strandIdxOfSpring.maxCoeff() < net.nrOfStrands,
                    "Partial spring must map to full spring, which must have "
                    "a lower index.");
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
  for (size_t link_idx = 0; link_idx < net.nrOfLinks; ++link_idx) {
    std::vector<size_t> thisLinksSprings = net.strandIndicesOfLink[link_idx];
    std::sort(thisLinksSprings.begin(), thisLinksSprings.end());
    auto last = std::unique(thisLinksSprings.begin(), thisLinksSprings.end());
    RUNTIME_EXP_IFN(last == thisLinksSprings.end(),
                    "Expect each link to only have one back-link to the "
                    "springs, found back-links " +
                    pylimer_tools::utils::join(thisLinksSprings.begin(),
                      thisLinksSprings.end(),
                      std::string("_")) +
                    " for link " + std::to_string(link_idx) + ".");
    for (size_t spring_idx : thisLinksSprings) {
      std::vector<size_t> thisSpringsLinks =
        net.linkIndicesOfStrand[spring_idx];
      RUNTIME_EXP_IFN(std::find(thisSpringsLinks.begin(),
                        thisSpringsLinks.end(),
                        link_idx) != thisSpringsLinks.end(),
                      "Spring must have a connection to the link, too. Did "
                      "not find link " +
                      std::to_string(link_idx) + " in spring " +
                      std::to_string(spring_idx) + ".");
    }
  }

  /**
   * Test the assumptions on slip-links
   */
  for (size_t linkIdx = 0; linkIdx < net.nrOfLinks;
       ++linkIdx) {
    if (net.linkIsEntanglement[linkIdx]) {
      RUNTIME_EXP_IFN(
        net.strandIndicesOfLink[linkIdx].size() == 2 ||
        net.strandIndicesOfLink[linkIdx].size() == 1,
        "Expect each slip-link to be involved in exactly one or two "
        "springs, "
        "got " +
        std::to_string(net.strandIndicesOfLink[linkIdx].size()) + ".");
    }
  }

  /**
   * Test the validitiy of springs and their mapping
   */
  Eigen::ArrayXi nrOfMentions = Eigen::ArrayXi::Zero(net.nrOfLinks);
  for (size_t i = 0; i < net.nrOfStrands; ++i) {
    RUNTIME_EXP_IFN(net.linkIndicesOfStrand[i].size() >= 2,
                    "Each spring requires at least two links, got " +
                    std::to_string(net.linkIndicesOfStrand[i].size()) +
                    " at i = " + std::to_string(i) + ".");
    RUNTIME_EXP_IFN(net.springContourLength[i] > 0,
                    "Unexpected spring contour length, got " +
                    std::to_string(net.springContourLength[i]) +
                    " for spring " + std::to_string(i) + ".");
    RUNTIME_EXP_IFN(
      net.springIndicesOfStrand[i].size() ==
      net.linkIndicesOfStrand[i].size() - 1,
      "Require a global index for each local one, got " +
      std::to_string(net.springIndicesOfStrand[i].size()) +
      " != " + std::to_string(net.linkIndicesOfStrand[i].size() - 1) +
      " for spring " + std::to_string(i) + ".");
    for (size_t partialIdx = 0;
         partialIdx < net.springIndicesOfStrand[i].size();
         ++partialIdx) {
      size_t partialSpringIdx = net.springIndicesOfStrand[i][partialIdx];
      size_t partner0 = net.linkIndicesOfStrand[i][partialIdx];
      size_t partner1 = net.linkIndicesOfStrand[i][partialIdx + 1];
      RUNTIME_EXP_IFN(
        ((net.springIndexA[partialSpringIdx] == partner0 &&
          net.springIndexB[partialSpringIdx] == partner1)),
        "Expect linkIndicesOfSprings and localToGlobalSpringIndex "
        "ordering to correspond. Got partner0 = " +
        std::to_string(partner0) + ", partner1 = " +
        std::to_string(partner1) + " vs. springs part indices " +
        std::to_string(net.springIndexA[partialSpringIdx]) + " and " +
        std::to_string(net.springIndexB[partialSpringIdx]) + " in spring " +
        std::to_string(i) + " (partial: " + std::to_string(partialSpringIdx) +
        ") with global indices " +
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
    for (size_t j = 0; j < links.size(); ++j) {
      size_t link_idx = links[j];
      nrOfMentions[link_idx] += 1;
      RUNTIME_EXP_IFN(net.linkIsEntanglement[link_idx] ==
                      ((j != 0) && (j != (links.size() - 1))),
                      "Cross-links must be first and last in a spring, "
                      "slip-links in-between. Found discrepancy at " +
                      std::to_string(j) + "/" + std::to_string(links.size()) +
                      " in spring " + std::to_string(i) + ".")
      std::vector<size_t> thisLinksSprings = net.strandIndicesOfLink[link_idx];
      RUNTIME_EXP_IFN(
        std::find(thisLinksSprings.begin(), thisLinksSprings.end(), i) !=
        thisLinksSprings.end(),
        "Link must have a connection to the spring, too. Did not find "
        "spring " +
        std::to_string(i) + " in link " + std::to_string(link_idx) + ".");
    }
  }
  for (size_t i = net.nrOfNodes; i < net.nrOfLinks; ++i) {
    RUNTIME_EXP_IFN(nrOfMentions[i] == 2,
                    "Expect each slip-link to be mentioned twice in the "
                    "links-of-springs mapping, but " +
                    std::to_string(i) + " was mentioned " +
                    std::to_string(nrOfMentions[i]) + " times.");
  }

  /**
   * Test the validity of partial springs and their mapping
   */
  for (size_t i = 0; i < net.nrOfSprings; i++) {
    const size_t partialEndA = net.springIndexA[i];
    const size_t partialEndB = net.springIndexB[i];
    RUNTIME_EXP_IFN(partialEndA < net.nrOfLinks,
                    "Cannot have a spring (" + std::to_string(i) +
                    ") part larger " + std::to_string(partialEndA) +
                    " than the nr of links (" +
                    std::to_string(net.nrOfLinks) + ").");
    RUNTIME_EXP_IFN(partialEndB < net.nrOfLinks,
                    "Cannot have a spring (" + std::to_string(i) +
                    ") part larger " + std::to_string(partialEndB) +
                    " than the nr of links (" +
                    std::to_string(net.nrOfLinks) + ").");
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