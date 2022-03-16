#include "Universe.h"
#include "../utils/GraphUtils.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
#include "Box.h"

extern "C" {
#include <igraph/igraph.h>
}
#include <cassert>

#include <algorithm>
#include <iterator> // for back_inserter
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace pylimer_tools {
namespace entities {
Universe::Universe(Box box) {
  /* turn on attribute handling: TODO: move to some main() function  */
  igraph_set_attribute_table(&igraph_cattribute_table);
  this->box = box;

  // igraph_vector_t gtypes, vtypes, etypes;
  // igraph_strvector_t gnames, vnames, enames;

  // igraph_vector_init(&gtypes, 0);
  // igraph_vector_init(&vtypes, 0);
  // igraph_vector_init(&etypes, 0);
  // igraph_strvector_init(&gnames, 0);
  // igraph_strvector_init(&vnames, 0);
  // igraph_strvector_init(&enames, 0);

  // start setting properties
  igraph_empty(&this->graph, 0, IGRAPH_UNDIRECTED);

  //
  // igraph_cattribute_list(&this->graph, &gnames, &gtypes, &vnames, &vtypes,
  //                        &enames, &etypes);

  // // not sure if the above is really needed as we can destroy the vectors
  // here already without problems
  // /* Destroy */
  // igraph_vector_destroy(&gtypes);
  // igraph_vector_destroy(&vtypes);
  // igraph_vector_destroy(&etypes);
  // igraph_strvector_destroy(&gnames);
  // igraph_strvector_destroy(&vnames);
  // igraph_strvector_destroy(&enames);
}

Universe::Universe(const double Lx, const double Ly, const double Lz)
    : Universe(Box(Lx, Ly, Lz)) {}

// 1. destructor (to destroy the graph)
Universe::~Universe() {
  // in addition to basic fields being deleted, we need to clean up the graph
  // as is done in parent
  igraph_destroy(&this->graph);
};

// 2. copy constructor
Universe::Universe(const Universe &src) {
  this->timestep = src.timestep;
  this->NAtoms = src.NAtoms;
  this->NBonds = src.NBonds;
  this->angleFrom = src.angleFrom;
  this->angleTo = src.angleTo;
  this->angleVia = src.angleVia;
  // using copy assignement operators ourselfes
  this->box = src.box;
  this->atomIdToVectorIdx = src.atomIdToVectorIdx;
  this->weightPerType = src.weightPerType;
  igraph_copy(&this->graph, &src.graph);
};

// 3. copy assignment operator
Universe &Universe::operator=(Universe src) {
  std::swap(this->timestep, src.timestep);
  std::swap(this->NAtoms, src.NAtoms);
  std::swap(this->NBonds, src.NBonds);
  std::swap(this->angleFrom, src.angleFrom);
  std::swap(this->angleTo, src.angleTo);
  std::swap(this->angleVia, src.angleVia);
  std::swap(this->box, src.box);
  std::swap(this->graph, src.graph);
  std::swap(this->atomIdToVectorIdx, src.atomIdToVectorIdx);
  std::swap(this->weightPerType, src.weightPerType);

  return *this;
};

// other functions
void Universe::addAtoms(std::vector<long int> newIds, std::vector<int> newTypes,
                        std::vector<double> newX, std::vector<double> newY,
                        std::vector<double> newZ, std::vector<int> newNx,
                        std::vector<int> newNy, std::vector<int> newNz) {
  this->addAtoms(newIds.size(), newIds, newTypes, newX, newY, newZ, newNx,
                 newNy, newNz);
}

void Universe::addAtoms(const size_t NNewAtoms, std::vector<long int> newIds,
                        std::vector<int> newTypes, std::vector<double> newX,
                        std::vector<double> newY, std::vector<double> newZ,
                        std::vector<int> newNx, std::vector<int> newNy,
                        std::vector<int> newNz) {
  if (newTypes.size() != NNewAtoms || newIds.size() != newTypes.size() ||
      newX.size() != newNx.size() || newY.size() != newNy.size() ||
      newZ.size() != newNz.size() || newX.size() != newY.size() ||
      NNewAtoms != newZ.size()) {
    throw std::invalid_argument("All atom inputs must have the same size.");
  }
  // actually add the vertices
  if (igraph_add_vertices(&this->graph, NNewAtoms, 0)) {
    throw std::runtime_error("Failed to add new atoms to graph.");
  }
  this->atomIdToVectorIdx.reserve(this->NAtoms + NNewAtoms);
  // do map for easy access afterwards
  for (size_t i = 0; i < NNewAtoms; ++i) {
    if (this->atomIdToVectorIdx.contains(newIds[i])) {
      throw std::invalid_argument("Atom with id " + std::to_string(newIds[i]) +
                                  " already exists");
    }
    this->atomIdToVectorIdx.emplace(newIds[i], this->NAtoms + i);
  }
  // append attributes
  // it is empirically more efficient to do it this split up way,
  // though there might be even more efficient intermediate splits
  if (this->NAtoms == 0) {
    // NOTE: using the same vector over an over might be bad for performance?
    igraph_vector_t valueVec;
    igraph_vector_init(&valueVec, NNewAtoms);
    pylimer_tools::utils::StdVectorToIgraphVectorT(newIds, &valueVec);
    igraph_cattribute_VAN_setv(&this->graph, "id", &valueVec);
    pylimer_tools::utils::StdVectorToIgraphVectorT(newX, &valueVec);
    igraph_cattribute_VAN_setv(&this->graph, "x", &valueVec);
    pylimer_tools::utils::StdVectorToIgraphVectorT(newY, &valueVec);
    igraph_cattribute_VAN_setv(&this->graph, "y", &valueVec);
    pylimer_tools::utils::StdVectorToIgraphVectorT(newZ, &valueVec);
    igraph_cattribute_VAN_setv(&this->graph, "z", &valueVec);
    pylimer_tools::utils::StdVectorToIgraphVectorT(newTypes, &valueVec);
    igraph_cattribute_VAN_setv(&this->graph, "type", &valueVec);
    pylimer_tools::utils::StdVectorToIgraphVectorT(newNx, &valueVec);
    igraph_cattribute_VAN_setv(&this->graph, "nx", &valueVec);
    pylimer_tools::utils::StdVectorToIgraphVectorT(newNy, &valueVec);
    igraph_cattribute_VAN_setv(&this->graph, "ny", &valueVec);
    pylimer_tools::utils::StdVectorToIgraphVectorT(newNz, &valueVec);
    igraph_cattribute_VAN_setv(&this->graph, "nz", &valueVec);
    igraph_vector_destroy(&valueVec);
  } else {
    for (size_t i = 0; i < NNewAtoms; ++i) {
      igraph_cattribute_VAN_set(&this->graph, "id", this->NAtoms + i,
                                newIds[i]);
      igraph_cattribute_VAN_set(&this->graph, "x", this->NAtoms + i, newX[i]);
      igraph_cattribute_VAN_set(&this->graph, "y", this->NAtoms + i, newY[i]);
      igraph_cattribute_VAN_set(&this->graph, "z", this->NAtoms + i, newZ[i]);
      igraph_cattribute_VAN_set(&this->graph, "type", this->NAtoms + i,
                                newTypes[i]);
      igraph_cattribute_VAN_set(&this->graph, "nx", this->NAtoms + i, newNx[i]);
      igraph_cattribute_VAN_set(&this->graph, "ny", this->NAtoms + i, newNy[i]);
      igraph_cattribute_VAN_set(&this->graph, "nz", this->NAtoms + i, newNz[i]);
    }
  }
  // this->NAtoms += NNewAtoms;
  this->NAtoms = igraph_vcount(&this->graph);
}

void Universe::addBonds(std::vector<long int> from, std::vector<long int> to) {
  this->addBonds(from.size(), from, to);
}
void Universe::addBonds(const size_t NNewBonds, std::vector<long int> from,
                        std::vector<long int> to) {
  this->addBonds(NNewBonds, from, to, std::vector<int>());
}

void Universe::addAngles(std::vector<long int> from, std::vector<long int> via,
                         std::vector<long int> to) {
  if (from.size() != to.size() || from.size() != via.size()) {
    throw std::invalid_argument("All angle inputs must have the same size.");
  }

  this->angleFrom.insert(std::end(this->angleFrom), std::begin(from),
                         std::end(from));
  this->angleVia.insert(std::end(this->angleVia), std::begin(via),
                        std::end(via));
  this->angleTo.insert(std::end(this->angleTo), std::begin(to), std::end(to));
}

void Universe::addBonds(const size_t NNewBonds, std::vector<long int> from,
                        std::vector<long int> to, std::vector<int> bondTypes,
                        const bool ignoreNonExistentAtoms,
                        const bool simplify) {
  if (from.size() != to.size() || from.size() != NNewBonds) {
    throw std::invalid_argument("All bond inputs must have the same size.");
  }
  std::vector<long int> newEdgesVector =
      pylimer_tools::utils::interleave(from, to);
  size_t edgesSize = newEdgesVector.size();
  // translate from atomId to VertexIdx
  igraph_vector_t newEdges;
  size_t actualNrOfBondsAdded = 0;
  igraph_vector_init(&newEdges, edgesSize);
  int innerIndex = 0;
  for (size_t i = 1; i < edgesSize; i += 2) {
    if (this->atomIdToVectorIdx.contains(newEdgesVector[i - 1]) &&
        this->atomIdToVectorIdx.contains(newEdgesVector[i])) {
      igraph_vector_set(&newEdges, innerIndex,
                        this->atomIdToVectorIdx.at(newEdgesVector[i - 1]));
      innerIndex += 1;
      igraph_vector_set(&newEdges, innerIndex,
                        this->atomIdToVectorIdx.at(newEdgesVector[i]));
      innerIndex += 1;
      actualNrOfBondsAdded += 1;
    } else if (!ignoreNonExistentAtoms) {
      throw std::invalid_argument("Bond with atom with id " +
                                  std::to_string(newEdgesVector[i]) +
                                  " impossible as atom is not added yet.");
    }
  }
  assert(innerIndex == 2 * actualNrOfBondsAdded);
  igraph_vector_resize(&newEdges, 2 * actualNrOfBondsAdded);
  // add the new edges
  if (igraph_add_edges(&this->graph, &newEdges, 0)) {
    throw std::runtime_error("Failed to add edges to graph.");
  }
  igraph_vector_destroy(&newEdges);
  if (actualNrOfBondsAdded > 0) {
    // add attributes
    // if (bondTypes.size() == NNewBonds && this->NBonds ==
    // igraph_ecount(&this->graph) - NNewBonds)
    // {
    //   for (size_t i = 0; i < NNewBonds; ++i)
    //   {
    //     // append attributes
    //     igraph_cattribute_EAN_set(&this->graph, "type", this->NBonds + i,
    //     bondTypes[i]);
    //   }
    // }
    // else: too risky to add bond attributes
    // simplify graph
    // this->NBonds += NNewBonds;
    if (simplify) {
      this->simplify();
    } else {
      this->NBonds = igraph_ecount(&this->graph);
    }
  }
};

void Universe::simplify() {
  igraph_attribute_combination_t comb;
  igraph_attribute_combination_init(&comb);
  igraph_simplify(&this->graph, /*multiple=*/1, /*loops=*/1, &comb);
  igraph_attribute_combination_destroy(&comb);
  this->NBonds = igraph_ecount(&this->graph);
}

/**
 * @brief Count how many of each atom type there are in the universe
 *
 * @return std::map<int, int>
 */
std::map<int, int> Universe::countAtomTypes() const {
  std::vector<int> atomTypes = this->getAtomTypes();
  std::map<int, int> result;
  for (int atomType : atomTypes) {
    result[atomType] += 1;
  }
  return result;
};

/**
 * @brief Set the masses of the atoms in this universe
 *
 * @param weightPerType the weight per type
 */
void Universe::setMasses(std::map<int, double> weightPerType) {
  this->weightPerType = weightPerType;
}

std::map<int, double> Universe::getMasses() { return this->weightPerType; };

/**
 * @brief Get the standalone components of the network
 *
 * @return std::vector<Molecule>
 */
std::vector<Molecule> Universe::getClusters() const {
  std::vector<Molecule> molecules;
  if (this->getNrOfAtoms() == 0) {
    return molecules;
  }

  // split the copy into the separate components
  igraph_vector_ptr_t components;
  igraph_vector_ptr_init(&components, 0);
  if (igraph_decompose(&graph, &components, IGRAPH_WEAK, -1, 0)) {
    throw std::runtime_error("Failed to decompose graph.");
  }
  size_t NComponents = igraph_vector_ptr_size(&components);
  // std::cout << NComponents << " clusters found." << std::endl;
  molecules.reserve(NComponents);
  for (size_t i = 0; i < NComponents; ++i) {
    // make the molecule the owner of the graph
    igraph_t *g = (igraph_t *)VECTOR(components)[i];

    if (igraph_vcount(g)) {
      molecules.push_back(Molecule(&this->box, g, MoleculeType::UNDEFINED,
                                   this->weightPerType));
    }
  }
  igraph_decompose_destroy(&components);
  igraph_vector_ptr_destroy(&components);
  return molecules;
}

/**
 * @brief Decompose this universe into chains/molecules by splitting them into
 * clusters
 *
 * @param atomTypeToOmit The atom type to remove to get more clusters
 * @return std::vector<Molecule>
 */
std::vector<Molecule> Universe::getMolecules(const int atomTypeToOmit) const {
  std::vector<Molecule> molecules;
  if (this->getNrOfAtoms() == 0) {
    return molecules;
  }
  // make a copy to remove crosslinkers from
  igraph_t graphWithoutCrosslinkers;
  if (igraph_copy(&graphWithoutCrosslinkers, &this->graph)) {
    throw std::runtime_error("Failed to copy graph.");
  }

  // select vertices of type
  std::vector<long int> indicesToRemove =
      this->getIndicesOfType(atomTypeToOmit);
  std::sort(indicesToRemove.rbegin(), indicesToRemove.rend());
  if (indicesToRemove.size() > 0) {
    igraph_vs_t verticesToRemove = this->getVerticesByIndices(indicesToRemove);

    // remove elements of type
    if (igraph_delete_vertices(&graphWithoutCrosslinkers, verticesToRemove)) {
      throw std::runtime_error("Failed to delete crosslinkers from graph.");
    }

    igraph_vs_destroy(&verticesToRemove);
  }

  // split the copy into the separate components
  igraph_vector_ptr_t components;
  igraph_vector_ptr_init(&components, this->getNrOfAtoms());
  if (igraph_decompose(&graphWithoutCrosslinkers, &components, IGRAPH_WEAK, -1,
                       0)) {
    throw std::runtime_error("Failed to decompose graph.");
  }
  size_t NComponents = igraph_vector_ptr_size(&components);
  // std::cout << NComponents << " molecules found. Removed " <<
  // indicesToRemove.size()
  //           << " vertices. Size now: " <<
  //           igraph_vcount(&graphWithoutCrosslinkers) << " atoms with " <<
  //           igraph_ecount(&graphWithoutCrosslinkers) << " bonds." <<
  //           std::endl;
  molecules.reserve(NComponents);
  for (size_t i = 0; i < NComponents; ++i) {
    // make the molecule the owner of the graph
    igraph_t *g = (igraph_t *)VECTOR(components)[i];

    if (igraph_vcount(g)) {
      molecules.push_back(Molecule(&this->box, g, MoleculeType::UNDEFINED,
                                   this->weightPerType));
    }
  }
  igraph_decompose_destroy(&components);
  igraph_vector_ptr_destroy(&components);
  igraph_destroy(&graphWithoutCrosslinkers);
  return molecules;
}

/**
 * @brief Get a vertex selector to find all vertices with a certain type
 *
 * @param type the type to select
 * @return igraph_vs_t
 */
igraph_vs_t Universe::getVerticesOfType(const int type) const {
  std::vector<long int> indices = this->getIndicesOfType(type);
  return this->getVerticesByIndices(indices);
}

/**
 * @brief Get a vertex selector to find all vertices with a certain index
 *
 * @param indices the vertex indices to select
 * @return igraph_vs_t
 */
igraph_vs_t
Universe::getVerticesByIndices(std::vector<long int> indices) const {
  igraph_vector_t indicesToSelect;
  igraph_vector_init(&indicesToSelect, indices.size());
  pylimer_tools::utils::StdVectorToIgraphVectorT(indices, &indicesToSelect);
  igraph_vs_t result;
  if (igraph_vs_vector_copy(&result, &indicesToSelect)) {
    throw std::runtime_error("Failed to select vertices");
  }
  igraph_vector_destroy(&indicesToSelect);
  return result;
}

/**
 * @brief Get the vertex indices of atoms with a certain type
 *
 * @param type the type to select
 * @return std::vector<long int>
 */
std::vector<long int> Universe::getIndicesOfType(const int type) const {
  std::vector<long int> indices;
  if (this->getNrOfAtoms() == 0) {
    return indices;
  }

  igraph_vector_t types;
  igraph_vector_init(&types, this->getNrOfAtoms());
  VANV(&this->graph, "type", &types);
  for (int i = 0; i < this->NAtoms; ++i) {
    if (VECTOR(types)[i] == type) {
      indices.push_back(i);
    }
  }
  igraph_vector_destroy(&types);
  return indices;
}

/**
 * @brief Decompose the network into clusters, re-adding the atoms omitted to
 * get more clusters
 *
 * @param crosslinkerType the type of the atoms to omit and re-add
 * @return std::vector<Molecule>
 */
std::vector<Molecule>
Universe::getChainsWithCrosslinker(const int crosslinkerType) const {
  std::vector<Molecule> molecules;
  if (this->getNrOfAtoms() == 0) {
    return molecules;
  }
  // make a copy to remove crosslinkers from
  igraph_t graphWithoutCrosslinkers;
  if (igraph_copy(&graphWithoutCrosslinkers, &this->graph)) {
    throw std::runtime_error("Failed to copy graph.");
  }
  // select vertices of type
  std::vector<long int> indicesToRemove =
      this->getIndicesOfType(crosslinkerType);
  std::sort(indicesToRemove.rbegin(), indicesToRemove.rend());
  if (indicesToRemove.size() > 0) {
    igraph_vs_t verticesToRemove = this->getVerticesByIndices(indicesToRemove);

    // remove elements of type
    if (igraph_delete_vertices(&graphWithoutCrosslinkers, verticesToRemove)) {
      throw std::runtime_error("Failed to delete crosslinkers from graph.");
    }

    igraph_vs_destroy(&verticesToRemove);
  }

  // split the copy into the separate
  igraph_vector_ptr_t components;
  igraph_vector_ptr_init(&components, 3);
  if (igraph_decompose(&graphWithoutCrosslinkers, &components, IGRAPH_STRONG,
                       -1, 0)) {
    throw std::runtime_error("Failed to decompose graph.");
  }
  size_t NComponents = igraph_vector_ptr_size(&components);
  molecules.reserve(NComponents);
  for (size_t i = 0; i < NComponents; ++i) {
    // loop the chains to add the crosslinkers back
    igraph_t *chain = (igraph_t *)VECTOR(components)[i];
    int moleculeLengthBefore = igraph_vcount(chain);
    // also select ones of degree 0 for dangling atoms
    std::vector<long int> endNodeIndices =
        pylimer_tools::utils::getVerticesWithDegree(chain, {{0, 1}});
    igraph_vector_t endNodeSelectorVector;
    igraph_vector_init(&endNodeSelectorVector, endNodeIndices.size());
    pylimer_tools::utils::StdVectorToIgraphVectorT(endNodeIndices,
                                                   &endNodeSelectorVector);
    MoleculeType molType = MoleculeType::UNDEFINED;
    bool isLoop = false;

    if (moleculeLengthBefore >= 1) {
      // no care about single atoms for now (though they are included)
      igraph_vit_t endNodeVit;
      igraph_vit_create(chain, igraph_vss_vector(&endNodeSelectorVector),
                        &endNodeVit);
      // collect atoms to add, since adding them would invalidate the iterator
      std::vector<long int> atomsToAdd;
      std::vector<std::vector<long int>> bondsToAdd;
      // loop end nodes
      while (!IGRAPH_VIT_END(endNodeVit)) {
        long int newEndNodeVertexId = (long int)IGRAPH_VIT_GET(endNodeVit);
        long int oldEndNodeId = (long int)VAN(chain, "id", newEndNodeVertexId);
        long int originalEndNodeVertexId =
            this->findVertexIdForProperty("id", oldEndNodeId);
        igraph_vector_t neighbours;
        igraph_vector_init(&neighbours, 0);

        if (igraph_neighbors(&graph, &neighbours, originalEndNodeVertexId,
                             IGRAPH_ALL)) {
          throw std::runtime_error("Failed to get neighbours in graph");
        }

        std::vector<long int> neighborsVec;
        pylimer_tools::utils::igraphVectorTToStdVector(&neighbours,
                                                       neighborsVec);

        // loop neighbours
        for (long int originalNeighbourId : neighborsVec) {
          int originalNeighbourType =
              igraph_cattribute_VAN(&graph, "type", originalNeighbourId);

          if (originalNeighbourType == crosslinkerType) {
            // found a crosslinker neighbour
            long int originalNeighbourAtomId =
                igraph_cattribute_VAN(&graph, "id", originalNeighbourId);
            atomsToAdd.push_back(originalNeighbourId);
            bondsToAdd.push_back({{newEndNodeVertexId, originalNeighbourId}});
          }
        }

        if (atomsToAdd.size() == 2 && atomsToAdd[0] == atomsToAdd[1]) {
          isLoop = true;
          // we only want to add it once -> remove
          atomsToAdd.pop_back();
        }

        IGRAPH_VIT_NEXT(endNodeVit);
        igraph_vector_destroy(&neighbours);
      } // loop end nodes

      std::unordered_map<long int, long int> newAtomsMap;
      // actually add the atoms...
      for (auto atomToAddOriginalId : atomsToAdd) {
        igraph_add_vertices(chain, 1, 0);
        long int newCrosslinkerVertexIdx = igraph_vcount(chain) - 1;
        newAtomsMap.insert_or_assign(atomToAddOriginalId,
                                     newCrosslinkerVertexIdx);
        // additional loop check
        long int originalNeighbourAtomId =
            igraph_cattribute_VAN(&graph, "id", atomToAddOriginalId);
        if (pylimer_tools::utils::graphHasVertexWithProperty(
                chain, "id", originalNeighbourAtomId)) {
          isLoop = true;
        }

        // including all attributes
        for (auto property : {"id", "type", "x", "y", "z", "nx", "ny", "nz"}) {
          SETVAN(chain, property, newCrosslinkerVertexIdx,
                 VAN(&graph, property, atomToAddOriginalId));
        }
      }
      // ...and bonds
      for (auto bond : bondsToAdd) {
        igraph_add_edge(chain, bond[0], newAtomsMap[bond[1]]);
      }
      igraph_vit_destroy(&endNodeVit);
    } // if molecule length
    igraph_vector_destroy(&endNodeSelectorVector);
    // decide on molecule type
    int newMoleculeLength = igraph_vcount(chain);
    if (newMoleculeLength == moleculeLengthBefore) {
      molType = MoleculeType::FREE_CHAIN;
    } else if (newMoleculeLength == moleculeLengthBefore + 1) {
      molType = MoleculeType::DANGLING_CHAIN;
    } else if (newMoleculeLength == moleculeLengthBefore + 2) {
      molType = MoleculeType::NETWORK_STRAND;
    }
    if (isLoop) {
      molType = MoleculeType::PRIMARY_LOOP;
    }

    // finally, create the molecule/chain
    molecules.push_back(
        Molecule(&this->box, chain, molType, this->weightPerType));
  }
  igraph_decompose_destroy(&components);
  igraph_vector_ptr_destroy(&components);
  igraph_destroy(&graphWithoutCrosslinkers);

  return molecules;
}

/**
 * @brief Find the loops in the network
 *
 * NOTE: there are exponentially many paths between two vertices of a graph,
 * and you may run out of memory when using this function, if your graph is
 * lattice-like.
 *
 * @param crosslinkerType
 * @param maxLength
 * @return std::map<int, std::vector<std::vector<Atom>>>
 */
std::map<int, std::vector<std::vector<Atom>>>
Universe::findLoops(const int crosslinkerType, const int maxLength,
                    bool skipSelfLoops) const {
  // NOTE: there are exponentially many paths between two vertices of a graph,
  // and you may run out of memory when using this function, if your graph is
  // lattice-like.
  std::map<int, std::vector<std::vector<Atom>>> results;

  std::vector<long int> startingCrosslinkers =
      this->getIndicesOfType(crosslinkerType);
  std::unordered_set<int> processedPathsKeys;

  // note: this algorithm is not particularly efficient
  // it is of the order of O(n*n!)
  for (long int startingCrosslinkerVertexId : startingCrosslinkers) {
    // ideally, we would only select the neighbouring *cross-linkers* here to
    // reduce the overhead. but well: let's leave that to the user with
    // #getNetworkOfCrosslinker
    igraph_vector_t neighbours;
    igraph_vector_init(&neighbours, 0);

    if (igraph_neighbors(&graph, &neighbours, startingCrosslinkerVertexId,
                         IGRAPH_ALL)) {
      throw std::runtime_error("Failed to get neighbours in graph");
    }

    // loop neighbours
    igraph_vector_int_t paths;
    igraph_vector_int_init(&paths, 0);
    // for each neighbour, we search the simple paths
    if (igraph_get_all_simple_paths(
            &this->graph, &paths, startingCrosslinkerVertexId,
            igraph_vss_vector(&neighbours), maxLength, IGRAPH_ALL)) {
      throw std::runtime_error("Failed to get simple paths in graph");
    }

    igraph_vector_destroy(&neighbours);
    // translate the paths we found
    std::vector<Atom> currentPath;
    int currentFunctionality = 0;
    long int currentPathKey = 0;
    size_t n = igraph_vector_int_size(&paths);
    for (size_t i = 0; i < n; ++i) {
      const long int currentVal = igraph_vector_int_e(&paths, i);
      if (currentVal == -1) {
        // skip self-loops and duplicates
        if ((!skipSelfLoops || currentPath.size() > 3) &&
            !processedPathsKeys.contains(currentPathKey)) {
          results[currentFunctionality].push_back(currentPath);
          processedPathsKeys.insert(currentPathKey);
        }
        currentPath.clear();
        currentFunctionality = 0;
        currentPathKey = 0;
      } else {
        Atom newAtom = this->getAtomByVertexIdx(currentVal);
        currentPathKey = currentPathKey xor currentVal; // compute hash
        currentPath.push_back(newAtom);
        if (newAtom.getType() == crosslinkerType) {
          currentFunctionality += 1;
        }
      }
    }
    igraph_vector_int_destroy(&paths);

    // additional check for self-loops
    if (!skipSelfLoops &&
        !processedPathsKeys.contains(0 xor startingCrosslinkerVertexId)) {
      std::vector<long int> crosslinkersBonds =
          this->getVertexIdxsConnectedTo(startingCrosslinkerVertexId);
      if (std::find(crosslinkersBonds.begin(), crosslinkersBonds.end(),
                    startingCrosslinkerVertexId) != crosslinkersBonds.end()) {

        currentPath.clear();
        currentPath.push_back(
            this->getAtomByVertexIdx(startingCrosslinkerVertexId));
        results[1].push_back(currentPath);
        currentPath.clear();
      }
    }
  }

  return results;
};

/**
 * @brief Check whether the universe contains a loop that crosses the periodic
 * boundaries an odd times
 *
 * NOTE: there are exponentially many paths between two vertices of a graph,
 * and you may run out of memory when using this function, if your graph is
 * lattice-like.
 *
 * @param crosslinkerType
 * @param maxLength the maximum length of the loop ()
 * @return true
 * @return false
 */
bool Universe::hasInfiniteStrand(const int crosslinkerType,
                                 const int maxLength) const {
  // NOTE: there are exponentially many paths between two vertices of a graph,
  // and you may run out of memory when using this function, if your graph is
  // lattice-like.
  std::map<int, std::vector<std::vector<Atom>>> results;

  std::vector<long int> startingCrosslinkers =
      this->getIndicesOfType(crosslinkerType);
  std::unordered_set<int> processedPathsKeys;

  // note: this algorithm is not particularly efficient
  // it is of the order of O(n*n!)
  for (long int startingCrosslinkerVertexId : startingCrosslinkers) {
    // select all neighbouring atoms as possible directions for the loop
    igraph_vector_t neighbours;
    igraph_vector_init(&neighbours, 0);

    if (igraph_neighbors(&graph, &neighbours, startingCrosslinkerVertexId,
                         IGRAPH_ALL)) {
      throw std::runtime_error("Failed to get neighbours in graph");
    }

    // loop neighbours
    igraph_vector_int_t paths;
    igraph_vector_int_init(&paths, 0);
    // for each neighbour, we search the simple paths
    if (igraph_get_all_simple_paths(
            &this->graph, &paths, startingCrosslinkerVertexId,
            igraph_vss_vector(&neighbours), maxLength, IGRAPH_ALL)) {
      throw std::runtime_error("Failed to get simple paths in graph");
    }

    igraph_vector_destroy(&neighbours);
    // translate the paths we found
    std::vector<Atom> currentPath;
    int nrOfTraversalsX = 0;
    int nrOfTraversalsY = 0;
    int nrOfTraversalsZ = 0;
    size_t n = igraph_vector_int_size(&paths);
    for (int i = 0; i < n; ++i) {
      const long int currentVal = igraph_vector_int_e(&paths, i);
      if (currentVal == -1) {
        // finished a loop. Check.
        // we have an infinite loop if the box boundary was passed in one
        // direction only NOTE: this neglects infinite networks (≠ infinite
        // loops) such as ones caused by entanglement between images
        if (nrOfTraversalsX != 0 || nrOfTraversalsY != 0 ||
            nrOfTraversalsZ != 0) {
          igraph_vector_int_destroy(&paths);
          return true;
        }
        // then reset
        currentPath.clear();
        nrOfTraversalsX = 0;
        nrOfTraversalsY = 0;
        nrOfTraversalsZ = 0;
      } else {
        Atom newAtom = this->getAtomByVertexIdx(currentVal);
        if (!currentPath.empty()) {
          Atom lastAtom = currentPath.back();
          double dx = newAtom.getX() - lastAtom.getX();
          nrOfTraversalsX += (dx) > 0.5 * (this->box.getLx())
                                 ? 1
                                 : (dx < -0.5 * (this->box.getLx()) ? -1 : 0);
          double dy = newAtom.getY() - lastAtom.getY();
          nrOfTraversalsY += (dx) > 0.5 * (this->box.getLy())
                                 ? 1
                                 : (dy < -0.5 * (this->box.getLy()) ? -1 : 0);
          double dz = newAtom.getZ() - lastAtom.getZ();
          nrOfTraversalsZ += (dz) > 0.5 * (this->box.getLz())
                                 ? 1
                                 : (dz < -0.5 * (this->box.getLz()) ? -1 : 0);
        }
        currentPath.push_back(newAtom);
      }
    }
    igraph_vector_int_destroy(&paths);
  }

  return false;
}

/**
 * @brief Determine the maximum functionality per atom type
 *
 * @return std::map<int, int>
 */
std::map<int, int> Universe::determineFunctionalityPerType() const {
  std::map<int, int> result;
  igraph_vector_t degrees;
  if (igraph_vector_init(&degrees, 0)) {
    throw std::runtime_error("Failed to instantiate result vector.");
  }
  igraph_vs_t allVertexIds;
  igraph_vs_all(&allVertexIds);
  // complexity: O(|v|*d)
  if (igraph_degree(&this->graph, &degrees, allVertexIds, IGRAPH_ALL, false)) {
    throw std::runtime_error("Failed to determine degree of vertices");
  }

  std::vector<int> types = this->getPropertyValues<int>("type");
  std::set<int> uniqueTypes(types.begin(), types.end());
  // make sure the keys are (re)set, for every type
  for (int type : uniqueTypes) {
    result[type] = 0;
  }

  // complexity: O(|V|)
  igraph_vit_t vit;
  igraph_vit_create(&graph, allVertexIds, &vit);
  while (!IGRAPH_VIT_END(vit)) {
    long int vertexId = static_cast<long int>(IGRAPH_VIT_GET(vit));
    result[types[vertexId]] = std::max((int)igraph_vector_e(&degrees, vertexId),
                                       result[types[vertexId]]);
    IGRAPH_VIT_NEXT(vit);
  }
  igraph_vit_destroy(&vit);
  igraph_vs_destroy(&allVertexIds);
  igraph_vector_destroy(&degrees);

  return result;
}

/**
 * @brief determine the effective functionality of each atom type
 *
 * @return std::map<int, double>
 */
std::map<int, double> Universe::determineEffectiveFunctionalityPerType() const {
  std::map<int, double> result;
  igraph_vector_t degrees;
  if (igraph_vector_init(&degrees, 0)) {
    throw std::runtime_error("Failed to instantiate result vector.");
  }
  igraph_vs_t allVertexIds;
  igraph_vs_all(&allVertexIds);
  // complexity: O(|v|*d)
  if (igraph_degree(&this->graph, &degrees, allVertexIds, IGRAPH_ALL, false)) {
    throw std::runtime_error("Failed to determine degree of vertices");
  }

  std::vector<int> types = this->getPropertyValues<int>("type");
  std::set<int> uniqueTypes(types.begin(), types.end());
  // make sure the keys are (re)set, for every type
  for (int type : uniqueTypes) {
    result[type] = 0;
  }

  // complexity: O(|V|)
  igraph_vit_t vit;
  igraph_vit_create(&graph, allVertexIds, &vit);
  while (!IGRAPH_VIT_END(vit)) {
    long int vertexId = static_cast<long int>(IGRAPH_VIT_GET(vit));
    result[types[vertexId]] += igraph_vector_e(&degrees, vertexId);
    IGRAPH_VIT_NEXT(vit);
  }
  igraph_vit_destroy(&vit);
  igraph_vs_destroy(&allVertexIds);
  igraph_vector_destroy(&degrees);

  std::map<int, int> typeCounts = this->countAtomTypes();
  for (const auto typePair : typeCounts) {
    result[typePair.first] /= typePair.second;
  }

  return result;
}

/**
 * @brief Compute the weight fractions of each atom type in the network.
 *
 * @return std::map<int, double> $\\vec{W_i}$ (dict): using the type i as a key,
 * this dict contains the weight fractions ($\\frac{W_i}{W_{tot}}$)
 */
std::map<int, double> Universe::computeWeightFractions() const {
  std::map<int, double> partialMasses;
  if (this->getNrOfAtoms() == 0) {
    return partialMasses;
  }

  std::vector<int> types = this->getPropertyValues<int>("type");
  double totalMass = 0.0;
  for (int type : types) {
    totalMass += this->weightPerType.at(type);
    partialMasses.try_emplace(type, 0.0);
    partialMasses[type] += this->weightPerType.at(type);
  }

  if (totalMass == 0.0) {
    return partialMasses;
  }

  //  loop to turn partial masses into weight fractions
  for (const auto &partialMassPair : partialMasses) {
    partialMasses[partialMassPair.first] = partialMassPair.second / totalMass;
  }

  return partialMasses;
}

/**
 * @brief Get an atom id by its vertex index
 *
 * @param vertexId
 * @return long int
 */
long int Universe::getAtomIdByIdx(const int vertexId) const {
  return VAN(&this->graph, "id", vertexId);
}

/**
 * @brief Get a vertex index by the atom id
 *
 * @param atomId
 * @return long int
 */
long int Universe::getIdxByAtomId(const int atomId) const {
  if (!this->atomIdToVectorIdx.contains(atomId)) {
    throw std::invalid_argument(
        "Universe cannot return idx of atom id: atom with this id (" +
        std::to_string(atomId) + ") does not exist");
  }
  return this->atomIdToVectorIdx.at(atomId);
}

/**
 * @brief Get an atom by its id
 *
 * @param atomId
 * @return Atom
 */
Atom Universe::getAtom(const int atomId) const {
  return this->getAtomByVertexIdx(this->getIdxByAtomId(atomId));
}

/**
 * @brief Get all atoms in this universe
 *
 * @return std::vector<Atom>
 */
std::vector<Atom> Universe::getAtoms() {
  std::vector<Atom> atoms;
  atoms.reserve(this->getNrOfAtoms());
  igraph_vit_t vit;
  igraph_vit_create(&graph, igraph_vss_all(), &vit);
  while (!IGRAPH_VIT_END(vit)) {
    long int vertexId = static_cast<long int>(IGRAPH_VIT_GET(vit));
    atoms.push_back(this->getAtomByVertexIdx(vertexId));
    IGRAPH_VIT_NEXT(vit);
  }
  igraph_vit_destroy(&vit);
  return atoms;
}

/**
 * @brief Get all angles stored in this universe
 *
 * @return std::map<std::string, std::vector<long int>>
 */
std::map<std::string, std::vector<long int>> Universe::getAngles() const {
  std::map<std::string, std::vector<long int>> results;
  results.insert_or_assign("angle_from", this->angleFrom);
  results.insert_or_assign("angle_to", this->angleTo);
  results.insert_or_assign("angle_via", this->angleVia);

  return results;
}

/**
 * @brief Find all angles that appear in this universe
 *
 * @return std::map<std::string, std::vector<long int>>
 */
std::map<std::string, std::vector<long int>> Universe::detectAngles() const {
  std::set<int> anglesFound;
  std::vector<long int> angleFromFound;
  std::vector<long int> angleToFound;
  std::vector<long int> angleViaFound;
  // query all atoms
  igraph_vit_t vit;
  igraph_vit_create(&this->graph, igraph_vss_all(), &vit);
  while (!IGRAPH_VIT_END(vit)) {
    const long int vertexIdx = static_cast<long int>(IGRAPH_VIT_GET(vit));
    // find the connected atoms
    std::vector<long int> connections =
        this->getVertexIdxsConnectedTo(vertexIdx);
    // with 1 or 0 connections, there is no angle
    if (connections.size() < 2) {
      IGRAPH_VIT_NEXT(vit);
      continue;
    }
    // loop the connections to find angles
    for (size_t connectionI = 0; connectionI < connections.size();
         ++connectionI) {
      for (size_t connectionJ = connectionI + 1;
           connectionJ < connections.size(); ++connectionJ) {
        int angleKey =
            vertexIdx xor connections[connectionI] xor connections[connectionJ];
        if (!anglesFound.contains(angleKey)) {
          angleFromFound.push_back(
              this->getAtomIdByIdx(connections[connectionI]));
          angleViaFound.push_back(this->getAtomIdByIdx(vertexIdx));
          angleToFound.push_back(
              this->getAtomIdByIdx(connections[connectionJ]));
          anglesFound.insert(angleKey);
        }
      }
    }
    IGRAPH_VIT_NEXT(vit);
  }
  igraph_vit_destroy(&vit);

  std::map<std::string, std::vector<long int>> results;
  results.insert_or_assign("angle_from", angleFromFound);
  results.insert_or_assign("angle_to", angleToFound);
  results.insert_or_assign("angle_via", angleViaFound);

  return results;
}

/**
 * @brief Reduce the network to cross-linkers only.
 * Includes self-loops (from primary ones).
 * You can remove them using #simplify
 *
 * @param crosslinkerType the atom type of the cross-linker beads
 * @return Universe
 */
Universe Universe::getNetworkOfCrosslinker(const int crosslinkerType) const {
  // How this works:
  // 1. find all crosslinkers
  // 2. from each crosslinker, walk in all directions
  // 3. if the walk reaches another crosslinker, we found a
  //    crosslinker-crosslinker connection.
  //    To reduce duplicates, we only take the ones where we started from a
  //    crosslinker with a smaller (or equal, for self-/primary-loops) vertex
  //    index
  Universe newUniverse =
      Universe(this->box.getLx(), this->box.getLy(), this->box.getLz());
  std::vector<long int> bondFrom;
  std::vector<long int> bondTo;
  std::vector<long int> crosslinkers = this->getIndicesOfType(crosslinkerType);
  for (long int crosslinker : crosslinkers) {
    std::vector<long int> connections =
        this->getVertexIdxsConnectedTo(crosslinker);
    for (long int connection : connections) {
      long int currentCenter = connection;
      long int lastCenter = crosslinker;
      std::vector<long int> subConnections =
          this->getVertexIdxsConnectedTo(currentCenter);
      while (subConnections.size() > 0) {
        if (this->getPropertyValue<int>("type", currentCenter) ==
            crosslinkerType) {
          // found cross-linker
          if (currentCenter >= crosslinker) {
            bondFrom.push_back(
                this->getPropertyValue<int>("id", currentCenter));
            bondTo.push_back(this->getPropertyValue<int>("id", crosslinker));
          }
          break;
        }

        if (subConnections.size() == 1) {
          break;
        }
        // we assume a functionality of 2 for ordinary strands
        assert(subConnections.size() == 2);
        int subConnectionDirection = (subConnections[0] == lastCenter) ? 1 : 0;
        lastCenter = currentCenter;
        currentCenter = subConnections[subConnectionDirection];
        subConnections = this->getVertexIdxsConnectedTo(currentCenter);
      }
    }
  }

  std::vector<int> zeros =
      pylimer_tools::utils::initializeWithValue(crosslinkers.size(), 0);

  newUniverse.addAtoms(crosslinkers.size(),
                       this->getPropertyValues<long int>("id", crosslinkers),
                       this->getPropertyValues<int>("type", crosslinkers),
                       this->getPropertyValues<double>("x", crosslinkers),
                       this->getPropertyValues<double>("y", crosslinkers),
                       this->getPropertyValues<double>("z", crosslinkers),
                       zeros, zeros, zeros);
  newUniverse.addBonds(
      bondFrom.size(), bondFrom, bondTo,
      pylimer_tools::utils::initializeWithValue(bondFrom.size(), 0), false,
      false);

  return newUniverse;
};

/**
 * @brief Get the number of angles stored in this universe.
 *
 * @return const int
 */
size_t Universe::getNrOfAngles() const {
  assert(this->angleFrom.size() == this->angleTo.size());
  assert(this->angleFrom.size() == this->angleVia.size());
  return this->angleFrom.size();
}

/**
 * @brief Find the vertex id where the vertex has a certain property
 *
 * @tparam IN the type of the property you want the vertex to have
 * @param propertyName the name of the property
 * @param propertyValue the objective value
 * @return long int the vertex index, -1 if not found.
 */
template <typename IN>
long int Universe::findVertexIdForProperty(const char *propertyName,
                                           IN propertyValue) const {
  igraph_vector_t allValues;
  igraph_vector_init(&allValues, this->getNrOfAtoms());
  VANV(&this->graph, propertyName, &allValues);
  for (int i = 0; i < this->NAtoms; ++i) {
    if (VECTOR(allValues)[i] == propertyValue) {
      igraph_vector_destroy(&allValues);
      return i;
    }
  }
  igraph_vector_destroy(&allValues);
  return -1;
}

/**
 * @brief Compute the x distance for all bonds passed in
 *
 * @param bondFrom
 * @param bondTo
 * @return std::vector<double>
 */
std::vector<double> Universe::computeDxs(const std::vector<int> bondFrom,
                                         const std::vector<int> bondTo) {
  return this->computeDs(bondFrom, bondTo, "x", this->box.getLx());
};

/**
 * @brief Compute the y distance for all bonds passed in
 *
 * @param bondFrom
 * @param bondTo
 * @return std::vector<double>
 */
std::vector<double> Universe::computeDys(const std::vector<int> bondFrom,
                                         const std::vector<int> bondTo) {
  return this->computeDs(bondFrom, bondTo, "y", this->box.getLy());
};

/**
 * @brief Compute the z distance for all bonds passed in
 *
 * @param bondFrom
 * @param bondTo
 * @return std::vector<double>
 */
std::vector<double> Universe::computeDzs(const std::vector<int> bondFrom,
                                         const std::vector<int> bondTo) {
  return this->computeDs(bondFrom, bondTo, "z", this->box.getLz());
};

/**
 * @brief Compute the distance for all bonds passed in in a certain direction
 *
 * @param bondFrom
 * @param bondTo
 * @param direction
 * @param boxLimit
 * @return std::vector<double>
 */
std::vector<double> Universe::computeDs(const std::vector<int> bondFrom,
                                        const std::vector<int> bondTo,
                                        std::string direction,
                                        double boxLimit) const {
  if (bondFrom.size() != bondTo.size()) {
    throw std::invalid_argument(
        "bond from and bond to must have the same size");
  }

  size_t nBonds = bondFrom.size();

  igraph_vector_t vertexIdFrom;
  igraph_vector_init(&vertexIdFrom, nBonds);
  igraph_vector_t vertexIdTo;
  igraph_vector_init(&vertexIdTo, nBonds);

  for (size_t i = 0; i < nBonds; ++i) {
    igraph_vector_set(&vertexIdFrom, i,
                      this->atomIdToVectorIdx.at(bondFrom[i]));
    igraph_vector_set(&vertexIdTo, i, this->atomIdToVectorIdx.at(bondTo[i]));
  }

  igraph_vector_t dValuesFrom;
  igraph_vector_init(&dValuesFrom, nBonds);
  igraph_vector_t dValuesTo;
  igraph_vector_init(&dValuesTo, nBonds);

  std::string property = direction;
  igraph_cattribute_VANV(&this->graph, property.c_str(),
                         igraph_vss_vector(&vertexIdFrom), &dValuesFrom);
  igraph_cattribute_VANV(&this->graph, property.c_str(),
                         igraph_vss_vector(&vertexIdTo), &dValuesTo);

  igraph_vector_destroy(&vertexIdFrom);
  igraph_vector_destroy(&vertexIdTo);

  std::vector<double> results;
  results.reserve(nBonds);

  for (int i = 0; i < nBonds; ++i) {
    double currentD = (double)igraph_vector_e(&dValuesTo, i) -
                      (double)igraph_vector_e(&dValuesFrom, i);
    while (std::fabs(currentD) > 0.5 * boxLimit) {
      if (currentD < 0.0) {
        currentD += boxLimit;
      } else {
        currentD -= boxLimit;
      }
    }
    results.push_back(currentD);
  }

  igraph_vector_destroy(&dValuesFrom);
  igraph_vector_destroy(&dValuesTo);

  return results;
};

/**
 * @brief Get the mean number of beads between beads with the passed type
 *
 * @param junctionType
 * @return double
 */
double Universe::getMeanStrandLength(int junctionType) {
  std::vector<Molecule> molecules = this->getMolecules(junctionType);

  double multiplier = 1.0 / molecules.size();
  double meanStrandLength = 0;

  for (Molecule molecule : molecules) {
    meanStrandLength += molecule.getLength() * multiplier;
  }
  return meanStrandLength;
}

/**
 * @brief Get the mean end to end distance
 *
 * Does not take loops etc. into account.
 *
 * @param junctionType
 * @return double
 */
double Universe::computeMeanEndToEndDistance(int junctionType) {
  std::vector<Molecule> molecules =
      this->getChainsWithCrosslinker(junctionType);

  double meanEndToEndDistance = 0;
  int validMolecules = 0;

  for (Molecule molecule : molecules) {
    double dist = molecule.computeEndToEndDistance();
    if (dist > 0.0) {
      meanEndToEndDistance += dist;
      validMolecules += 1;
    }
  }

  return meanEndToEndDistance / (double)validMolecules;
}

/**
 * @brief Compute the mean of all bond lengths
 *
 * @return double
 */
double Universe::computeMeanBondLength() {
  double length = 0.0;
  if (this->getNrOfBonds() == 0) {
    return length;
  }
  // construct iterator
  igraph_eit_t bondIterator;
  if (igraph_eit_create(&this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID),
                        &bondIterator)) {
    throw std::runtime_error("Cannot create iterator to loop bonds");
  }

  while (!IGRAPH_EIT_END(bondIterator)) {
    long int edgeId = static_cast<long int>(IGRAPH_EIT_GET(bondIterator));
    int bondFrom;
    int bondTo;
    igraph_edge(&this->graph, edgeId, &bondFrom, &bondTo);
    // TODO: this is more intensive than needed
    // check whether the compiler optimizes this or not
    Atom atom1 = this->getAtomByVertexIdx(bondFrom);
    Atom atom2 = this->getAtomByVertexIdx(bondTo);
    length += atom1.distanceTo(atom2, &this->box);
    IGRAPH_EIT_NEXT(bondIterator);
  }

  igraph_eit_destroy(&bondIterator);
  return length / (double)this->getNrOfBonds();
}

/**
 * @brief check whether the internal counts are equal to those of the graph
 *
 * @throws std::runtime_error if the validation fails
 * @return true
 */
bool Universe::validate() {
  if (this->getNrOfAtoms() != igraph_vcount(&this->graph)) {
    throw std::runtime_error(
        "Validation failed: " + std::to_string(this->getNrOfAtoms()) +
        " atoms for " + std::to_string(igraph_vcount(&this->graph)) +
        " vertices.");
  }
  if (this->getNrOfBonds() != igraph_ecount(&this->graph)) {
    throw std::runtime_error(
        "Validation failed: " + std::to_string(this->getNrOfBonds()) +
        " bonds for " + std::to_string(igraph_ecount(&this->graph)) +
        " edges.");
  }
  return true;
}

/**
 * @brief Compute the volume of the underlying box
 *
 * @return double
 */
double Universe::getVolume() { return this->box.getVolume(); }

size_t Universe::getNrOfAtoms() const { return this->NAtoms; }

size_t Universe::getNrOfBonds() const { return this->NBonds; }

void Universe::setBox(Box passedBox) { this->box = passedBox; }

void Universe::setBoxLengths(const double Lx, const double Ly,
                             const double Lz) {
  this->setBox(Box(Lx, Ly, Lz));
}

Box Universe::getBox() { return this->box; }
} // namespace entities
} // namespace pylimer_tools
