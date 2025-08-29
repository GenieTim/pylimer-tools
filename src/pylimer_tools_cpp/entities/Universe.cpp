#include "Universe.h"
#include "../topo/TopologyCalc.h"
#include "../utils/BoolUtils.h"
#include "../utils/Counter.h"
#include "../utils/GraphUtils.h"
#include "../utils/StringUtils.h"
#include "../utils/VectorUtils.h"
#include "Box.h"
#include "NeighbourList.h"

#include <Eigen/Dense>

extern "C"
{
#include <igraph/igraph.h>
}
#include <algorithm>
#include <cassert>
#include <iterator> // for back_inserter
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace pylimer_tools::entities {

igraph_error_t
count_found_cycle(const igraph_vector_int_t* vertices,
                  const igraph_vector_int_t* edges,
                  void* arg)
{
  pylimer_tools::utils::Counter<int>* loopLengths =
    static_cast<pylimer_tools::utils::Counter<int>*>(arg);
  loopLengths->increment(igraph_vector_int_size(vertices));
  return IGRAPH_SUCCESS;
};

Universe::Universe(const Box& box)
{
  /* turn on attribute handling: TODO: move to some main() function  */
  if (!igraph_has_attribute_table()) {
    igraph_set_attribute_table(&igraph_cattribute_table);
  }
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
  : Universe(Box(Lx, Ly, Lz))
{
}

// 1. destructor (to destroy the graph)
Universe::~Universe()
{
  // in addition to basic fields being deleted, we need to clean up the graph
  // as is done in parent
  igraph_destroy(&this->graph);
};

// 2. copy constructor
Universe::Universe(const Universe& src)
{
  this->timestep = src.timestep;
  this->NAtoms = src.NAtoms;
  this->NBonds = src.NBonds;
  // angles
  this->angleFrom = src.angleFrom;
  this->angleVia = src.angleVia;
  this->angleTo = src.angleTo;
  this->angleType = src.angleType;

  // dihedral angles
  this->dihedralAngleFrom = src.dihedralAngleFrom;
  this->dihedralAngleVia1 = src.dihedralAngleVia1;
  this->dihedralAngleVia2 = src.dihedralAngleVia2;
  this->dihedralAngleTo = src.dihedralAngleTo;
  this->dihedralAngleType = src.dihedralAngleType;

  // using copy assignement operators ourselfes
  this->box = src.box;
  this->atomIdToVertexIdx = src.atomIdToVertexIdx;
  this->atomsHaveCustomAttributes = src.atomsHaveCustomAttributes;
  this->massPerType = src.massPerType;
  igraph_copy(&this->graph, &src.graph);
};

// 3. copy assignment operator
Universe&
Universe::operator=(Universe src)
{
  std::swap(this->timestep, src.timestep);
  std::swap(this->NAtoms, src.NAtoms);
  std::swap(this->NBonds, src.NBonds);

  // angles
  std::swap(this->angleFrom, src.angleFrom);
  std::swap(this->angleVia, src.angleVia);
  std::swap(this->angleTo, src.angleTo);
  std::swap(this->angleType, src.angleType);

  // dihedrals
  std::swap(this->dihedralAngleFrom, src.dihedralAngleFrom);
  std::swap(this->dihedralAngleVia1, src.dihedralAngleVia1);
  std::swap(this->dihedralAngleVia2, src.dihedralAngleVia2);
  std::swap(this->dihedralAngleTo, src.dihedralAngleTo);
  std::swap(this->dihedralAngleType, src.dihedralAngleType);

  //
  std::swap(this->box, src.box);
  std::swap(this->graph, src.graph);
  std::swap(this->atomsHaveCustomAttributes, src.atomsHaveCustomAttributes);
  std::swap(this->atomIdToVertexIdx, src.atomIdToVertexIdx);
  std::swap(this->massPerType, src.massPerType);

  return *this;
};

// equality operator
bool
Universe::operator==(const Universe& other) const
{
  // Compare basic properties
  if (this->timestep != other.timestep || this->NAtoms != other.NAtoms ||
      this->NBonds != other.NBonds || this->box != other.box ||
      this->atomIdToVertexIdx != other.atomIdToVertexIdx ||
      this->atomsHaveCustomAttributes != other.atomsHaveCustomAttributes ||
      this->massPerType != other.massPerType) {
    return false;
  }

  // Compare angles
  if (this->angleFrom != other.angleFrom || this->angleVia != other.angleVia ||
      this->angleTo != other.angleTo || this->angleType != other.angleType) {
    return false;
  }

  // Compare dihedral angles
  if (this->dihedralAngleFrom != other.dihedralAngleFrom ||
      this->dihedralAngleVia1 != other.dihedralAngleVia1 ||
      this->dihedralAngleVia2 != other.dihedralAngleVia2 ||
      this->dihedralAngleTo != other.dihedralAngleTo ||
      this->dihedralAngleType != other.dihedralAngleType) {
    return false;
  }

  // Compare graph structure
  // For a complete equality check, we need to compare all vertex and edge
  // attributes

  // start with comparing the connectivity
  igraph_t difference;
  igraph_difference(&difference, &this->graph, &other.graph);
  const igraph_integer_t diff_count = igraph_ecount(&difference);
  igraph_destroy(&difference);
  if (diff_count != 0) {
    return false;
  }

  // Compare vertex attributes
  std::vector<Atom> atoms = this->getAtoms();
  std::vector<Atom> otherAtoms = other.getAtoms();
  for (size_t i = 0; i < atoms.size(); i++) {
    if (atoms[i] != otherAtoms[i]) {
      return false;
    }
  }

  // TODO: Compare edge attributes / check whether igraph_difference
  // should take care of it
  return true;
}

bool
Universe::operator!=(const Universe& other) const
{
  return !(*this == other);
}

void
Universe::initializeFromGraph(const igraph_t* ingraph)
{
  igraph_destroy(&this->graph);
  igraph_copy(&this->graph, ingraph);
  this->NAtoms = igraph_vcount(&this->graph);
  this->NBonds = igraph_ecount(&this->graph);
  // load the ids
  igraph_vector_t allIds;
  igraph_vector_init(&allIds, this->NAtoms);
  VANV(&this->graph, "id", &allIds);
  if (igraph_cattribute_VANV(&this->graph, "id", igraph_vss_all(), &allIds)) {
    throw std::runtime_error(
      "Universes's graph's attribute id is not accessible.");
  };
  std::vector<int> ids;
  pylimer_tools::utils::igraphVectorTToStdVector(&allIds, ids);
  igraph_vector_destroy(&allIds);
  if (ids.size() == 0 && this->NAtoms > 0) {
    throw std::runtime_error(
      "Universes's graph's attribute id was not queried.");
  }
  this->atomIdToVertexIdx.reserve(ids.size());
  for (int i = 0; i < ids.size(); ++i) {
    this->atomIdToVertexIdx[ids[i]] = i;
  }
  // detect whether the graph has more than the standard attributes
  this->atomsHaveCustomAttributes = this->checkIfAtomsHaveCustomAttributes();
}

void
Universe::addAtoms(const std::vector<long int>& newIds,
                   const std::vector<int>& newTypes,
                   const std::vector<double>& newX,
                   const std::vector<double>& newY,
                   const std::vector<double>& newZ,
                   const std::vector<int>& newNx,
                   const std::vector<int>& newNy,
                   const std::vector<int>& newNz)
{
  std::unordered_map<std::string, std::vector<double>> additionalData;
  this->addAtoms(
    newIds, newTypes, newX, newY, newZ, newNx, newNy, newNz, additionalData);
}

void
Universe::addAtoms(
  const std::vector<long int>& newIds,
  const std::vector<int>& newTypes,
  const std::vector<double>& newX,
  const std::vector<double>& newY,
  const std::vector<double>& newZ,
  const std::vector<int>& newNx,
  const std::vector<int>& newNy,
  const std::vector<int>& newNz,
  const std::unordered_map<std::string, std::vector<double>>& additionalData)
{
  size_t NNewAtoms = newNy.size();
  INVALIDARG_EXP_IFN(all_equal<size_t>(8,
                                       newTypes.size(),
                                       newIds.size(),
                                       newX.size(),
                                       newNx.size(),
                                       newY.size(),
                                       newNy.size(),
                                       newZ.size(),
                                       newNz.size()),
                     "All atom inputs must have the same size.");

  if (!additionalData.empty()) {
    for (const auto& [key, value] : additionalData) {
      INVALIDARG_EXP_IFN(NNewAtoms == value.size(),
                         "Key " + key +
                           " in additional atom data has not the same length "
                           "as the other atom properties.");
    }
  }
  // actually add the vertices
  REQUIRE_IGRAPH_SUCCESS(igraph_add_vertices(&this->graph, NNewAtoms, 0));
  this->atomIdToVertexIdx.reserve(this->NAtoms + NNewAtoms);
  // do map for easy access afterwards
  // simultaneously, check that the input ids are unique
  for (size_t i = 0; i < NNewAtoms; ++i) {
    bool wasAdded =
      (this->atomIdToVertexIdx.emplace(newIds[i], this->NAtoms + i)).second;
    if (!wasAdded) {
      // remove the added atoms again
      for (size_t j = 0; j < i; ++j) {
        this->atomIdToVertexIdx.erase(newIds[j]);
      }
      throw std::invalid_argument("Atom ids must be unique; at least " +
                                  std::to_string(newIds[i]) +
                                  " is found twice.");
    }
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
    if (!additionalData.empty()) {
      this->atomsHaveCustomAttributes = true;
      for (const auto& [key, value] : additionalData) {
        pylimer_tools::utils::StdVectorToIgraphVectorT(value, &valueVec);
        igraph_cattribute_VAN_setv(&this->graph, key.c_str(), &valueVec);
      }
    }
    igraph_vector_destroy(&valueVec);
  } else {
    for (size_t i = 0; i < NNewAtoms; ++i) {
      igraph_cattribute_VAN_set(
        &this->graph, "id", this->NAtoms + i, newIds[i]);
      igraph_cattribute_VAN_set(&this->graph, "x", this->NAtoms + i, newX[i]);
      igraph_cattribute_VAN_set(&this->graph, "y", this->NAtoms + i, newY[i]);
      igraph_cattribute_VAN_set(&this->graph, "z", this->NAtoms + i, newZ[i]);
      igraph_cattribute_VAN_set(
        &this->graph, "type", this->NAtoms + i, newTypes[i]);
      igraph_cattribute_VAN_set(&this->graph, "nx", this->NAtoms + i, newNx[i]);
      igraph_cattribute_VAN_set(&this->graph, "ny", this->NAtoms + i, newNy[i]);
      igraph_cattribute_VAN_set(&this->graph, "nz", this->NAtoms + i, newNz[i]);
      if (!additionalData.empty()) {
        this->atomsHaveCustomAttributes = true;
        for (const auto& [key, value] : additionalData) {
          igraph_cattribute_VAN_set(
            &this->graph, key.c_str(), this->NAtoms + i, value[i]);
        }
      }
    }
  }
  // this->NAtoms += NNewAtoms;
  this->NAtoms = igraph_vcount(&this->graph);
}

void
Universe::replaceAtom(const long int id, const Atom& replacement)
{
  const long int vertexIdx = this->getIdxByAtomId(id);
  if (replacement.getId() != id) {
    throw std::invalid_argument("The replacement atom's id must be the same "
                                "as the one of the atom to be replaced.");
  }
  // this->atomIdToVertexIdx[replacement.getId()] = vertexIdx;
  igraph_cattribute_VAN_set(&this->graph, "x", vertexIdx, replacement.getX());
  igraph_cattribute_VAN_set(&this->graph, "y", vertexIdx, replacement.getY());
  igraph_cattribute_VAN_set(&this->graph, "z", vertexIdx, replacement.getZ());
  igraph_cattribute_VAN_set(&this->graph, "nx", vertexIdx, replacement.getNX());
  igraph_cattribute_VAN_set(&this->graph, "ny", vertexIdx, replacement.getNY());
  igraph_cattribute_VAN_set(&this->graph, "nz", vertexIdx, replacement.getNZ());
  igraph_cattribute_VAN_set(
    &this->graph, "type", vertexIdx, replacement.getType());
  for (const auto& [key, value] : replacement.getExtraData()) {
    igraph_cattribute_VAN_set(&this->graph, key.c_str(), vertexIdx, value);
    this->atomsHaveCustomAttributes = true;
  }
}

void
Universe::replaceAtomType(const long int id, const int newType)
{
  const long int vertexIdx = this->getIdxByAtomId(id);
  igraph_cattribute_VAN_set(&this->graph, "type", vertexIdx, newType);
}

void
Universe::resampleVelocities(const double mean,
                             const double variance,
                             std::string seed,
                             const bool is2d)
{
  // initialize randomness
  std::mt19937 e2;
  if (seed.empty()) {
    std::random_device rd;
    e2 = std::mt19937(rd());
  } else {
    std::seed_seq seed2(seed.begin(), seed.end());
    e2 = std::mt19937(seed2);
  }

  std::normal_distribution<double> dist(mean, variance);

  this->atomsHaveCustomAttributes = true;
  for (igraph_integer_t i = 0; i < this->NAtoms; ++i) {
    igraph_cattribute_VAN_set(&this->graph, "vx", i, dist(e2));
    igraph_cattribute_VAN_set(&this->graph, "vy", i, dist(e2));
    if (!is2d) {
      igraph_cattribute_VAN_set(&this->graph, "vz", i, dist(e2));
    }
  }
}

void
Universe::removeAtoms(const std::vector<long int>& ids)
{
  igraph_vector_int_t vertexIds;
  igraph_vector_int_init(&vertexIds, ids.size());
  for (size_t i = 0; i < ids.size(); ++i) {
    igraph_vector_int_set(&vertexIds, i, this->getIdxByAtomId(ids[i]));
  }
  igraph_delete_vertices(&this->graph, igraph_vss_vector(&vertexIds));
  igraph_vector_int_destroy(&vertexIds);

  // now, we need to update the id-atomId map
  this->atomIdToVertexIdx.clear();

  igraph_vector_t atomIds;
  igraph_vector_init(&atomIds, igraph_vcount(&this->graph));
  igraph_cattribute_VANV(&this->graph, "id", igraph_vss_all(), &atomIds);

  for (size_t i = 0; i < igraph_vector_size(&atomIds); i++) {
    long int atomId = castToIgraphInt(igraph_vector_get(&atomIds, i));
    this->atomIdToVertexIdx.emplace(atomId, i);
  }

  igraph_vector_destroy(&atomIds);

  assert(igraph_vcount(&this->graph) == this->NAtoms - ids.size());

  this->NAtoms = igraph_vcount(&this->graph);
  this->NBonds = igraph_ecount(&this->graph);
}

void
Universe::removeBonds(const std::vector<long int>& atomIdsFrom,
                      const std::vector<long int>& atomIdsTo)
{
  INVALIDARG_EXP_IFN(atomIdsFrom.size() == atomIdsTo.size(),
                     "Vertex ids from and to must have the same length.");
  for (size_t i = 0; i < atomIdsFrom.size(); ++i) {
    std::vector<igraph_integer_t> edgeIds = this->getEdgeIdsFromTo(
      this->getIdxByAtomId(atomIdsFrom[i]), this->getIdxByAtomId(atomIdsTo[i]));
    if (edgeIds.size() > 1) {
      igraph_vector_int_t edgeIdsV;
      igraph_vector_int_init(&edgeIdsV, edgeIds.size());
      pylimer_tools::utils::StdVectorToIgraphVectorT(edgeIds, &edgeIdsV);
      igraph_delete_edges(&this->graph, igraph_ess_vector(&edgeIdsV));
      igraph_vector_int_destroy(&edgeIdsV);
    } else if (edgeIds.size() == 1) {
      assert(edgeIds[0] < this->getNrOfBonds());
      igraph_delete_edges(&this->graph, igraph_ess_1(edgeIds[0]));
    }
  }

  this->NBonds = igraph_ecount(&this->graph);
}

void
Universe::removeBondsOfType(const int bondType)
{
  RUNTIME_EXP_IFN(
    igraph_cattribute_has_attr(&this->graph, IGRAPH_ATTRIBUTE_EDGE, "type"),
    "The graph does not have any bond types associated.");
  // load types
  igraph_vector_t typesVec;
  igraph_vector_init(&typesVec, this->getNrOfBonds());
  REQUIRE_IGRAPH_SUCCESS(igraph_cattribute_EANV(
    &this->graph, "type", igraph_ess_all(IGRAPH_EDGEORDER_ID), &typesVec));
  // enumerate the bonds to delete
  igraph_vector_int_t edges_to_remove;
  igraph_vector_int_init(&edges_to_remove, 0);
  for (size_t i = 0; i < igraph_vector_size(&typesVec); ++i) {
    int currentBondType = static_cast<int>(igraph_vector_get(&typesVec, i));
    if (currentBondType == bondType) {
      igraph_vector_int_push_back(&edges_to_remove, i);
    }
  }
  igraph_vector_destroy(&typesVec);

  // actually remove the bonds
  if (igraph_vector_int_size(&edges_to_remove) > 0) {
    igraph_delete_edges(&this->graph, igraph_ess_vector(&edges_to_remove));
  }
  igraph_vector_int_destroy(&edges_to_remove);
  this->NBonds = igraph_ecount(&this->graph);
}

void
Universe::addBonds(const std::vector<long int>& from,
                   const std::vector<long int>& to)
{
  this->addBonds(from.size(), from, to);
}

void
Universe::addBonds(const size_t NNewBonds,
                   const std::vector<long int>& from,
                   const std::vector<long int>& to)
{
  this->addBonds(NNewBonds, from, to, std::vector<int>());
}

void
Universe::addBonds(const std::vector<long int>& from,
                   const std::vector<long int>& to,
                   const std::vector<int>& types)
{
  this->addBonds(from.size(), from, to, types);
}

void
Universe::addBonds(const size_t NNewBonds,
                   const std::vector<long int>& from,
                   const std::vector<long int>& to,
                   const std::vector<int>& bondTypes,
                   const bool ignoreNonExistentAtoms,
                   const bool simplify)
{
  INVALIDARG_EXP_IFN(from.size() == to.size() && from.size() == NNewBonds &&
                       (bondTypes.size() == NNewBonds || bondTypes.size() == 0),
                     "All bond inputs must have the same size. Got " +
                       std::to_string(from.size()) + " atoms from, " +
                       std::to_string(to.size()) + " to, and " +
                       std::to_string(bondTypes.size()) + " types, alleged " +
                       std::to_string(NNewBonds) + ".");

  std::vector<long int> newEdgesVector =
    pylimer_tools::utils::interleave(from, to);
  size_t edgesSize = newEdgesVector.size();
  assert(edgesSize == NNewBonds * 2);
  // translate from atomId to VertexIdx
  igraph_vector_int_t newEdges;
  size_t actualNrOfBondsAdded = 0;
  igraph_vector_int_init(&newEdges, edgesSize);
  int innerIndex = 0;
  for (size_t i = 1; i < edgesSize; i += 2) {
    try {
      // two at once to throw in case one end cannot be resolved
      size_t vertexFrom = this->atomIdToVertexIdx.at(newEdgesVector[i - 1]);
      size_t vertexTo = this->atomIdToVertexIdx.at(newEdgesVector[i]);
      igraph_vector_int_set(&newEdges, innerIndex, vertexFrom);
      innerIndex += 1;
      igraph_vector_int_set(&newEdges, innerIndex, vertexTo);
      innerIndex += 1;
      actualNrOfBondsAdded += 1;
    } catch ([[maybe_unused]] std::out_of_range& ex) {
      if (!ignoreNonExistentAtoms) {
        igraph_vector_int_destroy(&newEdges);
        throw std::invalid_argument(
          "Bond with atom with id " + std::to_string(newEdgesVector[i]) +
          " and " + std::to_string(newEdgesVector[i - 1]) +
          " impossible as atom is not added yet.");
      }
    }
  }
  assert(innerIndex == 2 * actualNrOfBondsAdded);
  if (!ignoreNonExistentAtoms) {
    assert(actualNrOfBondsAdded == NNewBonds);
  }
  igraph_vector_int_resize(&newEdges, 2 * actualNrOfBondsAdded);
  // add the new edges
  REQUIRE_IGRAPH_SUCCESS(igraph_add_edges(&this->graph, &newEdges, 0));
  igraph_vector_int_destroy(&newEdges);
  if (actualNrOfBondsAdded > 0) {
    // add attributes
    if (bondTypes.size() == NNewBonds &&
        this->NBonds == (igraph_ecount(&this->graph) - NNewBonds)) {
      if (this->NBonds == 0) {
        // fast track
        igraph_vector_t types_igraph_vec;
        igraph_vector_init(&types_igraph_vec, NNewBonds);
        pylimer_tools::utils::StdVectorToIgraphVectorT(bondTypes,
                                                       &types_igraph_vec);
        REQUIRE_IGRAPH_SUCCESS(
          igraph_cattribute_EAN_setv(&this->graph, "type", &types_igraph_vec));
        igraph_vector_destroy(&types_igraph_vec);
      } else {
        for (size_t i = 0; i < NNewBonds; ++i) {
          // append attributes
          REQUIRE_IGRAPH_SUCCESS(igraph_cattribute_EAN_set(
            &this->graph, "type", this->NBonds + i, bondTypes[i]));
        }
      }
    }
    // else: too risky to add bond attributes
    // simplify graph
    // this->NBonds += NNewBonds;
    if (simplify) {
      this->simplify();
    } else {
      this->NBonds = igraph_ecount(&this->graph);
    }
  }
}

/**
 * @brief Add additional angles to the universe
 *
 * @param from
 * @param via
 * @param to
 * @param types
 */
void
Universe::addAngles(const std::vector<long int>& from,
                    const std::vector<long int>& via,
                    const std::vector<long int>& to,
                    const std::vector<int>& types)
{
  if (!all_equal<size_t>(4, from.size(), to.size(), via.size(), types.size())) {
    throw std::invalid_argument("All angle inputs must have the same size.");
  }

  this->angleFrom.insert(
    std::end(this->angleFrom), std::begin(from), std::end(from));
  this->angleVia.insert(
    std::end(this->angleVia), std::begin(via), std::end(via));
  this->angleTo.insert(std::end(this->angleTo), std::begin(to), std::end(to));
  this->angleType.insert(
    std::end(this->angleType), std::begin(types), std::end(types));

  RUNTIME_EXP_IFN(all_equal<size_t>(4,
                                    this->angleFrom.size(),
                                    this->angleTo.size(),
                                    this->angleVia.size(),
                                    this->angleType.size()),
                  "Angles' state is inconsistent.");
}

void
Universe::removeAllAngles()
{
  this->angleFrom.clear();
  this->angleVia.clear();
  this->angleTo.clear();
  this->angleType.clear();
}

/**
 * @brief Add additional angles to the universe
 *
 * @param from
 * @param via1
 * @param via2
 * @param to
 * @param types
 */
void
Universe::addDihedralAngles(const std::vector<long int>& from,
                            const std::vector<long int>& via1,
                            const std::vector<long int>& via2,
                            const std::vector<long int>& to,
                            const std::vector<int>& types)
{
  INVALIDARG_EXP_IFN(
    all_equal<size_t>(
      5, from.size(), to.size(), via1.size(), via2.size(), types.size()),
    "All dihedral angle inputs must have the same size.");

  this->dihedralAngleFrom.insert(
    std::end(this->dihedralAngleFrom), std::begin(from), std::end(from));
  this->dihedralAngleVia1.insert(
    std::end(this->dihedralAngleVia1), std::begin(via1), std::end(via1));
  this->dihedralAngleVia2.insert(
    std::end(this->dihedralAngleVia2), std::begin(via2), std::end(via2));
  this->dihedralAngleTo.insert(
    std::end(this->dihedralAngleTo), std::begin(to), std::end(to));
  this->dihedralAngleType.insert(
    std::end(this->dihedralAngleType), std::begin(types), std::end(types));
}

void
Universe::removeAllDihedralAngles()
{
  this->dihedralAngleFrom.clear();
  this->dihedralAngleVia1.clear();
  this->dihedralAngleVia2.clear();
  this->dihedralAngleTo.clear();
  this->dihedralAngleType.clear();
}

void
Universe::inferCoordinates(const int crosslinkerType)
{
  std::vector<Molecule> molecules = this->getMolecules(crosslinkerType);
  for (const Molecule& molecule : molecules) {
    std::vector<Atom> linedUpAtoms =
      molecule.getAtomsLinedUp(crosslinkerType, true, false);
    for (const Atom& atom : linedUpAtoms) {
      this->replaceAtom(atom.getId(), atom);
    }
  }
};

/**
 * @brief Simplify the underlying graph by removing double edges etc.
 *
 */
void
Universe::simplify()
{
  igraph_attribute_combination_t comb;
  igraph_attribute_combination_init(&comb);
  // how to combine two edges and their attributes
  // currently, only the type attribute exists.
  // let's take the mean.
  REQUIRE_IGRAPH_SUCCESS(igraph_attribute_combination_add(
    &comb, NULL, IGRAPH_ATTRIBUTE_COMBINE_MEAN, NULL));
  igraph_simplify(&this->graph, /*multiple=*/1, /*loops=*/1, &comb);
  igraph_attribute_combination_destroy(&comb);
  this->NBonds = igraph_ecount(&this->graph);
}

/**
 * @brief Count how many of each atom type there are in the universe
 *
 * @return std::map<int, int>
 */
std::map<int, int>
Universe::countAtomTypes() const
{
  std::vector<int> atomTypes = this->getAtomTypes();
  std::map<int, int> result;
  for (int atomType : atomTypes) {
    result[atomType] += 1;
  }
  return result;
}

/**
 * @brief Count the number of atoms within a certain distance.
 *
 * @param distances
 * @param unwrapped
 * @return std::vector<size_t>
 */
std::vector<size_t>
Universe::countAtomsInSkinDistance(const std::vector<double>& distances,
                                   const bool unwrapped) const
{
  if (distances.size() <= 1) {
    return std::vector<size_t>();
  }
  std::vector<size_t> result =
    pylimer_tools::utils::initializeWithValue<size_t>(distances.size() - 1, 0);

  // first, validate the input distances
  INVALIDARG_EXP_IFN(distances[0] >= 0.0, "Distances must be positive.");
  for (size_t i = 1; i < distances.size(); ++i) {
    if (distances[i] <= distances[i - 1]) {
      throw std::invalid_argument(
        "Distances must be increasing and unique, got " +
        std::to_string(distances[i]) + " and " +
        std::to_string(distances[i - 1]) + " at i = " + std::to_string(i) +
        " and i-1.");
    }
  }

  // then, start processing using a neighbour list
  const std::vector<Atom> atoms = this->getAtoms();
  NeighbourList neighbourList = NeighbourList(
    atoms, this->getBox(), pylimer_tools::utils::last<double>(distances));

  for (size_t i = 1; i < distances.size(); ++i) {
    for (const Atom& a : atoms) {
      std::vector<Atom> closeAtoms = neighbourList.getAtomsCloseTo(
        a, distances[i], distances[i - 1], unwrapped, true);
      // std::cout << closeAtoms.size() << " atoms close to " << a.getId()
      //           << " between " << distances[i] << " and " << distances[i -
      //           1] << std::endl;
      result[i - 1] += closeAtoms.size();
    }
  }

  return result;
}

/**
 * @brief Set the masses of the atoms in this universe
 *
 * @param atomMassPerType the weight per type
 */
void
Universe::setMasses(const std::map<int, double>& atomMassPerType)
{
  this->massPerType = atomMassPerType;
}

void
Universe::setMassForType(const int atomType, const double mass)
{
  this->massPerType[atomType] = mass;
}

std::map<int, double>
Universe::getMasses()
{
  return this->massPerType;
};

/**
 * @brief Get the standalone components of the network
 *
 * @return std::vector<Universe>
 */
std::vector<Universe>
Universe::getClusters() const
{
  std::vector<Universe> clusters;
  if (this->getNrOfAtoms() == 0) {
    return clusters;
  }

  // split the copy into the separate components
  igraph_graph_list_t components;
  igraph_graph_list_init(&components, 0);
  REQUIRE_IGRAPH_SUCCESS(
    igraph_decompose(&graph, &components, IGRAPH_WEAK, -1, 0));
  size_t NComponents = igraph_graph_list_size(&components);
  // std::cout << NComponents << " clusters found." << std::endl;
  clusters.reserve(NComponents);
  for (size_t i = 0; i < NComponents; ++i) {
    // make the molecule the owner of the graph
    igraph_t* g = igraph_graph_list_get_ptr(&components, i);

    if (igraph_vcount(g)) {
      Universe newUniverse = Universe(this->box);
      newUniverse.initializeFromGraph(g);
      newUniverse.setMasses(this->massPerType);
      // TODO: add angles, dihedrals etc.
      clusters.push_back(newUniverse);
    }

    igraph_destroy(g);
  }
  igraph_graph_list_destroy(&components);
  return clusters;
}

StrandAffiliationInfo
Universe::getStrandAffiliation(const int crosslinkerType) const
{
  // reserve results
  StrandAffiliationInfo result;
  result.strandIdOfVertex = pylimer_tools::utils::initializeWithValue<long int>(
    this->getNrOfAtoms(), -1);
  result.indexOfVertexInStrand =
    pylimer_tools::utils::initializeWithValue<long int>(this->getNrOfAtoms(),
                                                        -1);

  // find possible starting points for the depth-first search
  std::vector<int> vertexAtomTypes = this->getPropertyValues<int>("type");
  std::vector<int> vertexDegree = this->getVertexDegrees();
  assert(vertexAtomTypes.size() == vertexDegree.size());
  assert(vertexAtomTypes.size() == this->getNrOfAtoms());
  std::vector<bool> vertexIsStartingPoint =
    pylimer_tools::utils::initializeWithValue<bool>(this->getNrOfAtoms(),
                                                    false);
  // assemble the starting points
  std::vector<igraph_integer_t> startingPoints;
  for (size_t i = 0; i < this->getNrOfAtoms(); ++i) {
    if (vertexAtomTypes[i] == crosslinkerType || vertexDegree[i] != 2) {
      vertexIsStartingPoint[i] = true;
      startingPoints.push_back(i);
    }
  }
  // variables to count
  long int currentStrandId = 0;
  long int currentIndexOfVertexInStrand = 0;

  // do the depth-first search
  igraph_adjlist_t adjlist;
  igraph_adjlist_init(
    &this->graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE);
  std::vector<bool> visited = pylimer_tools::utils::initializeWithValue<bool>(
    this->getNrOfAtoms(), false);
  std::stack<long int> verticesToVisit;
  for (size_t idxInStartingPoints = 0;
       idxInStartingPoints < startingPoints.size();
       ++idxInStartingPoints) {
    verticesToVisit.push(startingPoints[idxInStartingPoints]);
  }
  while (!verticesToVisit.empty()) {
    const long int currentVertex = verticesToVisit.top();
    visited[currentVertex] = true;
    // TODO: this fails for primary ends, since they will be their own strand.
    // while we want them as starting points, they should only increase the
    // strandId if actually used as starting points
    if (vertexIsStartingPoint[currentVertex]) {
      currentStrandId += 1;
      currentIndexOfVertexInStrand = 0;
    }
    // by using these starting points,
    // we can guarantee that the search will always be along a strand
    if (vertexDegree[currentVertex] <= 2 &&
        vertexAtomTypes[currentVertex] != crosslinkerType) {
      result.strandIdOfVertex[currentVertex] = currentStrandId;
      result.indexOfVertexInStrand[currentVertex] =
        currentIndexOfVertexInStrand;
    }
    // add unvisited neighbours to stack to visit
    size_t nUnvisitedNeighbours = 0;
    const igraph_vector_int_t* neighbors =
      igraph_adjlist_get(&adjlist, currentVertex);
    for (long int j = 0; j < igraph_vector_int_size(neighbors); ++j) {
      if (!visited[igraph_vector_int_get(neighbors, j)]) {
        verticesToVisit.push(igraph_vector_int_get(neighbors, j));
        nUnvisitedNeighbours += 1;
      }
    }
    if (nUnvisitedNeighbours == 0) {
      verticesToVisit.pop();
    }
  }

  igraph_adjlist_destroy(&adjlist);

  return result;
}

/**
 * @brief Decompose this universe into chains/molecules by splitting them into
 * clusters
 *
 * @param atomTypeToOmit The atom type to remove to get more clusters
 * @return std::vector<Molecule>
 */
std::vector<Molecule>
Universe::getMolecules(const int atomTypeToOmit) const
{
  std::vector<Molecule> molecules;
  if (this->getNrOfAtoms() == 0) {
    return molecules;
  }
  // make a copy to remove crossLinkers from
  igraph_t graphWithoutCrosslinkers;
  REQUIRE_IGRAPH_SUCCESS(igraph_copy(&graphWithoutCrosslinkers, &this->graph));

  // select vertices of type
  std::vector<igraph_integer_t> indicesToRemove =
    this->getIndicesWithAttribute("type", atomTypeToOmit);
  std::sort(indicesToRemove.rbegin(), indicesToRemove.rend());
  if (indicesToRemove.size() > 0) {
    igraph_vs_t verticesToRemove = this->getVerticesByIndices(indicesToRemove);

    // remove elements of type
    REQUIRE_IGRAPH_SUCCESS(
      igraph_delete_vertices(&graphWithoutCrosslinkers, verticesToRemove));

    igraph_vs_destroy(&verticesToRemove);
  }

  // split the copy into the separate components
  igraph_graph_list_t components;
  igraph_graph_list_init(&components, 1);
  REQUIRE_IGRAPH_SUCCESS(igraph_decompose(
    &graphWithoutCrosslinkers, &components, IGRAPH_WEAK, -1, 0));
  size_t NComponents = igraph_graph_list_size(&components);
  // std::cout << NComponents << " molecules found. Removed " <<
  // indicesToRemove.size()
  //           << " vertices. Size now: " <<
  //           igraph_vcount(&graphWithoutCrosslinkers) << " atoms with " <<
  //           igraph_ecount(&graphWithoutCrosslinkers) << " bonds." <<
  //           std::endl;
  molecules.reserve(NComponents);
  for (size_t i = 0; i < NComponents; ++i) {
    // make the molecule the owner of the graph
    igraph_t* g = igraph_graph_list_get_ptr(&components, i);

    if (igraph_vcount(g)) {
      molecules.push_back(
        Molecule(this->box, g, MoleculeType::UNDEFINED, this->massPerType));
    }
  }
  igraph_graph_list_destroy(&components);
  igraph_destroy(&graphWithoutCrosslinkers);
  return molecules;
}

/**
 * @brief Get a vertex selector to find all vertices with a certain type
 *
 * @param type the type to select
 * @return igraph_vs_t
 */
igraph_vs_t
Universe::getVerticesOfType(const int type) const
{
  std::vector<igraph_integer_t> indices =
    this->getIndicesWithAttribute("type", type);
  return this->getVerticesByIndices(indices);
}

/**
 * @brief Get a vertex selector to find all vertices with a certain index
 *
 * @param indices the vertex indices to select
 * @return igraph_vs_t
 */
igraph_vs_t
Universe::getVerticesByIndices(std::vector<igraph_integer_t> indices) const
{
  igraph_vector_int_t indicesToSelect;
  igraph_vector_int_init(&indicesToSelect, indices.size());
  pylimer_tools::utils::StdVectorToIgraphVectorT(indices, &indicesToSelect);
  igraph_vs_t result;
  REQUIRE_IGRAPH_SUCCESS(igraph_vs_vector_copy(&result, &indicesToSelect));
  igraph_vector_int_destroy(&indicesToSelect);
  return result;
}

/**
 * @brief Decompose the network into clusters, re-adding the atoms omitted to
 * get more clusters
 *
 * @param crossLinkerType the type of the atoms to omit and re-add
 * @return std::vector<Molecule>
 */
std::vector<Molecule>
Universe::getChainsWithCrosslinker(const int crossLinkerType) const
{
  std::vector<Molecule> molecules;
  if (this->getNrOfAtoms() == 0) {
    return molecules;
  }
  std::vector<bool> vertexIsJunction =
    pylimer_tools::utils::initializeWithValue(this->getNrOfAtoms(), false);
  // make a copy to remove junctions from
  igraph_t graphWithoutJunctions;
  REQUIRE_IGRAPH_SUCCESS(igraph_copy(&graphWithoutJunctions, &this->graph));
  // select vertices of junctions type
  std::vector<int> vertexDegrees = this->getVertexDegrees();
  std::vector<int> vertexTypes = this->getPropertyValues<int>("type");
  std::vector<igraph_integer_t> indicesToRemove;
  for (long int vIdx = 0; vIdx < this->getNrOfAtoms(); ++vIdx) {
    if (vertexDegrees[vIdx] > 2 || vertexTypes[vIdx] == crossLinkerType) {
      vertexIsJunction[vIdx] = true;
      indicesToRemove.push_back(vIdx);
    }
  }
  std::sort(indicesToRemove.rbegin(), indicesToRemove.rend());
  if (indicesToRemove.size() > 0) {
    igraph_vs_t verticesToRemove = this->getVerticesByIndices(indicesToRemove);

    // remove elements of type
    REQUIRE_IGRAPH_SUCCESS(
      igraph_delete_vertices(&graphWithoutJunctions, verticesToRemove));

    igraph_vs_destroy(&verticesToRemove);
    RUNTIME_EXP_IFN(igraph_vcount(&graphWithoutJunctions) ==
                      this->getNrOfAtoms() - indicesToRemove.size(),
                    "Expected all junctions to be removed from the graph.");
  }

  // load the properties that will need to be copied
  std::pair<std::vector<std::string>, std::vector<std::string>>
    vertexAndEdgeProperties = this->getVertexAndEdgePropertyNames();

  // split the copy into the separate components
  igraph_graph_list_t components;
  igraph_graph_list_init(&components, 1);
  if (igraph_decompose(
        &graphWithoutJunctions, &components, IGRAPH_STRONG, -1, 0)) {
    throw std::runtime_error("Failed to decompose graph.");
  }
  size_t NComponents = igraph_graph_list_size(&components);
  molecules.reserve(NComponents);
  for (size_t i = 0; i < NComponents; ++i) {
    // loop the chains to add the crossLinkers back
    igraph_t* chain = igraph_graph_list_get_ptr(&components, i);
    int moleculeLengthBefore = igraph_vcount(chain);
    // also select ones of degree 0 for dangling atoms
    std::vector<long int> endNodeIndices =
      this->getVerticesWithDegree(chain, { { 0, 1 } });
    igraph_vector_int_t endNodeSelectorVector;
    igraph_vector_int_init(&endNodeSelectorVector, endNodeIndices.size());
    pylimer_tools::utils::StdVectorToIgraphVectorT(endNodeIndices,
                                                   &endNodeSelectorVector);
    MoleculeType molType = MoleculeType::UNDEFINED;
    bool isLoop = false;

    if (moleculeLengthBefore >= 1) {
      // no care about single atoms for now (though they are included)
      igraph_vit_t endNodeVit;
      igraph_vit_create(
        chain, igraph_vss_vector(&endNodeSelectorVector), &endNodeVit);
      // this check is in case the primary loop is "free"
      isLoop = IGRAPH_VIT_SIZE(endNodeVit) == 0 &&
               (moleculeLengthBefore > 1 || igraph_ecount(chain) > 0);
      // collect atoms to add, since adding them would invalidate the iterator
      std::vector<long int> atomsToAdd;
      std::vector<std::vector<long int>> bondsToAdd;
      // loop end nodes
      while (!IGRAPH_VIT_END(endNodeVit)) {
        long int newEndNodeVertexId = (long int)IGRAPH_VIT_GET(endNodeVit);
        long int oldEndNodeId =
          igraphRealToInt<long int>(VAN(chain, "id", newEndNodeVertexId));
        long int originalEndNodeVertexId =
          this->atomIdToVertexIdx.at(oldEndNodeId);
        // this->findVertexIdForProperty("id", oldEndNodeId);
        igraph_vector_int_t neighbors;
        igraph_vector_int_init(&neighbors, 0);

        REQUIRE_IGRAPH_SUCCESS(igraph_neighbors(
          &graph, &neighbors, originalEndNodeVertexId, IGRAPH_ALL));

        // loop neighbors of this end node
        for (igraph_integer_t neighIdx = 0;
             neighIdx < igraph_vector_int_size(&neighbors);
             ++neighIdx) {
          igraph_integer_t neighbourOriginalId =
            igraph_vector_int_get(&neighbors, neighIdx);

          if (vertexIsJunction[neighbourOriginalId]) {
            // found a crosslinker neighbour
            long int originalNeighbourAtomId =
              igraph_cattribute_VAN(&graph, "id", neighbourOriginalId);
            atomsToAdd.push_back(neighbourOriginalId);
            bondsToAdd.push_back({ {
              originalEndNodeVertexId,
              newEndNodeVertexId,
              static_cast<long int>(neighbourOriginalId),
            } });
          }
        }

        IGRAPH_VIT_NEXT(endNodeVit);
        igraph_vector_int_destroy(&neighbors);
      } // loop end nodes

      if (atomsToAdd.size() == 2 && atomsToAdd[0] == atomsToAdd[1]) {
        isLoop = true;
        // we only want to add it once -> remove
        atomsToAdd.pop_back();
      }

      std::unordered_map<long int, long int> newAtomsMap;
      // actually add the atoms...
      for (auto atomToAddOriginalIdx : atomsToAdd) {
        REQUIRE_IGRAPH_SUCCESS(igraph_add_vertices(chain, 1, nullptr));
        long int newCrosslinkerVertexIdx = igraph_vcount(chain) - 1;
        newAtomsMap.insert_or_assign(atomToAddOriginalIdx,
                                     newCrosslinkerVertexIdx);

        // deprecated additional loop check
        long int originalNeighbourAtomId = igraphRealToInt<long int>(
          igraph_cattribute_VAN(&graph, "id", atomToAddOriginalIdx));
        assert(originalNeighbourAtomId == 0 ||
               !pylimer_tools::utils::graphHasVertexWithProperty<long int>(
                 chain, "id", originalNeighbourAtomId));

        // including all attributes
        pylimer_tools::utils::copyVertexProperties(
          &this->graph,
          atomToAddOriginalIdx,
          chain,
          newCrosslinkerVertexIdx,
          vertexAndEdgeProperties.first);
      }
      // ...and bonds
      for (auto bond : bondsToAdd) {
        REQUIRE_IGRAPH_SUCCESS(
          igraph_add_edge(chain, bond[1], newAtomsMap.at(bond[2])));
        // also copy the bond attributes
        std::vector<igraph_integer_t> oldEIds =
          this->getEdgeIdsFromTo(bond[0], bond[2]);
        RUNTIME_EXP_IFN(oldEIds.size() > 0,
                        "Expected at least one edge between the same atoms.");
        pylimer_tools::utils::copyEdgeProperties(
          &this->graph,
          oldEIds[0],
          chain,
          igraph_ecount(chain) - 1,
          vertexAndEdgeProperties.second);
      }
      igraph_vit_destroy(&endNodeVit);
    } // if molecule length
    igraph_vector_int_destroy(&endNodeSelectorVector);
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
    molecules.emplace_back(this->box, chain, molType, this->massPerType);
  }
  // what was neglected so far: springs between junctions!
  for (long int junctionIdx : indicesToRemove) {
    std::vector<igraph_integer_t> connections =
      this->getVertexIdxsConnectedTo(junctionIdx);
    std::vector<igraph_integer_t> edgeIds =
      this->getIncidentEdgeIds(junctionIdx);
    assert(connections.size() == edgeIds.size());
    igraph_attribute_combination_t comb;
    igraph_attribute_combination_init(&comb);
    // how to combine two edges and their attributes
    // currently, only the type attribute exists.
    // let's take the mean.
    REQUIRE_IGRAPH_SUCCESS(igraph_attribute_combination_add(
      &comb, nullptr, IGRAPH_ATTRIBUTE_COMBINE_MEAN, nullptr));
    bool hasPrimary = false;
    for (size_t i = 0; i < connections.size(); ++i) {
      long int connectedVertexIdx = connections[i];
      if (connectedVertexIdx == junctionIdx) {
        hasPrimary = true;
      }
      if (connectedVertexIdx > junctionIdx &&
          vertexIsJunction[connectedVertexIdx]) {
        // new Molecule from just these two junctions
        MoleculeType molType = MoleculeType::UNDEFINED;
        if (this->getVertexDegree(junctionIdx) == 1 &&
            this->getVertexDegree(connectedVertexIdx) == 1) {
          molType = MoleculeType::FREE_CHAIN;
        } else if (this->getVertexDegree(junctionIdx) > 1 &&
                   this->getVertexDegree(connectedVertexIdx) > 1) {
          molType = MoleculeType::NETWORK_STRAND;
        } else if (this->getVertexDegree(junctionIdx) > 1 ||
                   this->getVertexDegree(connectedVertexIdx) > 1) {
          molType = MoleculeType::DANGLING_CHAIN;
        }

        // could also use igraph_induced_subgraph, but performance is very
        // bad for large structures, given how this is more or less O(1)
        igraph_t chain;
        REQUIRE_IGRAPH_SUCCESS(igraph_empty(&chain, 2, IGRAPH_UNDIRECTED));
        pylimer_tools::utils::copyVertexProperties(
          &this->graph, junctionIdx, &chain, 0, vertexAndEdgeProperties.first);
        pylimer_tools::utils::copyVertexProperties(
          &this->graph,
          connectedVertexIdx,
          &chain,
          1,
          vertexAndEdgeProperties.first);
        REQUIRE_IGRAPH_SUCCESS(igraph_add_edge(&chain, 0, 1));
        pylimer_tools::utils::copyEdgeProperties(
          &this->graph, edgeIds[i], &chain, 0, vertexAndEdgeProperties.second);
        molecules.emplace_back(this->box, &chain, molType, this->massPerType);
        igraph_destroy(&chain);
      }
    }

    if (hasPrimary) {
      // new Molecule from just these two junctions
      std::vector<igraph_integer_t> primaryLoopEdgeIdsVec =
        this->getEdgeIdsFromTo(junctionIdx, junctionIdx);
      std::set<long int> primaryLoopEdgeIds(primaryLoopEdgeIdsVec.begin(),
                                            primaryLoopEdgeIdsVec.end());
      size_t nAdditionalChains = primaryLoopEdgeIds.size();
      assert(nAdditionalChains >= 1);

      MoleculeType molType = MoleculeType::PRIMARY_LOOP;

      for (long int edgeId : primaryLoopEdgeIds) {
        // could also use igraph_induced_subgraph, but performance is very
        // bad for large structures, given how this is more or less O(1)
        igraph_t chain;
        REQUIRE_IGRAPH_SUCCESS(igraph_empty(&chain, 1, IGRAPH_UNDIRECTED));
        pylimer_tools::utils::copyVertexProperties(
          &this->graph, junctionIdx, &chain, 0, vertexAndEdgeProperties.first);
        REQUIRE_IGRAPH_SUCCESS(igraph_add_edge(&chain, 0, 0));
        pylimer_tools::utils::copyEdgeProperties(
          &this->graph, edgeId, &chain, 0, vertexAndEdgeProperties.second);
        molecules.emplace_back(this->box, &chain, molType, this->massPerType);
        igraph_destroy(&chain);
      }
    }

    igraph_attribute_combination_destroy(&comb);
  }

  igraph_graph_list_destroy(&components);
  igraph_destroy(&graphWithoutJunctions);

  return molecules;
};

/**
 * Identify vertices associated with dangling and free chains.
 *
 * Note: There are special cases that fail here, in particular,
 * if you have a primary or secondary loops without internal beads,
 * sandwiched between two dangling chains.
 *
 * @return a list of molecule types,
 * for each vertex in the universe the anticipated type
 * of the corresponding molecule, if the type is either dangling or free chain.
 */
std::vector<pylimer_tools::entities::MoleculeType>
Universe::identifyObviouslyDanglingAtoms(const bool distinguishFree) const
{
  std::vector<pylimer_tools::entities::MoleculeType> result =
    pylimer_tools::utils::initializeWithValue(
      this->getNrOfAtoms(), pylimer_tools::entities::MoleculeType::UNDEFINED);

  // it works like this:
  // we find the definitive ends of dangling chains, namely f = 0 and 1 atoms.
  // from these ends, we walk along the chain, marking each encountered atom
  // as visited and dangling, iff it's not a junction.
  // it it's a junction, we require all but one of the neighbours to be
  // dangling, in order to mark the junction as dangling as well, and progress
  // with the walk down the remaining chain.
  std::vector<bool> vertexVisited =
    pylimer_tools::utils::initializeWithValue(this->getNrOfAtoms(), false);
  std::vector<bool> edgeVisited =
    pylimer_tools::utils::initializeWithValue(this->getNrOfBonds(), false);
  std::vector<int> vertexDegrees = this->getVertexDegrees();
  igraph_lazy_adjlist_t adjlist;
  REQUIRE_IGRAPH_SUCCESS(igraph_lazy_adjlist_init(
    &this->graph, &adjlist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, IGRAPH_MULTIPLE));
  igraph_vector_int_t* neighbors = nullptr;
  igraph_lazy_inclist_t inclist;
  REQUIRE_IGRAPH_SUCCESS(igraph_lazy_inclist_init(
    &this->graph, &inclist, IGRAPH_ALL, IGRAPH_LOOPS_TWICE));
  igraph_vector_int_t* neighbourEdges = nullptr;

  for (igraph_integer_t startingVertexIdx = 0;
       startingVertexIdx < this->getNrOfAtoms();
       ++startingVertexIdx) {
    if (vertexDegrees[startingVertexIdx] > 1) {
      continue;
    }

    // free bead; we don't mark it as visited, but could in principle
    if (vertexDegrees[startingVertexIdx] == 0) {
      result[startingVertexIdx] =
        pylimer_tools::entities::MoleculeType::FREE_CHAIN;
      continue;
    }

    // finally an f = 1 bead we can start walking from, down the chain
    igraph_integer_t currentVertexIdx = startingVertexIdx;
    // actually start walking down the chain
    do {
      assert(!vertexVisited[currentVertexIdx]);
      vertexVisited[currentVertexIdx] = true;
      neighbors = igraph_lazy_adjlist_get(&adjlist, currentVertexIdx);
      neighbourEdges = igraph_lazy_inclist_get(&inclist, currentVertexIdx);
      assert(vertexDegrees[currentVertexIdx] ==
             igraph_vector_int_size(neighbors));
      assert(vertexDegrees[currentVertexIdx] ==
             igraph_vector_int_size(neighbourEdges));
      size_t nVisitedNeighbours = 0;
      size_t nVisitedNeighbourEdges = 0;
      size_t nDanglingNeighbours = 0;
      igraph_integer_t unvisitedNeighbour = -1;
      igraph_integer_t unvisitedNeighbourEdge = -1;
      // check the neighbours
      for (igraph_integer_t neighbourIdx = 0;
           neighbourIdx < igraph_vector_int_size(neighbors);
           ++neighbourIdx) {
        igraph_integer_t neighbourVertex =
          igraph_vector_int_get(neighbors, neighbourIdx);
        igraph_integer_t neighbourEdge =
          igraph_vector_int_get(neighbourEdges, neighbourIdx);

#ifndef NDEBUG
        // validate that the inclist & adjlist are actually compliant
        // by checking that such an edge exists in the graph
        igraph_integer_t a = 0;
        igraph_integer_t b = 0;
        igraph_edge(&graph, neighbourEdge, &a, &b);
        assert(a == currentVertexIdx || b == currentVertexIdx);
        assert(a == neighbourVertex || b == neighbourVertex);
#endif

        nVisitedNeighbours += vertexVisited[neighbourVertex];
        nVisitedNeighbourEdges += edgeVisited[neighbourEdge];
        nDanglingNeighbours +=
          (result[neighbourVertex] ==
           pylimer_tools::entities::MoleculeType::DANGLING_CHAIN);
        // we have to check both edges and vertices, since we could have,
        // for example, a dangling chain connected to a primary loop,
        // connected to a dangling chain, in which case we do not want to mark
        // the primary loop as dangling
        if (!vertexVisited[neighbourVertex] && !edgeVisited[neighbourEdge]) {
          unvisitedNeighbour = neighbourVertex;
          unvisitedNeighbourEdge = neighbourEdge;
        }
      }
      // if we found mostly dangling neighbour, we mark this vertex as dangling
      // as well
      if (nDanglingNeighbours >= vertexDegrees[currentVertexIdx] - 1) {
        result[currentVertexIdx] =
          pylimer_tools::entities::MoleculeType::DANGLING_CHAIN;
      }
      // then, find where to continue the walk
      if ((nVisitedNeighbourEdges == vertexDegrees[currentVertexIdx] - 1) &&
          (nVisitedNeighbours == vertexDegrees[currentVertexIdx] - 1)) {
        assert(unvisitedNeighbour != -1);
        assert(unvisitedNeighbourEdge != -1);
        currentVertexIdx = unvisitedNeighbour;
        edgeVisited[unvisitedNeighbourEdge] = true;
      } else {
        // we can't walk further from this vertex,
        // we need to go back to a start
        // however, we may revisit this vertex later, if needed
        vertexVisited[currentVertexIdx] = false;
        break;
      }
    } while (true);
  }

  // then, when we have found all obviously dangling atoms,
  // some of them might be free chains instead.
  // We can find those by checking whether for a cluster of vertices
  // all vertices are considered dangling.
  if (distinguishFree) {
    // reset visited flags
    edgeVisited =
      pylimer_tools::utils::initializeWithValue(this->getNrOfBonds(), false);
    vertexVisited =
      pylimer_tools::utils::initializeWithValue(this->getNrOfAtoms(), false);
    std::vector<igraph_integer_t> verticesOfCluster;
    verticesOfCluster.reserve(0.01 * this->getNrOfAtoms());

    for (igraph_integer_t startingVertexIdx = 0;
         startingVertexIdx < this->getNrOfAtoms();
         ++startingVertexIdx) {
      if (vertexVisited[startingVertexIdx]) {
        continue;
      }
      if (vertexDegrees[startingVertexIdx] == 1) {
        bool isFree = true;
        // start searching all connected vertices
        std::stack<igraph_integer_t> vertexSearchStack;
        vertexSearchStack.push(startingVertexIdx);
        do {
          igraph_integer_t currentVertex = vertexSearchStack.top();
          vertexSearchStack.pop();
          vertexVisited[currentVertex] = true;
          verticesOfCluster.push_back(currentVertex);
          const igraph_vector_int_t* neighbours =
            igraph_lazy_adjlist_get(&adjlist, currentVertex);
          for (igraph_integer_t neighbourIdx = 0;
               neighbourIdx < igraph_vector_int_size(neighbours);
               ++neighbourIdx) {
            igraph_integer_t neighbourVertex =
              igraph_vector_int_get(neighbours, neighbourIdx);
            if (!vertexVisited[neighbourVertex]) {
              vertexSearchStack.push(neighbourVertex);
            }
            if (result[neighbourVertex] !=
                pylimer_tools::entities::MoleculeType::DANGLING_CHAIN) {
              isFree = false;
              break;
            }
          }
        } while (!vertexSearchStack.empty() && isFree);

        // propagate this change from dangling to free chains for all connected
        // vertices
        if (isFree) {
          for (const igraph_integer_t vertex : verticesOfCluster) {
            result[vertex] = pylimer_tools::entities::MoleculeType::FREE_CHAIN;
          }
          verticesOfCluster.clear();
        }
      }
    }
  }

  igraph_lazy_adjlist_destroy(&adjlist);
  igraph_lazy_inclist_destroy(&inclist);

  return result;
}

/**
 * @brief Detect loops (cycles) in the graph
 *
 * @param crossLinkerType
 * @param maxLength
 * @param skipSelfLoops
 * @param edges
 * @return std::vector<std::vector<long int>> the list of vertices in each
 * loop
 */
std::vector<std::vector<igraph_integer_t>>
Universe::findLoops(const int crossLinkerType,
                    const int maxLength,
                    bool skipSelfLoops,
                    std::vector<std::vector<igraph_integer_t>>* edges) const
{
  // NOTE: There are exponentially many paths between two vertices of a graph,
  // and you may run out of memory when using this function, if your graph is
  // lattice-like.

  // TODO: shrink the graph by collapsing all 2-functional vertices?

  igraph_vector_int_list_t v_results;
  igraph_vector_int_list_init(&v_results, 0);
  igraph_vector_int_list_t e_results;
  igraph_vector_int_list_init(&e_results, 0);

  REQUIRE_IGRAPH_SUCCESS(igraph_simple_cycles(
    &this->graph, &v_results, &e_results, IGRAPH_ALL, -1, maxLength));

  std::vector<std::vector<igraph_integer_t>> results;

  for (igraph_integer_t listIdx = 0;
       listIdx < igraph_vector_int_list_size(&v_results);
       ++listIdx) {
    std::vector<igraph_integer_t> currentPath;
    for (igraph_integer_t i = 0;
         i < igraph_vector_int_size(
               igraph_vector_int_list_get_ptr(&v_results, listIdx));
         ++i) {
      currentPath.push_back(igraph_vector_int_get(
        igraph_vector_int_list_get_ptr(&v_results, listIdx), i));
    }
    results.push_back(currentPath);
  }

  if (edges != nullptr) {
    edges->clear();
    for (igraph_integer_t listIdx = 0;
         listIdx < igraph_vector_int_list_size(&e_results);
         ++listIdx) {
      std::vector<igraph_integer_t> currentPath;
      for (igraph_integer_t i = 0;
           i < igraph_vector_int_size(
                 igraph_vector_int_list_get_ptr(&e_results, listIdx));
           ++i) {
        currentPath.push_back(igraph_vector_int_get(
          igraph_vector_int_list_get_ptr(&e_results, listIdx), i));
      }
      edges->push_back(currentPath);
    }
  }

  igraph_vector_int_list_destroy(&v_results);
  igraph_vector_int_list_destroy(&e_results);

  return results;
}

/**
 * @brief Detect loops (cycles) in the graph and count how long they are
 *
 * @param maxLength
 * @return std::map<int, int> the number of vertices per loop
 */
std::unordered_map<int, int>
Universe::countLoopLengths(const int maxLength) const
{
  // NOTE: There are exponentially many paths between two vertices of a graph,
  // and you may run out of memory when using this function, if your graph is
  // lattice-like.

  pylimer_tools::utils::Counter<int> loopLengths;

  // TODO: shrink the graph by collapsing all 2-functional vertices?
  REQUIRE_IGRAPH_SUCCESS(igraph_simple_cycles_callback(
    &this->graph, IGRAPH_ALL, -1, maxLength, &count_found_cycle, &loopLengths));

  return loopLengths.asMap();
}

/**
 * @brief Find the loops in the network
 *
 * NOTE: There are exponentially many paths between two vertices of a graph,
 * and you may run out of memory when using this function, if your graph is
 * lattice-like.
 *
 * @param crossLinkerType
 * @param maxLength
 * @return std::map<int, std::vector<std::vector<Atom>>>
 */
std::map<int, std::vector<std::vector<Atom>>>
Universe::findLoopsOfAtoms(const int crossLinkerType,
                           const int maxLength,
                           const bool skipSelfLoops) const
{
  // NOTE: There are exponentially many paths between two vertices of a graph,
  // and you may run out of memory when using this function, if your graph is
  // lattice-like.
  std::map<int, std::vector<std::vector<Atom>>> results;

  std::vector<std::vector<igraph_integer_t>> loops =
    this->findLoops(crossLinkerType, maxLength, skipSelfLoops);

  for (size_t i = 0; i < loops.size(); ++i) {
    std::vector<igraph_integer_t> loop = loops[i];
    std::vector<Atom> currentPath;
    currentPath.reserve(loop.size());
    int currentFunctionality = 0;
    for (size_t j = 0; j < loop.size(); ++j) {
      Atom newAtom = this->getAtomByVertexIdx(loop[j]);
      if (newAtom.getType() == crossLinkerType) {
        currentFunctionality += 1;
      }
      currentPath.push_back(newAtom);
    }
    results[currentFunctionality].push_back(currentPath);
  }

  return results;
};

/**
 * @brief Find the loops in the network starting with one connection
 *
 * NOTE: There are exponentially many paths between two vertices of a graph,
 * and you may run out of memory when using this function, if your graph is
 * lattice-like.
 *
 * @param loopStart
 * @param loopStep1
 * @param maxLength
 * @return std::vector<Atom>
 */
std::vector<Atom>
Universe::findMinimalOrderLoopFrom(const long int loopStart,
                                   const long int loopStep1,
                                   const int maxLength,
                                   const bool skipSelfLoops) const
{
  long int startingCrosslinkerVertexId = this->getIdxByAtomId(loopStart);
  long int nextStepVertexId = this->getIdxByAtomId(loopStep1);

  std::vector<Atom> minimalPath;

  // First, check the two simplest cases
  // check for self-loops
  if (!skipSelfLoops && loopStart == loopStep1) {
    std::vector<igraph_integer_t> crossLinkersBonds =
      this->getVertexIdxsConnectedTo(startingCrosslinkerVertexId);
    if (std::find(crossLinkersBonds.begin(),
                  crossLinkersBonds.end(),
                  startingCrosslinkerVertexId) != crossLinkersBonds.end()) {
      minimalPath.push_back(
        this->getAtomByVertexIdx(startingCrosslinkerVertexId));
      return minimalPath;
    }
  }

  // check for second order loops
  if (this->getEdgeIdsFromTo(startingCrosslinkerVertexId, nextStepVertexId)
        .size() > 1) {
    minimalPath.push_back(
      this->getAtomByVertexIdx(startingCrosslinkerVertexId));
    minimalPath.push_back(this->getAtomByVertexIdx(nextStepVertexId));
    return minimalPath;
  }

  // NOTE: There are exponentially many paths between two vertices of a graph,
  // and you may run out of memory when using this function, if your graph is
  // lattice-like.
  // note: this algorithm is not particularly efficient
  // it is of the order of O(n*n!)
  std::unordered_set<int> processedPathsKeys;

  // loop neighbours
  int currentMaxLength = 4;
  igraph_vector_int_t paths;
  igraph_vector_int_init(&paths, 0);
  while (igraph_vector_int_size(&paths) <= 4 &&
         (maxLength < 0 || currentMaxLength <= maxLength)) {
    // for this specified neighbour, we search the simple paths
    // (multiple times as we, for memory issue prevention, increase the )
    REQUIRE_IGRAPH_SUCCESS(
      igraph_get_all_simple_paths(&this->graph,
                                  &paths,
                                  startingCrosslinkerVertexId,
                                  igraph_vss_1(nextStepVertexId),
                                  currentMaxLength,
                                  IGRAPH_ALL));
    if (currentMaxLength == maxLength || currentMaxLength < 1) {
      break;
    }
    currentMaxLength *= 2;
    if (maxLength > 0 && currentMaxLength > maxLength) {
      currentMaxLength = maxLength;
    }
    if (maxLength < 0 && currentMaxLength > 64) {
      currentMaxLength = -1;
    }
  }

  // translate the paths we found
  std::vector<Atom> currentPath;
  long int currentPathKey = 0;
  size_t n = igraph_vector_int_size(&paths);
  for (size_t i = 0; i < n; ++i) {
    const long int currentVal = igraph_vector_int_get(&paths, i);
    if (currentVal == -1) {
      // skip self-loops and duplicates
      bool allowedSelfLoop = (!skipSelfLoops || currentPath.size() > 2);
      bool secondOrderLoopValid =
        (currentPath.size() != 2 ||
         this->getEdgeIdsFromTo(startingCrosslinkerVertexId, nextStepVertexId)
             .size() > 1);
      bool alreadyExistingLoop =
        pylimer_tools::utils::map_has_key(processedPathsKeys, currentPathKey);
      if (allowedSelfLoop && !alreadyExistingLoop && secondOrderLoopValid) {
        bool pathIsNewMinimal = (currentPath.size() <= minimalPath.size() &&
                                 currentPath.size() > 0) ||
                                minimalPath.size() <= 1;
        if (pathIsNewMinimal) {
          minimalPath = currentPath;
        }
        processedPathsKeys.insert(currentPathKey);
      }
      currentPath.clear();
      currentPathKey = 0;
    } else {
      Atom newAtom = this->getAtomByVertexIdx(currentVal);
      currentPathKey = currentPathKey xor currentVal; // compute hash
      currentPath.push_back(newAtom);
    }
  }

  igraph_vector_int_destroy(&paths);

  return minimalPath;
};

/**
 * @brief Check whether the universe contains a loop that crosses the periodic
 * boundaries an odd times
 *
 * NOTE: There are exponentially many paths between two vertices of a graph,
 * and you may run out of memory when using this function, if your graph is
 * lattice-like.
 *
 * @param crossLinkerType
 * @param maxLength the maximum length of the loop ()
 * @return true
 * @return false
 */
bool
Universe::hasInfiniteStrand(const int crossLinkerType,
                            const int maxLength) const
{
  // NOTE: There are exponentially many paths between two vertices of a graph,
  // and you may run out of memory when using this function, if your graph is
  // lattice-like.
  std::vector<igraph_integer_t> startingCrosslinkers =
    this->getIndicesWithAttribute("type", crossLinkerType);

  // note: this algorithm is not particularly efficient
  // it is of the order of O(n*n!)
  for (long int startingCrosslinkerVertexId : startingCrosslinkers) {
    // select all neighbouring atoms as possible directions for the loop
    igraph_vector_int_t neighbours;
    igraph_vector_int_init(&neighbours, 0);

    REQUIRE_IGRAPH_SUCCESS(igraph_neighbors(
      &graph, &neighbours, startingCrosslinkerVertexId, IGRAPH_ALL));

    // loop neighbours
    igraph_vector_int_t paths;
    igraph_vector_int_init(&paths, 0);
    // for each neighbour, we search the simple paths
    REQUIRE_IGRAPH_SUCCESS(
      igraph_get_all_simple_paths(&this->graph,
                                  &paths,
                                  startingCrosslinkerVertexId,
                                  igraph_vss_vector(&neighbours),
                                  maxLength,
                                  IGRAPH_ALL));

    igraph_vector_int_destroy(&neighbours);
    // translate the paths we found
    std::vector<Atom> currentPath;
    int nrOfTraversalsX = 0;
    int nrOfTraversalsY = 0;
    int nrOfTraversalsZ = 0;
    size_t n = igraph_vector_int_size(&paths);
    for (int i = 0; i < n; ++i) {
      const long int currentVal = igraph_vector_int_get(&paths, i);
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
          nrOfTraversalsX += ((dx) > 0.5 * (this->box.getLx()))
                               ? 1
                               : (dx < -0.5 * (this->box.getLx()) ? -1 : 0);
          double dy = newAtom.getY() - lastAtom.getY();
          nrOfTraversalsY += ((dx) > 0.5 * (this->box.getLy()))
                               ? 1
                               : (dy < -0.5 * (this->box.getLy()) ? -1 : 0);
          double dz = newAtom.getZ() - lastAtom.getZ();
          nrOfTraversalsZ += ((dz) > 0.5 * (this->box.getLz()))
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
std::map<int, int>
Universe::determineFunctionalityPerType() const
{
  std::map<int, int> result;
  igraph_vector_int_t degrees;
  REQUIRE_IGRAPH_SUCCESS(igraph_vector_int_init(&degrees, 0));
  igraph_vs_t allVertexIds;
  igraph_vs_all(&allVertexIds);
  // complexity: O(|v|*d)
  REQUIRE_IGRAPH_SUCCESS(
    igraph_degree(&this->graph, &degrees, allVertexIds, IGRAPH_ALL, false));

  std::vector<int> types = this->getPropertyValues<int>("type");
  std::set<int> uniqueTypes(types.begin(), types.end());
  // make sure the keys are (re)set, for every type
  for (int type : uniqueTypes) {
    result[type] = 0;
  }
  for (const auto& [type, mass] : this->massPerType) {
    result[type] = 0;
  }

  // complexity: O(|V|)
  igraph_vit_t vit;
  REQUIRE_IGRAPH_SUCCESS(igraph_vit_create(&graph, allVertexIds, &vit));
  while (!IGRAPH_VIT_END(vit)) {
    long int vertexId = static_cast<long int>(IGRAPH_VIT_GET(vit));
    result[types[vertexId]] =
      std::max(static_cast<int>(igraph_vector_int_get(&degrees, vertexId)),
               result[types[vertexId]]);
    IGRAPH_VIT_NEXT(vit);
  }
  igraph_vit_destroy(&vit);
  igraph_vs_destroy(&allVertexIds);
  igraph_vector_int_destroy(&degrees);

  return result;
}

/**
 * @brief determine the effective functionality of each atom type
 *
 * @return std::map<int, double>
 */
std::map<int, double>
Universe::determineEffectiveFunctionalityPerType() const
{
  std::map<int, double> result;
  igraph_vector_int_t degrees;
  REQUIRE_IGRAPH_SUCCESS(igraph_vector_int_init(&degrees, 0));

  igraph_vs_t allVertexIds;
  igraph_vs_all(&allVertexIds);
  // complexity: O(|v|*d)
  REQUIRE_IGRAPH_SUCCESS(
    igraph_degree(&this->graph, &degrees, allVertexIds, IGRAPH_ALL, false));

  std::vector<int> types = this->getPropertyValues<int>("type");
  std::set<int> uniqueTypes(types.begin(), types.end());
  // make sure the keys are (re)set, for every type
  for (int type : uniqueTypes) {
    result[type] = 0;
  }

  // complexity: O(|V|)
  igraph_vit_t vit;
  REQUIRE_IGRAPH_SUCCESS(igraph_vit_create(&graph, allVertexIds, &vit));
  while (!IGRAPH_VIT_END(vit)) {
    long int vertexId = static_cast<long int>(IGRAPH_VIT_GET(vit));
    result[types[vertexId]] += igraph_vector_int_get(&degrees, vertexId);
    IGRAPH_VIT_NEXT(vit);
  }
  igraph_vit_destroy(&vit);
  igraph_vs_destroy(&allVertexIds);
  igraph_vector_int_destroy(&degrees);

  std::map<int, int> typeCounts = this->countAtomTypes();
  for (const auto typePair : typeCounts) {
    result[typePair.first] /= typePair.second;
  }

  return result;
}

/**
 * @brief Compute the weight fractions of each atom type in the network.
 *
 * @return std::map<int, double> $\\vec{W_i}$ (dict): using the type i as a
 * key, this dict contains the weight fractions ($\\frac{W_i}{W_{tot}}$)
 */
std::map<int, double>
Universe::computeWeightFractions() const
{
  std::map<int, double> partialMasses;
  if (this->getNrOfAtoms() == 0) {
    return partialMasses;
  }

  std::map<int, int> numberPerType = this->countAtomTypes();

  // if we do not have any masses stored, we return the "general" fractions
  if (this->massPerType.empty()) {
    // Cast the int_map to double_map
    for (const auto& [key, value] : numberPerType) {
      partialMasses[key] = static_cast<double>(value) / this->getNrOfAtoms();
    }
    return partialMasses;
  }

  double totalMass = 0.0;

  // make sure we have a record for all types
  for (const auto& [type, mass] : this->massPerType) {
    partialMasses[type] = 0.;
  }

  for (const auto& [type, count] : numberPerType) {
    totalMass += this->massPerType.at(type) * count;
    partialMasses[type] += this->massPerType.at(type) * count;
  }

  if (totalMass == 0.0) {
    return partialMasses;
  }

  //  loop to turn partial masses into weight fractions
  for (const auto& partialMassPair : partialMasses) {
    partialMasses[partialMassPair.first] = partialMassPair.second / totalMass;
  }

  return partialMasses;
}

/**
 * @brief Compute the weight fraction of clusters which contain one of the
 * atoms given
 *
 * @param atomIds
 * @return double
 */
double
Universe::computeWeightFractionOfClustersAssociatedWith(
  std::vector<long int> atomIds) const
{
  double totalMass = 0.0;
  double partialMass = 0.0;

  std::vector<pylimer_tools::entities::Universe> clusters = this->getClusters();
  for (const pylimer_tools::entities::Universe& cluster : clusters) {
    double clusterMass = cluster.computeTotalMass();
    totalMass += clusterMass;
    if (std::ranges::any_of(atomIds, [cluster](const long int atomId) {
          return cluster.containsAtomWithId(atomId);
        })) {
      partialMass += clusterMass;
    }
  }

  if (totalMass > 0.0) {
    return partialMass / totalMass;
  } else {
    return totalMass;
  }
}

/**
 * @brief Get an atom id by its vertex index
 *
 * @param vertexId
 * @return long int
 */
long int
Universe::getAtomIdByIdx(const igraph_integer_t vertexId) const
{
  return igraphRealToInt<long int>(VAN(&this->graph, "id", vertexId));
}

/**
 * @brief Get a vertex index by the atom id
 *
 * @param atomId
 * @return long int
 */
igraph_integer_t
Universe::getIdxByAtomId(const long int atomId) const
{
  // if (!pylimer_tools::utils::map_has_key(this->atomIdToVertexIdx, atomId))
  // {
  //   throw std::invalid_argument(
  //     "Universe cannot return idx of atom id: atom with this id (" +
  //     std::to_string(atomId) + ") does not exist");
  // }
  return this->atomIdToVertexIdx.at(atomId);
}

/**
 * @brief Check whether the given atom id is present in this universe
 *
 * @param atomId
 * @return true
 * @return false
 */
bool
Universe::containsAtomWithId(const int atomId) const
{
  return pylimer_tools::utils::map_has_key(this->atomIdToVertexIdx, atomId);
}

/**
 * @brief Check whether the given atom is present in this universe
 *
 * @param atom
 * @return true
 * @return false
 */
bool
Universe::containsAtom(const Atom& atom) const
{
  return this->containsAtomWithId(atom.getId()) &&
         (this->getAtom(atom.getId()) == atom);
}

/**
 * @brief Get an atom by its id
 *
 * @param atomId
 * @return Atom
 */
Atom
Universe::getAtom(const int atomId) const
{
  return this->getAtomByVertexIdx(this->getIdxByAtomId(atomId));
}

/**
 * @brief Get all atoms in this universe
 *
 * @return std::vector<Atom>
 */
std::vector<Atom>
Universe::getAtoms() const
{
  std::vector<Atom> atoms;
  atoms.reserve(this->getNrOfAtoms());
  igraph_vit_t vit;
  REQUIRE_IGRAPH_SUCCESS(igraph_vit_create(&graph, igraph_vss_all(), &vit));
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
std::map<std::string, std::vector<long int>>
Universe::getAngles() const
{
  RUNTIME_EXP_IFN(all_equal<size_t>(4,
                                    this->angleFrom.size(),
                                    this->angleTo.size(),
                                    this->angleVia.size(),
                                    this->angleType.size()),
                  "Angles' state is inconsistent.");
  std::map<std::string, std::vector<long int>> results;
  results.insert_or_assign("angle_from", this->angleFrom);
  results.insert_or_assign("angle_to", this->angleTo);
  results.insert_or_assign("angle_via", this->angleVia);
  std::vector<long int> angleTypes(this->angleType.begin(),
                                   this->angleType.end());
  results.insert_or_assign("angle_type", angleTypes);

  return results;
}

/**
 * @brief Compute the magnitude of the angles stored in this universe
 *
 * @return std::vector<double>
 */
std::vector<double>
Universe::computeAngles() const
{
  RUNTIME_EXP_IFN(all_equal<size_t>(4,
                                    this->angleFrom.size(),
                                    this->angleTo.size(),
                                    this->angleVia.size(),
                                    this->angleType.size()),
                  "Angles' state is inconsistent.");
  std::vector<double> results;
  results.reserve(this->getNrOfAngles());
  for (size_t i = 0; i < this->angleFrom.size(); ++i) {
    Eigen::Vector3d vec1 =
      this->getAtom(this->angleVia[i])
        .vectorTo(this->getAtom(this->angleTo[i]), this->getBox());
    Eigen::Vector3d vec2 =
      this->getAtom(this->angleVia[i])
        .vectorTo(this->getAtom(this->angleFrom[i]), this->getBox());
    results.push_back(std::acos(vec1.dot(vec2) / (vec1.norm() * vec2.norm())));
  }
  return results;
}

/**
 * @brief Get all dihedral angles stored in this universe
 *
 * @return std::map<std::string, std::vector<long int>>
 */
std::map<std::string, std::vector<long int>>
Universe::getDihedralAngles() const
{
  RUNTIME_EXP_IFN(all_equal<size_t>(5,
                                    this->dihedralAngleFrom.size(),
                                    this->dihedralAngleTo.size(),
                                    this->dihedralAngleVia1.size(),
                                    this->dihedralAngleVia2.size(),
                                    this->dihedralAngleType.size()),
                  "Dihedral angles' state is inconsistent.");
  std::map<std::string, std::vector<long int>> results;
  results.insert_or_assign("dihedral_angle_from", this->dihedralAngleFrom);
  results.insert_or_assign("dihedral_angle_via1", this->dihedralAngleVia1);
  results.insert_or_assign("dihedral_angle_via2", this->dihedralAngleVia2);
  results.insert_or_assign("dihedral_angle_to", this->dihedralAngleTo);
  std::vector<long int> angleTypes(this->dihedralAngleType.begin(),
                                   this->dihedralAngleType.end());
  results.insert_or_assign("dihedral_angle_type", angleTypes);

  return results;
}

/**
 * @brief Find all angles that appear in this universe
 *
 * @return std::map<std::string, std::vector<long int>>
 */
std::map<std::string, std::vector<long int>>
Universe::detectAngles() const
{
  std::vector<long int> angleFromFound;
  std::vector<long int> angleToFound;
  std::vector<long int> angleViaFound;
  std::vector<long int> angleTypeFound;
  // query all atoms
  igraph_vit_t vit;
  REQUIRE_IGRAPH_SUCCESS(
    igraph_vit_create(&this->graph, igraph_vss_all(), &vit));
  while (!IGRAPH_VIT_END(vit)) {
    const igraph_integer_t vertexIdx = IGRAPH_VIT_GET(vit);
    // find the connected atoms
    std::vector<igraph_integer_t> connections =
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
           connectionJ < connections.size();
           ++connectionJ) {
        // TODO: check again that we do not get duplicates this way
        int atomTypeFrom =
          this->getPropertyValue<int>("type", connections[connectionI]);
        int atomTypeVia = this->getPropertyValue<int>("type", vertexIdx);
        int atomTypeTo =
          this->getPropertyValue<int>("type", connections[connectionJ]);
        angleFromFound.push_back(
          this->getAtomIdByIdx(connections[connectionI]));
        angleViaFound.push_back(this->getAtomIdByIdx(vertexIdx));
        angleToFound.push_back(this->getAtomIdByIdx(connections[connectionJ]));
        angleTypeFound.push_back(
          this->hashAngleType(atomTypeFrom, atomTypeVia, atomTypeTo));
      }
    }
    IGRAPH_VIT_NEXT(vit);
  }
  igraph_vit_destroy(&vit);

  std::map<std::string, std::vector<long int>> results;
  results.insert_or_assign("angle_from", angleFromFound);
  results.insert_or_assign("angle_to", angleToFound);
  results.insert_or_assign("angle_via", angleViaFound);
  results.insert_or_assign("angle_type", angleTypeFound);

  return results;
}

/**
 * @brief Find all dihedral angles that appear in this universe
 *
 * @return std::map<std::string, std::vector<long int>>
 */
std::map<std::string, std::vector<long int>>
Universe::detectDihedralAngles() const
{
  std::vector<long int> angleFromFound;
  std::vector<long int> angleVia1Found;
  std::vector<long int> angleVia2Found;
  std::vector<long int> angleToFound;
  std::vector<long int> angleTypeFound;
  // query all atoms
  igraph_vit_t vit;
  REQUIRE_IGRAPH_SUCCESS(
    igraph_vit_create(&this->graph, igraph_vss_all(), &vit));
  while (!IGRAPH_VIT_END(vit)) {
    const long int vertexIdx = static_cast<long int>(IGRAPH_VIT_GET(vit));
    // find the connected atoms
    igraph_vector_int_t dihedral_paths_v;
    igraph_vector_int_init(&dihedral_paths_v, 4);
    REQUIRE_IGRAPH_SUCCESS(igraph_get_all_simple_paths(&this->graph,
                                                       &dihedral_paths_v,
                                                       vertexIdx,
                                                       igraph_vss_all(),
                                                       5,
                                                       IGRAPH_ALL));

    std::vector<std::vector<int>> dihedral_sets;
    // std::cout << "Found paths: " <<
    // igraph_vector_int_size(&dihedral_paths_v)
    //           << std::endl;
    dihedral_sets.reserve(igraph_vector_int_size(&dihedral_paths_v) / 4);
    std::vector<int> current_dihedral_path;
    current_dihedral_path.reserve(4);
    for (size_t i = 0; i < igraph_vector_int_size(&dihedral_paths_v); i++) {
      const int current_vertex = igraph_vector_int_get(&dihedral_paths_v, i);
      if (current_vertex == -1) {
        if (current_dihedral_path.size() == 4) {
          // std::cout << "Found dihedral path of length 4" << std::endl;
          dihedral_sets.push_back(current_dihedral_path);
        }
        current_dihedral_path.clear();
      } else {
        current_dihedral_path.push_back(current_vertex);
      }
    }
    if (current_dihedral_path.size() == 4) {
      dihedral_sets.push_back(current_dihedral_path);
    }

    igraph_vector_int_destroy(&dihedral_paths_v);

    // std::cout << "Found dihedral paths: " << dihedral_sets.size() <<
    // std::endl;

    for (std::vector<int> dihedral_path : dihedral_sets) {
      RUNTIME_EXP_IFN(dihedral_path.size() == 4,
                      "Expected 4 neighbors, got " +
                        std::to_string(dihedral_path.size()) +
                        " when detecting dihedrals.");
      long int vertexIdxFrom = dihedral_path[0];
      long int vertexIdxVia1 = dihedral_path[1];
      long int vertexIdxVia2 = dihedral_path[2];
      long int vertexIdxTo = dihedral_path[3];
      RUNTIME_EXP_IFN(
        vertexIdxFrom == vertexIdx || vertexIdxTo == vertexIdx,
        "Expected the original vertex to be found in the neighbourhood.");
      if (vertexIdxFrom < vertexIdxTo) {
        // we list every dihedral twice otherwise
        continue;
        // std::swap(vertexIdxFrom, vertexIdxTo);
        // std::swap(vertexIdxVia1, vertexIdxVia2);
      }

      int atomTypeFrom = this->getPropertyValue<int>("type", vertexIdxFrom);
      int atomTypeVia1 = this->getPropertyValue<int>("type", vertexIdxVia1);
      int atomTypeVia2 = this->getPropertyValue<int>("type", vertexIdxVia2);
      int atomTypeTo = this->getPropertyValue<int>("type", vertexIdxTo);

      angleFromFound.push_back(this->getAtomIdByIdx(vertexIdxFrom));
      angleVia1Found.push_back(this->getAtomIdByIdx(vertexIdxVia1));
      angleVia2Found.push_back(this->getAtomIdByIdx(vertexIdxVia2));
      angleToFound.push_back(this->getAtomIdByIdx(vertexIdxTo));
      angleTypeFound.push_back(this->hashDihedralAngleType(
        atomTypeFrom, atomTypeVia1, atomTypeVia2, atomTypeTo));
    }
    IGRAPH_VIT_NEXT(vit);
  }
  igraph_vit_destroy(&vit);

  std::map<std::string, std::vector<long int>> results;
  results.insert_or_assign("dihedral_angle_from", angleFromFound);
  results.insert_or_assign("dihedral_angle_to", angleToFound);
  results.insert_or_assign("dihedral_angle_via1", angleVia1Found);
  results.insert_or_assign("dihedral_angle_via2", angleVia2Found);
  results.insert_or_assign("dihedral_angle_type", angleTypeFound);

  return results;
}

/**
 * @brief Reduce the network to crosslinkers only.
 * Includes self-loops (from primary ones).
 * You can remove them using #simplify
 *
 * @param crossLinkerType the atom type of the crosslinker beads
 * @return Universe
 */
Universe
Universe::getNetworkOfCrosslinker(const int crossLinkerType) const
{
  // TODO: prevent duplicate of self-loops
  // How this works:
  // 1. find all crossLinkers
  // 2. from each crossLinker, walk in all directions
  // 3. if the walk reaches another crossLinker, we found a
  //    crossLinker-crossLinker connection.
  //    To reduce duplicates, we only take the ones where we started from a
  //    crossLinker with a smaller (or equal, for self-/primary-loops)
  //    vertex index
  Universe newUniverse =
    Universe(this->box.getLx(), this->box.getLy(), this->box.getLz());
  std::vector<long int> bondFrom;
  std::vector<long int> bondTo;
  std::vector<igraph_integer_t> crossLinkers =
    this->getIndicesWithAttribute("type", crossLinkerType);
  for (long int crossLinker : crossLinkers) {
    std::vector<igraph_integer_t> connections =
      this->getVertexIdxsConnectedTo(crossLinker);
    for (long int connection : connections) {
      long int currentCenter = connection;
      long int lastCenter = crossLinker;
      std::vector<igraph_integer_t> subConnections =
        this->getVertexIdxsConnectedTo(currentCenter);
      while (subConnections.size() > 0) {
        if (this->getPropertyValue<int>("type", currentCenter) ==
            crossLinkerType) {
          // found crosslinker
          if (currentCenter >= crossLinker) {
            bondFrom.push_back(
              this->getPropertyValue<int>("id", currentCenter));
            bondTo.push_back(this->getPropertyValue<int>("id", crossLinker));
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

  assert(bondTo.size() == bondFrom.size());
  // with the algorithm above, self-loops are counted twice.
  // let's just remove the second (and/or fourth) one where needed
  // NOTE: some assumptions are made here that could be problematic;
  // for example, that there are not more than 1 self-loops per crosslink
  // in the beginning
  std::vector<size_t> indicesToRemove;
  std::map<int, int> nrOfSelfLoops;
  for (size_t i = 0; i < bondTo.size(); ++i) {
    if (bondTo[i] == bondFrom[i]) {
      if (!pylimer_tools::utils::map_has_key(nrOfSelfLoops, bondTo[i])) {
        nrOfSelfLoops.emplace(bondTo[i], 0);
      }
      nrOfSelfLoops[bondTo[i]] += 1;
      if (nrOfSelfLoops.at(bondTo[i]) % 2 == 0) {
        indicesToRemove.push_back(i);
      }
    }
  }
  // actually remove them. Reverse in order to not mess with the indices
  if (indicesToRemove.size() > 0) {
    // CAUTION: do not use size_t here!
    for (int i = indicesToRemove.size() - 1; i >= 0; --i) {
      bondTo.erase(bondTo.begin() + indicesToRemove[i]);
      bondFrom.erase(bondFrom.begin() + indicesToRemove[i]);
    }
  }

  std::vector<int> zeros =
    pylimer_tools::utils::initializeWithValue(crossLinkers.size(), 0);

  newUniverse.addAtoms(this->getPropertyValues<long int>("id", crossLinkers),
                       this->getPropertyValues<int>("type", crossLinkers),
                       this->getPropertyValues<double>("x", crossLinkers),
                       this->getPropertyValues<double>("y", crossLinkers),
                       this->getPropertyValues<double>("z", crossLinkers),
                       zeros,
                       zeros,
                       zeros);
  newUniverse.addBonds(
    bondFrom.size(),
    bondFrom,
    bondTo,
    pylimer_tools::utils::initializeWithValue(bondFrom.size(), 0),
    false,
    false);

  return newUniverse;
};

/**
 * @brief Combine vertices along a bond, for a certain bond type.
 *
 * This function replaces two vertices that are connected by an edge of the type
 * `bondType` with a single vertex.
 *
 * @param bondType
 * @return Universe
 */
Universe
Universe::contractVerticesAlongBondType(const int bondType) const
{
  Universe result = Universe(*this);

  bool lastRoundDidRemove = true;

  std::vector<std::string> edgeProperties =
    this->getVertexAndEdgePropertyNames().second;

  igraph_vector_int_t edgesToCopy;
  igraph_vector_int_init(&edgesToCopy, 0);
  while (lastRoundDidRemove) {
    lastRoundDidRemove = false;

    // Collect all operations to batch them
    std::vector<igraph_integer_t> edgesToAdd;
    std::vector<std::pair<igraph_integer_t, igraph_integer_t>>
      edgePropertyCopyPairs; // (old_edge, new_edge_offset)
    std::vector<igraph_integer_t> edgesToDelete;
    std::vector<igraph_integer_t> verticesToDelete;
    std::vector<bool> vertexAffectedThisIteration =
      pylimer_tools::utils::initializeWithValue<bool>(
        igraph_vcount(&result.graph), false);

    for (igraph_integer_t edgeIdx = 0; edgeIdx < igraph_ecount(&result.graph);
         ++edgeIdx) {
      if (result.getEdgePropertyValue<int>("type", edgeIdx) == bondType) {
        lastRoundDidRemove = true;

        igraph_integer_t vertex1OfEdge;
        igraph_integer_t vertex2OfEdge;
        REQUIRE_IGRAPH_SUCCESS(
          igraph_edge(&result.graph, edgeIdx, &vertex1OfEdge, &vertex2OfEdge));

        if (vertexAffectedThisIteration[vertex1OfEdge] ||
            vertexAffectedThisIteration[vertex2OfEdge]) {
          // This vertex was already processed in this round, skip it
          continue;
        }

        vertexAffectedThisIteration[vertex1OfEdge] = true;
        vertexAffectedThisIteration[vertex2OfEdge] = true;

        if (vertex1OfEdge == vertex2OfEdge) {
          // Self-loop case
          edgesToDelete.push_back(edgeIdx);
        } else {
          // Get incident edges of vertex1 to redirect to vertex2
          REQUIRE_IGRAPH_SUCCESS(igraph_incident(
            &result.graph, &edgesToCopy, vertex1OfEdge, IGRAPH_ALL));

          for (size_t i = 0; i < igraph_vector_int_size(&edgesToCopy); ++i) {
            igraph_integer_t edgeToCopy = VECTOR(edgesToCopy)[i];
            if (edgeToCopy != edgeIdx) {
              // Add new edge from vertex2 to the other end
              igraph_integer_t otherVertex =
                IGRAPH_OTHER(&result.graph, edgeToCopy, vertex1OfEdge);
              edgesToAdd.push_back(vertex2OfEdge);
              edgesToAdd.push_back(otherVertex);

              vertexAffectedThisIteration[otherVertex] = true;

              // Record which edge properties to copy later
              edgePropertyCopyPairs.emplace_back(edgeToCopy,
                                                 (edgesToAdd.size() / 2) - 1);
            }
          }

          edgesToDelete.push_back(edgeIdx);
          verticesToDelete.push_back(vertex1OfEdge);
        }
      }
    }

    if (!lastRoundDidRemove) {
      break;
    }

    // Batch add all new edges at once
    if (!edgesToAdd.empty()) {
      igraph_vector_int_t newEdgesVec;
      igraph_vector_int_init(&newEdgesVec, edgesToAdd.size());
      pylimer_tools::utils::StdVectorToIgraphVectorT(edgesToAdd, &newEdgesVec);

      igraph_integer_t edgeCountBefore = igraph_ecount(&result.graph);
      REQUIRE_IGRAPH_SUCCESS(
        igraph_add_edges(&result.graph, &newEdgesVec, nullptr));

      // Copy edge properties for all new edges
      for (const auto& copyPair : edgePropertyCopyPairs) {
        igraph_integer_t newEdgeIdx = edgeCountBefore + copyPair.second;
        pylimer_tools::utils::copyEdgeProperties(&result.graph,
                                                 copyPair.first,
                                                 &result.graph,
                                                 newEdgeIdx,
                                                 edgeProperties);
      }

      igraph_vector_int_destroy(&newEdgesVec);
    }

    // Batch delete vertices (this also removes incident edges)
    if (!verticesToDelete.empty()) {
      // Sort and remove duplicates
      std::sort(verticesToDelete.begin(), verticesToDelete.end());
      verticesToDelete.erase(
        std::unique(verticesToDelete.begin(), verticesToDelete.end()),
        verticesToDelete.end());

      igraph_vector_int_t verticesToDeleteVec;
      igraph_vector_int_init(&verticesToDeleteVec, verticesToDelete.size());
      for (size_t i = 0; i < verticesToDelete.size(); ++i) {
        VECTOR(verticesToDeleteVec)[i] = verticesToDelete[i];
      }

      REQUIRE_IGRAPH_SUCCESS(igraph_delete_vertices(
        &result.graph, igraph_vss_vector(&verticesToDeleteVec)));

      igraph_vector_int_destroy(&verticesToDeleteVec);
    }
  }

  igraph_vector_int_destroy(&edgesToCopy);
  result.NAtoms = igraph_vcount(&result.graph);
  result.NBonds = igraph_ecount(&result.graph);
  result.resetAtomIdMapping();
  result.removeAllAngles();
  result.removeAllDihedralAngles();

  return result;
}

/**
 * @brief Get the number of angles stored in this universe.
 *
 * @return const int
 */
size_t
Universe::getNrOfAngles() const
{
  RUNTIME_EXP_IFN(all_equal<size_t>(4,
                                    this->angleFrom.size(),
                                    this->angleTo.size(),
                                    this->angleVia.size(),
                                    this->angleType.size()),
                  "Invalid internal state: The angle info is not consistent");
  return this->angleFrom.size();
}

/**
 * @brief Get the number of angles stored in this universe.
 *
 * @return const int
 */
size_t
Universe::getNrOfDihedralAngles() const
{
  RUNTIME_EXP_IFN(
    all_equal<size_t>(5,
                      this->dihedralAngleFrom.size(),
                      this->dihedralAngleTo.size(),
                      this->dihedralAngleVia1.size(),
                      this->dihedralAngleVia2.size(),
                      this->dihedralAngleType.size()),
    "Invalid internal state: The dihedral info is not consistent");
  return this->dihedralAngleFrom.size();
}

/**
 * @brief Find the vertex id where the vertex has a certain property
 *
 * @tparam IN the type of the property you want the vertex to have
 * @param propertyName the name of the property
 * @param propertyValue the objective value
 * @return long int the vertex index, -1 if not found.
 */
template<typename IN>
long int
Universe::findVertexIdForProperty(const char* propertyName,
                                  IN propertyValue) const
{
  igraph_vector_t allValues;
  REQUIRE_IGRAPH_SUCCESS(igraph_vector_init(&allValues, this->getNrOfAtoms()));
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
 * @brief An algorithm to determine whether two loops are entangled
 *
 * @param vertexIndicesLoop1
 * @param vertexIndicesLoop2
 * @return true
 * @return false
 */
std::vector<LoopIntersectionInfo>
Universe::findLoopEntanglements(
  const std::vector<igraph_integer_t>& vertexIndicesLoop1,
  const std::vector<igraph_integer_t>& vertexIndicesLoop2,
  const std::vector<igraph_integer_t>& edgeIndicesLoop1,
  const std::vector<igraph_integer_t>& edgeIndicesLoop2) const
{
  std::vector<LoopIntersectionInfo> results;
  Eigen::Vector3d helperNode = Eigen::Vector3d::Zero();
  double sizeDenominator = 1.0 / static_cast<double>(vertexIndicesLoop1.size());
  INVALIDARG_EXP_IFN(
    (edgeIndicesLoop1.size() == vertexIndicesLoop1.size()) ||
      (edgeIndicesLoop1.size() == 0),
    "Every loop must consist of equal number of edges and vertices");
  INVALIDARG_EXP_IFN(
    (edgeIndicesLoop2.size() == vertexIndicesLoop2.size()) ||
      (edgeIndicesLoop2.size() == 0),
    "Every loop must consist of equal number of edges and vertices");
  for (long int i = 0; i < vertexIndicesLoop1.size(); ++i) {
    // find PBC corrected mean position
    Eigen::Vector3d vertexIdxPos = this->getPositionVectorForVertex(i);
    Eigen::Vector3d distance = vertexIdxPos - helperNode;
    this->box.handlePBC(distance);
    helperNode += ((distance).array() * sizeDenominator).matrix();
  }
  int intersections = 0;
  // now that we have the helper node, we can check all edges of loop2, how
  // often they intersect the triangles
  for (long int i = 0; i < vertexIndicesLoop1.size(); ++i) {
    if (i == 0) {
      // random estimate
      results.reserve((vertexIndicesLoop1.size() + vertexIndicesLoop2.size()) /
                      3);
    }
    if (i == 1) {
      // better estimate
      results.reserve(results.size() * (results.size() - 1));
    }
    Eigen::Vector3d vertex0 =
      this->getPositionVectorForVertex(vertexIndicesLoop1[i]);
    long int vertex1Index =
      (i == 0) ? (vertexIndicesLoop1.size() - 1) : (i - 1);
    Eigen::Vector3d vertex1 =
      this->getPositionVectorForVertex(vertexIndicesLoop1[vertex1Index]);
    // the triangle is now spawned by vertex0, vertex1 and the helperNode
    for (size_t j = 0; j < vertexIndicesLoop2.size(); ++j) {
      Eigen::Vector3d rayOrigin =
        this->getPositionVectorForVertex(vertexIndicesLoop2[j]);
      long int directionIdx = (j == 0) ? vertexIndicesLoop2.size() - 1 : j - 1;
      Eigen::Vector3d rayTarget =
        this->getPositionVectorForVertex(vertexIndicesLoop2[directionIdx]);
      Eigen::Vector3d intersectionPoint;
      if (pylimer_tools::topo::segmentIntersectsTriangle(
            rayOrigin,
            rayTarget,
            vertex0,
            vertex1,
            helperNode,
            intersectionPoint,
            [&](Eigen::Vector3d vec) {
              this->box.handlePBC<Eigen::Vector3d>(vec);
              return vec;
            })) {
        std::cout << "Intersection found at indices " << i << ", " << j
                  << " from to " << rayOrigin[0] << ", " << rayOrigin[1] << ", "
                  << rayOrigin[2] << "; " << rayTarget[0] << ", "
                  << rayTarget[1] << ", " << rayTarget[2] << " in triangle "
                  << vertex0[0] << ", " << vertex0[1] << ", " << vertex0[2]
                  << "; " << vertex1[0] << ", " << vertex1[1] << ", "
                  << vertex1[2] << "; " << helperNode[0] << ", "
                  << helperNode[1] << ", " << helperNode[2] << "." << std::endl;
        intersections += 1;
        // TODO: maybe do something with the intersection point?
        LoopIntersectionInfo result;
        std::vector<Atom> involvedAtoms;
        involvedAtoms.reserve(4);
        involvedAtoms.push_back(
          this->getAtomByVertexIdx(vertexIndicesLoop1[i]));
        involvedAtoms.push_back(
          this->getAtomByVertexIdx(vertexIndicesLoop1[vertex1Index]));
        involvedAtoms.push_back(
          this->getAtomByVertexIdx(vertexIndicesLoop2[j]));
        involvedAtoms.push_back(
          this->getAtomByVertexIdx(vertexIndicesLoop2[directionIdx]));
        result.involvedAtoms = involvedAtoms;
        if (edgeIndicesLoop1.size() > 0) {
          result.edge1 = edgeIndicesLoop1[vertex1Index];
        } else {
          result.edge1 = -1;
        }
        if (edgeIndicesLoop2.size() > 0) {
          result.edge2 = edgeIndicesLoop2[directionIdx];
        } else {
          result.edge2 = -1;
        }
        result.intersectionPoint = intersectionPoint;
        // TODO: check that these are always in the same direction
        Eigen::Vector3d triangleNormal =
          (helperNode - vertex0).cross(helperNode - vertex1);
        result.direction = triangleNormal.dot(rayTarget - rayOrigin);
        results.push_back(result);
      }
    }
  }
  // we have an entanglement, iff the triangle spawned is crossed an odd
  // number of times
  // return intersections > 0 && intersections % 2 == 0;
  return results;
};

std::vector<std::pair<size_t, size_t>>
Universe::interpolateEdges(const int crossLinkerType,
                           const double interpolationFactor) const
{
  INVALIDARG_EXP_IFN(interpolationFactor >= 0.0,
                     "Negative interpolation factor does not make sense");
  std::vector<long int> vertexIdToNewIdx =
    pylimer_tools::utils::initializeWithValue<long int>(this->getNrOfAtoms(),
                                                        -1);
  size_t currentMaxIdx = 0;
  std::vector<std::pair<size_t, size_t>> results;
  std::vector<Molecule> chainsToInterpolate =
    this->getChainsWithCrosslinker(crossLinkerType);
  for (const auto& chain : chainsToInterpolate) {
    std::vector<Atom> ends = chain.getChainEnds(crossLinkerType, true);
    RUNTIME_EXP_IFN(ends.size() == 2, "Chain must have exactly two ends");
    size_t end0 = this->getIdxByAtomId(ends[0].getId());
    size_t end1 = this->getIdxByAtomId(ends[1].getId());
    if (vertexIdToNewIdx[end0] == -1) {
      vertexIdToNewIdx[end0] = currentMaxIdx;
      currentMaxIdx += 1;
    }
    size_t previousId = vertexIdToNewIdx[end0];
    size_t numExtraBonds = static_cast<size_t>(std::round(
      static_cast<double>(chain.getNrOfBonds()) * interpolationFactor));
    for (size_t i = 1; i < numExtraBonds; ++i) {
      results.emplace_back(previousId, currentMaxIdx);
      previousId = currentMaxIdx;
      currentMaxIdx += 1;
    }
    if (vertexIdToNewIdx[end1] == -1) {
      vertexIdToNewIdx[end1] = currentMaxIdx;
      currentMaxIdx += 1;
    }
    results.emplace_back(vertexIdToNewIdx[end1], previousId);
  }

  return results;
}

/**
 * @brief Compute the x distance for all bonds passed in
 *
 * @param bondFrom
 * @param bondTo
 * @return std::vector<double>
 */
std::vector<double>
Universe::computeDxs(const std::vector<long int>& bondFrom,
                     const std::vector<long int>& bondTo) const
{
  return this->computeDs(bondFrom, bondTo, "x", this->box.getLx());
};

/**
 * @brief Compute the y distance for all bonds passed in
 *
 * @param bondFrom
 * @param bondTo
 * @return std::vector<double>
 */
std::vector<double>
Universe::computeDys(const std::vector<long int>& bondFrom,
                     const std::vector<long int>& bondTo) const
{
  return this->computeDs(bondFrom, bondTo, "y", this->box.getLy());
};

/**
 * @brief Compute the z distance for all bonds passed in
 *
 * @param bondFrom
 * @param bondTo
 * @return std::vector<double>
 */
std::vector<double>
Universe::computeDzs(const std::vector<long int>& bondFrom,
                     const std::vector<long int>& bondTo) const
{
  return this->computeDs(bondFrom, bondTo, "z", this->box.getLz());
};

/**
 * @brief Compute the distance for all bonds passed in in a certain
 * direction
 *
 * @param bondFrom
 * @param bondTo
 * @param direction
 * @param boxLimit
 * @return std::vector<double>
 */
std::vector<double>
Universe::computeDs(const std::vector<long int>& bondFrom,
                    const std::vector<long int>& bondTo,
                    const std::string& direction,
                    const double boxLimit) const
{
  INVALIDARG_EXP_IFN(bondFrom.size() == bondTo.size(),
                     "bond from and bond to must have the same size");

  size_t nBonds = bondFrom.size();

  igraph_vector_int_t vertexIdFrom;
  REQUIRE_IGRAPH_SUCCESS(igraph_vector_int_init(&vertexIdFrom, nBonds));
  igraph_vector_int_t vertexIdTo;
  REQUIRE_IGRAPH_SUCCESS(igraph_vector_int_init(&vertexIdTo, nBonds));

  for (size_t i = 0; i < nBonds; ++i) {
    igraph_vector_int_set(
      &vertexIdFrom, i, this->atomIdToVertexIdx.at(bondFrom[i]));
    igraph_vector_int_set(
      &vertexIdTo, i, this->atomIdToVertexIdx.at(bondTo[i]));
  }

  igraph_vector_t dValuesFrom;
  REQUIRE_IGRAPH_SUCCESS(igraph_vector_init(&dValuesFrom, nBonds));
  igraph_vector_t dValuesTo;
  REQUIRE_IGRAPH_SUCCESS(igraph_vector_init(&dValuesTo, nBonds));

  std::string property = direction;
  REQUIRE_IGRAPH_SUCCESS(
    igraph_cattribute_VANV(&this->graph,
                           property.c_str(),
                           igraph_vss_vector(&vertexIdFrom),
                           &dValuesFrom));
  REQUIRE_IGRAPH_SUCCESS(igraph_cattribute_VANV(&this->graph,
                                                property.c_str(),
                                                igraph_vss_vector(&vertexIdTo),
                                                &dValuesTo));

  igraph_vector_int_destroy(&vertexIdFrom);
  igraph_vector_int_destroy(&vertexIdTo);

  std::vector<double> results;
  results.reserve(nBonds);

  for (int i = 0; i < nBonds; ++i) {
    double currentD = (double)igraph_vector_get(&dValuesTo, i) -
                      (double)igraph_vector_get(&dValuesFrom, i);
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

double
Universe::computeTemperature(const int dimensions, const double kb) const
{
  INVALIDARG_EXP_IFN(this->vertexPropertyExists("vx") &&
                       this->vertexPropertyExists("vy") &&
                       this->vertexPropertyExists("vz"),
                     "Velocities are not present in this universe. Cannot "
                     "compute temperature.");
  std::vector<Atom> atoms = this->getAtoms();
  double kineticEnergy = 0.;
  for (const Atom& atom : atoms) {
    Eigen::Vector3d velocity;
    velocity << atom.getProperty("vx"), atom.getProperty("vy"),
      atom.getProperty("vz");
    kineticEnergy +=
      0.5 * this->massPerType.at(atom.getType()) * velocity.squaredNorm();
  }

  return kineticEnergy /
         ((static_cast<double>(dimensions) / 2.) * (atoms.size()) * kb);
}

/**
 * @brief Get the mean number of beads between beads with the passed type
 *
 * @param crossLinkerType
 * @return double
 */
double
Universe::getMeanStrandLength(const int crossLinkerType) const
{
  std::vector<Molecule> molecules = this->getMolecules(crossLinkerType);

  double multiplier = 1.0 / static_cast<double>(molecules.size());
  double meanStrandLength = std::accumulate(
    molecules.begin(),
    molecules.end(),
    0.0,
    [multiplier](const double val, const Molecule molecule) {
      return val + static_cast<double>(molecule.getLength()) * multiplier;
    });

  return meanStrandLength;
}

/**
 * @brief Get the mean end to end distance
 *
 * Does not take loops into account as a contributor to the mean.
 * Returns 0 for systems without any qualifying strands.
 *
 * @param crossLinkerType
 * @return double
 */
std::vector<double>
Universe::computeEndToEndDistances(const int crossLinkerType,
                                   bool implyImageFlags) const
{
  std::vector<Molecule> molecules =
    this->getChainsWithCrosslinker(crossLinkerType);

  std::vector<double> distances;
  distances.reserve(molecules.size());
  std::transform(
    molecules.begin(),
    molecules.end(),
    std::back_inserter(distances),
    [implyImageFlags](Molecule molecule) {
      return implyImageFlags
               ? molecule.computeEndToEndDistanceWithDerivedImageFlags()
               : molecule.computeEndToEndDistance();
    });

  return distances;
}

/**
 * @brief Get the mean end to end distance
 *
 * Does not take loops into account as a contributor to the mean.
 * Returns 0 for systems without any qualifying strands.
 *
 * @param crossLinkerType
 * @return double
 */
double
Universe::computeMeanEndToEndDistance(const int crossLinkerType,
                                      const bool implyImageFlags) const
{
  std::vector<Molecule> molecules =
    this->getChainsWithCrosslinker(crossLinkerType);

  double meanEndToEndDistance = 0.0;
  int validMolecules = 0;

  for (Molecule molecule : molecules) {
    double dist = implyImageFlags
                    ? molecule.computeEndToEndDistanceWithDerivedImageFlags()
                    : molecule.computeEndToEndDistance();
    if (dist > 0.0) {
      meanEndToEndDistance += dist;
      validMolecules += 1;
    }
  }

  if (validMolecules == 0) {
    return 0.0;
  }

  return meanEndToEndDistance / static_cast<double>(validMolecules);
}

/**
 * @brief Get the mean end to end distance
 *
 * Does not take loops into account as a contributor to the mean.
 * Returns 0 for systems without any qualifying strands.
 *
 * @param crossLinkerType
 * @return double
 */
double
Universe::computeMeanSquareEndToEndDistance(
  const int crossLinkerType,
  const bool onlyThoseWithTwoCrosslinkers,
  const bool implyImageFlags) const
{
  std::vector<Molecule> molecules =
    this->getChainsWithCrosslinker(crossLinkerType);

  double meanEndToEndDistance = 0.0;
  int validMolecules = 0;

  for (Molecule molecule : molecules) {
    if (!onlyThoseWithTwoCrosslinkers ||
        molecule.getAtomsOfType(crossLinkerType).size() == 2) {
      double dist = implyImageFlags
                      ? molecule.computeEndToEndDistanceWithDerivedImageFlags()
                      : molecule.computeEndToEndDistance();
      if (dist > 0.0) {
        meanEndToEndDistance += dist * dist;
        validMolecules += 1;
      }
    }
  }

  if (validMolecules == 0) {
    return 0.0;
  }

  return meanEndToEndDistance / static_cast<double>(validMolecules);
}

/**
 * @brief Compute the total mass of this universe
 *
 * @return double
 */
double
Universe::computeTotalMass() const
{
  return this->computeTotalMassWithMasses(this->massPerType);
}

double
Universe::computeTotalMassWithMasses(
  std::map<int, double> massPerTypeToUse) const
{
  std::vector<int> types = this->getAtomTypes();
  double weight = 0.0;
  for (int i = 0; i < types.size(); i++) {
    weight += massPerTypeToUse.at(types[i]);
  }
  return weight;
}

/**
 * @brief Compute the mean of all bond lengths
 *
 * @return double
 */
double
Universe::computeMeanBondLength() const
{
  double length = 0.0;
  if (this->getNrOfBonds() == 0) {
    return length;
  }
  // construct iterator
  igraph_eit_t bondIterator;
  REQUIRE_IGRAPH_SUCCESS(igraph_eit_create(
    &this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &bondIterator));

  while (!IGRAPH_EIT_END(bondIterator)) {
    long int edgeId = static_cast<long int>(IGRAPH_EIT_GET(bondIterator));
    igraph_integer_t bondFrom;
    igraph_integer_t bondTo;
    REQUIRE_IGRAPH_SUCCESS(
      igraph_edge(&this->graph, edgeId, &bondFrom, &bondTo));
    // TODO: this is more intensive than needed
    // check whether the compiler optimizes this or not
    Atom atom1 = this->getAtomByVertexIdx(bondFrom);
    Atom atom2 = this->getAtomByVertexIdx(bondTo);
    length += atom1.distanceTo(atom2, this->box);
    IGRAPH_EIT_NEXT(bondIterator);
  }

  igraph_eit_destroy(&bondIterator);
  return length / static_cast<double>(this->getNrOfBonds());
}
/**
 * @brief Compute the weight average molecular weight
 *
 * @source https://www.pslc.ws/macrog/average.htm
 *
 * @param crossLinkerType
 * @return double
 */
double
Universe::computeWeightAverageMolecularWeight(const int crossLinkerType) const
{
  double weightAverage = 0.0;

  std::vector<Molecule> molecules = this->getMolecules(crossLinkerType);
  std::map<int, double> massPerTypeToUse(this->massPerType);
  massPerTypeToUse[crossLinkerType] = 0.0;
  double totalMass = this->computeTotalMassWithMasses(massPerTypeToUse);
  double massDivisor = 1.0 / totalMass;
  for (Molecule molecule : molecules) {
    double moleculeMass = molecule.computeTotalMass();
    weightAverage += moleculeMass * moleculeMass * massDivisor;
  }

  return weightAverage;
};

/**
 * @brief Compute the number average molecular weight
 *
 * @source https://www.pslc.ws/macrog/average.htm
 *
 * @param crossLinkerType
 * @return double
 */
double
Universe::computeNumberAverageMolecularWeight(const int crossLinkerType) const
{
  std::vector<Molecule> molecules = this->getMolecules(crossLinkerType);
  std::map<int, double> massPerTypeToUse(this->massPerType);
  massPerTypeToUse[crossLinkerType] = 0.0;
  double totalMass = this->computeTotalMassWithMasses(massPerTypeToUse);

  return totalMass / static_cast<double>(molecules.size());
};

/**
 * @brief Compute the polydispersity index (weight average molecular weight
 * over number average molecular weight)
 *
 * @source https://www.pslc.ws/macrog/average.htm
 *
 * @param crossLinkerType
 * @return double
 */
double
Universe::computePolydispersityIndex(const int crossLinkerType) const
{
  // TODO: check assembly whether the double getMolecules() and
  // computeTotalMass() is cancelled out or not
  return this->computeWeightAverageMolecularWeight(crossLinkerType) /
         this->computeNumberAverageMolecularWeight(crossLinkerType);
}

void
Universe::resetAtomIdMapping()
{
  this->atomIdToVertexIdx.clear();
  std::vector<int> atomIds = this->getPropertyValues<int>("id");
  assert(atomIds.size() == igraph_vcount(&this->graph));
  for (igraph_integer_t i = 0; i < atomIds.size(); i++) {
    this->atomIdToVertexIdx.emplace(atomIds[i], i);
  }
}

/**
 * @brief check whether the internal counts are equal to those of the graph
 *
 * @throws std::runtime_error if the validation fails
 * @return true
 */
bool
Universe::validate() const
{
  RUNTIME_EXP_IFN(this->getNrOfAtoms() == igraph_vcount(&this->graph),
                  "Validation failed: " + std::to_string(this->getNrOfAtoms()) +
                    " atoms for " +
                    std::to_string(igraph_vcount(&this->graph)) + " vertices.");

  RUNTIME_EXP_IFN(this->getNrOfBonds() == igraph_ecount(&this->graph),
                  "Validation failed: " + std::to_string(this->getNrOfBonds()) +
                    " bonds for " +
                    std::to_string(igraph_ecount(&this->graph)) + " edges.");

  return true;
}

/**
 * @brief Compute the volume of the underlying box
 *
 * @return double
 */
double
Universe::getVolume() const
{
  return this->box.getVolume();
}

size_t
Universe::getNrOfAtoms() const
{
  return this->NAtoms;
  // return this->getNrOfVertices();
}

size_t
Universe::getNrOfBonds() const
{
  return this->NBonds;
  // return this->getNrOfEdges();
}

void
Universe::setBox(const Box& passedBox, const bool rescaleAtomCoordinates)
{
  if (rescaleAtomCoordinates) {
    double scalingFactorX = passedBox.getLx() / this->box.getLx();
    double offsetX = passedBox.getLowX() - scalingFactorX * this->box.getLowX();
    igraph_vector_t xValueVec;
    REQUIRE_IGRAPH_SUCCESS(
      igraph_vector_init(&xValueVec, this->getNrOfAtoms()));
    REQUIRE_IGRAPH_SUCCESS(
      igraph_cattribute_VANV(&this->graph, "x", igraph_vss_all(), &xValueVec));

    double scalingFactorY = passedBox.getLy() / this->box.getLy();
    double offsetY = passedBox.getLowY() - scalingFactorY * this->box.getLowY();
    igraph_vector_t yValueVec;
    REQUIRE_IGRAPH_SUCCESS(
      igraph_vector_init(&yValueVec, this->getNrOfAtoms()));
    REQUIRE_IGRAPH_SUCCESS(
      igraph_cattribute_VANV(&this->graph, "y", igraph_vss_all(), &yValueVec));

    double scalingFactorZ = passedBox.getLz() / this->box.getLz();
    double offsetZ = passedBox.getLowZ() - scalingFactorZ * this->box.getLowZ();
    igraph_vector_t zValueVec;
    REQUIRE_IGRAPH_SUCCESS(
      igraph_vector_init(&zValueVec, this->getNrOfAtoms()));
    REQUIRE_IGRAPH_SUCCESS(
      igraph_cattribute_VANV(&this->graph, "z", igraph_vss_all(), &zValueVec));

    for (size_t i = 0; i < this->getNrOfAtoms(); ++i) {
      igraph_vector_set(&xValueVec,
                        i,
                        igraph_vector_get(&xValueVec, i) * scalingFactorX +
                          offsetX);
      igraph_vector_set(&yValueVec,
                        i,
                        igraph_vector_get(&yValueVec, i) * scalingFactorY +
                          offsetY);
      igraph_vector_set(&zValueVec,
                        i,
                        igraph_vector_get(&zValueVec, i) * scalingFactorZ +
                          offsetZ);
    }

    REQUIRE_IGRAPH_SUCCESS(
      igraph_cattribute_VAN_setv(&this->graph, "x", &xValueVec));
    REQUIRE_IGRAPH_SUCCESS(
      igraph_cattribute_VAN_setv(&this->graph, "y", &yValueVec));
    REQUIRE_IGRAPH_SUCCESS(
      igraph_cattribute_VAN_setv(&this->graph, "z", &zValueVec));

    igraph_vector_destroy(&xValueVec);
    igraph_vector_destroy(&yValueVec);
    igraph_vector_destroy(&zValueVec);
  }
  this->box = passedBox;
}

void
Universe::setBoxLengths(const double Lx,
                        const double Ly,
                        const double Lz,
                        bool rescaleAtomCoordinates)
{
  this->setBox(Box(Lx, Ly, Lz));
}

Box
Universe::getBox() const
{
  return this->box;
}
} // namespace pylimer_tools::entities
