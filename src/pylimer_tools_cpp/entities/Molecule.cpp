#include "Molecule.h"
#include "../utils/GraphUtils.h"
#include "../utils/StringUtils.h"
#include "Atom.h"
#include <set>
extern "C" {
#include <igraph/igraph.h>
}
#include <iostream>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

namespace pylimer_tools {
namespace entities {

Molecule::Molecule(const Box *parent, const igraph_t *ingraph,
                   MoleculeType type, std::map<int, double> weightPerType) {
  this->parent = parent;
  this->initializeFromGraph(ingraph);
  this->typeOfThisMolecule = type;
  this->weightPerType = weightPerType;
};

void Molecule::initializeFromGraph(const igraph_t *ingraph) {
  igraph_copy(&this->graph, ingraph);
  this->size = igraph_vcount(&this->graph);
  // construct a key for this molecule: a concatenation of all ids in this
  // molecule
  if (!igraph_cattribute_has_attr(&this->graph, IGRAPH_ATTRIBUTE_VERTEX,
                                  "id")) {
    throw std::runtime_error("Molecule's graph does not have attribute id");
  }
  igraph_vector_t allIds;
  igraph_vector_init(&allIds, this->size);
  VANV(&this->graph, "id", &allIds);
  if (igraph_cattribute_VANV(&this->graph, "id", igraph_vss_all(), &allIds)) {
    throw std::runtime_error(
        "Molecule's graph's attribute id is not accessible.");
  };
  std::vector<int> ids;
  pylimer_tools::utils::igraphVectorTToStdVector(&allIds, ids);
  if (ids.size() == 0 && this->size > 0) {
    throw std::runtime_error(
        "Molecule's graph's attribute id was not queried.");
  }
  this->atomIdToVectorIdx.reserve(ids.size());
  for (int i = 0; i < ids.size(); ++i) {
    this->atomIdToVectorIdx[ids[i]] = i;
  }
  std::sort(ids.begin(), ids.end());
  this->key =
      pylimer_tools::utils::join(ids.begin(), ids.end(), std::string("-"));
  igraph_vector_destroy(&allIds);
};

// rule of three:
// 1. destructor (to destroy the graph)
Molecule::~Molecule() {
  // in addition to basic fields being deleted, we need to clean up the graph
  // as is done in parent
  igraph_destroy(&this->graph);
};
// 2. copy constructor
Molecule::Molecule(const Molecule &src)
    : Molecule(src.parent, &src.graph, src.typeOfThisMolecule,
               src.weightPerType){};
// 3. copy assignment operator
Molecule &Molecule::operator=(Molecule src) {
  std::swap(this->parent, src.parent);
  std::swap(this->typeOfThisMolecule, src.typeOfThisMolecule);
  std::swap(this->size, src.size);
  std::swap(this->key, src.key);
  std::swap(this->_boxNoUse, src._boxNoUse);
  std::swap(this->weightPerType, src.weightPerType);
  std::swap(this->graph, src.graph);

  return *this;
};

double Molecule::computeEndToEndDistance() {
  if (this->size < 2) {
    return 0.0;
  }

  std::vector<Atom> endNodes = this->getAtomsOfDegree(1);

  double distance = -1.0; // TODO: find a nice default for "no end to end"

  // we only compute an end-to-end distance if we have exactly two ends.
  // this is clearly not optimal, but at least unambiguous
  if (endNodes.size() == 2) {
    // TODO: this is more intensive than needed
    // check whether the compiler optimizes this or not
    Atom atom1 = endNodes[0];
    Atom atom2 = endNodes[1];
    distance = atom1.distanceTo(atom2, this->parent);
  }
  return distance;
}
/**
 * @brief compute the weight of this molecule
 *
 * @return double the total weight
 */
double Molecule::computeWeight() {
  std::vector<int> presentTypes = this->getPropertyValues<int>("type");
  double totalWeight = 0.0;
  for (int type : presentTypes) {
    totalWeight += this->weightPerType[type];
  }
  return totalWeight;
}


long int Molecule::getAtomIdByIdx(const int vertexId) const {
  return VAN(&this->graph, "id", vertexId);
};

long int Molecule::getIdxByAtomId(const int atomId) const {
  if (!this->atomIdToVectorIdx.contains(atomId)) {
    throw std::invalid_argument("Molecule cannot return vertex idx of this atom: an atom with this id (" + std::to_string(atomId) +
                                ") does not exist");
  }
  return this->atomIdToVectorIdx.at(atomId);
};

/**
 * @brief Get the nr of atoms in the molecule
 *
 * @return int
 */
int Molecule::getLength() const { return this->size; };

/**
 * @brief Get the nr of atoms in the molecule
 *
 * @return int
 */
int Molecule::getNrOfAtoms() const { return this->size; }

/**
 * @brief Get the nr of bonds in the molecule
 *
 * @return int
 */
int Molecule::getNrOfBonds() const { return igraph_ecount(&this->graph); }

/**
 * @brief Get the type of the molecule
 *
 * @return MoleculeType
 */
MoleculeType Molecule::getType() { return this->typeOfThisMolecule; };

const Box *Molecule::getBox() const { return this->parent; }

double Molecule::computeRadiusOfGyration() {
  double meanX = 0.0, meanY = 0.0, meanZ = 0.0;
  // would be faster to just query the attributes.
  // But the OOP interface is just too tempting
  // as long as there are no external additional performance demands
  std::vector<Atom> allAtoms = this->getAtoms();
  double multiplier = 1 / allAtoms.size();

#pragma omp parallel for reduction(+ : meanX, meanY, meanZ)
  for (Atom a : allAtoms) {
    meanX += multiplier * a.getUnwrappedX(this->parent);
    meanY += multiplier * a.getUnwrappedY(this->parent);
    meanZ += multiplier * a.getUnwrappedZ(this->parent);
  }

  Atom virtualCenterAtom = Atom(0, 0, meanX, meanY, meanZ, 0, 0, 0);

  double Rg2 = 0.0;
  for (Atom a : allAtoms) {
    double dist = a.distanceTo(virtualCenterAtom, this->parent);
    Rg2 += dist * dist;
  }

  return Rg2 * multiplier;
}

std::string Molecule::getKey() const { return this->key; }

std::vector<Atom> Molecule::getAtoms() {
  std::vector<Atom> results;
  size_t nrOfAtoms = this->getNrOfAtoms();
  results.reserve(nrOfAtoms);

  // #pragma omp declare reduction (merge : std::vector<Atom> :
  // omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
  // #pragma omp parallel for reduction(merge: results)
  for (size_t i = 0; i < nrOfAtoms; ++i) {
    results.push_back(this->getAtomByVertexIdx(i));
  }

  return results;
};

std::vector<Atom> Molecule::getAtomsLinedUp() {
  std::vector<Atom> results;
  size_t nrOfAtoms = this->getNrOfAtoms();
  results.reserve(nrOfAtoms);

  long int vertexIdToStartWith = 0;
  std::vector<long int> ends =
      pylimer_tools::utils::getVerticesWithDegree(&this->graph, 1);
  if (ends.size() > 0) {
    vertexIdToStartWith = ends[0];
  }

  std::vector<long int> connections =
      this->getVertexIdxsConnectedTo(vertexIdToStartWith);
  results.push_back(this->getAtomByVertexIdx(vertexIdToStartWith));
  for (long int connection : connections) {
    long int currentCenter = connection;
      results.push_back(this->getAtomByVertexIdx(currentCenter));
    long int lastCenter = vertexIdToStartWith;
    std::vector<long int> subConnections =
        this->getVertexIdxsConnectedTo(currentCenter);
    while (subConnections.size() > 0) {
      if (subConnections.size() == 1) {
        break;
      }
      // we assume a functionality of 2 for ordinary strands
      if (subConnections.size() != 2) {
        throw std::runtime_error(
            "Failed to align all atoms on one strand, as a functionality of " +
            std::to_string(subConnections.size()) +
            " was found and 1 or 2 expected.");
      }
      int subConnectionDirection = (subConnections[0] == lastCenter) ? 1 : 0;
      if (subConnectionDirection == vertexIdToStartWith) {
        break;
      }
      lastCenter = currentCenter;
      currentCenter = subConnections[subConnectionDirection];
      results.push_back(this->getAtomByVertexIdx(currentCenter));
      subConnections = this->getVertexIdxsConnectedTo(currentCenter);
    }
  }

  if (results.size() != this->getNrOfAtoms()) {
    throw std::runtime_error(
        "Failed to align all atoms on one strand: Lined up " +
        std::to_string(results.size()) + " instead of " +
        std::to_string(this->getNrOfAtoms()) + " atoms.");
  }
  return results;
};

} // namespace entities
} // namespace pylimer_tools
