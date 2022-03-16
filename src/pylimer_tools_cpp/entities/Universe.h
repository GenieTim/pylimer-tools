#ifndef UNIVERSE_H
#define UNIVERSE_H

extern "C" {
#include <igraph/igraph.h>
}
#include "Atom.h"
#include "AtomGraphParent.h"
#include "Molecule.h"
#include <map>
#include <unordered_map>
#include <vector>

namespace pylimer_tools {
namespace entities {

class Universe : public AtomGraphParent {
public:
  Universe(const double Lx, const double Ly, const double Lz);
  Universe(Box box);

  // rule of three:
  // 1. destructor (to destroy the graph)
  ~Universe();
  // 2. copy constructor
  Universe(const Universe &src);
  // 3. copy assignment operator
  Universe &operator=(Universe src);

  // initilaization/setters
  void setBoxLengths(const double Lx, const double Ly, const double Lz);
  void addAtoms(std::vector<long int> ids, std::vector<int> types,
                std::vector<double> x, std::vector<double> y,
                std::vector<double> z, std::vector<int> nx, std::vector<int> ny,
                std::vector<int> nz);
  void addAtoms(const size_t NNewAtoms, std::vector<long int> ids,
                std::vector<int> types, std::vector<double> x,
                std::vector<double> y, std::vector<double> z,
                std::vector<int> nx, std::vector<int> ny, std::vector<int> nz);
  void addBonds(std::vector<long int> from, std::vector<long int> to);
  void addBonds(const size_t NNewBonds, std::vector<long int> from,
                std::vector<long int> to);
  void addBonds(const size_t NNewBonds, std::vector<long int> from,
                std::vector<long int> to, std::vector<int> bondTypes,
                const bool ignoreNonExistentAtoms = false,
                const bool simplify = true);
  void addAngles(std::vector<long int> from, std::vector<long int> via,
                 std::vector<long int> to);
  void setMasses(std::map<int, double> weightPerType);
  void setBox(Box box);
  void setTimestep(long int timestep) { this->timestep = timestep; };
  void simplify();

  // getters
  Atom getAtom(const int atomId) const;
  std::vector<Atom> getAtoms();
  // std::map<std::string, std::vector<long int>> getBonds() const;
  std::map<std::string, std::vector<long int>> getAngles() const;
  std::vector<Molecule> getClusters() const;
  std::vector<Molecule> getMolecules(const int atomTypeToOmit = -1) const;
  std::vector<Molecule>
  getChainsWithCrosslinker(const int crosslinkerType) const;
  Universe getNetworkOfCrosslinker(const int crosslinkerType) const;
  // TODO: find & implement a better return type, e.g. std::vector<Molecule>
  std::map<int, std::vector<std::vector<Atom>>>
  findLoops(const int crosslinkerType, const int maxLength = -1, bool skipSelfLoops = false) const;
  bool hasInfiniteStrand(const int crosslinkerType, const int maxLength = -1) const;
  std::vector<int> getAtomTypes() const {
    return this->getPropertyValues<int>("type");
  }
  std::map<int, int> countAtomTypes() const;
  template <typename IN>
  long int findVertexIdForProperty(const char *propertyName,
                                   IN propertyValue) const;
  Box getBox();
  double getVolume();
  size_t getNrOfAtoms() const;
  size_t getNrOfBonds() const;
  size_t getNrOfAngles() const;
  std::map<int, double> getMasses();
  long int getTimestep() { return this->timestep; };
  long int getAtomIdByIdx(const int vertexId) const;
  long int getIdxByAtomId(const int atomId) const;

  // operators
  Atom operator[](size_t index) const { return this->getAtom(index); }

  // computations
  std::map<std::string, std::vector<long int>> detectAngles() const;
  std::map<int, int> determineFunctionalityPerType() const;
  std::map<int, double> determineEffectiveFunctionalityPerType() const;
  std::map<int, double> computeWeightFractions() const;
  std::vector<double> computeDxs(const std::vector<int> bondFrom,
                                 const std::vector<int> bondTo);
  std::vector<double> computeDys(const std::vector<int> bondFrom,
                                 const std::vector<int> bondTo);
  std::vector<double> computeDzs(const std::vector<int> bondFrom,
                                 const std::vector<int> bondTo);
  std::vector<double> computeBondLengths() {
    return AtomGraphParent::computeBondLengths(&this->box);
  };
  double getMeanStrandLength(int junctionType);
  double computeMeanEndToEndDistance(int junctionType);
  double computeMeanBondLength();
  bool validate();

protected:
  // properties of the universe
  long int timestep;
  size_t NAtoms = 0;
  size_t NBonds = 0;
  Box box;
  // connectivity
  // igraph_t graph;
  std::unordered_map<int, int> atomIdToVectorIdx;
  // extra info
  std::vector<long int> angleFrom;
  std::vector<long int> angleTo;
  std::vector<long int> angleVia;

  // type's properties
  std::map<int, double>
      weightPerType; // a dictionary with key: type, and value: weight per atom
                     // of this atom type.

  // internal functions
  igraph_vs_t getVerticesOfType(const int type) const;
  std::vector<long int> getIndicesOfType(const int type) const;
  igraph_vs_t getVerticesByIndices(std::vector<long int> indices) const;
  std::vector<double> computeDs(const std::vector<int> bondFrom,
                                const std::vector<int> bondTo,
                                std::string direction, double boxLimit) const;
};
} // namespace entities
} // namespace pylimer_tools

#endif
