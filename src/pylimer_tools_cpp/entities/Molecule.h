#ifndef MOLECULE_H
#define MOLECULE_H

extern "C" {
#include <igraph/igraph.h>
}
#include "Atom.h"
#include "AtomGraphParent.h"
#include "Box.h"
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace pylimer_tools {
namespace entities {

enum MoleculeType {
  UNDEFINED,
  NETWORK_STRAND,
  PRIMARY_LOOP,
  DANGLING_CHAIN,
  FREE_CHAIN
};

class Molecule : public AtomGraphParent {
public:
  Molecule(const Box *parent, const igraph_t *graph, MoleculeType type,
           std::map<int, double> weightPerType);

  // rule of three:
  // 1. destructor (to destroy the graph)
  ~Molecule();
  // 2. copy constructor
  Molecule(const Molecule &src);
  // 3. copy assignment operator
  Molecule &operator=(Molecule src);
  // getters
  int getLength() const;
  MoleculeType getType();
  std::vector<Atom> getAtoms();
  // std::map<std::string, std::vector<long int>> getBonds() const;
  std::vector<Atom> getAtomsLinedUp();
  int getNrOfBonds() const;
  int getNrOfAtoms() const;
  const Box *getBox() const;
  std::string getKey() const;
  std::vector<int> getAtomTypes() {
    return this->getPropertyValues<int>("type");
  }
  long int getAtomIdByIdx(const int vertexId) const;
  long int getIdxByAtomId(const int atomId) const;

  // computations
  double computeEndToEndDistance();
  double computeRadiusOfGyration();
  double computeWeight();
  std::vector<double> computeBondLengths() {
    return AtomGraphParent::computeBondLengths(this->parent);
  };

  // operators
  Atom operator[](size_t index) const {
    return this->getAtomByVertexIdx(index);
  }

private:
  Box _boxNoUse;
  const Box *parent;
  MoleculeType typeOfThisMolecule;
  int size;
  std::string key;
  std::map<int, double> weightPerType;
  std::unordered_map<int, int> atomIdToVectorIdx;

  void initializeFromGraph(const igraph_t *ingraph);
};
} // namespace entities
} // namespace pylimer_tools

#endif
