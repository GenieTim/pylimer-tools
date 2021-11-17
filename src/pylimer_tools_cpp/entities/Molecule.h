#ifndef MOLECULE_H
#define MOLECULE_H

extern "C" {
#include <igraph/igraph.h>
}
#include "Box.h"
#include "Atom.h"
#include <vector>
#include <string>

namespace pylimer_tools
{
  namespace entities
  {

    enum MoleculeType
    {
      UNDEFINED,
      NETWORK_STRAND,
      PRIMARY_LOOP,
      DANGLING_CHAIN,
      FREE_CHAIN
    };

    class Molecule
    {
    public:
      Molecule(Box *parent, igraph_t *graph, MoleculeType type);
      double computeEndToEndDistance();
      double computeRadiusOfGyration();
      std::vector<double> computeBondLengths();
      int getLength();
      MoleculeType getType();
      Atom getAtomForVertexId(long int vertexIdx);
      std::vector<Atom> getAtoms();
      std::vector<Atom> getAtomsWithType(const int atomType);
      int getNrOfBonds();
      int getNrOfAtoms();
      Box *getBox();
      std::string getKey();
      template <typename OUT>
      std::vector<OUT> getPropertyValues(const char *propertyName);
      std::vector<int> getAtomTypes() { return this->getPropertyValues<int>("type"); }

    private:
      Box *parent;
      MoleculeType typeOfThisMolecule;
      igraph_t *graph;
      int size;
      std::string key;
    };
  }
}

#endif
