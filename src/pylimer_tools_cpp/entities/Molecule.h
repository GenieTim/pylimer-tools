#ifndef MOLECULE_H
#define MOLECULE_H

extern "C"
{
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
      // getters
      int getLength() const;
      MoleculeType getType();
      Atom getAtomForVertexId(long int vertexIdx) const;
      std::vector<Atom> getAtoms();
      std::vector<Atom> getAtomsWithType(const int atomType);
      std::vector<Atom> getAtomsOfDegree(const int degree);
      int getNrOfBonds() const;
      int getNrOfAtoms() const;
      Box *getBox();
      std::string getKey();
      template <typename OUT>
      std::vector<OUT> getPropertyValues(const char *propertyName);
      std::vector<int> getAtomTypes() { return this->getPropertyValues<int>("type"); }

      // computations
      double computeEndToEndDistance();
      double computeRadiusOfGyration();
      std::vector<double> computeBondLengths();

      // operators
      Atom operator[](size_t index) const
      {
        return this->getAtomForVertexId(index);
      }
 
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
