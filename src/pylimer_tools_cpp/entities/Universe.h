#ifndef UNIVERSE_H
#define UNIVERSE_H

extern "C"
{
#include <igraph/igraph.h>
}
#include <vector>
#include <map>
#include "Molecule.h"
#include "Atom.h"

namespace pylimer_tools
{
  namespace entities
  {

    class Universe
    {
    public:
      Universe(const double Lx, const double Ly, const double Lz);
      // rule of three:
      // 1. destructor (to destroy the graph)
      ~Universe();
      // 2. copy constructor
      Universe(const Universe &src);
      // 3. copy assignment operator
      Universe &operator=(Universe src);

      // initilaization/setters
      void setBoxLengths(const double Lx, const double Ly, const double Lz);
      void addAtoms(const size_t NNewAtoms, std::vector<long int> ids, std::vector<int> types, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<int> nx, std::vector<int> ny, std::vector<int> nz);
      void addBonds(const size_t NNewBonds, std::vector<long int> from, std::vector<long int> to);
      void addBonds(const size_t NNewBonds, std::vector<long int> from, std::vector<long int> to, std::vector<int> bondTypes);
      void setMasses(std::map<int, double> weightPerType);
      void setBox(Box box);
      void setTimestep(long int timestep) { this->timestep = timestep; };

      // getters
      Atom getAtom(const int atomId) const;
      std::vector<Atom> getAtomsWithType(const int atomType);
      std::vector<Atom> getAtoms();
      Atom getAtomByIdx(const int vertexIdx) const;
      std::vector<Molecule> getMolecules(const int atomTypeToOmit = -1);
      std::vector<Molecule> getChainsWithCrosslinker(const int crosslinkerType);
      template <typename OUT>
      std::vector<OUT> getPropertyValues(const char *propertyName);
      std::vector<int> getAtomTypes() { return this->getPropertyValues<int>("type"); }
      template <typename IN>
      long int findVertexIdForProperty(const char *propertyName, IN propertyValue);
      Box getBox();
      std::map<std::string, std::vector<long int>> getBonds();
      double getVolume();
      const int getNrOfAtoms() const;
      const int getNrOfBonds() const;
      std::map<int, double> getMasses();
      long int getTimestep() { return this->timestep; };
      int getNrOfBondsOfAtom(const long int atomId);
      int getNrOfBondsOfVertex(const long int vertexId);

      // operators
      Atom operator[](size_t index) const { return this->getAtom(index); }

      // computations
      std::map<int, int> determineFunctionalityPerType();
      std::map<int, double> computeWeightFractions();
      double getMeanStrandLength(int junctionType);
      bool validate();

    protected:
      // properties of the universe
      long int timestep;
      int NAtoms = 0;
      int NBonds = 0;
      Box box;
      // connectivity
      igraph_t graph;
      std::map<int, int> atomIdToVectorIdx;

      // type's properties
      std::map<int, double> weightPerType; // a dictionary with key: type, and value: weight per atom of this atom type.

      // internal functions
      igraph_vs_t
      getVerticesOfType(const int type);
      std::vector<long int> getIndicesOfType(const int type);
      igraph_vs_t getVerticesByIndices(std::vector<long int> indices);
    };
  }
}

#endif
