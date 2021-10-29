#ifndef UNIVERSE_H
#define UNIVERSE_H

#include <igraph/igraph.h>
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
      void setBoxLengths(const double Lx, const double Ly, const double Lz);
      void addAtoms(const int NNewAtoms, std::vector<long int> ids, std::vector<int> types, std::vector<double> x, std::vector<double> y, std::vector<double> z, std::vector<int> nx, std::vector<int> ny, std::vector<int> nz);
      void addBonds(const int NNewBonds, std::vector<long int> from, std::vector<long int> to);
      std::vector<Molecule> getMolecules(const int atomTypeToOmit);
      std::vector<Molecule> getChainsWithCrosslinker(const int crosslinkerType);
      std::map<int, int> determineFunctionalityPerType();
      std::map<int, double> computeWeightFractions(std::map<int, double> weightPerType);
      template <typename OUT> 
      std::vector<OUT> getPropertyValues(const char *propertyName);
      template <typename IN>
      long int findVertexIdForProperty(const char *propertyName, IN propertyValue);
      Atom getAtom(const int atomId);
      std::vector<Atom> getAtomsWithType(const int atomType);
      Atom getAtomByIdx(const int vertexIdx);
      double getMeanStrandLength(int junctionType);
      void setBox(Box box);
      Box getBox();
      double getVolume();
      int getNrOfAtoms();
      int getNrOfBonds();

    protected:
      // properties of the box
      int NAtoms = 0;
      int NBonds = 0;
      Box box;
      // connectivity
      igraph_t graph;
      std::map<int, int> atomIdToVectorIdx;

      igraph_vs_t getVerticesOfType(const int type);
      std::vector<long int> getIndicesOfType(const int type);
      igraph_vs_t getVerticesByIndices(std::vector<long int> indices);
    };
  }
}

#endif
