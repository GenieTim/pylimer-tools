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
      Atom getAtom(const int atomId);
      Atom *getAtomsWithType(const int atomType);
      void setBox(Box box);
      Box getBox();
      double getVolume();
      int getNrOfAtoms();
    };
  }
}

#endif
