#ifndef MOLECULE_H
#define MOLECULE_H

#include <igraph/igraph.h>
#include "Box.h"
#include "Atom.h"
#include "Universe.h"

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
      Molecule(Universe *parent, igraph_t *graph, MoleculeType type);
      Molecule *decomposeFurther(int atomTypeToOmit = 0);
      double computeEndToEndDistance();
      std::vector<double> computeBondLengths();
      int getLength();
      MoleculeType getType();
    };
  }
}

#endif
