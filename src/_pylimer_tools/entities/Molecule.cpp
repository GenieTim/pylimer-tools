#include "Molecule.h"
#include "Atom.h"
#include "GraphUtils.h"

namespace pylimer_tools
{
  namespace entities
  {
    class Molecule
    {
      Molecule(Universe *parent, igraph_t *graph, MoleculeType type)
      {
        this->parent = parent;
        this->graph = graph;
        this->size = igraph_vcount(graph);
        this->typeOfThisMolecule = type;
      };

      double computeEndToEndDistance()
      {
        if (this->size < 2)
        {
          return NULL;
        }
        igraph_vs_t endNodes = pylimer_tools::utils::getVerticesWithDegree(this->graph, 1);
        igraph_vit_t vit;
        igraph_vit_create(graph, endNodes, &vit);

        double distance;
        Box box = this->parent->getBox();

        // we only compute an end-to-end distance if we have exactly two ends.
        // this is clearly not optimal, but at least unambiguous
        if (IGRAPH_VIT_SIZE(vit) == 2)
        {
          long int vertexId1 = (long int)IGRAPH_VIT_GET(vit);
          IGRAPH_VIT_NEXT(vit);
          long int vertexId2 = (long int)IGRAPH_VIT_GET(vit);
          // TODO: this is more intensive than needed
          // check whether the compiler optimizes this or not
          Atom atom1 = this->getAtomForVertexId(vertexId1);
          Atom atom2 = this->getAtomForVertexId(vertexId2);
          distance = atom1.distanceTo(atom2, box);
        }
        igraph_vit_destroy(&vit);
        return distance;
      }

      std::vector<double> computeBondLengths()
      {
        Box box = this->parent->getBox();
        std::vector<double> lengths;
        lengths.reserve(this->getNumBonds());
        // construct iterator
        igraph_eit_t bondIterator;
        if (igraph_eit_create(this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &bondIterator))
        {
          throw std::runtime_error("Cannot create iterator to loop bonds");
        }

        while (!IGRAPH_EIT_END(bondIterator))
        {
          long int edgeId = (long int)IGRAPH_EIT_GET(bondIterator);
          int bondFrom;
          int bondTo;
          igraph_edge(this->graph, edgeId, &bondFrom, &bondTo);
          // TODO: this is more intensive than needed
          // check whether the compiler optimizes this or not
          Atom atom1 = this->getAtomForVertexId(bondFrom);
          Atom atom2 = this->getAtomForVertexId(bondTo);
          lengths.push_back(atom1.distanceTo(atom2, box));
          IGRAPH_EIT_NEXT(bondIterator);
        }

        igraph_eit_destroy(&bondIterator);
        return lengths;
      }

      Atom getAtomForVertexId(long int vertexIdx)
      {
        return Atom(VAN(this->graph, "id", vertexIdx), VAN(this->graph, "type", vertexIdx), VAN(this->graph, "x", vertexIdx), VAN(this->graph, "y", vertexIdx), VAN(this->graph, "z", vertexIdx),
                    VAN(this->graph, "nx", vertexIdx), VAN(this->graph, "ny", vertexIdx), VAN(this->graph, "nz", vertexIdx));
      }

      int getLength()
      {
        return this->size;
      };

      int getNumBonds() { return igraph_ecount(this->graph); }

      MoleculeType getType()
      {
        return this->typeOfThisMolecule;
      };

      Universe *parent;
      MoleculeType typeOfThisMolecule;
      igraph_t *graph;
      int size;
    };
  }
}
