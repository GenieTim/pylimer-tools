#include "Molecule.h"
#include "Atom.h"
#include "../utils/GraphUtils.h"
#include "../utils/StringUtil.h"
#include <igraph/igraph.h>
#include <iostream>

namespace pylimer_tools
{
  namespace entities
  {
    Molecule::Molecule(Box *parent, igraph_t *graph, MoleculeType type)
    {
      this->parent = parent;
      this->graph = graph;
      this->size = igraph_vcount(graph);
      this->typeOfThisMolecule = type;

      // construct a key for this molecule: a concatenation of all ids in this molecule
      if (!igraph_cattribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX, "id"))
      {
        throw std::runtime_error("Molecule's graph does not have attribute id");
      }
      igraph_vector_t allIds;
      igraph_vector_init(&allIds, this->size);
      VANV(this->graph, "id", &allIds);
      if (igraph_cattribute_VANV(this->graph, "id", igraph_vss_all(), &allIds))
      {
        throw std::runtime_error("Molecule's graph's attribute id is not accessible.");
      };
      std::vector<int> ids;
      pylimer_tools::utils::igraphVectorTToStdVector(&allIds, ids);
      if (ids.size() == 0 && this->size > 0) {
        throw std::runtime_error("Molecule's graph's attribute id was not queried.");
      }
      std::sort(ids.begin(), ids.end());
      this->key = pylimer_tools::utils::join(ids.begin(), ids.end(), std::string("-"));
      igraph_vector_destroy(&allIds);
    };

    double Molecule::computeEndToEndDistance()
    {
      if (this->size < 2)
      {
        return 0.0;
      }
      igraph_vs_t endNodes = pylimer_tools::utils::getVerticesWithDegree(this->graph, 1);
      igraph_vit_t vit;
      igraph_vit_create(graph, endNodes, &vit);

      double distance = -1.0; // TODO: find a nice default for "no end to end"
      Box *box = this->getBox();

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

    std::vector<double> Molecule::computeBondLengths()
    {
      Box *box = this->getBox();
      std::vector<double> lengths;
      lengths.reserve(this->getNrOfBonds());
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

    Atom Molecule::getAtomForVertexId(long int vertexIdx)
    {
      return Atom(VAN(this->graph, "id", vertexIdx), VAN(this->graph, "type", vertexIdx), VAN(this->graph, "x", vertexIdx), VAN(this->graph, "y", vertexIdx), VAN(this->graph, "z", vertexIdx),
                  VAN(this->graph, "nx", vertexIdx), VAN(this->graph, "ny", vertexIdx), VAN(this->graph, "nz", vertexIdx));
    }

    int Molecule::getLength()
    {
      return this->size;
    };

    int Molecule::getNrOfAtoms() { return this->size; }

    int Molecule::getNrOfBonds() { return igraph_ecount(this->graph); }

    MoleculeType Molecule::getType()
    {
      return this->typeOfThisMolecule;
    };

    Box *Molecule::getBox() { return this->parent; }

    template <typename OUT>
    std::vector<OUT> Molecule::getPropertyValues(const char *propertyName)
    {
      igraph_vector_t allValues;
      igraph_vector_init(&allValues, this->getNrOfAtoms());
      VANV(this->graph, propertyName, &allValues);
      std::vector<OUT> results;
      pylimer_tools::utils::igraphVectorTToStdVector(&allValues, results);
      igraph_vector_destroy(&allValues);
      return results;
    }

    std::string Molecule::getKey() { return this->key; }

    std::vector<Atom> Molecule::getAtomsWithType(const int atomType)
    {
      std::vector<int> types = this->getPropertyValues<int>("type");
      std::vector<Atom> results;

      for (size_t i = 0; i < types.size(); ++i)
      {
        if (types[i] == atomType)
        {
          results.push_back(this->getAtomForVertexId(i));
        }
      }

      return results;
    };
  }
}
