#include "Molecule.h"
#include "Atom.h"
#include "../utils/GraphUtils.h"
#include "../utils/StringUtil.h"
#include <igraph/igraph.h>
#include <iostream>
#ifdef OPENMP_FOUND
#include <omp.h>
#endif

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
      if (ids.size() == 0 && this->size > 0)
      {
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
      std::vector<long int> endNodeIndices = pylimer_tools::utils::getVerticesWithDegree(this->graph, 1);
      igraph_vector_t endNodeSelectorVector;
      igraph_vector_init(&endNodeSelectorVector, endNodeIndices.size());
      pylimer_tools::utils::StdVectorToIgraphVectorT(endNodeIndices, &endNodeSelectorVector);
      igraph_vit_t vit;
      igraph_vit_create(graph, igraph_vss_vector(&endNodeSelectorVector), &vit);

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
      if (this->getNrOfBonds() == 0)
      {
        return lengths;
      }
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
      std::vector<OUT> results;
      if (this->getNrOfAtoms() == 0) {
        return results;
      }
      igraph_vector_t allValues;
      igraph_vector_init(&allValues, this->getNrOfAtoms());
      if (igraph_cattribute_VANV(this->graph, propertyName, igraph_vss_all(), &allValues))
      {
        throw std::runtime_error("Failed to query properties of molecule.");
      }
      pylimer_tools::utils::igraphVectorTToStdVector(&allValues, results);
      igraph_vector_destroy(&allValues);
      return results;
    }

    double Molecule::computeRadiusOfGyration()
    {
      double meanX, meanY, meanZ;
      // would be faster to just query the attributes.
      // But the OOP interface is just too tempting
      // as long as there are no external additional performance demands
      std::vector<Atom> allAtoms = this->getAtoms();
      double multiplier = 1 / allAtoms.size();

#pragma omp parallel for reduction(+ \
                                   : meanX, meanY, meanZ)
      for (Atom a : allAtoms)
      {
        meanX += multiplier * a.getUnwrappedX(this->parent);
        meanY += multiplier * a.getUnwrappedY(this->parent);
        meanZ += multiplier * a.getUnwrappedZ(this->parent);
      }

      Atom virtualCenterAtom = Atom(0, 0, meanX, meanY, meanZ, 0, 0, 0);

      double Rg2 = 0.0;
      for (Atom a : allAtoms)
      {
        double dist = a.distanceTo(virtualCenterAtom, this->parent);
        Rg2 += dist * dist;
      }

      return Rg2 * multiplier;
    }

    std::string Molecule::getKey() { return this->key; }

    std::vector<Atom> Molecule::getAtoms()
    {
      std::vector<Atom> results;
      size_t nrOfAtoms = this->getNrOfAtoms();
      results.reserve(nrOfAtoms);

      // #pragma omp declare reduction (merge : std::vector<Atom> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
      // #pragma omp parallel for reduction(merge: results)
      for (size_t i = 0; i < nrOfAtoms; ++i)
      {
        results.push_back(this->getAtomForVertexId(i));
      }

      return results;
    };

    std::vector<Atom> Molecule::getAtomsWithType(const int atomType)
    {
      std::vector<Atom> results;
      const std::vector<int> types = this->getPropertyValues<int>("type");
      size_t nrOfTypes = types.size();

      // #pragma omp declare reduction (merge : std::vector<Atom> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
      // #pragma omp parallel for reduction(merge: results)
      for (size_t i = 0; i < nrOfTypes; ++i)
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
