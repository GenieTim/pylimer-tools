#include "Universe.h"
#include "../utils/VectorUtils.h"
#include "../utils/GraphUtils.h"
#include "Box.h"

#include <vector>
#include <map>
#include <iterator> // for back_inserter

namespace pylimer_tools
{
  namespace entities
  {

    Universe::Universe(const double Lx, const double Ly, const double Lz)
    {
      box = Box(Lx, Ly, Lz);
      /* turn on attribute handling: TODO: move to some main() function  */
      igraph_set_attribute_table(&igraph_cattribute_table);
      igraph_vector_t gtypes, vtypes, etypes;
      igraph_strvector_t gnames, vnames, enames;
      igraph_vector_t vec;
      igraph_strvector_t svec;

      igraph_vector_init(&gtypes, 0);
      igraph_vector_init(&vtypes, 0);
      igraph_vector_init(&etypes, 0);
      igraph_strvector_init(&gnames, 0);
      igraph_strvector_init(&vnames, 0);
      igraph_strvector_init(&enames, 0);

      // start setting properties
      igraph_empty(&this->graph, 0, IGRAPH_UNDIRECTED);

      //
      igraph_cattribute_list(&this->graph, &gnames, &gtypes, &vnames, &vtypes,
                             &enames, &etypes);
    }

    void Universe::addAtoms(const int NNewAtoms, std::vector<long int> newIds, std::vector<int> newTypes, std::vector<double> newX, std::vector<double> newY, std::vector<double> newZ, std::vector<int> newNx, std::vector<int> newNy, std::vector<int> newNz)
    {
      if (newTypes.size() != NNewAtoms || newIds.size() != newTypes.size() || newX.size() != newNx.size() || newY.size() != newNy.size() || newZ.size() != newNz.size() || newX.size() != newY.size() || NNewAtoms != newZ.size())
      {
        throw std::invalid_argument("All inputs must have the same size.");
      }
      // actually add the vertices
      if (igraph_add_vertices(&this->graph, NNewAtoms, 0))
      {
        throw std::runtime_error("Failed to add new atoms to graph.");
      }
      // do map for easy access afterwards
      for (int i = 0; i < NNewAtoms; ++i)
      {
        this->atomIdToVectorIdx[newIds[i]] = this->NAtoms + i;
        // append attributes
        igraph_cattribute_VAN_set(&this->graph, "id", this->NAtoms + i, newIds[i]);
        igraph_cattribute_VAN_set(&this->graph, "x", this->NAtoms + i, newX[i]);
        igraph_cattribute_VAN_set(&this->graph, "y", this->NAtoms + i, newY[i]);
        igraph_cattribute_VAN_set(&this->graph, "z", this->NAtoms + i, newZ[i]);
        igraph_cattribute_VAN_set(&this->graph, "type", this->NAtoms + i, newTypes[i]);
        igraph_cattribute_VAN_set(&this->graph, "nx", this->NAtoms + i, newNx[i]);
        igraph_cattribute_VAN_set(&this->graph, "ny", this->NAtoms + i, newNy[i]);
        igraph_cattribute_VAN_set(&this->graph, "nz", this->NAtoms + i, newNz[i]);
      }
      this->NAtoms += NNewAtoms;
    }

    void Universe::addBonds(const int NNewBonds, std::vector<long int> from, std::vector<long int> to)
    {
      if (from.size() != to.size() || from.size() != NNewBonds)
      {
        throw std::invalid_argument("All inputs must have the same size");
      }
      std::vector<long int> newEdgesVector = pylimer_tools::utils::interleave(from, to);
      // translate from atomId to VertexIdx
      for (int i = 0; i < newEdgesVector.size(); ++i)
      {
        newEdgesVector[i] = this->atomIdToVectorIdx[i];
      }
      // add the new edges
      igraph_vector_t newEdges;
      pylimer_tools::utils::StdVectorToIgraphVectorT(newEdgesVector, &newEdges);
      if (igraph_add_edges(&this->graph, &newEdges, 0))
      {
        throw std::runtime_error("Failed to add edges to graph.");
      }
      this->NBonds += NNewBonds;
    }

    std::vector<Molecule> Universe::getMolecules(const int atomTypeToOmit)
    {
      std::vector<Molecule> molecules;
      if (this->getNrOfAtoms() == 0)
      {
        return molecules;
      }
      // make a copy to remove crosslinkers from
      igraph_t graphWithoutCrosslinkers;
      if (igraph_copy(&graphWithoutCrosslinkers, &this->graph))
      {
        throw std::runtime_error("Failed to copy graph.");
      }
      // select vertices of type
      std::vector<long int> indicesToRemove = this->getIndicesOfType(atomTypeToOmit);
      igraph_vs_t verticesToRemove = this->getVerticesByIndices(indicesToRemove);
      // remove elements of type
      if (igraph_delete_vertices(&graphWithoutCrosslinkers, verticesToRemove))
      {
        throw std::runtime_error("Failed to delete crosslinkers from graph.");
      }
      // split the copy into the separate components
      igraph_vector_ptr_t components;
      igraph_vector_ptr_init(&components, 3);
      int NComponents = igraph_decompose(&graphWithoutCrosslinkers, &components, IGRAPH_STRONG, -1, 0);
      molecules.reserve(NComponents);
      for (int i = 0; i < NComponents; ++i)
      {
        igraph_t *g = (igraph_t *)VECTOR(components)[i];

        molecules.push_back(Molecule(&this->box, g, MoleculeType::UNDEFINED));
      }
      igraph_decompose_destroy(&components);
      igraph_destroy(&graphWithoutCrosslinkers);
      return molecules;
    }

    igraph_vs_t Universe::getVerticesOfType(const int type)
    {
      std::vector<long int> indices = this->getIndicesOfType(type);
      return this->getVerticesByIndices(indices);
    }

    igraph_vs_t Universe::getVerticesByIndices(std::vector<long int> indices)
    {
      igraph_vector_t indicesToSelect;
      igraph_vector_init(&indicesToSelect, indices.size());
      pylimer_tools::utils::StdVectorToIgraphVectorT(indices, &indicesToSelect);
      return igraph_vss_vector(&indicesToSelect);
    }

    std::vector<long int> Universe::getIndicesOfType(const int type)
    {
      igraph_vector_t types;
      igraph_vector_init(&types, this->getNrOfAtoms());
      VANV(&this->graph, "type", &types);
      std::vector<long int> indices;
      for (int i = 0; i < this->NAtoms; ++i)
      {
        if (VECTOR(types)[i] == type)
        {
          indices.push_back(i);
        }
      }
      return indices;
    }

    std::vector<Molecule> Universe::getChainsWithCrosslinker(const int crosslinkerType)
    {
      std::vector<Molecule> molecules;
      if (this->getNrOfAtoms() == 0)
      {
        return molecules;
      }
      // make a copy to remove crosslinkers from
      igraph_t graphWithoutCrosslinkers;
      if (igraph_copy(&graphWithoutCrosslinkers, &this->graph))
      {
        throw std::runtime_error("Failed to copy graph.");
      }
      // select vertices of type
      igraph_vs_t verticesToRemove = this->getVerticesOfType(crosslinkerType);
      // remove elements of type
      if (igraph_delete_vertices(&graphWithoutCrosslinkers, verticesToRemove))
      {
        throw std::runtime_error("Failed to delete crosslinkers from graph.");
      }
      // split the copy into the separate
      igraph_vector_ptr_t components;
      igraph_vector_ptr_init(&components, 3);
      int NComponents = igraph_decompose(&graphWithoutCrosslinkers, &components, IGRAPH_STRONG, -1, 0);
      molecules.reserve(NComponents);
      for (int i = 0; i < NComponents; ++i)
      {
        // loop the chains to add the crosslinkers back
        igraph_t *chain = (igraph_t *)VECTOR(components)[i];
        int moleculeLengthBefore = igraph_vcount(chain);
        igraph_vs_t endNodes = pylimer_tools::utils::getVerticesWithDegree(chain, 1);
        MoleculeType molType = MoleculeType::UNDEFINED;
        bool isLoop = false;

        if (moleculeLengthBefore > 1)
        {
          // no care about single atoms for now (though they are included)
          igraph_vit_t endNodeVit;
          igraph_vit_create(chain, endNodes, &endNodeVit);
          // loop end nodes
          while (!IGRAPH_VIT_END(endNodeVit))
          {
            long int newVertexId = (long int)IGRAPH_VIT_GET(endNodeVit);
            long int originalVertexId = this->findVertexIdForProperty("id", VAN(&graph, "id", newVertexId));
            igraph_vs_t neighbours;
            if (igraph_vs_adj(&neighbours, originalVertexId, IGRAPH_ALL))
            {
              throw std::runtime_error("Failed to get neighbours in graph");
            }

            igraph_vit_t originalNeighbourVit;
            igraph_vit_create(&graph, neighbours, &originalNeighbourVit);

            // loop neighbours
            while (!IGRAPH_VIT_END(originalNeighbourVit))
            {
              long int neighbourId = (long int)IGRAPH_VIT_GET(originalNeighbourVit);
              if (igraph_cattribute_VAN(&graph, "type", neighbourId) == crosslinkerType)
              {
                // check if this crosslinker exists in the current chain to find loops
                igraph_vs_t neighbourInChain;
                if (igraph_vs_1(&neighbourInChain, neighbourId))
                {
                  isLoop = true;
                }
                else
                {
                  // check whether the neighbour has been found in the chain
                  igraph_vit_t oneNeighbourVit;
                  igraph_vit_create(chain, neighbourInChain, &oneNeighbourVit);
                  if (IGRAPH_VIT_SIZE(oneNeighbourVit) == 1)
                  {
                    isLoop = true;
                  }
                  else
                  {
                    // add crosslinker back to chain
                    igraph_add_vertices(chain, 1, 0);
                    long int newCrosslinkerVertexIdx = igraph_vcount(chain) - 1;
                    // including bond, of course
                    igraph_add_edge(chain, newVertexId, newCrosslinkerVertexIdx);
                    // copy all attributes
                    for (auto property : {"id", "type", "x", "y", "z", "nx", "ny", "nz"})
                    {
                      SETVAN(chain, property, newCrosslinkerVertexIdx, VAN(&graph, property, neighbourId));
                    }
                  }
                }
                igraph_vs_destroy(&neighbourInChain);
              }

              IGRAPH_VIT_NEXT(originalNeighbourVit);
            }

            igraph_vit_destroy(&originalNeighbourVit);
            igraph_vs_destroy(&neighbours);
            IGRAPH_VIT_NEXT(endNodeVit);
          }
          igraph_vit_destroy(&endNodeVit);
        }
        igraph_vs_destroy(&endNodes);
        // decide on molecule type
        int newMoleculeLength = igraph_vcount(chain);
        if (newMoleculeLength == moleculeLengthBefore)
        {
          molType = MoleculeType::FREE_CHAIN;
        }
        else if (newMoleculeLength == moleculeLengthBefore + 1)
        {
          molType = MoleculeType::DANGLING_CHAIN;
        }
        else if (newMoleculeLength == moleculeLengthBefore + 2)
        {
          molType = MoleculeType::NETWORK_STRAND;
        }
        if (isLoop)
        {
          molType = MoleculeType::PRIMARY_LOOP;
        }

        // finally, create the molecule/chain
        molecules.push_back(Molecule(&this->box, chain, molType));
      }
      igraph_decompose_destroy(&components);
      igraph_destroy(&graphWithoutCrosslinkers);

      return molecules;
    }

    std::map<int, int> Universe::determineFunctionalityPerType()
    {
      std::map<int, int> result;
      igraph_vector_t degrees;
      if (igraph_vector_init(&degrees, 0))
      {
        throw std::runtime_error("Failed to instantiate result vector.");
      }
      igraph_vs_t allVertexIds;
      igraph_vs_all(&allVertexIds);
      // complexity: O(|v|*d)
      if (igraph_degree(&this->graph, &degrees, allVertexIds, IGRAPH_ALL, false))
      {
        throw std::runtime_error("Failed to determine degree of vertices");
      }

      igraph_vector_t typesVec;
      igraph_vector_init(&typesVec, this->getNrOfAtoms());
      VANV(&this->graph, "type", &typesVec);
      std::vector<long int> types;
      pylimer_tools::utils::igraphVectorTToStdVector(&typesVec, types);
      igraph_vector_destroy(&typesVec);
      std::vector<long int> uniqueTypes;
      copy(std::begin(types), std::end(types), std::back_inserter(uniqueTypes));
      auto uniqueTypesIter = std::unique(std::begin(types), std::end(types));
      uniqueTypes.erase(uniqueTypesIter, uniqueTypes.end());
      // make sure the keys are (re)set, for every type
      for (int type : uniqueTypes)
      {
        result[type] = 0;
      }

      // complexity: O(|V|)
      igraph_vit_t vit;
      igraph_vit_create(&graph, allVertexIds, &vit);
      while (!IGRAPH_VIT_END(vit))
      {
        long int vertexId = (long int)IGRAPH_VIT_GET(vit);
        result[types[vertexId]] = std::max((int)igraph_vector_e(&degrees, vertexId), result[types[vertexId]]);
        IGRAPH_VIT_NEXT(vit);
      }
      igraph_vit_destroy(&vit);
      igraph_vs_destroy(&allVertexIds);

      return result;
    }

    Atom Universe::getAtom(const int atomId)
    {
      int atomIdx = this->atomIdToVectorIdx[atomId];
      return this->getAtomByIdx(atomIdx);
    }

    Atom Universe::getAtomByIdx(const int vertexIdx)
    {
      return Atom(VAN(&this->graph, "id", vertexIdx), VAN(&this->graph, "type", vertexIdx), VAN(&this->graph, "x", vertexIdx), VAN(&this->graph, "y", vertexIdx), VAN(&this->graph, "z", vertexIdx),
                  VAN(&this->graph, "nx", vertexIdx), VAN(&this->graph, "ny", vertexIdx), VAN(&this->graph, "nz", vertexIdx));
    }

    std::vector<Atom> Universe::getAtomsWithType(const int atomType)
    {
      std::vector<Atom> atoms;
      auto indicesWithType = this->getIndicesOfType(atomType);
      atoms.reserve(indicesWithType.size());
      for (auto idx : indicesWithType)
      {
        atoms.push_back(this->getAtomByIdx(idx));
      }

      return atoms;
    }

    template <typename IN>
    long int Universe::findVertexIdForProperty(const char *propertyName, IN propertyValue)
    {
      igraph_vector_t allValues;
      igraph_vector_init(&allValues, this->getNrOfAtoms());
      VANV(&this->graph, propertyName, &allValues);
      for (int i = 0; i < this->NAtoms; ++i)
      {
        if (VECTOR(allValues)[i] == propertyValue)
        {
          return i;
        }
      }
      return -1;
    }

    double Universe::getVolume()
    {
      return this->box.getVolume();
    }

    int Universe::getNrOfAtoms()
    {
      return this->NAtoms;
    }

    void Universe::setBox(Box box)
    {
      this->box = box;
    }

    Box Universe::getBox() { return this->box; }
  }
}
