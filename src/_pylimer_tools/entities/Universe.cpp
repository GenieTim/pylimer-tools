#include "Universe.h"
#include "vector_utils.h"
#include "graph_utils.h"

#include <vector>
#include <map>
#include <iterator> // for back_inserter

namespace pylimer_tools
{
  namespace entities
  {

    class Universe
    {
      Universe(const double Lx, const double Ly, const double Lz)
      {
        this->setBoxLengths(Lx, Ly, Lz);
        igraph_empty(&this->graph, 0, IGRAPH_UNDIRECTED);
      }

      void setBoxLengths(const double Lx, const double Ly, const double Lz)
      {
        this->Lx = Lx;
        this->Ly = Ly;
        this->Lz = Lz;
      }

      void addAtoms(const int NNewAtoms, std::vector<int> newIds, std::vector<int> newTypes, std::vector<double> newX, std::vector<double> newY, std::vector<double> newZ, std::vector<double> newNx, std::vector<double> newNy, std::vector<double> newNz)
      {
        if (types.size() != NNewAtoms || x.size() != nx.size() || y.size() != ny.size() || z.size() != nz.size())
        {
          throw std::invalid_argument("All inputs must have the same size.");
        }
        // actually add the vertices
        if (igraph_add_vertices(&this->graph, NNewAtoms, 0))
        {
          throw std::runtime_error("Failed to add new atoms to graph.");
        }
        // append attributes
        ids.insert(std::end(ids), std::begin(newIds), std::end(newIds));
        types.insert(std::end(types), std::begin(newTypes), std::end(newTypes));
        x.insert(std::end(x), std::begin(newX), std::end(newX));
        y.insert(std::end(y), std::begin(newY), std::end(newY));
        z.insert(std::end(z), std::begin(newZ), std::end(newZ));
        nx.insert(std::end(nx), std::begin(newNx), std::end(newNx));
        ny.insert(std::end(ny), std::begin(newNy), std::end(newNy));
        nz.insert(std::end(nz), std::begin(newNz), std::end(newNz));
        // do map for easy access afterwards
        for (int i = 0; i < NNewAtoms; ++i)
        {
          this->atomIdToVectorIdx[ids[i]] = this->NAtoms + i;
        }
        this->NAtoms += NNewAtoms;
      }

      void addBonds(const int NNewBonds, std::vector<int> from, std::vector<int> to)
      {
        if (from.size() != to.size() || from.size() != NNewBonds)
        {
          throw std::invalid_argument("All inputs must have the same size");
        }
        std::vector<long int> newEdgesVector(from.size() + to.size());
        pylimer_tools::utils::interleave(from.begin(), from.end(), to.begin(), to.end(), newEdgesVector);
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

      std::vector<Molecule> getMolecules(const int atomTypeToOmit)
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
        igraph_vs_t verticesToRemove = this->getVerticesOfType(atomTypeToOmit);
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
          igraph_t *g = (igraph_t *)VECTOR(components)[i];
          molecules.push_back(Molecule(g, MoleculeType::UNDEFINED));
        }
        igraph_decompose_destroy(&components);
        igraph_destroy(&graphWithoutCrosslinkers);
        return molecules;
      }

      igraph_vs_t getVerticesOfType(const int type)
      {
        std::vector<long int> indices;
        for (int i = 0; i < this->NAtoms; ++i)
        {
          if (this->types[i] == type)
          {
            indices.push_back(i);
          }
        }
        igraph_vector_t indicesToSelect;
        pylimer_tools::utils::StdVectorToIgraphVectorT(indices, &indicesToSelect);
        return igraph_vss_vector(&indicesToSelect);
      }

      std::vector<Molecule> getChainsWithCrosslinker(const int crosslinkerType)
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
            igraph_vit_create(&graph, endNodes, &endNodeVit);
            // loop end nodes
            while (!IGRAPH_VIT_END(endNodeVit))
            {
              long int vertexId = (long int)IGRAPH_VIT_GET(endNodeVit);
              igraph_vs_t neighbours;
              if (igraph_vs_adj(&neighbours, vertexId, IGRAPH_ALL))
              {
                throw std::runtime_error("Failed to get neighbours in graph");
              }

              igraph_vit_t neighbourVit;
              igraph_vit_create(&graph, neighbours, &neighbourVit);

              // loop neighbours
              while (!IGRAPH_VIT_END(neighbourVit))
              {
                long int neighbourId = (long int)IGRAPH_VIT_GET(neighbourVit);
                if (this->types[neighbourId] == crosslinkerType) {
                  // check if this crosslinker exists in the current chain to find loops
                  igraph_vs_t neighbourInChain;
                  if (igraph_vs_1(&neighbourInChain, neighbourId)) {
                    isLoop = true;
                  } else {
                    // check whether the neighbour has been found
                    if (IGRAPH_VIT_SIZE(oneNeighbourVit) == 1) {
                      isLoop = true;
                    } else {
                      
                    }
                  }
                  igraph_vs_destroy(&neighbourInChain);
                }

                IGRAPH_VIT_NEXT(neighbourVit);
              }

              igraph_vit_destroy(&neighbourVit);
              igraph_vs_destroy(&neighbours);
              IGRAPH_VIT_NEXT(endNodeVit);
            }
            igraph_vit_destroy(&endNodeVit);
          }
          igraph_vs_destroy(&endNodes);
          // TODO: add crosslinkers back to chain/molecule

          molecules.push_back(Molecule(chain, molType));
        }
        igraph_decompose_destroy(&components);
        igraph_destroy(&graphWithoutCrosslinkers);
      }

      std::map<int, int> determineFunctionalityPerType()
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

        std::vector<int> uniqueTypes;
        copy(std::begin(this->types), std::end(this->types), std::back_inserter(uniqueTypes));
        auto uniqueTypesIter = std::unique(this->types.begin(), this->types.end());
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
          result[this->types[vertexId]] = std::max((int)igraph_vector_e(&degrees, vertexId), result[this->types[vertexId]]);
          IGRAPH_VIT_NEXT(vit);
        }
        igraph_vit_destroy(&vit);
        igraph_vs_destroy(&allVertexIds);

        return result;
      }

      Atom getAtom(const int atomId)
      {
        int atomIdx = this->atomIdToVectorIdx[atomId];
        return this->getAtomByIdx(atomIdx);
      }

      Atom getAtomByIdx(const int vertexIdx)
      {
        return Atom(this->ids[vertexIdx], this->types[vertexIdx], this->x[vertexIdx], this->y[vertexIdx], this->z[vertexIdx]);
      }

      std::vector<Atom> getAtomsWithType(const int atomType)
      {
        std::vector<Atom> atoms;
        for (int i = 0; i < this->NAtoms; ++i)
        {
          if (this->types[i] == atomType)
          {
            atoms.push_back(this->getAtomByIdx(i));
          }
        }
        return atoms;
      }

      double getVolume()
      {
        return this->Lx * this->Ly * this->Lz;
      }

      int getNrOfAtoms()
      {
        return this->NAtoms;
      }

      // properties of the box
      double Lx, Ly, Lz;
      int NAtoms = 0;
      int NBonds = 0;
      // connectivity
      igraph_t graph;
      std::map<int, int> atomIdToVectorIdx;
      // properties of the atoms
      std::vector<long int> ids;
      std::vector<int> types;
      std::vector<double> x, y, z;
      std::vector<double> nx, ny, nz;
    };
  }
}
