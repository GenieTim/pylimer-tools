#include "Universe.h"
#include "../utils/VectorUtils.h"
#include "../utils/GraphUtils.h"
#include "Box.h"
extern "C" {
#include <igraph/igraph.h>
}

#include <vector>
#include <set>
#include <map>
#include <iterator> // for back_inserter
#include <algorithm>

namespace pylimer_tools
{
  namespace entities
  {

    Universe::Universe(const double Lx, const double Ly, const double Lz)
    {
      igraph_set_attribute_table(&igraph_cattribute_table);
      box = Box(Lx, Ly, Lz);

      /* turn on attribute handling: TODO: move to some main() function  */
      // igraph_set_attribute_table(&igraph_cattribute_table);
      igraph_vector_t gtypes, vtypes, etypes;
      igraph_strvector_t gnames, vnames, enames;

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

    void Universe::addAtoms(const size_t NNewAtoms, std::vector<long int> newIds, std::vector<int> newTypes, std::vector<double> newX, std::vector<double> newY, std::vector<double> newZ, std::vector<int> newNx, std::vector<int> newNy, std::vector<int> newNz)
    {
      if (newTypes.size() != NNewAtoms || newIds.size() != newTypes.size() || newX.size() != newNx.size() || newY.size() != newNy.size() || newZ.size() != newNz.size() || newX.size() != newY.size() || NNewAtoms != newZ.size())
      {
        throw std::invalid_argument("All atom inputs must have the same size.");
      }
      // actually add the vertices
      if (igraph_add_vertices(&this->graph, NNewAtoms, 0))
      {
        throw std::runtime_error("Failed to add new atoms to graph.");
      }
      // do map for easy access afterwards
      for (size_t i = 0; i < NNewAtoms; ++i)
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
      // this->NAtoms += NNewAtoms;
      this->NAtoms = igraph_vcount(&this->graph);
    }

    void Universe::addBonds(const size_t NNewBonds, std::vector<long int> from, std::vector<long int> to)
    {
      if (from.size() != to.size() || from.size() != NNewBonds)
      {
        throw std::invalid_argument("All bond inputs must have the same size.");
      }
      std::vector<long int> newEdgesVector = pylimer_tools::utils::interleave(from, to);
      // translate from atomId to VertexIdx
      size_t edgesSize = newEdgesVector.size();
      for (size_t i = 0; i < edgesSize; ++i)
      {
        newEdgesVector[i] = this->atomIdToVectorIdx[newEdgesVector[i]];
      }
      // add the new edges
      igraph_vector_t newEdges;
      igraph_vector_init(&newEdges, edgesSize);
      pylimer_tools::utils::StdVectorToIgraphVectorT(newEdgesVector, &newEdges);
      if (igraph_add_edges(&this->graph, &newEdges, 0))
      {
        throw std::runtime_error("Failed to add edges to graph.");
      }
      igraph_attribute_combination_t comb;
      igraph_attribute_combination_init(&comb);
      igraph_simplify(&this->graph, /*multiple=*/1, /*loops=*/1, &comb);
      igraph_attribute_combination_destroy(&comb);
      // this->NBonds += NNewBonds;
      this->NBonds = igraph_ecount(&this->graph);
    }

    void Universe::setMasses(std::map<int, double> weightPerType)
    {
      this->weightPerType = weightPerType;
    }

    std::map<int, double> Universe::getMasses()
    {
      return this->weightPerType;
    };

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
      std::sort(indicesToRemove.rbegin(), indicesToRemove.rend());
      if (indicesToRemove.size() > 0)
      {
        igraph_vs_t verticesToRemove = this->getVerticesByIndices(indicesToRemove);

        // remove elements of type
        if (igraph_delete_vertices(&graphWithoutCrosslinkers, verticesToRemove))
        {
          throw std::runtime_error("Failed to delete crosslinkers from graph.");
        }

        igraph_vs_destroy(&verticesToRemove);
      }

      // split the copy into the separate components
      igraph_vector_ptr_t components;
      igraph_vector_ptr_init(&components, this->getNrOfAtoms());
      if (igraph_decompose(&graphWithoutCrosslinkers, &components, IGRAPH_WEAK, -1, 0))
      {
        throw std::runtime_error("Failed to decompose graph.");
      }
      size_t NComponents = igraph_vector_ptr_size(&components);
      // std::cout << NComponents << " molecules found. Removed " << indicesToRemove.size()
      //           << " vertices. Size now: " << igraph_vcount(&graphWithoutCrosslinkers) << " atoms with " << igraph_ecount(&graphWithoutCrosslinkers) << " bonds." << std::endl;
      molecules.reserve(NComponents);
      for (size_t i = 0; i < NComponents; ++i)
      {
        igraph_t *g = (igraph_t *)VECTOR(components)[i];

        if (igraph_vcount(g))
        {
          molecules.push_back(Molecule(&this->box, g, MoleculeType::UNDEFINED));
        }
      }
      // igraph_decompose_destroy(&components);
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
      igraph_vs_t result;
      if (igraph_vs_vector_copy(&result, &indicesToSelect))
      {
        throw std::runtime_error("Failed to select vertices");
      }
      return result;
    }

    std::vector<long int> Universe::getIndicesOfType(const int type)
    {
      std::vector<long int> indices;
      if (this->getNrOfAtoms() == 0)
      {
        return indices;
      }

      igraph_vector_t types;
      igraph_vector_init(&types, this->getNrOfAtoms());
      VANV(&this->graph, "type", &types);
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
      std::vector<long int> indicesToRemove = this->getIndicesOfType(crosslinkerType);
      std::sort(indicesToRemove.rbegin(), indicesToRemove.rend());
      if (indicesToRemove.size() > 0)
      {
        igraph_vs_t verticesToRemove = this->getVerticesByIndices(indicesToRemove);

        // remove elements of type
        if (igraph_delete_vertices(&graphWithoutCrosslinkers, verticesToRemove))
        {
          throw std::runtime_error("Failed to delete crosslinkers from graph.");
        }

        igraph_vs_destroy(&verticesToRemove);
      }

      // split the copy into the separate
      igraph_vector_ptr_t components;
      igraph_vector_ptr_init(&components, 3);
      if (igraph_decompose(&graphWithoutCrosslinkers, &components, IGRAPH_STRONG, -1, 0))
      {
        throw std::runtime_error("Failed to decompose graph.");
      }
      size_t NComponents = igraph_vector_ptr_size(&components);
      molecules.reserve(NComponents);
      for (size_t i = 0; i < NComponents; ++i)
      {
        // loop the chains to add the crosslinkers back
        igraph_t *chain = (igraph_t *)VECTOR(components)[i];
        int moleculeLengthBefore = igraph_vcount(chain);
        std::vector<long int> endNodeIndices = pylimer_tools::utils::getVerticesWithDegree(chain, 1);
        igraph_vector_t endNodeSelectorVector;
        igraph_vector_init(&endNodeSelectorVector, endNodeIndices.size());
        pylimer_tools::utils::StdVectorToIgraphVectorT(endNodeIndices, &endNodeSelectorVector);
        MoleculeType molType = MoleculeType::UNDEFINED;
        bool isLoop = false;

        if (moleculeLengthBefore >= 1)
        {
          // no care about single atoms for now (though they are included)
          igraph_vit_t endNodeVit;
          igraph_vit_create(chain, igraph_vss_vector(&endNodeSelectorVector), &endNodeVit);
          // collect atoms to add, since adding them would invalidate the iterator
          std::vector<long int> atomsToAdd;
          std::vector<std::vector<long int>> bondsToAdd;
          // loop end nodes
          while (!IGRAPH_VIT_END(endNodeVit))
          {
            long int newEndNodeVertexId = (long int)IGRAPH_VIT_GET(endNodeVit);
            long int originalEndNodeVertexId = this->findVertexIdForProperty("id", VAN(&graph, "id", newEndNodeVertexId));
            igraph_vs_t neighbours;
            if (igraph_vs_adj(&neighbours, originalEndNodeVertexId, IGRAPH_ALL))
            {
              throw std::runtime_error("Failed to get neighbours in graph");
            }

            igraph_vit_t originalNeighbourVit;
            igraph_vit_create(&graph, neighbours, &originalNeighbourVit);

            // loop neighbours
            while (!IGRAPH_VIT_END(originalNeighbourVit))
            {
              long int originalNeighbourId = (long int)IGRAPH_VIT_GET(originalNeighbourVit);

              if (igraph_cattribute_VAN(&graph, "type", originalNeighbourId) == crosslinkerType)
              {
                // found a crosslinker neighbour
                long int originalNeighbourAtomId = igraph_cattribute_VAN(&graph, "id", originalNeighbourId);
                atomsToAdd.push_back(originalNeighbourId);
                bondsToAdd.push_back({{newEndNodeVertexId, originalNeighbourId}});
              }

              IGRAPH_VIT_NEXT(originalNeighbourVit);
            }

            if (atomsToAdd.size() == 2 && atomsToAdd[0] == atomsToAdd[1])
            {
              isLoop = true;
              // we only want to add it once -> remove
              atomsToAdd.pop_back();
            }

            IGRAPH_VIT_NEXT(endNodeVit);
            igraph_vit_destroy(&originalNeighbourVit);
            igraph_vs_destroy(&neighbours);
          } // loop end nodes

          std::map<long int, long int> newAtomsMap;
          // actually add the atoms...
          for (auto atomToAddOriginalId : atomsToAdd)
          {
            igraph_add_vertices(chain, 1, 0);
            long int newCrosslinkerVertexIdx = igraph_vcount(chain) - 1;
            newAtomsMap.insert_or_assign(atomToAddOriginalId, newCrosslinkerVertexIdx);
            // additional loop check
            long int originalNeighbourAtomId = igraph_cattribute_VAN(&graph, "id", atomToAddOriginalId);
            if (pylimer_tools::utils::graphHasVertexWithProperty(chain, "id", originalNeighbourAtomId))
            {
              isLoop = true;
            }

            // including all attributes
            for (auto property : {"id", "type", "x", "y", "z", "nx", "ny", "nz"})
            {
              SETVAN(chain, property, newCrosslinkerVertexIdx, VAN(&graph, property, atomToAddOriginalId));
            }
          }
          // ...and bonds
          for (auto bond : bondsToAdd)
          {
            igraph_add_edge(chain, bond[0], newAtomsMap[bond[2]]);
          }
          igraph_vit_destroy(&endNodeVit);
        } // if molecule length
        igraph_vector_destroy(&endNodeSelectorVector);
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
      // igraph_decompose_destroy(&components);
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

      std::vector<int> types = this->getPropertyValues<int>("type");
      std::set<int> uniqueTypes(types.begin(), types.end());
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

    /*
    Compute the weight fractions of each atom type in the network.

    Returns:
      - $\\vec{W_i}$ (dict): using the type i as a key, this dict contains the weight fractions ($\\frac{W_i}{W_{tot}}$)
    */
    std::map<int, double> Universe::computeWeightFractions()
    {
      std::map<int, double> partialMasses;
      if (this->getNrOfAtoms() == 0)
      {
        return partialMasses;
      }

      std::vector<int> types = this->getPropertyValues<int>("type");
      double totalMass = 0.0;
      for (int type : types)
      {
        totalMass += this->weightPerType[type];
        partialMasses.try_emplace(type, 0.0);
        partialMasses[type] += weightPerType[type];
      }

      if (totalMass == 0.0)
      {
        return partialMasses;
      }

      //  loop to turn partial masses into weight fractions
      for (const auto &partialMassPair : partialMasses)
      {
        partialMasses[partialMassPair.first] = partialMassPair.second / totalMass;
      }

      return partialMasses;
    }

    Atom Universe::getAtom(const int atomId)
    {
      if (!this->atomIdToVectorIdx.contains(atomId))
      {
        throw std::invalid_argument("Atom with this id (" + std::to_string(atomId) + ") does not exist");
      }
      int atomIdx = this->atomIdToVectorIdx[atomId];
      return this->getAtomByIdx(atomIdx);
    }

    Atom Universe::getAtomByIdx(const int vertexIdx)
    {
      if (vertexIdx > this->getNrOfAtoms())
      {
        throw std::invalid_argument("Atom with this vertex id (" + std::to_string(vertexIdx) + ") does not exist");
      }
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

    template <typename OUT>
    std::vector<OUT> Universe::getPropertyValues(const char *propertyName)
    {
      igraph_vector_t allValues;
      igraph_vector_init(&allValues, this->getNrOfAtoms());
      VANV(&this->graph, propertyName, &allValues);
      std::vector<OUT> results;
      pylimer_tools::utils::igraphVectorTToStdVector(&allValues, results);
      igraph_vector_destroy(&allValues);
      return results;
    }

    double Universe::getMeanStrandLength(int junctionType)
    {
      std::vector<Molecule> molecules = this->getMolecules(junctionType);

      double multiplier = 1.0 / molecules.size();
      double meanStrandLength = 0;

      for (Molecule molecule : molecules)
      {
        meanStrandLength += molecule.getLength() * multiplier;
      }
      return meanStrandLength;
    }

    int Universe::getNrOfBondsOfAtom(const long int atomId)
    {
      return this->getNrOfBondsOfVertex(this->atomIdToVectorIdx[atomId]);
    }

    int Universe::getNrOfBondsOfVertex(const long int vertexId)
    {
      igraph_vector_t results;
      igraph_vector_init(&results, 1);
      if (igraph_degree(&this->graph, &results, igraph_vss_1(vertexId), IGRAPH_ALL, false /** don't count loops */))
      {
        throw std::runtime_error("Failed to query degree.");
      }
      return VECTOR(results)[0];
    }

    bool Universe::validate()
    {
      if (this->getNrOfAtoms() != igraph_vcount(&this->graph))
      {
        throw std::runtime_error("Validation failed: " + std::to_string(this->getNrOfAtoms()) + " atoms for " + std::to_string(igraph_vcount(&this->graph)) + " vertices.");
      }
      if (this->getNrOfBonds() != igraph_ecount(&this->graph))
      {
        throw std::runtime_error("Validation failed: " + std::to_string(this->getNrOfBonds()) + " bonds for " + std::to_string(igraph_ecount(&this->graph)) + " edges.");
      }
      return true;
    }

    double Universe::getVolume()
    {
      return this->box.getVolume();
    }

    const int Universe::getNrOfAtoms()
    {
      return this->NAtoms;
    }

    const int Universe::getNrOfBonds()
    {
      return this->NBonds;
    }

    void Universe::setBox(Box box)
    {
      this->box = box;
    }

    void Universe::setBoxLengths(const double Lx, const double Ly, const double Lz)
    {
      this->setBox(Box(Lx, Ly, Lz));
    }

    Box Universe::getBox() { return this->box; }
  }
}
