#include "Universe.h"
#include "Box.h"
#include "../utils/GraphUtils.h"
#include "../utils/StringUtil.h"
#include "../utils/VectorUtils.h"
extern "C"
{
#include <igraph/igraph.h>
}
#include <cassert>

#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <iterator> // for back_inserter
#include <algorithm>
#include <unordered_set>

namespace pylimer_tools
{
  namespace entities
  {

    Universe::Universe(const double Lx, const double Ly, const double Lz)
    {
      /* turn on attribute handling: TODO: move to some main() function  */
      igraph_set_attribute_table(&igraph_cattribute_table);
      box = Box(Lx, Ly, Lz);

      // igraph_vector_t gtypes, vtypes, etypes;
      // igraph_strvector_t gnames, vnames, enames;

      // igraph_vector_init(&gtypes, 0);
      // igraph_vector_init(&vtypes, 0);
      // igraph_vector_init(&etypes, 0);
      // igraph_strvector_init(&gnames, 0);
      // igraph_strvector_init(&vnames, 0);
      // igraph_strvector_init(&enames, 0);

      // start setting properties
      igraph_empty(&this->graph, 0, IGRAPH_UNDIRECTED);

      //
      // igraph_cattribute_list(&this->graph, &gnames, &gtypes, &vnames, &vtypes,
      //                        &enames, &etypes);

      // // not sure if the above is really needed as we can destroy the vectors here already without problems
      // /* Destroy */
      // igraph_vector_destroy(&gtypes);
      // igraph_vector_destroy(&vtypes);
      // igraph_vector_destroy(&etypes);
      // igraph_strvector_destroy(&gnames);
      // igraph_strvector_destroy(&vnames);
      // igraph_strvector_destroy(&enames);
    }

    // 1. destructor (to destroy the graph)
    Universe::~Universe()
    {
      // in addition to basic fields being deleted, we need to clean up the graph
      igraph_destroy(&this->graph);
    };

    // 2. copy constructor
    Universe::Universe(const Universe &src)
    {
      igraph_copy(&this->graph, &src.graph);
      this->timestep = src.timestep;
      this->NAtoms = src.NAtoms;
      this->NBonds = src.NBonds;
      // using copy assignement operators ourselfes
      this->box = src.box;
      this->atomIdToVectorIdx = src.atomIdToVectorIdx;
      this->weightPerType = src.weightPerType;
    };

    // 3. copy assignment operator
    Universe &Universe::operator=(Universe src)
    {
      std::swap(this->timestep, src.timestep);
      std::swap(this->NAtoms, src.NAtoms);
      std::swap(this->NBonds, src.NBonds);
      std::swap(this->box, src.box);
      std::swap(this->graph, src.graph);
      std::swap(this->atomIdToVectorIdx, src.atomIdToVectorIdx);
      std::swap(this->weightPerType, src.weightPerType);

      return *this;
    };

    // other functions

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
      this->atomIdToVectorIdx.reserve(this->NAtoms + NNewAtoms);
      // do map for easy access afterwards
      for (size_t i = 0; i < NNewAtoms; ++i)
      {
        if (this->atomIdToVectorIdx.contains(newIds[i]))
        {
          throw std::invalid_argument("Atom with id " + std::to_string(newIds[i]) + " already exists");
        }
        this->atomIdToVectorIdx.emplace(newIds[i], this->NAtoms + i);
      }
      // append attributes
      // it is empirically more efficient to do it this split up way
      if (this->NAtoms == 0)
      {
        // NOTE: using the same vector over an over might be bad for performance?
        igraph_vector_t valueVec;
        igraph_vector_init(&valueVec, NNewAtoms);
        pylimer_tools::utils::StdVectorToIgraphVectorT(newIds, &valueVec);
        igraph_cattribute_VAN_setv(&this->graph, "id", &valueVec);
        pylimer_tools::utils::StdVectorToIgraphVectorT(newX, &valueVec);
        igraph_cattribute_VAN_setv(&this->graph, "x", &valueVec);
        pylimer_tools::utils::StdVectorToIgraphVectorT(newY, &valueVec);
        igraph_cattribute_VAN_setv(&this->graph, "y", &valueVec);
        pylimer_tools::utils::StdVectorToIgraphVectorT(newZ, &valueVec);
        igraph_cattribute_VAN_setv(&this->graph, "z", &valueVec);
        pylimer_tools::utils::StdVectorToIgraphVectorT(newTypes, &valueVec);
        igraph_cattribute_VAN_setv(&this->graph, "type", &valueVec);
        pylimer_tools::utils::StdVectorToIgraphVectorT(newNx, &valueVec);
        igraph_cattribute_VAN_setv(&this->graph, "nx", &valueVec);
        pylimer_tools::utils::StdVectorToIgraphVectorT(newNy, &valueVec);
        igraph_cattribute_VAN_setv(&this->graph, "ny", &valueVec);
        pylimer_tools::utils::StdVectorToIgraphVectorT(newNz, &valueVec);
        igraph_cattribute_VAN_setv(&this->graph, "nz", &valueVec);
        igraph_vector_destroy(&valueVec);
      }
      else
      {
        for (size_t i = 0; i < NNewAtoms; ++i)
        {
          igraph_cattribute_VAN_set(&this->graph, "id", this->NAtoms + i, newIds[i]);
          igraph_cattribute_VAN_set(&this->graph, "x", this->NAtoms + i, newX[i]);
          igraph_cattribute_VAN_set(&this->graph, "y", this->NAtoms + i, newY[i]);
          igraph_cattribute_VAN_set(&this->graph, "z", this->NAtoms + i, newZ[i]);
          igraph_cattribute_VAN_set(&this->graph, "type", this->NAtoms + i, newTypes[i]);
          igraph_cattribute_VAN_set(&this->graph, "nx", this->NAtoms + i, newNx[i]);
          igraph_cattribute_VAN_set(&this->graph, "ny", this->NAtoms + i, newNy[i]);
          igraph_cattribute_VAN_set(&this->graph, "nz", this->NAtoms + i, newNz[i]);
        }
      }
      // this->NAtoms += NNewAtoms;
      this->NAtoms = igraph_vcount(&this->graph);
    }

    void Universe::addBonds(const size_t NNewBonds, std::vector<long int> from, std::vector<long int> to)
    {
      this->addBonds(NNewBonds, from, to, std::vector<int>());
    }

    void Universe::addBonds(const size_t NNewBonds, std::vector<long int> from, std::vector<long int> to, std::vector<int> bondTypes)
    {

      this->addBonds(NNewBonds, from, to, bondTypes, false);
    }

    void Universe::addBonds(const size_t NNewBonds, std::vector<long int> from, std::vector<long int> to, std::vector<int> bondTypes, const bool ignoreNonExistentAtoms)
    {
      if (from.size() != to.size() || from.size() != NNewBonds)
      {
        throw std::invalid_argument("All bond inputs must have the same size.");
      }
      std::vector<long int> newEdgesVector = pylimer_tools::utils::interleave(from, to);
      size_t edgesSize = newEdgesVector.size();
      // translate from atomId to VertexIdx
      igraph_vector_t newEdges;
      size_t actualNrOfBondsAdded = 0;
      igraph_vector_init(&newEdges, edgesSize);
      for (size_t i = 0; i < edgesSize; ++i)
      {
        if (this->atomIdToVectorIdx.contains(newEdgesVector[i]))
        {
          igraph_vector_set(&newEdges, i, this->atomIdToVectorIdx.at(newEdgesVector[i]));
          actualNrOfBondsAdded += 1;
        }
        else if (!ignoreNonExistentAtoms)
        {
          throw std::invalid_argument("Bond with atom with id " + std::to_string(newEdgesVector[i]) + " impossible as atom is not added yet.");
        }
      }
      igraph_vector_resize(&newEdges, actualNrOfBondsAdded);
      // add the new edges
      if (igraph_add_edges(&this->graph, &newEdges, 0))
      {
        throw std::runtime_error("Failed to add edges to graph.");
      }
      igraph_vector_destroy(&newEdges);
      if (actualNrOfBondsAdded > 0)
      {
        // add attributes
        // if (bondTypes.size() == NNewBonds && this->NBonds == igraph_ecount(&this->graph) - NNewBonds)
        // {
        //   for (size_t i = 0; i < NNewBonds; ++i)
        //   {
        //     // append attributes
        //     igraph_cattribute_EAN_set(&this->graph, "type", this->NBonds + i, bondTypes[i]);
        //   }
        // }
        // else: too risky to add bond attributes
        // simplify graph
        igraph_attribute_combination_t comb;
        igraph_attribute_combination_init(&comb);
        igraph_simplify(&this->graph, /*multiple=*/1, /*loops=*/1, &comb);
        igraph_attribute_combination_destroy(&comb);
        // this->NBonds += NNewBonds;
        this->NBonds = igraph_ecount(&this->graph);
      }
    };

    void Universe::setMasses(std::map<int, double> weightPerType)
    {
      this->weightPerType = weightPerType;
    }

    std::map<int, double> Universe::getMasses()
    {
      return this->weightPerType;
    };

    /**
     * @brief Get the standalone components of the network
     * 
     * @return std::vector<Molecule> 
     */
    std::vector<Molecule> Universe::getClusters()
    {
      std::vector<Molecule> molecules;
      if (this->getNrOfAtoms() == 0)
      {
        return molecules;
      }

      // split the copy into the separate components
      igraph_vector_ptr_t components;
      igraph_vector_ptr_init(&components, 0);
      if (igraph_decompose(&graph, &components, IGRAPH_WEAK, -1, 0))
      {
        throw std::runtime_error("Failed to decompose graph.");
      }
      size_t NComponents = igraph_vector_ptr_size(&components);
      // std::cout << NComponents << " clusters found." << std::endl;
      molecules.reserve(NComponents);
      for (size_t i = 0; i < NComponents; ++i)
      {
        // make the molecule the owner of the graph
        igraph_t *g = (igraph_t *)VECTOR(components)[i];

        if (igraph_vcount(g))
        {
          molecules.push_back(Molecule(&this->box, g, MoleculeType::UNDEFINED, this->weightPerType));
        }
      }
      igraph_decompose_destroy(&components);
      igraph_vector_ptr_destroy(&components);
      return molecules;
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
        // make the molecule the owner of the graph
        igraph_t *g = (igraph_t *)VECTOR(components)[i];

        if (igraph_vcount(g))
        {
          molecules.push_back(Molecule(&this->box, g, MoleculeType::UNDEFINED, this->weightPerType));
        }
      }
      igraph_decompose_destroy(&components);
      igraph_vector_ptr_destroy(&components);
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
      igraph_vector_destroy(&indicesToSelect);
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
      igraph_vector_destroy(&types);
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
        // also select ones of degree 0 for dangling atoms
        std::vector<long int> endNodeIndices = pylimer_tools::utils::getVerticesWithDegree(chain, {{0, 1}});
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
            long int oldEndNodeId = (long int)VAN(chain, "id", newEndNodeVertexId);
            long int originalEndNodeVertexId = this->findVertexIdForProperty("id", oldEndNodeId);
            igraph_vector_t neighbours;
            igraph_vector_init(&neighbours, 0);

            if (igraph_neighbors(&graph, &neighbours, originalEndNodeVertexId, IGRAPH_ALL))
            {
              throw std::runtime_error("Failed to get neighbours in graph");
            }

            std::vector<long int> neighborsVec;
            pylimer_tools::utils::igraphVectorTToStdVector(&neighbours, neighborsVec);

            // loop neighbours
            for (long int originalNeighbourId : neighborsVec)
            {
              int originalNeighbourType = igraph_cattribute_VAN(&graph, "type", originalNeighbourId);

              if (originalNeighbourType == crosslinkerType)
              {
                // found a crosslinker neighbour
                long int originalNeighbourAtomId = igraph_cattribute_VAN(&graph, "id", originalNeighbourId);
                atomsToAdd.push_back(originalNeighbourId);
                bondsToAdd.push_back({{newEndNodeVertexId, originalNeighbourId}});
              }
            }

            if (atomsToAdd.size() == 2 && atomsToAdd[0] == atomsToAdd[1])
            {
              isLoop = true;
              // we only want to add it once -> remove
              atomsToAdd.pop_back();
            }

            IGRAPH_VIT_NEXT(endNodeVit);
            igraph_vector_destroy(&neighbours);
          } // loop end nodes

          std::unordered_map<long int, long int> newAtomsMap;
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
            igraph_add_edge(chain, bond[0], newAtomsMap[bond[1]]);
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
        molecules.push_back(Molecule(&this->box, chain, molType, this->weightPerType));
      }
      igraph_decompose_destroy(&components);
      igraph_vector_ptr_destroy(&components);
      igraph_destroy(&graphWithoutCrosslinkers);

      return molecules;
    }

    std::map<int, std::vector<std::vector<Atom>>> Universe::findLoops(const int crosslinkerType, const int maxLength)
    {
      // NOTE: there are exponentially many paths between two vertices of a graph,
      // and you may run out of memory when using this function, if your graph is lattice-like.
      std::map<int, std::vector<std::vector<Atom>>> results;

      std::vector<long int> startingCrosslinkers = this->getIndicesOfType(crosslinkerType);
      std::unordered_set<int> processedPathsKeys;

      // note: this algorithm is not particularly efficient
      // it is of the order of O(n*n!)
      for (long int startingCrosslinkerVertexId : startingCrosslinkers)
      {
        // ideally, we would only select the neighbouring *crosslinkers* here to reduce the overhead.
        // but well.
        igraph_vector_t neighbours;
        igraph_vector_init(&neighbours, 0);

        if (igraph_neighbors(&graph, &neighbours, startingCrosslinkerVertexId, IGRAPH_ALL))
        {
          throw std::runtime_error("Failed to get neighbours in graph");
        }

        // loop neighbours
        igraph_vector_int_t paths;
        igraph_vector_int_init(&paths, 0);
        // for each neighbour, we search the simple paths
        if (igraph_get_all_simple_paths(&this->graph, &paths, startingCrosslinkerVertexId, igraph_vss_vector(&neighbours), maxLength, IGRAPH_ALL))
        {
          throw std::runtime_error("Failed to get simple paths in graph");
        }

        igraph_vector_destroy(&neighbours);
        // translate the paths we found
        std::vector<Atom> currentPath;
        int currentFunctionality = 0;
        long int currentPathKey = 0;
        size_t n = igraph_vector_int_size(&paths);
        for (int i = 0; i < n; ++i)
        {
          const long int currentVal = igraph_vector_int_e(&paths, i);
          if (currentVal == -1)
          {
            // skip self-loops and duplicates
            if (currentPath.size() > 3 && !processedPathsKeys.contains(currentPathKey))
            {
              results[currentFunctionality].push_back(currentPath);
              processedPathsKeys.insert(currentPathKey);
            }
            currentPath.clear();
            currentFunctionality = 0;
            currentPathKey = 0;
          }
          else
          {
            Atom newAtom = this->getAtomByIdx(currentVal);
            currentPathKey = currentPathKey xor currentVal; // compute hash
            currentPath.push_back(newAtom);
            if (newAtom.getType() == crosslinkerType)
            {
              currentFunctionality += 1;
            }
          }
        }
        igraph_vector_int_destroy(&paths);
      }

      return results;
    };

    bool Universe::hasInfiniteStrand(const int crosslinkerType, const int maxLength)
    {
      // NOTE: there are exponentially many paths between two vertices of a graph,
      // and you may run out of memory when using this function, if your graph is lattice-like.
      std::map<int, std::vector<std::vector<Atom>>> results;

      std::vector<long int> startingCrosslinkers = this->getIndicesOfType(crosslinkerType);
      std::unordered_set<int> processedPathsKeys;

      // note: this algorithm is not particularly efficient
      // it is of the order of O(n*n!)
      for (long int startingCrosslinkerVertexId : startingCrosslinkers)
      {
        // select all neighbouring atoms as possible directions for the loop
        igraph_vector_t neighbours;
        igraph_vector_init(&neighbours, 0);

        if (igraph_neighbors(&graph, &neighbours, startingCrosslinkerVertexId, IGRAPH_ALL))
        {
          throw std::runtime_error("Failed to get neighbours in graph");
        }

        // loop neighbours
        igraph_vector_int_t paths;
        igraph_vector_int_init(&paths, 0);
        // for each neighbour, we search the simple paths
        if (igraph_get_all_simple_paths(&this->graph, &paths, startingCrosslinkerVertexId, igraph_vss_vector(&neighbours), maxLength, IGRAPH_ALL))
        {
          throw std::runtime_error("Failed to get simple paths in graph");
        }

        igraph_vector_destroy(&neighbours);
        // translate the paths we found
        std::vector<Atom> currentPath;
        int nrOfTraversalsX = 0;
        int nrOfTraversalsY = 0;
        int nrOfTraversalsZ = 0;
        size_t n = igraph_vector_int_size(&paths);
        for (int i = 0; i < n; ++i)
        {
          const long int currentVal = igraph_vector_int_e(&paths, i);
          if (currentVal == -1)
          {
            // finished a loop. Check.
            // we have an infinite loop if the box boundary was passed in one direction only
            // NOTE: this neglects infinite networks (≠ infinite loops) such as ones caused by 
            // entanglement between images
            if (nrOfTraversalsX != 0 || nrOfTraversalsY != 0 || nrOfTraversalsZ != 0)
            {
              igraph_vector_int_destroy(&paths);
              return true;
            }
            // then reset
            currentPath.clear();
            nrOfTraversalsX = 0;
            nrOfTraversalsY = 0;
            nrOfTraversalsZ = 0;
          }
          else
          {
            Atom newAtom = this->getAtomByIdx(currentVal);
            if (!currentPath.empty())
            {
              Atom lastAtom = currentPath.back();
              double dx = newAtom.getX() - lastAtom.getX();
              nrOfTraversalsX += (dx) > 0.5 * (this->box.getLx()) ? 1 : (dx < -0.5 * (this->box.getLx()) ? -1 : 0);
              double dy = newAtom.getY() - lastAtom.getY();
              nrOfTraversalsY += (dx) > 0.5 * (this->box.getLy()) ? 1 : (dy < -0.5 * (this->box.getLy()) ? -1 : 0);
              double dz = newAtom.getZ() - lastAtom.getZ();
              nrOfTraversalsZ += (dz) > 0.5 * (this->box.getLz()) ? 1 : (dz < -0.5 * (this->box.getLz()) ? -1 : 0);
            }
            currentPath.push_back(newAtom);
          }
        }
        igraph_vector_int_destroy(&paths);
      }

      return false;
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
      igraph_vector_destroy(&degrees);

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

    Atom Universe::getAtom(const int atomId) const
    {
      if (!this->atomIdToVectorIdx.contains(atomId))
      {
        throw std::invalid_argument("Atom with this id (" + std::to_string(atomId) + ") does not exist");
      }
      int atomIdx = this->atomIdToVectorIdx.at(atomId);
      return this->getAtomByIdx(atomIdx);
    }

    Atom Universe::getAtomByIdx(const int vertexIdx) const
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

    std::vector<Atom> Universe::getAtoms()
    {
      std::vector<Atom> atoms;
      atoms.reserve(this->getNrOfAtoms());
      igraph_vit_t vit;
      igraph_vit_create(&graph, igraph_vss_all(), &vit);
      while (!IGRAPH_VIT_END(vit))
      {
        long int vertexId = (long int)IGRAPH_VIT_GET(vit);
        atoms.push_back(this->getAtomByIdx(vertexId));
        IGRAPH_VIT_NEXT(vit);
      }
      igraph_vit_destroy(&vit);
      return atoms;
    }

    std::map<std::string, std::vector<long int>> Universe::getBonds()
    {
      igraph_vector_t allEdges;
      igraph_vector_init(&allEdges, this->getNrOfBonds());
      if (igraph_edges(&this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &allEdges))
      {
        throw std::runtime_error("Failed to get all edges");
      }

      std::vector<long int> from;
      from.reserve(this->getNrOfBonds());
      std::vector<long int> to;
      to.reserve(this->getNrOfBonds());
      std::vector<long int> type;
      type.reserve(this->getNrOfBonds());

      for (long int i = 0; i < igraph_vector_size(&allEdges); i++)
      {
        if (i % 2 == 0)
        {
          from.push_back(igraph_vector_e(&allEdges, i));
        }
        else
        {
          to.push_back(igraph_vector_e(&allEdges, i));
        }
      }

      igraph_vector_destroy(&allEdges);

      // if (igraph_cattribute_has_attr(&this->graph, IGRAPH_ATTRIBUTE_EDGE, "type"))
      // {
      //   igraph_vector_t typesVec;
      //   igraph_vector_init(&typesVec, 0);
      //   igraph_cattribute_EANV(&this->graph, "type", igraph_ess_all(IGRAPH_EDGEORDER_ID), &typesVec);
      //   pylimer_tools::utils::igraphVectorTToStdVector(&typesVec, type);
      //   igraph_vector_destroy(&typesVec);
      // }
      // else
      {
        for (int i = 0; i < this->NBonds; ++i)
        {
          type.push_back(-1); // TODO: find a nice default
        }
      }

      std::map<std::string, std::vector<long int>> results;
      results.insert_or_assign("bond_from", from);
      results.insert_or_assign("bond_to", to);
      results.insert_or_assign("bond_type", type);

      return results;
    };

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
          igraph_vector_destroy(&allValues);
          return i;
        }
      }
      igraph_vector_destroy(&allValues);
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

    std::vector<double> Universe::computeDxs(const std::vector<int> bondFrom, const std::vector<int> bondTo)
    {
      return this->computeDs(bondFrom, bondTo, "x", this->box.getLx());
    };

    std::vector<double> Universe::computeDys(const std::vector<int> bondFrom, const std::vector<int> bondTo)
    {
      return this->computeDs(bondFrom, bondTo, "y", this->box.getLy());
    };

    std::vector<double> Universe::computeDzs(const std::vector<int> bondFrom, const std::vector<int> bondTo)
    {
      return this->computeDs(bondFrom, bondTo, "z", this->box.getLz());
    };

    std::vector<double> Universe::computeDs(const std::vector<int> bondFrom, const std::vector<int> bondTo, std::string direction, double boxLimit)
    {
      if (bondFrom.size() != bondTo.size())
      {
        throw std::invalid_argument("bond from and bond to must have the same size");
      }

      int nBonds = bondFrom.size();

      igraph_vector_t vertexIdFrom;
      igraph_vector_init(&vertexIdFrom, nBonds);
      igraph_vector_t vertexIdTo;
      igraph_vector_init(&vertexIdTo, nBonds);

      for (int i = 0; i < nBonds; ++i)
      {
        igraph_vector_set(&vertexIdFrom, i, this->atomIdToVectorIdx[bondFrom[i]]);
        igraph_vector_set(&vertexIdTo, i, this->atomIdToVectorIdx[bondTo[i]]);
      }

      igraph_vector_t dValuesFrom;
      igraph_vector_init(&dValuesFrom, nBonds);
      igraph_vector_t dValuesTo;
      igraph_vector_init(&dValuesTo, nBonds);

      std::string property = direction;
      igraph_cattribute_VANV(&this->graph, property.c_str(), igraph_vss_vector(&vertexIdFrom), &dValuesFrom);
      igraph_cattribute_VANV(&this->graph, property.c_str(), igraph_vss_vector(&vertexIdTo), &dValuesTo);

      igraph_vector_destroy(&vertexIdFrom);
      igraph_vector_destroy(&vertexIdTo);

      std::vector<double> results;
      results.reserve(nBonds);

      for (int i = 0; i < nBonds; ++i)
      {
        double currentD = (double)igraph_vector_e(&dValuesTo, i) - (double)igraph_vector_e(&dValuesFrom, i);
        while (std::fabs(currentD) > 0.5 * boxLimit)
        {
          if (currentD < 0.0)
          {
            currentD += boxLimit;
          }
          else
          {
            currentD -= boxLimit;
          }
        }
        results.push_back(currentD);
      }

      igraph_vector_destroy(&dValuesFrom);
      igraph_vector_destroy(&dValuesTo);

      return results;
    };

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
      int result = VECTOR(results)[0];
      igraph_vector_destroy(&results);
      return result;
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

    const int Universe::getNrOfAtoms() const
    {
      return this->NAtoms;
    }

    const int Universe::getNrOfBonds() const
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
