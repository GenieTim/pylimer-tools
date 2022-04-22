#ifndef ATOMGRAPHPARENT_H
#define ATOMGRAPHPARENT_H

extern "C"
{
#include <igraph/igraph.h>
}
#include "../utils/GraphUtils.h"
#include "../utils/StringUtils.h"
#include "Atom.h"
#include <algorithm>
#include <map>
#include <unordered_map>
#include <vector>

namespace pylimer_tools {
namespace entities {
  // abstract
  class AtomGraphParent
  {
  public:
    AtomGraphParent();
    // rule of three:
    // 1. destructor (to destroy the graph)
    virtual ~AtomGraphParent();
    // 2. copy constructor
    // AtomGraphParent(const AtomGraphParent &src) {
    //   igraph_copy(&this->graph, &src.graph);
    // };
    // 3. copy assignment operator
    // virtual AtomGraphParent &operator=(AtomGraphParent src) {
    //   std::swap(this->graph, src.graph);
    //   return *this;
    // };

    /**
     * @brief Get the vertex ids connected to a specified vertex Id
     *
     * @param vertexIdx the index of the vertex in the graph for which to get
     * the connected atoms
     * @return std::vector<long int>
     */
    std::vector<long int> getVertexIdxsConnectedTo(
      const long int vertexIdx) const;

    /**
     * @brief Get the Atoms Connected To an Atom specified by its vertex Id
     *
     * @param vertexIdx the index of the vertex in the graph for which to get
     * the connected atoms
     * @return std::vector<Atom>
     */
    std::vector<Atom> getAtomsConnectedTo(const long int vertexIdx) const;

    /**
     * @brief Get the number Of Atoms
     *
     * @return int
     */
    int getNrOfAtoms() const;

    /**
     * @brief Get the Nr Of Bonds
     *
     * @return int
     */
    int getNrOfBonds() const;

    /**
     * @brief Get all atoms of a certain type
     *
     * @param atomType the type to query for
     * @return std::vector<Atom>
     */
    std::vector<Atom> getAtomsWithType(const int atomType) const;

    /**
     * @brief Get the Atom Id By Idx object
     *
     * @param vertexId the index of the vertex
     * @return long int the atom's id
     */
    virtual long int getAtomIdByIdx(const int vertexId) const = 0;

    /**
     * @brief Get the vertex index by the Atom id
     *
     * @param atomId the id of the atom
     * @return long int the vertex index
     */
    virtual long int getIdxByAtomId(const int atomId) const = 0;

    /**
     * @brief Get an atom by its vertex id
     *
     * @param vertexIdx the id of the vertex on the graph
     * @return Atom
     */
    Atom getAtomByVertexIdx(const long int vertexIdx) const;

    /**
     * @brief Get the value of a property (attribute) of each and every vertex
     *
     * @tparam OUT
     * @param propertyName the name of the property to get
     * @return std::vector<OUT>
     */
    template<typename OUT>
    std::vector<OUT> getPropertyValues(const char* propertyName) const
    {
      std::vector<OUT> results;
      if (this->getNrOfAtoms() == 0) {
        return results;
      }
      igraph_vector_t allValues;
      igraph_vector_init(&allValues, this->getNrOfAtoms());
      if (igraph_cattribute_VANV(
            &this->graph, propertyName, igraph_vss_all(), &allValues)) {
        throw std::runtime_error("Failed to query properties of graph.");
      }
      pylimer_tools::utils::igraphVectorTToStdVector(&allValues, results);
      igraph_vector_destroy(&allValues);
      return results;
    }

    /**
     * @brief Get the value of a property (attribute) of certain vertices
     *
     * @tparam OUT
     * @param propertyName the name of the property to get
     * @param vertices the list of vertices to get the property for
     * @return std::vector<OUT>
     */
    template<typename OUT>
    std::vector<OUT> getPropertyValues(const char* propertyName,
                                       std::vector<long int> vertices) const
    {
      std::vector<OUT> results;
      if (vertices.size() == 0) {
        return results;
      }
      igraph_vector_t allValues;
      igraph_vector_init(&allValues, vertices.size());
      igraph_vector_t vertexIdxs;
      igraph_vector_init(&vertexIdxs, vertices.size());
      pylimer_tools::utils::StdVectorToIgraphVectorT(vertices, &vertexIdxs);
      if (igraph_cattribute_VANV(&this->graph,
                                 propertyName,
                                 igraph_vss_vector(&vertexIdxs),
                                 &allValues)) {
        throw std::runtime_error("Failed to query properties of graph.");
      }
      pylimer_tools::utils::igraphVectorTToStdVector(&allValues, results);
      igraph_vector_destroy(&allValues);
      igraph_vector_destroy(&vertexIdxs);
      return results;
    }

    /**
     * @brief Get the Property (attribute) of one vertex
     *
     * @tparam OUT
     * @param propertyName
     * @param vertexIdx
     * @return OUT
     */
    template<typename OUT>
    OUT getPropertyValue(const char* propertyName,
                         const long int vertexIdx) const
    {
      return igraph_cattribute_VAN(&this->graph, propertyName, vertexIdx);
    }

    /**
     * @brief Get all atoms with a certain number of bonds
     *
     * @param degree the number of bonds to search for
     * @return std::vector<Atom>
     */
    std::vector<Atom> getAtomsOfDegree(const int degree) const;

    /**
     * @brief compute the lengths of all bonds
     *
     * @return std::vector<double>
     */
    std::vector<double> computeBondLengths(const Box* box);

    /**
     * @brief Count the number of edges leading to/from one vertex
     *
     * @param vertexId
     * @return int
     */
    int computeFunctionalityForVertex(const long int vertexId);

    int computeFunctionalityForAtom(const long int atomId);

    /**
     * @brief Get all edges associated with this graph
     *
     * @return std::map<std::string, std::vector<long int>>
     */
    std::map<std::string, std::vector<long int>> getEdges() const;

    /**
     * @brief Get all bonds (edges) associated with this graph
     *
     * @return std::map<std::string, std::vector<long int>>
     */
    std::map<std::string, std::vector<long int>> getBonds() const;

  protected:
    igraph_t graph;

    igraph_vs_t getVerticesWithDegreeSelector(int degree) const;
    std::vector<long int> getVerticesWithDegree(int degree) const;
    std::vector<long int> getVerticesWithDegree(std::vector<int> ofDegrees) const;
    std::vector<long int> getVerticesWithDegree(const igraph_t* someGraph,
                                                std::vector<int> ofDegrees) const;
  };

} // namespace entities
} // namespace pylimer_tools

#endif
