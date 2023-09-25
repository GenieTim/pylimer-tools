#ifndef ATOMGRAPHPARENT_H
#define ATOMGRAPHPARENT_H

extern "C"
{
#include <igraph/igraph.h>
}
#include "../utils/GraphUtils.h"
#include "../utils/StringUtils.h"
#include "Atom.h"
#include <Eigen/Dense>
#include <algorithm>
#include <limits.h>
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

    std::vector<long int> getEdgeIdsFromTo(const long int vertexId1,
                                           const long int vertexId2) const;

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

    std::vector<Atom> getShortestPath(const long int vertexIdxFrom,
                                      const long int vertexIdxTo) const;

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
     * @brief convert the atom types involved in one angle into one long number
     *
     */
    long hashAngleType(int typeFrom, int typeVia, int typeTo) const;

    /**
     * @brief convert the atom types involved in one dihedral angle into one
     * long number
     *
     */
    long hashDihedralAngleType(int typeFrom,
                               int typeVia1,
                               int typeVia2,
                               int typeTo) const;

    /**
     * @brief Combine vertex indices to a hash, that respects the order of the
     * indices
     *
     * @param r a
     * @param count the number of vertex elements to combine
     * @param ... the vertex indices
     * @return unsigned long long
     */
    template<typename vertex_idx_type>
    unsigned long long hashVertexIndicesOrderRelevant(int r,
                                                      int count,
                                                      ...) const
    {
      unsigned long long hash = 0;
      int numVerticesTotal = this->getNrOfAtoms() + 1;
      // bits used: log(pow(numVerticesTotal, count))/log(2)
      // -> possible overflow?!?
      if (std::log(std::pow(r, count)) / std::log(2) >=
          sizeof(unsigned long long) * CHAR_BIT) {
        throw std::invalid_argument(
          "With this r and count, the hash will overflow.");
      }

      va_list args;
      va_start(args, count);

      for (int i = 0; i < count; i++) {
        const vertex_idx_type nextVertexIdx = va_arg(args, vertex_idx_type);
        hash *= numVerticesTotal;
        hash += nextVertexIdx;
      }

      va_end(args);

      return hash;
    }

    /**
     * @brief Combine vertex indices to a hash, that does not respect the order
     * of the indices
     *
     * @param count the number of vertex elements to combine
     * @param ... the vertex indices
     * @return unsigned long long
     */
    template<typename vertex_idx_type>
    unsigned long long hashVertexIndicesOrderIrrelevant(int count, ...) const
    {
      unsigned long long hash_product = 1;
      unsigned long long hash_sum = 0;
      unsigned int hash_xor = 0;

      va_list args;
      va_start(args, count);

      for (int i = 0; i < count; i++) {
        const vertex_idx_type nextVertexIdx = va_arg(args, vertex_idx_type);
        hash_product *= nextVertexIdx;
        hash_sum += nextVertexIdx;
        hash_xor ^= nextVertexIdx;
      }

      va_end(args);

      return hash_product + hash_sum + ((unsigned long long)hash_xor << 32);
    }

    /**
     * @brief Get all atoms of a certain type
     *
     * @param atomType the type to query for
     * @return std::vector<Atom>
     */
    std::vector<Atom> getAtomsOfType(const int atomType) const;

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
     * @brief Convert a list of vertex ids to a list of Atom instances
     *
     * @param vertexIds the list of vertex ids
     * @return std::vector<Atom>
     */
    std::vector<Atom> verticesToAtoms(
      const std::vector<long int>& vertexIds) const;

    bool vertexPropertyExists(const char* propertyName) const
    {
      return igraph_cattribute_has_attr(
        &this->graph, IGRAPH_ATTRIBUTE_VERTEX, propertyName);
    }

    /**
     * @brief Count how often a certain value appears in the vertex properties
     *
     * For example, to count the number of cross-links, one could use
     * this->countPropertyValue<int>("type", crosslinkerType)
     *
     * @tparam IN
     * @param propertyName
     * @param targetValue
     * @return int
     */
    template<typename IN>
    int countPropertyValue(const char* propertyName, IN targetValue) const
    {
      std::vector<IN> values = this->getPropertyValues<IN>(propertyName);
      int result = 0;
      for (IN value : values) {
        result += (value == targetValue);
      }
      return result;
    }

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

    template<typename IN>
    void setPropertyValue(const long int vertexId,
                          const char* propertyName,
                          IN value)
    {
      if (igraph_cattribute_VAN_set(
            &this->graph, propertyName, vertexId, value)) {
        throw std::runtime_error("Failed to set property value");
      }
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
    std::vector<OUT> getPropertyValues(
      const char* propertyName,
      const std::vector<long int>& vertices) const
    {
      std::vector<OUT> results;
      if (vertices.size() == 0) {
        return results;
      }
      igraph_vector_t allValues;
      igraph_vector_init(&allValues, vertices.size());
      igraph_vector_int_t vertexIdxs;
      igraph_vector_int_init(&vertexIdxs, vertices.size());
      pylimer_tools::utils::StdVectorToIgraphVectorT(vertices, &vertexIdxs);
      if (igraph_cattribute_VANV(&this->graph,
                                 propertyName,
                                 igraph_vss_vector(&vertexIdxs),
                                 &allValues)) {
        throw std::runtime_error("Failed to query properties of graph.");
      }
      pylimer_tools::utils::igraphVectorTToStdVector(&allValues, results);
      igraph_vector_destroy(&allValues);
      igraph_vector_int_destroy(&vertexIdxs);
      return results;
    }

    Eigen::VectorXd getUnwrappedVertexCoordinates(
      const igraph_vs_t selector,
      const pylimer_tools::entities::Box* box) const
    {
      igraph_vector_t xvalues;
      igraph_vector_init(&xvalues, 0);
      igraph_vector_t yvalues;
      igraph_vector_init(&yvalues, 0);
      igraph_vector_t zvalues;
      igraph_vector_init(&zvalues, 0);
      igraph_vector_t nxvalues;
      igraph_vector_init(&nxvalues, 0);
      igraph_vector_t nyvalues;
      igraph_vector_init(&nyvalues, 0);
      igraph_vector_t nzvalues;
      igraph_vector_init(&nzvalues, 0);

      if (igraph_cattribute_VANV(
            &this->graph, "x", igraph_vss_all(), &xvalues)) {
        throw std::runtime_error("Failed to query property x of graph.");
      }
      if (igraph_cattribute_VANV(
            &this->graph, "y", igraph_vss_all(), &yvalues)) {
        throw std::runtime_error("Failed to query property y of graph.");
      }
      if (igraph_cattribute_VANV(
            &this->graph, "z", igraph_vss_all(), &zvalues)) {
        throw std::runtime_error("Failed to query property z of graph.");
      }
      if (igraph_cattribute_VANV(
            &this->graph, "nx", igraph_vss_all(), &nxvalues)) {
        throw std::runtime_error("Failed to query property nx of graph.");
      }
      if (igraph_cattribute_VANV(
            &this->graph, "ny", igraph_vss_all(), &nyvalues)) {
        throw std::runtime_error("Failed to query property ny of graph.");
      }
      if (igraph_cattribute_VANV(
            &this->graph, "nz", igraph_vss_all(), &nzvalues)) {
        throw std::runtime_error("Failed to query property nz of graph.");
      }

      const size_t size = igraph_vector_size(&xvalues);
      Eigen::VectorXd results = Eigen::VectorXd::Zero(size * 3);
      for (size_t i = 0; i < size; i++) {
        results[3 * i] = igraph_vector_get(&xvalues, i) +
                         (box->getLx() * igraph_vector_get(&nxvalues, i));
        results[3 * i + 1] = igraph_vector_get(&yvalues, i) +
                             (box->getLy() * igraph_vector_get(&nyvalues, i));
        results[3 * i + 2] = igraph_vector_get(&zvalues, i) +
                             (box->getLz() * igraph_vector_get(&nzvalues, i));
      }

      igraph_vector_destroy(&xvalues);
      igraph_vector_destroy(&yvalues);
      igraph_vector_destroy(&zvalues);
      igraph_vector_destroy(&nxvalues);
      igraph_vector_destroy(&nyvalues);
      igraph_vector_destroy(&nzvalues);

      return results;
    }

    Eigen::VectorXd getUnwrappedVertexCoordinates(
      const pylimer_tools::entities::Box* box) const
    {
      return this->getUnwrappedVertexCoordinates(igraph_vss_all(), box);
    }

    Eigen::VectorXd getUnwrappedVertexCoordinates(
      igraph_vector_int_t& vertices,
      const pylimer_tools::entities::Box* box) const
    {
      return this->getUnwrappedVertexCoordinates(igraph_vss_vector(&vertices),
                                                 box);
    }

    Eigen::VectorXd getUnwrappedVertexCoordinates(
      std::vector<long int>& vertices,
      const pylimer_tools::entities::Box* box) const
    {
      igraph_vector_int_t vertices_v;
      igraph_vector_int_init(&vertices_v, vertices.size());
      pylimer_tools::utils::StdVectorToIgraphVectorT(vertices, &vertices_v);
      Eigen::VectorXd result =
        this->getUnwrappedVertexCoordinates(vertices_v, box);
      igraph_vector_int_destroy(&vertices_v);
      return result;
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

    /**
     * @brief Get the Assumed Vertex Coordinates
     * This means, the coordinates are derived ignoring the image flags,
     * assuming bonds between subsequent vertices,
     * which are assumed to be shorter than half the periodic box.
     *
     * @tparam OutVectorType
     * @param results
     * @param box
     * @param vertexIds
     */

    template<typename OutVectorType>
    OutVectorType getAssumedVertexCoordinates(
      OutVectorType& results,
      const Box* box,
      const std::vector<long int>& vertexIds) const
    {
      if (vertexIds.size() * 3 != results.size()) {
        throw std::invalid_argument(
          "The results must have size 3*the number of atoms to query.");
      }

      igraph_vector_int_t vertex_ids;
      igraph_vector_int_init(&vertex_ids, 0);
      pylimer_tools::utils::StdVectorToIgraphVectorT(vertexIds, &vertex_ids);

      Eigen::VectorXd coordinates =
        this->getUnwrappedVertexCoordinates(vertex_ids, box);
      for (size_t i = 0; i < coordinates.size(); ++i) {
      }

      igraph_vector_int_destroy(&vertex_ids);

      // take the distances
      Eigen::VectorXd distances =
        coordinates.segment(3, coordinates.size() - 3) -
        coordinates.segment(0, coordinates.size() - 3);

      // adjust them for the box
      box->handlePBC(distances);

      // and now:
      Eigen::Vector3d lastCoords = coordinates.segment(0, 3);
      results[0] = lastCoords[0];
      results[1] = lastCoords[1];
      results[2] = lastCoords[2];
      for (int i = 0; i < distances.size(); i += 3) {
        // for each next atom, we can use the shortest distance to the previous
        // in order to compensate/ignore the image flags while still enabling
        // larger end-to-end distances than the box size
        lastCoords += distances.segment(i, 3);
        results[i + 3] = lastCoords[0];
        results[i + 4] = lastCoords[1];
        results[i + 5] = lastCoords[2];
      }

      return results;
    }

    void writeGraphToFile(const std::string& filename) const;

  protected:
    igraph_t graph;
    igraph_vs_t getVerticesWithDegreeSelector(int degree) const;
    std::vector<long int> getVerticesWithDegree(int degree) const;
    std::vector<long int> getVerticesWithDegree(
      std::vector<int> ofDegrees) const;
    std::vector<long int> getVerticesWithDegree(
      const igraph_t* someGraph,
      std::vector<int> ofDegrees) const;
  };

} // namespace entities
} // namespace pylimer_tools

#endif
