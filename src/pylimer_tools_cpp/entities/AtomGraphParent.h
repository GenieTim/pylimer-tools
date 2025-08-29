#ifndef ATOMGRAPHPARENT_H
#define ATOMGRAPHPARENT_H

extern "C"
{
#include <igraph/igraph.h>
}

#include "../utils/CerealUtils.h"
#include "../utils/VectorUtils.h"
#include "Atom.h"
#include <Eigen/Dense>
#include <climits>
#include <functional>
#include <map>
#include <vector>

namespace pylimer_tools::entities {
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
   * @brief Get a mutable copy of the underlying igraph_t structure wrapped in
   * a unique_ptr
   *
   * This creates and returns a new igraph_t structure that is a copy of the
   * internal graph. The unique_ptr will automatically handle destruction of
   * both the igraph_t structure and the allocated memory when it goes out of
   * scope.
   *
   * @return igraph_t a copy of the graph
   */
  [[nodiscard]] igraph_t getCopyOfGraph() const
  {
    igraph_t graphCopy;
    igraph_copy(&graphCopy, &this->graph);
    return graphCopy;
  }

  [[nodiscard]] std::vector<igraph_integer_t> getEdgeIdsFromTo(
    igraph_integer_t vertexId1,
    igraph_integer_t vertexId2) const;

  [[nodiscard]] std::vector<igraph_integer_t> getIncidentEdgeIds(
    igraph_integer_t vertexId) const;

  /**
   * @brief Get the vertex ids connected to a specified vertex Id
   *
   * @param vertexIdx the index of the vertex in the graph for which to get
   * the connected atoms
   * @return std::vector<long int>
   */
  [[nodiscard]] std::vector<igraph_integer_t> getVertexIdxsConnectedTo(
    igraph_integer_t vertexIdx) const;

  /**
   * @brief Get the Atoms Connected To an Atom
   *
   * @param atom the atom for which to get the atoms connected to it
   * @return std::vector<Atom>
   */
  [[nodiscard]] std::vector<Atom> getConnectedAtoms(const Atom& atom) const;

  /**
   * @brief Get the Atoms Connected To an Atom specified by its vertex Id
   *
   * @param vertexIdx the index of the vertex in the graph for which to get
   * the connected atoms
   * @return std::vector<Atom>
   */
  [[nodiscard]] std::vector<Atom> getAtomsConnectedTo(long int vertexIdx) const;

  [[nodiscard]] std::vector<Atom> getShortestPath(long int vertexIdxFrom,
                                                  long int vertexIdxTo) const;

  /**
   * @brief Get the number of edges leading to/from one vertex
   *
   * @param vertexIdxFrom the vertex from which to start the path
   * @param vertexIdxTo the vertex to which to end the path
   * @param maxLength the maximum length of the path, -1 for no limit
   * @return the number of edges in the path, 0 if no path exists
   */
  [[nodiscard]] int getPathLength(long int vertexIdxFrom,
                                  long int vertexIdxTo,
                                  int maxLength = -1) const;

  /**
   * @brief Get the number Of Atoms
   *
   * @return int
   */
  [[nodiscard]] igraph_integer_t getNrOfVertices() const;

  /**
   * @brief Get the Nr Of Bonds
   *
   * @return int
   */
  [[nodiscard]] int getNrOfEdges() const;

  /**
   * @brief convert the atom types involved in one angle into one long number
   *
   */
  [[nodiscard]] long hashAngleType(int typeFrom, int typeVia, int typeTo) const;

  /**
   * @brief convert the atom types involved in one dihedral angle into one
   * long number
   *
   */
  [[nodiscard]] long hashDihedralAngleType(int typeFrom,
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
  [[nodiscard]] unsigned long long hashVertexIndicesOrderRelevant(int r,
                                                                  int count,
                                                                  ...) const
  {
    unsigned long long hash = 0;
    int numVerticesTotal = this->getNrOfVertices() + 1;
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
      hash *= static_cast<unsigned long long>(numVerticesTotal);
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
  [[nodiscard]] unsigned long long hashVertexIndicesOrderIrrelevant(int count,
                                                                    ...) const
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

    return hash_product + hash_sum +
           (static_cast<unsigned long long>(hash_xor) << 32);
  }

  /**
   * @brief Get all atoms of a certain type
   *
   * @param atomType the type to query for
   * @return std::vector<Atom>
   */
  [[nodiscard]] std::vector<Atom> getAtomsOfType(const int atomType) const;

  /**
   * @brief Get the Atom Id By Idx object
   *
   * @param vertexId the index of the vertex
   * @return long int the atom's id
   */
  [[nodiscard]] virtual long int getAtomIdByIdx(
    igraph_integer_t vertexId) const = 0;

  /**
   * @brief Get the vertex index by the Atom id
   *
   * @param atomId the id of the atom
   * @return long int the vertex index
   */
  virtual igraph_integer_t getIdxByAtomId(long int atomId) const = 0;

  /**
   * @brief Get an atom by its vertex id
   *
   * @param vertexIdx the id of the vertex on the graph
   * @return Atom
   */
  [[nodiscard]] Atom getAtomByVertexIdx(igraph_integer_t vertexIdx) const;

  [[nodiscard]] Atom getSimpleAtomByVertexIdx(igraph_integer_t vertexIdx) const;

  [[nodiscard]] Atom getComplexAtomByVertexIdx(const long int vertexIdx) const;

  /**
   * @brief Get the (wrapped) position vector for a vertex (ignoring image
   * flags)
   *
   * @param vertexId
   * @return Eigen::Vector3d
   */
  [[nodiscard]] Eigen::Vector3d getPositionVectorForVertex(
    const int vertexId) const;

  /**
   * @brief Get the unwrapped position vector for a vertex
   *
   * @param vertexId
   * @return Eigen::Vector3d
   */
  [[nodiscard]] Eigen::Vector3d getUnwrappedPositionVectorForVertex(
    const int vertexId,
    const pylimer_tools::entities::Box& box) const;

  /**
   * @brief Convert a list of vertex ids to a list of Atom instances
   *
   * @param vertexIds the list of vertex ids
   * @return std::vector<Atom>
   */
  [[nodiscard]] std::vector<Atom> verticesToAtoms(
    const std::vector<igraph_integer_t>& vertexIds) const;

  /**
   * @brief Check whether a vertex property exists
   *
   * @param propertyName the property to check
   * @return true if the attribute has been set at any point in time
   * @return false or not
   */
  [[nodiscard]] bool vertexPropertyExists(const char* propertyName) const;

  /**
   * @brief Count how often a certain value appears in the vertex properties
   *
   * for example, to count the number of crosslinkers, one could use
   * this->countPropertyValue<int>("type", crossLinkerType)
   *
   * @tparam IN
   * @param propertyName
   * @param targetValue
   * @return int
   */
  template<typename IN>
  [[nodiscard]] int countPropertyValue(const char* propertyName,
                                       IN targetValue) const
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
  [[nodiscard]] std::vector<OUT> getPropertyValues(
    const char* propertyName) const
  {
    std::vector<OUT> results;
    if (this->getNrOfVertices() == 0) {
      return results;
    }
    igraph_vector_t allValues;
    igraph_vector_init(&allValues, this->getNrOfVertices());
    if (igraph_cattribute_VANV(
          &this->graph, propertyName, igraph_vss_all(), &allValues)) {
      throw std::runtime_error("Failed to query properties of graph.");
    }
    pylimer_tools::utils::igraphVectorTToStdVector(&allValues, results);
    igraph_vector_destroy(&allValues);
    return results;
  }

  template<typename IN>
  void setPropertyValue(const igraph_integer_t vertexId,
                        const char* propertyName,
                        IN value)
  {
    INVALIDINDEX_EXP_IFN(vertexId < this->getNrOfVertices(),
                         "Index " + std::to_string(vertexId) +
                           " is out of bounds");
    if (igraph_cattribute_VAN_set(
          &this->graph, propertyName, vertexId, value)) {
      throw std::runtime_error("Failed to set property value");
    }
    if (!this->atomsHaveCustomAttributes) {
      this->atomsHaveCustomAttributes =
        this->checkIfAtomsHaveCustomAttributes();
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
    const std::vector<igraph_integer_t>& vertices) const
  {
    std::vector<OUT> results;
    if (vertices.size() == 0) {
      return results;
    }
    igraph_vector_t allValues;
    igraph_vector_init(&allValues,
                       static_cast<igraph_integer_t>(vertices.size()));
    igraph_vector_int_t vertexIdxs;
    igraph_vector_int_init(&vertexIdxs,
                           static_cast<igraph_integer_t>(vertices.size()));
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

  /**
   * @brief Get the unwrapped coordinates for a given set of vertices
   *
   * @param selector
   * @param box
   * @return Eigen::VectorXd
   */
  Eigen::VectorXd getUnwrappedVertexCoordinates(
    const igraph_vs_t selector,
    const pylimer_tools::entities::Box& box) const;

  Eigen::VectorXd getUnwrappedVertexCoordinates(
    const pylimer_tools::entities::Box& box) const;

  Eigen::VectorXd getUnwrappedVertexCoordinates(
    const igraph_vector_int_t& vertices,
    const pylimer_tools::entities::Box& box) const;

  Eigen::VectorXd getUnwrappedVertexCoordinates(
    const std::vector<igraph_integer_t>& vertices,
    const pylimer_tools::entities::Box& box) const;

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
                       const igraph_integer_t vertexIdx) const
  {
    return igraph_cattribute_VAN(&this->graph, propertyName, vertexIdx);
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
  OUT getEdgePropertyValue(const char* propertyName,
                           const long int vertexIdx) const
  {
    return igraph_cattribute_EAN(&this->graph, propertyName, vertexIdx);
  }

  /**
   * @brief Get the Vertex And Edge Property Names
   *
   * @return std::pair<std::vector<std::string>, std::vector<std::string>>
   * first the vertex property names, then the same for the edges
   */
  std::pair<std::vector<std::string>, std::vector<std::string>>
  getVertexAndEdgePropertyNames() const;

  /**
   * @brief Get all atoms with a certain number of bonds
   *
   * @param degree the number of bonds to search for
   * @return std::vector<Atom>
   */
  std::vector<Atom> getAtomsOfDegree(const int degree) const;

  /**
   * @brief Shorthand to query the degree of a vertex
   *
   * @param vertexId
   * @param loops
   * @return igraph_integer_t
   */
  igraph_integer_t getVertexDegree(const igraph_integer_t vertexId,
                                   const bool loops = true) const
  {
    igraph_integer_t degree;
    igraph_degree_1(&this->graph, &degree, vertexId, IGRAPH_ALL, loops);
    return degree;
  }

  std::vector<int> getVertexDegrees() const;

  /**
   * @brief compute the lengths of all bonds
   *
   * @return std::vector<double>
   */
  std::vector<double> computeBondLengths(const Box& box) const;

  double computeMeanSquaredBondLength(const Box& box) const;

  std::vector<Eigen::Vector3d> computeBondVectors(const Box& box) const;

  /**
   * @brief Count the number of edges leading to/from one vertex
   *
   * @param vertexId
   * @return int
   */
  int computeFunctionalityForVertex(const long int vertexId) const;

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

  Eigen::VectorXd getAssumedVertexCoordinates(const Box& box) const
  {
    igraph_vector_int_t order;
    igraph_vector_int_init(&order, this->getNrOfVertices());
    igraph_vector_int_t parents;
    igraph_vector_int_init(&parents, this->getNrOfVertices());
    igraph_dfs(&this->graph,
               0,
               IGRAPH_ALL,
               true,
               &order,
               nullptr,
               &parents,
               nullptr,
               nullptr,
               nullptr,
               nullptr);
    assert(igraph_vector_int_size(&order) == this->getNrOfVertices());
    assert(igraph_vector_int_size(&parents) == this->getNrOfVertices());

    Eigen::VectorXd coordinates = this->getUnwrappedVertexCoordinates(box);

    // use the parents and order to reconstruct the coordinates based on the
    // bonds
    for (igraph_integer_t i = 0; i < this->getNrOfVertices(); ++i) {
      if (igraph_vector_int_get(&order, i) >= 0 &&
          igraph_vector_int_get(&parents, i) >= 0) {
        Eigen::Vector3d distance =
          coordinates.segment<3>(3 * igraph_vector_int_get(&order, i), 3) -
          coordinates.segment<3>(3 * igraph_vector_int_get(&parents, i), 3);
        box.handlePBC(distance);
        coordinates.segment<3>(3 * igraph_vector_int_get(&order, i), 3) =
          coordinates.segment<3>(3 * igraph_vector_int_get(&parents, i), 3) +
          distance;
      }
    }

    igraph_vector_int_destroy(&order);
    igraph_vector_int_destroy(&parents);
    return coordinates;
  };

  /**
   * @brief Get the Assumed Vertex Coordinates
   * This means, the coordinates are derived ignoring the image flags,
   * assuming bonds between subsequent vertices,
   * which are assumed to be shorter than half the periodic box.
   *
   * @tparam OutVectorType
   * @param results
   * @param box
   * @param vertexIds the vertices lined up
   */
  template<typename OutVectorType>
  OutVectorType getAssumedVertexCoordinates(
    OutVectorType& results,
    const Box& box,
    const std::vector<igraph_integer_t>& vertexIds) const
  {
    if (vertexIds.size() * 3 != results.size()) {
      throw std::invalid_argument(
        "The results must have size 3*the number of atoms to query.");
    }

    if (vertexIds.empty()) {
      return results;
    }

    igraph_vector_int_t vertex_ids;
    igraph_vector_int_init(&vertex_ids, 0);
    pylimer_tools::utils::StdVectorToIgraphVectorT(vertexIds, &vertex_ids);

    Eigen::VectorXd coordinates =
      this->getUnwrappedVertexCoordinates(vertex_ids, box);

    igraph_vector_int_destroy(&vertex_ids);

    // take the distances
    Eigen::VectorXd distances = coordinates.segment(3, coordinates.size() - 3) -
                                coordinates.segment(0, coordinates.size() - 3);

    // adjust them for the box
    box.handlePBC(distances);

    // and now:
    Eigen::Vector3d lastCoords = coordinates.segment(0, 3);
    // put them into the box already
    box.handlePBC(lastCoords);
    results[0] = lastCoords[0];
    results[1] = lastCoords[1];
    results[2] = lastCoords[2];
    for (int i = 0; i < distances.size(); i += 3) {
      // TODO: rather than relying on the vertex ids to be subsequently
      // connected, it would be nice if we could support graph-like structures
      // as well (check if igraph_is_tree, then, start with root and do e.g.
      // recursively)

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

  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(graph);
  }

protected:
  igraph_t graph;
  bool atomsHaveCustomAttributes = false;

  std::vector<long int> getVerticesWithDegree(int degree) const;

  std::vector<long int> getVerticesWithDegree(
    std::function<bool(int)> selector) const;

  std::vector<long int> getVerticesWithDegree(
    const igraph_t* someGraph,
    std::function<bool(int)> selector) const;

  std::vector<long int> getVerticesWithDegree(const igraph_t* someGraph,
                                              std::vector<int> degrees) const;

  template<typename ValueT>
  std::vector<igraph_integer_t> getIndicesWithAttribute(
    const std::string& attributeName,
    ValueT value) const
  {
    if (this->getNrOfVertices() == 0) {
      return std::vector<igraph_integer_t>();
    }
    RUNTIME_EXP_IFN(
      igraph_cattribute_has_attr(
        &this->graph, IGRAPH_ATTRIBUTE_VERTEX, attributeName.c_str()),
      "The graph does not have any vertex attribute '" + attributeName + "'.");
    std::vector<igraph_integer_t> indices;
    igraph_vector_t attribute_vector;
    igraph_vector_init(&attribute_vector, this->getNrOfVertices());
    VANV(&this->graph, attributeName.c_str(), &attribute_vector);
    for (igraph_integer_t i = 0; i < igraph_vector_size(&attribute_vector);
         ++i) {
      if (static_cast<ValueT>(VECTOR(attribute_vector)[i]) == value) {
        indices.push_back(i);
      }
    }
    igraph_vector_destroy(&attribute_vector);
    return indices;
  };

  [[nodiscard]] bool checkIfAtomsHaveCustomAttributes() const;
};
} // namespace pylimer_tools::entities

#endif
