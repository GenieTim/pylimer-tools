extern "C"
{
#include <igraph/igraph.h>
}

#include "../utils/GraphUtils.h"
#include "../utils/StringUtils.h"
#include "Atom.h"
#include "AtomGraphParent.h"
#include <Eigen/Dense>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <vector>

namespace pylimer_tools::entities {
AtomGraphParent::AtomGraphParent() {}

// rule of three:
// 1. destructor (to destroy the graph)
AtomGraphParent::~AtomGraphParent()
{
  // in addition to basic fields being deleted, we need to clean up the
  // graph
  igraph_destroy(&this->graph);
}

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
std::vector<igraph_integer_t>
AtomGraphParent::getVertexIdxsConnectedTo(
  const igraph_integer_t vertexIdx) const
{
  std::vector<igraph_integer_t> incidentEdges =
    this->getIncidentEdgeIds(vertexIdx);
  std::vector<igraph_integer_t> results;
  results.reserve(incidentEdges.size());
  for (long int edgeId : incidentEdges) {
    igraph_integer_t vertex1OfEdge;
    igraph_integer_t vertex2OfEdge;
    igraph_edge(&this->graph, edgeId, &vertex1OfEdge, &vertex2OfEdge);

    if (vertex1OfEdge == vertexIdx) {
      results.push_back(vertex2OfEdge);
    } else {
      assert(vertex2OfEdge == vertexIdx);
      results.push_back(vertex1OfEdge);
    }
  }

  return results;
};

/**
 * @brief Get the Atoms Connected To an Atom
 *
 * @param atom the atom for which to get the atoms connected to it
 * @return std::vector<Atom>
 */
std::vector<Atom>
AtomGraphParent::getConnectedAtoms(const Atom& atom) const
{
  return this->getAtomsConnectedTo(this->getIdxByAtomId(atom.getId()));
};

/**
 * @brief Get the Atoms Connected To an Atom specified by its vertex Id
 *
 * @param vertexIdx the index of the vertex in the graph for which to get
 * the connected atoms
 * @return std::vector<Atom>
 */
std::vector<Atom>
AtomGraphParent::getAtomsConnectedTo(const long int vertexIdx) const
{
  std::vector<Atom> results;
  std::vector<igraph_integer_t> vertexIds =
    this->getVertexIdxsConnectedTo(vertexIdx);
  results.reserve(vertexIds.size());
  std::ranges::transform(vertexIds,
                         std::back_inserter(results),
                         [this](const long int vertexId) -> Atom {
                           return this->getAtomByVertexIdx(vertexId);
                         });
  return results;
};

/**
 * @brief Get the shortest sequence of atoms between two vertices
 *
 * @param vertexIdxFrom
 * @param vertexIdxTo
 * @return std::vector<Atom>
 */
std::vector<Atom>
AtomGraphParent::getShortestPath(const long int vertexIdxFrom,
                                 const long int vertexIdxTo) const
{
  std::vector<Atom> result;

  igraph_vector_int_t vertices;
  igraph_vector_int_init(&vertices, 0);

  igraph_get_shortest_path(
    &this->graph, &vertices, nullptr, vertexIdxFrom, vertexIdxTo, IGRAPH_ALL);

  result.reserve(igraph_vector_int_size(&vertices));
  for (igraph_integer_t i = 0; i < igraph_vector_int_size(&vertices); i++) {
    result.push_back(
      this->getAtomByVertexIdx(igraph_vector_int_get(&vertices, i)));
  }

  igraph_vector_int_destroy(&vertices);
  return result;
};

/**
 * @brief Get the number of edges leading to/from one vertex
 *
 * @param vertexIdxFrom the vertex from which to start the path
 * @param vertexIdxTo the vertex to which to end the path
 * @param maxLength the maximum length of the path, negative for no limit
 * @return the number of edges in the path, 0 if no path exists
 */
int
AtomGraphParent::getPathLength(const long int vertexIdxFrom,
                               const long int vertexIdxTo,
                               const int maxLength) const
{
  INVALIDINDEX_EXP_IFN(vertexIdxFrom >= 0 &&
                         vertexIdxFrom < this->getNrOfVertices(),
                       "Invalid vertex index from.");
  INVALIDINDEX_EXP_IFN(vertexIdxTo >= 0 &&
                         vertexIdxTo < this->getNrOfVertices(),
                       "Invalid vertex index target.");

  if (maxLength == 0) {
    return 0;
  }

  int foundAtDepth = 0;

  auto handler_lambda =
    [&foundAtDepth, vertexIdxTo, maxLength](const igraph_integer_t vid,
                                            const igraph_integer_t dist) {
      if (dist > maxLength && maxLength > 0) {
        return IGRAPH_STOP;
      }
      if (vid == vertexIdxTo) {
        foundAtDepth = dist;
        return IGRAPH_STOP;
      }
      return IGRAPH_SUCCESS;
    };

  igraph_bfs(
    &this->graph,
    vertexIdxFrom,
    nullptr,
    IGRAPH_ALL,
    false,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    [](const igraph_t* igraph,
       const igraph_integer_t vid,
       igraph_integer_t pred,
       igraph_integer_t succ,
       igraph_integer_t rank,
       const igraph_integer_t dist,
       void* f) {
      return (*static_cast<decltype(handler_lambda)*>(f))(vid, dist);
    },
    &handler_lambda);

  return foundAtDepth;
};

/**
 * @brief Get the edge ids of the edges between two vertices
 *
 * Main useage: check whether two vertices are connected twice
 *
 * @param vertexId1 id of the first vertex
 * @param vertexId2 id of the second vertex
 * @return std::vector<long int> the edge ids
 */
std::vector<igraph_integer_t>
AtomGraphParent::getEdgeIdsFromTo(const igraph_integer_t vertexId1,
                                  const igraph_integer_t vertexId2) const
{
  std::vector<igraph_integer_t> incidentEdges =
    this->getIncidentEdgeIds(vertexId1);
  std::vector<igraph_integer_t> results;
  results.reserve(incidentEdges.size());

  for (const long int edgeId : incidentEdges) {
    igraph_integer_t vertex1OfEdge;
    igraph_integer_t vertex2OfEdge;
    igraph_edge(&this->graph, edgeId, &vertex1OfEdge, &vertex2OfEdge);

    if ((vertex1OfEdge == vertexId1 && vertex2OfEdge == vertexId2) ||
        (vertex1OfEdge == vertexId2 && vertex2OfEdge == vertexId1)) {
      results.push_back(edgeId);
    }
  }

  return results;
}

std::vector<igraph_integer_t>
AtomGraphParent::getIncidentEdgeIds(const igraph_integer_t vertexId) const
{
  igraph_es_t edgeSelector;
  igraph_es_incident(&edgeSelector, vertexId, IGRAPH_ALL);
  igraph_eit_t iterator;
  igraph_eit_create(&this->graph, edgeSelector, &iterator);
  std::vector<igraph_integer_t> results;
  results.reserve(IGRAPH_EIT_SIZE(iterator));
  while (!IGRAPH_EIT_END(iterator)) {
    igraph_integer_t edgeId = IGRAPH_EIT_GET(iterator);
    results.push_back(edgeId);
    IGRAPH_EIT_NEXT(iterator);
  }
  igraph_eit_destroy(&iterator);
  igraph_es_destroy(&edgeSelector);

  return results;
};

/**
 * @brief Get the number Of Atoms
 *
 * @return int
 */
igraph_integer_t
AtomGraphParent::getNrOfVertices() const
{
  return igraph_vcount(&this->graph);
}

/**
 * @brief Get the Nr Of Bonds
 *
 * @return int
 */
int
AtomGraphParent::getNrOfEdges() const
{
  return igraph_ecount(&this->graph);
}

/**
 * @brief Get all atoms of a certain type
 *
 * @param atomType the type to query for
 * @return std::vector<Atom>
 */
std::vector<Atom>
AtomGraphParent::getAtomsOfType(const int atomType) const
{
  std::vector<Atom> results;
  std::vector<igraph_integer_t> indices =
    this->getIndicesWithAttribute<int>("type", atomType);

  for (igraph_integer_t i : indices) {
    results.push_back(this->getAtomByVertexIdx(i));
  }

  return results;
};

/**
 * @brief Get an atom by its vertex id
 *
 * @param vertexIdx the id of the vertex on the graph
 * @return Atom
 */
Atom
AtomGraphParent::getAtomByVertexIdx(const igraph_integer_t vertexIdx) const
{
  if (this->atomsHaveCustomAttributes) {
    return this->getComplexAtomByVertexIdx(vertexIdx);
  }
  return this->getSimpleAtomByVertexIdx(vertexIdx);
};

Atom
AtomGraphParent::getSimpleAtomByVertexIdx(
  const igraph_integer_t vertexIdx) const
{
  INVALIDINDEX_EXP_IFN(0 <= vertexIdx && vertexIdx < this->getNrOfVertices(),
                       "Atom with this vertex id (" +
                         std::to_string(vertexIdx) + ") does not exist");

  return Atom(
    static_cast<long int>(
      std::lround(igraph_cattribute_VAN(&this->graph, "id", vertexIdx))),
    static_cast<int>(
      std::rint(igraph_cattribute_VAN(&this->graph, "type", vertexIdx))),
    igraph_cattribute_VAN(&this->graph, "x", vertexIdx),
    igraph_cattribute_VAN(&this->graph, "y", vertexIdx),
    igraph_cattribute_VAN(&this->graph, "z", vertexIdx),
    static_cast<int>(
      std::rint(igraph_cattribute_VAN(&this->graph, "nx", vertexIdx))),
    static_cast<int>(
      std::rint(igraph_cattribute_VAN(&this->graph, "ny", vertexIdx))),
    static_cast<int>(
      std::rint(igraph_cattribute_VAN(&this->graph, "nz", vertexIdx))));
};

Atom
AtomGraphParent::getComplexAtomByVertexIdx(const long int vertexIdx) const
{
  if (vertexIdx > this->getNrOfVertices()) {
    throw std::invalid_argument("Atom with this vertex id (" +
                                std::to_string(vertexIdx) + ") does not exist");
  }

  igraph_strvector_t vnames;
  igraph_strvector_init(&vnames, 12);
  igraph_cattribute_list(
    &this->graph, nullptr, nullptr, &vnames, nullptr, nullptr, nullptr);

  // fetch all atom properties
  std::unordered_map<std::string, double> atomProperties;
  atomProperties.reserve(igraph_strvector_size(&vnames));
  for (igraph_integer_t i = 0; i < igraph_strvector_size(&vnames); ++i) {
    atomProperties[igraph_strvector_get(&vnames, i)] = igraph_cattribute_VAN(
      &this->graph, igraph_strvector_get(&vnames, i), vertexIdx);
  }

  igraph_strvector_destroy(&vnames);

  return Atom(atomProperties);
};

Eigen::Vector3d
AtomGraphParent::getPositionVectorForVertex(const int vertexId) const
{
  Eigen::Vector3d vertex = Eigen::Vector3d::Zero();
  vertex[0] = igraph_cattribute_VAN(&this->graph, "x", vertexId);
  vertex[1] = igraph_cattribute_VAN(&this->graph, "y", vertexId);
  vertex[2] = igraph_cattribute_VAN(&this->graph, "z", vertexId);
  return vertex;
}

Eigen::Vector3d
AtomGraphParent::getUnwrappedPositionVectorForVertex(
  const int vertexId,
  const pylimer_tools::entities::Box& box) const
{
  Eigen::Vector3d vertex = Eigen::Vector3d::Zero();
  vertex[0] = igraph_cattribute_VAN(&this->graph, "x", vertexId);
  vertex[1] = igraph_cattribute_VAN(&this->graph, "y", vertexId);
  vertex[2] = igraph_cattribute_VAN(&this->graph, "z", vertexId);

  Eigen::Array3d image = Eigen::Array3d::Zero();
  image[0] = igraph_cattribute_VAN(&this->graph, "nx", vertexId);
  image[1] = igraph_cattribute_VAN(&this->graph, "ny", vertexId);
  image[2] = igraph_cattribute_VAN(&this->graph, "nz", vertexId);

  return vertex + (box.getL() * image).matrix();
}

std::vector<Atom>
AtomGraphParent::verticesToAtoms(
  const std::vector<igraph_integer_t>& vertexIds) const
{
  std::vector<Atom> results;
  results.reserve(vertexIds.size());
  std::transform(vertexIds.begin(),
                 vertexIds.end(),
                 std::back_inserter(results),
                 [&](const long int vertexId) {
                   return this->getAtomByVertexIdx(vertexId);
                 });
  return results;
};

/**
 * @brief Check whether a vertex property exists
 *
 * @param propertyName the property to check
 * @return true if the attribute has been set at any point in time
 * @return false or not
 */
bool
AtomGraphParent::vertexPropertyExists(const char* propertyName) const
{
  return igraph_cattribute_has_attr(
    &this->graph, IGRAPH_ATTRIBUTE_VERTEX, propertyName);
};

/**
 * @brief Get the unwrapped coordinates for a given set of vertices
 *
 * @param selector
 * @param box
 * @return Eigen::VectorXd
 */
Eigen::VectorXd
AtomGraphParent::getUnwrappedVertexCoordinates(
  const igraph_vs_t selector,
  const pylimer_tools::entities::Box& box) const
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

  if (igraph_cattribute_VANV(&this->graph, "x", selector, &xvalues)) {
    throw std::runtime_error("Failed to query property x of graph.");
  }
  if (igraph_cattribute_VANV(&this->graph, "y", selector, &yvalues)) {
    throw std::runtime_error("Failed to query property y of graph.");
  }
  if (igraph_cattribute_VANV(&this->graph, "z", selector, &zvalues)) {
    throw std::runtime_error("Failed to query property z of graph.");
  }
  if (igraph_cattribute_VANV(&this->graph, "nx", selector, &nxvalues)) {
    throw std::runtime_error("Failed to query property nx of graph.");
  }
  if (igraph_cattribute_VANV(&this->graph, "ny", selector, &nyvalues)) {
    throw std::runtime_error("Failed to query property ny of graph.");
  }
  if (igraph_cattribute_VANV(&this->graph, "nz", selector, &nzvalues)) {
    throw std::runtime_error("Failed to query property nz of graph.");
  }

  const size_t size = igraph_vector_size(&xvalues);
  Eigen::VectorXd results = Eigen::VectorXd::Zero(size * 3);
  for (size_t i = 0; i < size; i++) {
    results[3 * i] = igraph_vector_get(&xvalues, i) +
                     (box.getLx() * igraph_vector_get(&nxvalues, i));
    results[3 * i + 1] = igraph_vector_get(&yvalues, i) +
                         (box.getLy() * igraph_vector_get(&nyvalues, i));
    results[3 * i + 2] = igraph_vector_get(&zvalues, i) +
                         (box.getLz() * igraph_vector_get(&nzvalues, i));
  }

  igraph_vector_destroy(&xvalues);
  igraph_vector_destroy(&yvalues);
  igraph_vector_destroy(&zvalues);
  igraph_vector_destroy(&nxvalues);
  igraph_vector_destroy(&nyvalues);
  igraph_vector_destroy(&nzvalues);

  return results;
}

Eigen::VectorXd
AtomGraphParent::getUnwrappedVertexCoordinates(
  const pylimer_tools::entities::Box& box) const
{
  return this->getUnwrappedVertexCoordinates(igraph_vss_all(), box);
}

Eigen::VectorXd
AtomGraphParent::getUnwrappedVertexCoordinates(
  const igraph_vector_int_t& vertices,
  const pylimer_tools::entities::Box& box) const
{
  return this->getUnwrappedVertexCoordinates(igraph_vss_vector(&vertices), box);
}

Eigen::VectorXd
AtomGraphParent::getUnwrappedVertexCoordinates(
  const std::vector<igraph_integer_t>& vertices,
  const pylimer_tools::entities::Box& box) const
{
  igraph_vector_int_t vertices_v;
  igraph_vector_int_init(&vertices_v, vertices.size());
  pylimer_tools::utils::StdVectorToIgraphVectorT(vertices, &vertices_v);
  Eigen::VectorXd result = this->getUnwrappedVertexCoordinates(vertices_v, box);
  igraph_vector_int_destroy(&vertices_v);
  assert(result.size() == vertices.size() * 3);
  return result;
}

/**
 * @brief Get the Vertex And Edge Property Names
 *
 * @return std::pair<std::vector<std::string>, std::vector<std::string>>
 * first the vertex property names, then the same for the edges
 */
std::pair<std::vector<std::string>, std::vector<std::string>>
AtomGraphParent::getVertexAndEdgePropertyNames() const
{
  igraph_strvector_t gnames;
  igraph_strvector_init(&gnames, 1);
  igraph_vector_int_t gtypes;
  igraph_vector_int_init(&gtypes, 1);
  igraph_strvector_t vnames;
  igraph_strvector_init(&vnames, 1);
  igraph_vector_int_t vtypes;
  igraph_vector_int_init(&vtypes, 1);
  igraph_strvector_t enames;
  igraph_strvector_init(&enames, 1);
  igraph_vector_int_t etypes;
  igraph_vector_int_init(&etypes, 1);
  igraph_cattribute_list(
    &this->graph, &gnames, &gtypes, &vnames, &vtypes, &enames, &etypes);

  std::vector<std::string> vertexPropertyNames;
  pylimer_tools::utils::igraphVectorTToStdVector(&vnames, vertexPropertyNames);
  std::vector<std::string> edgePropertyNames;
  pylimer_tools::utils::igraphVectorTToStdVector(&enames, edgePropertyNames);

  igraph_strvector_destroy(&gnames);
  igraph_strvector_destroy(&vnames);
  igraph_strvector_destroy(&enames);
  igraph_vector_int_destroy(&gtypes);
  igraph_vector_int_destroy(&vtypes);
  igraph_vector_int_destroy(&etypes);
  return std::make_pair(vertexPropertyNames, edgePropertyNames);
};

/**
 * @brief Get all atoms with a certain number of bonds
 *
 * @param degree the number of bonds to search for
 * @return std::vector<Atom>
 */
std::vector<Atom>
AtomGraphParent::getAtomsOfDegree(const int degree) const
{
  std::vector<long int> endNodeIndices = this->getVerticesWithDegree(degree);
  igraph_vector_int_t endNodeSelectorVector;
  igraph_vector_int_init(&endNodeSelectorVector, endNodeIndices.size());
  pylimer_tools::utils::StdVectorToIgraphVectorT(endNodeIndices,
                                                 &endNodeSelectorVector);
  igraph_vit_t vit;
  igraph_vit_create(
    &this->graph, igraph_vss_vector(&endNodeSelectorVector), &vit);

  std::vector<Atom> results;
  results.reserve(IGRAPH_VIT_SIZE(vit));
  while (!IGRAPH_VIT_END(vit)) {
    long int vertexId1 = static_cast<long int>(IGRAPH_VIT_GET(vit));
    Atom atom = this->getAtomByVertexIdx(vertexId1);
    results.push_back(atom);
    IGRAPH_VIT_NEXT(vit);
  }

  igraph_vector_int_destroy(&endNodeSelectorVector);
  igraph_vit_destroy(&vit);
  return results;
}

/**
 * @brief compute the lengths of all bonds
 *
 * @return std::vector<double>
 */
std::vector<double>
AtomGraphParent::computeBondLengths(const Box& box) const
{
  std::vector<double> lengths;
  lengths.reserve(this->getNrOfEdges());
  if (this->getNrOfEdges() == 0) {
    return lengths;
  }
  // construct iterator
  igraph_eit_t bondIterator;
  if (igraph_eit_create(
        &this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &bondIterator)) {
    throw std::runtime_error("Cannot create iterator to loop bonds");
  }

  while (!IGRAPH_EIT_END(bondIterator)) {
    long int edgeId = static_cast<long int>(IGRAPH_EIT_GET(bondIterator));
    Eigen::Vector3d distance =
      this->getPositionVectorForVertex(IGRAPH_TO(&this->graph, edgeId)) -
      this->getPositionVectorForVertex(IGRAPH_FROM(&this->graph, edgeId));
    box.handlePBC(distance);
    lengths.push_back(distance.norm());
    IGRAPH_EIT_NEXT(bondIterator);
  }

  igraph_eit_destroy(&bondIterator);
  return lengths;
}

/**
 * @brief compute the vectors of all bonds
 *
 * @return std::vector<Eigen::Vector3d>
 */
std::vector<Eigen::Vector3d>
AtomGraphParent::computeBondVectors(const Box& box) const
{
  std::vector<Eigen::Vector3d> lengths;
  lengths.reserve(this->getNrOfEdges());
  if (this->getNrOfEdges() == 0) {
    return lengths;
  }
  // construct iterator
  igraph_eit_t bondIterator;
  if (igraph_eit_create(
        &this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &bondIterator)) {
    throw std::runtime_error("Cannot create iterator to loop bonds");
  }

  while (!IGRAPH_EIT_END(bondIterator)) {
    long int edgeId = static_cast<long int>(IGRAPH_EIT_GET(bondIterator));
    Eigen::Vector3d distance =
      this->getPositionVectorForVertex(IGRAPH_TO(&this->graph, edgeId)) -
      this->getPositionVectorForVertex(IGRAPH_FROM(&this->graph, edgeId));
    lengths.push_back(distance);
    IGRAPH_EIT_NEXT(bondIterator);
  }

  igraph_eit_destroy(&bondIterator);
  return lengths;
}

double
AtomGraphParent::computeMeanSquaredBondLength(const Box& box) const
{
  double result = 0.0;

  // construct iterator
  igraph_eit_t bondIterator;
  if (igraph_eit_create(
        &this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &bondIterator)) {
    throw std::runtime_error("Cannot create iterator to loop bonds");
  }

  while (!IGRAPH_EIT_END(bondIterator)) {
    long int edgeId = static_cast<long int>(IGRAPH_EIT_GET(bondIterator));
    Eigen::Vector3d distance =
      this->getPositionVectorForVertex(IGRAPH_TO(&this->graph, edgeId)) -
      this->getPositionVectorForVertex(IGRAPH_FROM(&this->graph, edgeId));
    box.handlePBC(distance);
    result += distance.squaredNorm() /
              static_cast<double>(IGRAPH_EIT_SIZE(bondIterator));
    IGRAPH_EIT_NEXT(bondIterator);
  }

  igraph_eit_destroy(&bondIterator);
  return result;
}

/**
 * @brief Count the number of edges leading to/from one vertex
 *
 * @param vertexId
 * @return int
 */
int
AtomGraphParent::computeFunctionalityForVertex(const long int vertexId) const
{
  igraph_integer_t degree;
  if (igraph_degree_1(&this->graph, &degree, vertexId, IGRAPH_ALL, false)) {
    throw std::runtime_error("Failed to determine degree of vertex");
  }
  return degree;
}

int
AtomGraphParent::computeFunctionalityForAtom(const long int atomId)
{
  return this->computeFunctionalityForVertex(this->getIdxByAtomId(atomId));
}

/**
 * @brief convert the atom types involved in one angle into one long number
 *
 */
long
AtomGraphParent::hashAngleType(int typeFrom,
                               const int typeVia,
                               int typeTo) const
{
  if (typeFrom < 0 || typeTo < 0 || typeVia < 0) {
    throw std::invalid_argument(
      "AtomGraphParent::hashAngleType requires positve types");
  }
  if (typeFrom > 1000 || typeTo > 1000 || typeVia > 1000) {
    throw std::invalid_argument(
      "AtomGraphParent::hashAngleType requires types < 1000");
  }
  if (typeFrom < typeTo) {
    std::swap(typeFrom, typeTo);
  }
  long hash = (typeFrom << 20) | (typeVia << 10) | typeTo;
  return hash;
  // unhashing: typeFrom = hash >> 20; typeVia = (hash >> 10) & 0x3ff; typeTo
  // = (hash & 0x3ff);
};

/**
 * @brief convert the atom types involved in one angle into one long number
 *
 */
long
AtomGraphParent::hashDihedralAngleType(int typeFrom,
                                       int typeVia1,
                                       int typeVia2,
                                       int typeTo) const
{
  if (typeFrom < 0 || typeTo < 0 || typeVia1 < 0 || typeVia2 < 0) {
    throw std::invalid_argument(
      "AtomGraphParent::hashAngleType requires positve types");
  }
  if (typeFrom > 1000 || typeTo > 1000 || typeVia1 > 1000 || typeVia2 > 1000) {
    throw std::invalid_argument(
      "AtomGraphParent::hashAngleType requires types < 1000");
  }
  if (typeFrom < typeTo) {
    std::swap(typeFrom, typeTo);
    std::swap(typeVia1, typeVia2);
  }
  long hash = (typeFrom << 30) | (typeVia1 << 20) | (typeVia2 << 10) | typeTo;
  return hash;
};

/**
 * @brief Get all edges associated with this graph
 *
 * @return std::map<std::string, std::vector<long int>>
 */
std::map<std::string, std::vector<long int>>
AtomGraphParent::getEdges() const
{
  igraph_vector_int_t allEdges;
  igraph_vector_int_init(&allEdges, this->getNrOfEdges());
  if (igraph_edges(
        &this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &allEdges)) {
    throw std::runtime_error("Failed to get all edges");
  }

  std::vector<long int> from;
  from.reserve(this->getNrOfEdges());
  std::vector<long int> to;
  to.reserve(this->getNrOfEdges());
  std::vector<long int> type;
  type.reserve(this->getNrOfEdges());

  for (long int i = 0; i < igraph_vector_int_size(&allEdges); i++) {
    if (i % 2 == 0) {
      from.push_back(igraph_vector_int_get(&allEdges, i));
    } else {
      to.push_back(igraph_vector_int_get(&allEdges, i));
    }
  }

  igraph_vector_int_destroy(&allEdges);

  if (igraph_cattribute_has_attr(&this->graph, IGRAPH_ATTRIBUTE_EDGE, "type")) {
    igraph_vector_t typesVec;
    igraph_vector_init(&typesVec, 0);
    if (igraph_cattribute_EANV(&this->graph,
                               "type",
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &typesVec)) {
      throw std::runtime_error("Failed to fetch type attribute");
    }
    assert(igraph_vector_size(&typesVec) == this->getNrOfEdges());
    for (int i = 0; i < this->getNrOfEdges(); ++i) {
      type.push_back(
        igraphRealToInt<long int>(igraph_vector_get(&typesVec, i)));
    }
    igraph_vector_destroy(&typesVec);
  } else {
    for (int i = 0; i < this->getNrOfEdges(); ++i) {
      type.push_back(-1); // TODO: find a nice default
    }
  }

  std::map<std::string, std::vector<long int>> results;
  results.insert_or_assign("edge_from", from);
  results.insert_or_assign("edge_to", to);
  results.insert_or_assign("edge_type", type);

  return results;
};

/**
 * @brief Get all bonds (edges) associated with this graph
 *
 * @return std::map<std::string, std::vector<long int>>
 */
std::map<std::string, std::vector<long int>>
AtomGraphParent::getBonds() const
{
  std::map<std::string, std::vector<long int>> vertexResults = this->getEdges();

  std::vector<long int> newFrom;
  std::vector<long int> newTo;
  newFrom.reserve(this->getNrOfEdges());
  newTo.reserve(this->getNrOfEdges());

  std::vector<long int> oldFrom = vertexResults.at("edge_from");
  assert(oldFrom.size() == this->getNrOfEdges());
  std::vector<long int> oldTo = vertexResults.at("edge_to");
  assert(oldTo.size() == this->getNrOfEdges());

  for (int i = 0; i < this->getNrOfEdges(); ++i) {
    newFrom.push_back(this->getAtomIdByIdx(oldFrom[i]));
    newTo.push_back(this->getAtomIdByIdx(oldTo[i]));
  }

  assert(newFrom.size() == this->getNrOfEdges());
  assert(newTo.size() == this->getNrOfEdges());

  std::map<std::string, std::vector<long int>> results;
  results.insert_or_assign("bond_from", newFrom);
  results.insert_or_assign("bond_to", newTo);
  results.insert_or_assign("bond_type", vertexResults.at("edge_type"));

  return results;
}

std::vector<int>
AtomGraphParent::getVertexDegrees() const
{
  int graphSize = igraph_vcount(&this->graph);
  igraph_vector_int_t degrees;
  if (igraph_vector_int_init(&degrees, graphSize)) {
    throw std::runtime_error("Failed to instantiate result vector.");
  }
  // complexity: O(|v|*d)
  if (igraph_degree(&this->graph,
                    &degrees,
                    igraph_vss_all(),
                    IGRAPH_ALL,
                    IGRAPH_LOOPS_TWICE)) {
    throw std::runtime_error("Failed to determine degree of vertices");
  }

  std::vector<int> res;
  pylimer_tools::utils::igraphVectorTToStdVector(&degrees, res);
  igraph_vector_int_destroy(&degrees);

  return res;
};

std::vector<long int>
AtomGraphParent::getVerticesWithDegree(const igraph_t* someGraph,
                                       std::function<bool(int)> selector) const
{
  int graphSize = igraph_vcount(someGraph);
  igraph_vector_int_t degrees;
  if (igraph_vector_int_init(&degrees, graphSize)) {
    throw std::runtime_error("Failed to instantiate result vector.");
  }
  igraph_vs_t allVertexIds;
  igraph_vs_all(&allVertexIds);
  // complexity: O(|v|*d)
  if (igraph_degree(someGraph, &degrees, allVertexIds, IGRAPH_ALL, false)) {
    throw std::runtime_error("Failed to determine degree of vertices");
  }

  // NOTE: this is to omit the assumption, that the returned degree is
  // sequential for vertex 0, ..., |V|
  std::vector<long int> toSelect;
  igraph_vit_t vit;
  igraph_vit_create(someGraph, allVertexIds, &vit);
  while (!IGRAPH_VIT_END(vit)) {
    igraph_integer_t vertexId = IGRAPH_VIT_GET(vit);
    int currentDegree = igraph_vector_int_get(&degrees, vertexId);
    if (selector(currentDegree)) {
      toSelect.push_back(static_cast<long int>(vertexId));
    }
    IGRAPH_VIT_NEXT(vit);
  }
  igraph_vector_int_destroy(&degrees);
  igraph_vit_destroy(&vit);
  igraph_vs_destroy(&allVertexIds);

  return toSelect;
}

std::vector<long int>
AtomGraphParent::getVerticesWithDegree(const igraph_t* someGraph,
                                       std::vector<int> degrees) const
{
  return this->getVerticesWithDegree(
    someGraph, [degrees](const int currentDegree) {
      return std::find(degrees.begin(), degrees.end(), currentDegree) !=
             degrees.end();
    });
}

std::vector<long int>
AtomGraphParent::getVerticesWithDegree(std::function<bool(int)> selector) const
{
  return this->getVerticesWithDegree(&this->graph, selector);
}

std::vector<long int>
AtomGraphParent::getVerticesWithDegree(int degree) const
{
  return this->getVerticesWithDegree(
    [degree](const int currentDegree) { return currentDegree == degree; });
}

bool
AtomGraphParent::checkIfAtomsHaveCustomAttributes() const
{
  igraph_strvector_t vnames;
  igraph_strvector_init(&vnames, 8);
  igraph_cattribute_list(
    &this->graph, nullptr, nullptr, &vnames, nullptr, nullptr, nullptr);

  // default is 8 attributes: id, type, x, y, z, nx, ny, nz
  bool result = igraph_strvector_size(&vnames) > 8;
  igraph_strvector_destroy(&vnames);
  return result;
};

void
AtomGraphParent::writeGraphToFile(const std::string& filename) const
{
  FILE* fp = fopen(filename.c_str(), "w");
  igraph_write_graph_gml(&this->graph, fp, 0, NULL, "pylimer-tools");
  fclose(fp);
};
}
