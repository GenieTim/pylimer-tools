#ifndef ATOMGRAPHPARENT_H
#define ATOMGRAPHPARENT_H

extern "C" {
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
class AtomGraphParent {
public:
  AtomGraphParent() {}
  // rule of three:
  // 1. destructor (to destroy the graph)
  virtual ~AtomGraphParent() {
    // in addition to basic fields being deleted, we need to clean up the graph
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
   * @param vertexIdx the index of the vertex in the graph for which to get the
   * connected atoms
   * @return std::vector<long int>
   */
  std::vector<long int>
  getVertexIdxsConnectedTo(const long int vertexIdx) const {
    igraph_vs_t adjVs;
    if (igraph_vs_adj(&adjVs, vertexIdx, IGRAPH_ALL)) {
      throw std::runtime_error("Failed to find adjacent vertices of vertex.");
    }
    igraph_vit_t vit;
    igraph_vit_create(&this->graph, adjVs, &vit);
    std::vector<long int> results;
    results.reserve(IGRAPH_VIT_SIZE(vit));
    while (!IGRAPH_VIT_END(vit)) {
      results.push_back(static_cast<long int>(IGRAPH_VIT_GET(vit)));
      IGRAPH_VIT_NEXT(vit);
    }
    igraph_vs_destroy(&adjVs);
    igraph_vit_destroy(&vit);

    return results;
  }

  /**
   * @brief Get the Atoms Connected To an Atom specified by its vertex Id
   *
   * @param vertexIdx the index of the vertex in the graph for which to get the
   * connected atoms
   * @return std::vector<Atom>
   */
  std::vector<Atom> getAtomsConnectedTo(const long int vertexIdx) const {
    std::vector<Atom> results;
    std::vector<long int> vertexIds = this->getVertexIdxsConnectedTo(vertexIdx);
    results.reserve(vertexIds.size());
    for (long int vertexId : vertexIds) {
      results.push_back(this->getAtomByVertexIdx(vertexId));
    }
    return results;
  };

  /**
   * @brief Get the number Of Atoms
   *
   * @return int
   */
  int getNrOfAtoms() const { return igraph_vcount(&this->graph); }

  /**
   * @brief Get the Nr Of Bonds
   *
   * @return int
   */
  int getNrOfBonds() const { return igraph_ecount(&this->graph); }

  /**
   * @brief Get all atoms of a certain type
   *
   * @param atomType the type to query for
   * @return std::vector<Atom>
   */
  std::vector<Atom> getAtomsWithType(const int atomType) const {
    std::vector<Atom> results;
    const std::vector<int> types = this->getPropertyValues<int>("type");
    size_t nrOfTypes = types.size();

    // #pragma omp declare reduction (merge : std::vector<Atom> :
    // omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())) #pragma omp
    // parallel for reduction(merge: results)
    for (size_t i = 0; i < nrOfTypes; ++i) {
      if (types[i] == atomType) {
        results.push_back(this->getAtomByVertexIdx(i));
      }
    }

    return results;
  };

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
  Atom getAtomByVertexIdx(const long int vertexIdx) const {
    if (vertexIdx > this->getNrOfAtoms()) {
      throw std::invalid_argument("Atom with this vertex id (" +
                                  std::to_string(vertexIdx) +
                                  ") does not exist");
    }
    return Atom(
        VAN(&this->graph, "id", vertexIdx),
        VAN(&this->graph, "type", vertexIdx), VAN(&this->graph, "x", vertexIdx),
        VAN(&this->graph, "y", vertexIdx), VAN(&this->graph, "z", vertexIdx),
        VAN(&this->graph, "nx", vertexIdx), VAN(&this->graph, "ny", vertexIdx),
        VAN(&this->graph, "nz", vertexIdx));
  }

  /**
   * @brief Get the value of a property (attribute) of each and every vertex
   *
   * @tparam OUT
   * @param propertyName the name of the property to get
   * @return std::vector<OUT>
   */
  template <typename OUT>
  std::vector<OUT> getPropertyValues(const char *propertyName) const {
    std::vector<OUT> results;
    if (this->getNrOfAtoms() == 0) {
      return results;
    }
    igraph_vector_t allValues;
    igraph_vector_init(&allValues, this->getNrOfAtoms());
    if (igraph_cattribute_VANV(&this->graph, propertyName, igraph_vss_all(),
                               &allValues)) {
      throw std::runtime_error("Failed to query properties of molecule.");
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
  template <typename OUT>
  std::vector<OUT> getPropertyValues(const char *propertyName,
                                     std::vector<long int> vertices) const {
    std::vector<OUT> results;
    if (vertices.size() == 0) {
      return results;
    }
    igraph_vector_t allValues;
    igraph_vector_init(&allValues, vertices.size());
    igraph_vector_t vertexIdxs;
    igraph_vector_init(&vertexIdxs, vertices.size());
    pylimer_tools::utils::StdVectorToIgraphVectorT(vertices, &vertexIdxs);
    if (igraph_cattribute_VANV(&this->graph, propertyName,
                               igraph_vss_vector(&vertexIdxs), &allValues)) {
      throw std::runtime_error("Failed to query properties of molecule.");
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
  template <typename OUT>
  OUT getPropertyValue(const char *propertyName,
                       const long int vertexIdx) const {
    return igraph_cattribute_VAN(&this->graph, propertyName, vertexIdx);
  }

  /**
   * @brief Get all atoms with a certain number of bonds
   *
   * @param degree the number of bonds to search for
   * @return std::vector<Atom>
   */
  std::vector<Atom> getAtomsOfDegree(const int degree) const {
    std::vector<long int> endNodeIndices =
        pylimer_tools::utils::getVerticesWithDegree(&this->graph, degree);
    igraph_vector_t endNodeSelectorVector;
    igraph_vector_init(&endNodeSelectorVector, endNodeIndices.size());
    pylimer_tools::utils::StdVectorToIgraphVectorT(endNodeIndices,
                                                   &endNodeSelectorVector);
    igraph_vit_t vit;
    igraph_vit_create(&this->graph, igraph_vss_vector(&endNodeSelectorVector),
                      &vit);

    std::vector<Atom> results;
    results.reserve(IGRAPH_VIT_SIZE(vit));
    while (!IGRAPH_VIT_END(vit)) {
      long int vertexId1 = static_cast<long int>(IGRAPH_VIT_GET(vit));
      Atom atom = this->getAtomByVertexIdx(vertexId1);
      results.push_back(atom);
      IGRAPH_VIT_NEXT(vit);
    }

    igraph_vector_destroy(&endNodeSelectorVector);
    igraph_vit_destroy(&vit);
    return results;
  }

  /**
   * @brief compute the lengths of all bonds
   *
   * @return std::vector<double>
   */
  std::vector<double> computeBondLengths(const Box *box) {
    std::vector<double> lengths;
    lengths.reserve(this->getNrOfBonds());
    if (this->getNrOfBonds() == 0) {
      return lengths;
    }
    // construct iterator
    igraph_eit_t bondIterator;
    if (igraph_eit_create(&this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID),
                          &bondIterator)) {
      throw std::runtime_error("Cannot create iterator to loop bonds");
    }

    while (!IGRAPH_EIT_END(bondIterator)) {
      long int edgeId = static_cast<long int>(IGRAPH_EIT_GET(bondIterator));
      int bondFrom;
      int bondTo;
      igraph_edge(&this->graph, edgeId, &bondFrom, &bondTo);
      // TODO: this is more intensive than needed
      // check whether the compiler optimizes this or not
      Atom atom1 = this->getAtomByVertexIdx(bondFrom);
      Atom atom2 = this->getAtomByVertexIdx(bondTo);
      lengths.push_back(atom1.distanceTo(atom2, box));
      IGRAPH_EIT_NEXT(bondIterator);
    }

    igraph_eit_destroy(&bondIterator);
    return lengths;
  }

  /**
   * @brief Count the number of edges leading to/from one vertex
   *
   * @param vertexId
   * @return int
   */
  int computeFunctionalityForVertex(const long int vertexId) {
    igraph_vector_t degrees;
    if (igraph_vector_init(&degrees, 0)) {
      throw std::runtime_error("Failed to instantiate result vector.");
    }
    if (igraph_degree(&this->graph, &degrees, igraph_vss_1(vertexId),
                      IGRAPH_ALL, false)) {
      throw std::runtime_error("Failed to determine degree of vertex");
    }
    int result = igraph_vector_e(&degrees, 0);
    igraph_vector_destroy(&degrees);
    return result;
  }

  int computeFunctionalityForAtom(const long int atomId) {
    return this->computeFunctionalityForVertex(this->getIdxByAtomId(atomId));
  }

  /**
   * @brief Get all edges associated with this graph
   *
   * @return std::map<std::string, std::vector<long int>>
   */
  std::map<std::string, std::vector<long int>> getEdges() const {
    igraph_vector_t allEdges;
    igraph_vector_init(&allEdges, this->getNrOfBonds());
    if (igraph_edges(&this->graph, igraph_ess_all(IGRAPH_EDGEORDER_ID),
                     &allEdges)) {
      throw std::runtime_error("Failed to get all edges");
    }

    std::vector<long int> from;
    from.reserve(this->getNrOfBonds());
    std::vector<long int> to;
    to.reserve(this->getNrOfBonds());
    std::vector<long int> type;
    type.reserve(this->getNrOfBonds());

    for (long int i = 0; i < igraph_vector_size(&allEdges); i++) {
      if (i % 2 == 0) {
        from.push_back(igraph_vector_e(&allEdges, i));
      } else {
        to.push_back(igraph_vector_e(&allEdges, i));
      }
    }

    igraph_vector_destroy(&allEdges);

    // if (igraph_cattribute_has_attr(&this->graph, IGRAPH_ATTRIBUTE_EDGE,
    // "type"))
    // {
    //   igraph_vector_t typesVec;
    //   igraph_vector_init(&typesVec, 0);
    //   igraph_cattribute_EANV(&this->graph, "type",
    //   igraph_ess_all(IGRAPH_EDGEORDER_ID), &typesVec);
    //   pylimer_tools::utils::igraphVectorTToStdVector(&typesVec, type);
    //   igraph_vector_destroy(&typesVec);
    // }
    // else
    {
      for (size_t i = 0; i < this->getNrOfBonds(); ++i) {
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
  std::map<std::string, std::vector<long int>> getBonds() const {
    std::map<std::string, std::vector<long int>> vertexResults =
        this->getEdges();

    std::vector<long int> newFrom;
    std::vector<long int> newTo;
    newFrom.reserve(this->getNrOfBonds());
    newTo.reserve(this->getNrOfBonds());

    std::vector<long int> oldFrom = vertexResults.at("edge_from");
    assert(oldFrom.size() == this->getNrOfBonds());
    std::vector<long int> oldTo = vertexResults.at("edge_to");
    assert(oldTo.size() == this->getNrOfBonds());

    for (size_t i = 0; i < this->getNrOfBonds(); ++i) {
      newFrom.push_back(this->getAtomIdByIdx(oldFrom[i]));
      newTo.push_back(this->getAtomIdByIdx(oldTo[i]));
    }

    assert(newFrom.size() == this->getNrOfBonds());
    assert(newTo.size() == this->getNrOfBonds());

    std::map<std::string, std::vector<long int>> results;
    results.insert_or_assign("bond_from", newFrom);
    results.insert_or_assign("bond_to", newTo);
    results.insert_or_assign("bond_type", vertexResults.at("edge_type"));

    return results;
  };

protected:
  igraph_t graph;
};

} // namespace entities
} // namespace pylimer_tools

#endif
