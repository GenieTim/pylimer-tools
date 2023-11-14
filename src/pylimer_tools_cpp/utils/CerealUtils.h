#ifndef CEREAL_UTILS_H
#define CEREAL_UTILS_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>

#include "./ExtraEigenTypes.h"
#include "./VectorUtils.h"
#include <Eigen/Dense>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>
// #include <cereal/types/map.hpp>
// #include <cereal/types/set.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
extern "C"
{
#include <igraph/igraph.h>
}

namespace cereal {
////////////////////////////////////////////////////////////////
// serialization of Eigen objects

// if we can store binary data
template<class Archive, class Derived>
inline typename std::enable_if<
  traits::is_output_serializable<BinaryData<typename Derived::Scalar>,
                                 Archive>::value,
  void>::type
CEREAL_SAVE_FUNCTION_NAME(Archive& ar, Eigen::PlainObjectBase<Derived> const& m)
{
  typedef Eigen::PlainObjectBase<Derived> ArrT;
  if (ArrT::RowsAtCompileTime == Eigen::Dynamic)
    ar(m.rows());
  if (ArrT::ColsAtCompileTime == Eigen::Dynamic)
    ar(m.cols());
  ar(binary_data(m.data(), m.size() * sizeof(typename Derived::Scalar)));
}

template<class Archive, class Derived>
inline typename std::enable_if<
  traits::is_input_serializable<BinaryData<typename Derived::Scalar>,
                                Archive>::value,
  void>::type
CEREAL_LOAD_FUNCTION_NAME(Archive& ar, Eigen::PlainObjectBase<Derived>& m)
{
  typedef Eigen::PlainObjectBase<Derived> ArrT;
  Eigen::Index rows = ArrT::RowsAtCompileTime, cols = ArrT::ColsAtCompileTime;
  if (rows == Eigen::Dynamic)
    ar(rows);
  if (cols == Eigen::Dynamic)
    ar(cols);
  m.resize(rows, cols);
  ar(binary_data(
    m.data(),
    static_cast<std::size_t>(rows * cols * sizeof(typename Derived::Scalar))));
}

// if we cannot store binary data
template<class Archive, class Derived>
inline typename std::enable_if<
  !traits::is_output_serializable<BinaryData<typename Derived::Scalar>,
                                  Archive>::value,
  void>::type
CEREAL_SAVE_FUNCTION_NAME(Archive& ar, Eigen::PlainObjectBase<Derived> const& m)
{
  typedef Eigen::PlainObjectBase<Derived> ArrT;
  make_size_tag(m.rows());
  for (size_type i = 0; i < m.rows(); ++i) {
    make_size_tag(m.cols());
    for (size_type j = 0; j < m.cols(); ++j) {
      ar(m(i, j));
    }
  }
}

template<class Archive, class Derived>
inline typename std::enable_if<
  !traits::is_input_serializable<BinaryData<typename Derived::Scalar>,
                                 Archive>::value,
  void>::type
CEREAL_LOAD_FUNCTION_NAME(Archive& ar, Eigen::PlainObjectBase<Derived>& m)
{
  typedef Eigen::PlainObjectBase<Derived> ArrT;
  size_type rows;
  ar(make_size_tag(rows));
  size_type cols;
  ar(make_size_tag(cols));
  m.resize(rows, cols);
  for (size_type i = 0; i < rows; ++i) {
    if (i != 0) {
      ar(make_size_tag(cols));
    }
    for (size_type j = 0; j < cols; ++j) {
      ar(m(i, j));
    }
  }
}

////////////////////////////////////////////////////////////////
// serialization of igraph objects
template<class Archive>
inline void
CEREAL_SAVE_FUNCTION_NAME(Archive& ar, igraph_t const& graph)
{
  size_t numVertices = igraph_vcount(&graph);
  size_t numEdges = igraph_ecount(&graph);

  igraph_vector_int_t allEdges;
  igraph_vector_int_init(&allEdges, numEdges);
  if (igraph_edges(&graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &allEdges)) {
    throw std::runtime_error("Failed to get all edges");
  }
  std::vector<long int> edges;
  pylimer_tools::utils::igraphVectorTToStdVector(&allEdges, edges);
  igraph_vector_int_destroy(&allEdges);

  ar(numVertices);
  ar(numEdges);
  ar(edges);

  // after storing the edges, must also store the attributes
  // query them first
  igraph_strvector_t gnames;
  igraph_vector_int_t gtypes;
  igraph_strvector_t vnames;
  igraph_vector_int_t vtypes;
  igraph_strvector_t enames;
  igraph_vector_int_t etypes;
  igraph_cattribute_list(
    &graph, &gnames, &gtypes, &vnames, &vtypes, &enames, &etypes);

  if (igraph_strvector_size(&gnames) != 0) {
    throw std::runtime_error(
      "Graph attributes serialization not supported yet.");
  }

  // serizalize vertex attributes
  size_t numVertexAttributes = igraph_strvector_size(&vnames);
  ar(make_size_tag(numVertexAttributes));
  for (size_t i = 0; i < numVertexAttributes; i++) {
    const char* name = igraph_strvector_get(&vnames, i);
    ar(std::string(name));
    ar(igraph_vector_int_get(&vtypes, i));
    switch (igraph_vector_int_get(&vtypes, i)) {
      case IGRAPH_ATTRIBUTE_DEFAULT:
      case IGRAPH_ATTRIBUTE_NUMERIC: {
        igraph_vector_t results;
        igraph_vector_init(&results, numVertices);
        igraph_cattribute_VANV(
          &graph, igraph_strvector_get(&vnames, i), igraph_vss_all(), &results);
        std::vector<double> attributes;
        pylimer_tools::utils::igraphVectorTToStdVector(&results, attributes);
        ar(attributes);
        igraph_vector_destroy(&results);
      } break;
      case IGRAPH_ATTRIBUTE_STRING: {
        igraph_strvector_t strresults;
        igraph_strvector_init(&strresults, numVertices);
        igraph_cattribute_VASV(&graph,
                               igraph_strvector_get(&vnames, i),
                               igraph_vss_all(),
                               &strresults);
        std::vector<std::string> strattributes;
        pylimer_tools::utils::igraphVectorTToStdVector(&strresults,
                                                       strattributes);
        ar(strattributes);
        igraph_strvector_destroy(&strresults);
      } break;
      default:
        throw std::runtime_error(
          "This attribute type (" +
          std::to_string(igraph_vector_int_get(&vtypes, i)) +
          ") is not supported");
    }
  }

  // serizalize edge attributes
  size_t numEdgeAttributes = igraph_strvector_size(&enames);
  ar(make_size_tag(numEdgeAttributes));
  for (size_t i = 0; i < numEdgeAttributes; i++) {
    const char* name = igraph_strvector_get(&enames, i);
    ar(std::string(name));
    ar(igraph_vector_int_get(&etypes, i));
    switch (igraph_vector_int_get(&etypes, i)) {
      case IGRAPH_ATTRIBUTE_DEFAULT:
      case IGRAPH_ATTRIBUTE_NUMERIC: {
        igraph_vector_t results;
        igraph_vector_init(&results, numEdges);
        igraph_cattribute_EANV(&graph,
                               igraph_strvector_get(&enames, i),
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &results);
        std::vector<double> attributes;
        pylimer_tools::utils::igraphVectorTToStdVector(&results, attributes);
        ar(attributes);
        igraph_vector_destroy(&results);
      } break;
      case IGRAPH_ATTRIBUTE_STRING: {
        igraph_strvector_t strresults;
        igraph_strvector_init(&strresults, numEdges);
        igraph_cattribute_EASV(&graph,
                               igraph_strvector_get(&enames, i),
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &strresults);
        std::vector<std::string> strattributes;
        pylimer_tools::utils::igraphVectorTToStdVector(&strresults,
                                                       strattributes);
        ar(strattributes);
        igraph_strvector_destroy(&strresults);
      } break;
      default:
        throw std::runtime_error(
          "This attribute type (" +
          std::to_string(igraph_vector_int_get(&etypes, i)) +
          ") is not supported");
    }
  }

  igraph_strvector_destroy(&gnames);
  igraph_strvector_destroy(&enames);
  igraph_strvector_destroy(&vnames);

  igraph_vector_int_destroy(&gtypes);
  igraph_vector_int_destroy(&etypes);
  igraph_vector_int_destroy(&vtypes);
}

template<class Archive>
inline void
CEREAL_LOAD_FUNCTION_NAME(Archive& ar, igraph_t& graph)
{
  size_t numVertices;
  size_t numEdges;
  ar(numVertices);
  ar(numEdges);
  std::vector<long int> edges;
  ar(edges);
  igraph_vector_int_t allEdges;
  igraph_vector_int_init(&allEdges, numEdges);
  pylimer_tools::utils::StdVectorToIgraphVectorT(edges, &allEdges);

  igraph_add_edges(&graph, &allEdges, nullptr);

  // deserialize vertex attributes
  size_t numVertexAttributes;
  ar(make_size_tag(numVertexAttributes));
  for (size_t i = 0; i < numVertexAttributes; ++i) {
    std::string attributeName;
    ar(attributeName);
    int attributeType;
    ar(attributeType);
    switch (attributeType) {
      case IGRAPH_ATTRIBUTE_DEFAULT:
      case IGRAPH_ATTRIBUTE_NUMERIC: {
        std::vector<double> attributes;
        ar(attributes);
        igraph_vector_t results;
        igraph_vector_init(&results, 1);
        pylimer_tools::utils::StdVectorToIgraphVectorT(attributes, &results);
        igraph_cattribute_VAN_setv(&graph, attributeName.c_str(), &results);
        igraph_vector_destroy(&results);
      } break;
      case IGRAPH_ATTRIBUTE_STRING: {
        std::vector<std::string> strattributes;
        ar(strattributes);
        igraph_strvector_t strresults;
        igraph_strvector_init(&strresults, 1);
        pylimer_tools::utils::StdVectorToIgraphVectorT(strattributes,
                                                       &strresults);
        igraph_cattribute_VAS_setv(&graph, attributeName.c_str(), &strresults);
        igraph_strvector_destroy(&strresults);
      } break;
      default:
        throw std::runtime_error("This attribute type (" +
                                 std::to_string(attributeType) +
                                 ") is not supported");
    }
  }

  // and same for edge attributes
  size_t numEdgeAttributes;
  ar(make_size_tag(numEdgeAttributes));
  for (size_t i = 0; i < numEdgeAttributes; ++i) {
    std::string attributeName;
    ar(attributeName);
    int attributeType;
    ar(attributeType);
    switch (attributeType) {
      case IGRAPH_ATTRIBUTE_DEFAULT:
      case IGRAPH_ATTRIBUTE_NUMERIC: {
        std::vector<double> attributes;
        ar(attributes);
        igraph_vector_t results;
        igraph_vector_init(&results, 1);
        pylimer_tools::utils::StdVectorToIgraphVectorT(attributes, &results);
        igraph_cattribute_EAN_setv(&graph, attributeName.c_str(), &results);
        igraph_vector_destroy(&results);
      } break;
      case IGRAPH_ATTRIBUTE_STRING: {
        std::vector<std::string> strattributes;
        ar(strattributes);
        igraph_strvector_t strresults;
        igraph_strvector_init(&strresults, 1);
        pylimer_tools::utils::StdVectorToIgraphVectorT(strattributes,
                                                       &strresults);
        igraph_cattribute_EAS_setv(&graph, attributeName.c_str(), &strresults);
        igraph_strvector_destroy(&strresults);
      } break;
      default:
        throw std::runtime_error("This attribute type (" +
                                 std::to_string(attributeType) +
                                 ") is not supported");
    }
  }
}

}

namespace pylimer_tools {
namespace utils {

  template<typename T>
  void serializeToFile(T obj, std::string file)
  {
    std::ofstream os(file);
    cereal::BinaryOutputArchive oarchive(os);
    oarchive(obj);
  }

  template<typename T>
  void deserializeFromFile(T &obj, std::string file)
  {
    std::ifstream is(file);
    cereal::BinaryInputArchive iarchive(is);
    iarchive(obj);
  }
}
}

#endif
