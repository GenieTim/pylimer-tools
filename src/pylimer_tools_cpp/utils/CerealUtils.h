#pragma once
#ifdef CEREALIZABLE

#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <type_traits>

#include "./VectorUtils.h"
#include <Eigen/Dense>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
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
// MARK: randomness serialization
template<class Archive>
void
save(Archive& archive,
     const std::uniform_real_distribution<double>& distribution)
{
  archive(distribution.a(), distribution.b());
}

template<class Archive>
void
load(Archive& archive, std::uniform_real_distribution<double>& distribution)
{
  double a, b;
  archive(a, b);
  distribution = std::uniform_real_distribution<double>(a, b);
}

template<class Archive>
void
save(Archive& archive, const std::mt19937& mt)
{
  std::ostringstream oss;
  oss << mt;
  archive(oss.str());
}

template<class Archive>
void
load(Archive& archive, std::mt19937& mt)
{
  std::string str;
  archive(str);
  std::istringstream iss(str);
  iss >> mt;
}

////////////////////////////////////////////////////////////////
// MARK: Eigen objects serialization of

// if we can store binary data
template<class Archive, class Derived>
inline typename std::enable_if<
  traits::is_output_serializable<BinaryData<typename Derived::Scalar>,
                                 Archive>::value,
  void>::type
CEREAL_SAVE_FUNCTION_NAME(Archive& ar, Eigen::PlainObjectBase<Derived> const& m)
{
  typedef Eigen::PlainObjectBase<Derived> ArrT;
  if (ArrT::RowsAtCompileTime == Eigen::Dynamic) {
    ar(m.rows());
  }
  if (ArrT::ColsAtCompileTime == Eigen::Dynamic) {
    ar(m.cols());
  }
  ar(binary_data(m.data(),
                 static_cast<std::size_t>(m.size()) *
                   sizeof(typename Derived::Scalar)));
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
  if (rows == Eigen::Dynamic) {
    ar(rows);
  }
  if (cols == Eigen::Dynamic) {
    ar(cols);
  }
  m.resize(rows, cols);
  ar(binary_data(m.data(),
                 static_cast<std::size_t>(static_cast<std::size_t>(rows) *
                                          static_cast<std::size_t>(cols) *
                                          sizeof(typename Derived::Scalar))));
}

// if we cannot store binary data
template<class Archive, class Derived>
inline typename std::enable_if<
  !traits::is_output_serializable<BinaryData<typename Derived::Scalar>,
                                  Archive>::value,
  void>::type
CEREAL_SAVE_FUNCTION_NAME(Archive& ar, Eigen::PlainObjectBase<Derived> const& m)
{
  ar(m.rows());
  ar(m.cols());
  
  for (Eigen::Index i = 0; i < m.rows(); ++i) {
    for (Eigen::Index j = 0; j < m.cols(); ++j) {
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
  Eigen::Index rows;
  ar(rows);
  Eigen::Index cols;
  ar(cols);
  m.resize(static_cast<Eigen::Index>(rows), static_cast<Eigen::Index>(cols));
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(rows); ++i) {
    for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(cols); ++j) {
      ar(m(i, j));
    }
  }
}

////////////////////////////////////////////////////////////////
// MARK: igraph objects serialization

// a vector
template<class Archive>
inline void
CEREAL_SAVE_FUNCTION_NAME(Archive& ar, igraph_vector_int_t const& vec)
{
  const igraph_integer_t n = igraph_vector_int_size(&vec);
  ar(make_size_tag(static_cast<size_type>(n)));
  for (igraph_integer_t i = 0; i < n; ++i) {
    ar(igraph_vector_int_get(&vec, i));
  }
}

template<class Archive>
inline void
CEREAL_LOAD_FUNCTION_NAME(Archive& ar, igraph_vector_int_t& vec)
{
  size_type n;
  ar(make_size_tag(n));
  igraph_vector_int_resize(&vec, static_cast<igraph_integer_t>(n));
  for (igraph_integer_t i = 0; i < static_cast<igraph_integer_t>(n); ++i) {
    long int val;
    ar(val);
    igraph_vector_int_set(&vec, i, val);
  }
}

template<class Archive>
inline void
CEREAL_SAVE_FUNCTION_NAME(Archive& ar, igraph_vector_t const& vec)
{
  const igraph_integer_t n = igraph_vector_size(&vec);
  ar(make_size_tag(static_cast<size_type>(n)));
  for (igraph_integer_t i = 0; i < n; ++i) {
    ar(igraph_vector_get(&vec, i));
  }
}

template<class Archive>
inline void
CEREAL_LOAD_FUNCTION_NAME(Archive& ar, igraph_vector_t& vec)
{
  size_type n;
  ar(make_size_tag(n));
  igraph_vector_resize(&vec, static_cast<igraph_integer_t>(n));
  for (igraph_integer_t i = 0; i < static_cast<igraph_integer_t>(n); ++i) {
    double val;
    ar(val);
    igraph_vector_set(&vec, i, val);
  }
}

template<class Archive>
inline void
CEREAL_SAVE_FUNCTION_NAME(Archive& ar, igraph_strvector_t const& vec)
{
  const igraph_integer_t n = igraph_strvector_size(&vec);
  ar(make_size_tag(static_cast<size_type>(n)));
  for (igraph_integer_t i = 0; i < n; ++i) {
    std::string val = igraph_strvector_get(&vec, i);
    ar(val);
  }
}

template<class Archive>
inline void
CEREAL_LOAD_FUNCTION_NAME(Archive& ar, igraph_strvector_t& vec)
{
  size_type n;
  ar(make_size_tag(n));
  igraph_strvector_resize(&vec, static_cast<igraph_integer_t>(n));
  std::string val;
  val.reserve(50);
  for (igraph_integer_t i = 0; i < static_cast<igraph_integer_t>(n); ++i) {
    val.clear();
    ar(val);
    igraph_strvector_set(&vec, i, val.c_str());
  }
}

// the graph
template<class Archive>
inline void
CEREAL_SAVE_FUNCTION_NAME(Archive& ar,
                          igraph_t const& graph,
                          std::uint32_t const version)
{
  const igraph_integer_t numVertices = igraph_vcount(&graph);
  ar(make_nvp("num_vertices", static_cast<size_t>(numVertices)));
  const igraph_integer_t numEdges = igraph_ecount(&graph);
  ar(make_nvp("num_edges", static_cast<size_t>(numEdges)));

  igraph_vector_int_t allEdges;
  igraph_vector_int_init(&allEdges, numEdges);
  if (igraph_edges(&graph, igraph_ess_all(IGRAPH_EDGEORDER_ID), &allEdges)) {
    throw std::runtime_error("Failed to get all edges");
  }

  ar(make_nvp("edges", allEdges));
  igraph_vector_int_destroy(&allEdges);

  // after storing the edges, must also store the attributes
  // query them first
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
    &graph, &gnames, &gtypes, &vnames, &vtypes, &enames, &etypes);

  if (igraph_strvector_size(&gnames) != 0) {
    throw std::runtime_error(
      "Graph attributes serialization not supported yet.");
  }

  // serizalize vertex attributes
  const igraph_integer_t numVertexAttributes = igraph_strvector_size(&vnames);
  assert(igraph_strvector_size(&vnames) == igraph_vector_int_size(&vtypes));
  ar(make_nvp("vertex_attr_names", vnames));
  ar(make_nvp("vertex_attr_types", vtypes));
  //   ar(make_size_tag(numVertexAttributes));
  for (igraph_integer_t i = 0; i < numVertexAttributes; i++) {
    const char* name = igraph_strvector_get(&vnames, i);
    std::string namestr = std::string(name);
    switch (igraph_vector_int_get(&vtypes, i)) {
      // case IGRAPH_ATTRIBUTE_DEFAULT:
      case IGRAPH_ATTRIBUTE_NUMERIC: {
        igraph_vector_t results;
        igraph_vector_init(&results, numVertices);
        igraph_cattribute_VANV(
          &graph, igraph_strvector_get(&vnames, i), igraph_vss_all(), &results);
        ar(make_nvp("vertex_attr_" + namestr, results));
        igraph_vector_destroy(&results);
      } break;
      case IGRAPH_ATTRIBUTE_STRING: {
        igraph_strvector_t strresults;
        igraph_strvector_init(&strresults, numVertices);
        igraph_cattribute_VASV(&graph,
                               igraph_strvector_get(&vnames, i),
                               igraph_vss_all(),
                               &strresults);
        ar(make_nvp("vertex_attr_" + namestr, strresults));
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
  const igraph_integer_t numEdgeAttributes = igraph_strvector_size(&enames);
  assert(igraph_strvector_size(&enames) == igraph_vector_int_size(&etypes));
  ar(make_nvp("edge_attr_names", enames));
  ar(make_nvp("edge_attr_types", etypes));
  //   ar(make_size_tag(numEdgeAttributes * 3));
  for (igraph_integer_t i = 0; i < numEdgeAttributes; i++) {
    const char* name = igraph_strvector_get(&enames, i);
    std::string namestr = std::string(name);
    switch (igraph_vector_int_get(&etypes, i)) {
      // case IGRAPH_ATTRIBUTE_DEFAULT:
      case IGRAPH_ATTRIBUTE_NUMERIC: {
        igraph_vector_t results;
        igraph_vector_init(&results, numEdges);
        igraph_cattribute_EANV(&graph,
                               igraph_strvector_get(&enames, i),
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &results);
        ar(make_nvp("edge_attr_" + namestr, results));
        igraph_vector_destroy(&results);
      } break;
      case IGRAPH_ATTRIBUTE_STRING: {
        igraph_strvector_t strresults;
        igraph_strvector_init(&strresults, numEdges);
        igraph_cattribute_EASV(&graph,
                               igraph_strvector_get(&enames, i),
                               igraph_ess_all(IGRAPH_EDGEORDER_ID),
                               &strresults);
        ar(make_nvp("edge_attr_" + namestr, strresults));
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
CEREAL_LOAD_FUNCTION_NAME(Archive& ar,
                          igraph_t& graph,
                          std::uint32_t const version)
{
  size_t numVertices;
  ar(make_nvp("num_vertices", numVertices));
  size_t numEdges;
  ar(make_nvp("num_edges", numEdges));
  igraph_vector_int_t allEdges;
  igraph_vector_int_init(&allEdges, static_cast<igraph_integer_t>(numEdges));
  ar(make_nvp("edges", allEdges));

  igraph_add_vertices(
    &graph, static_cast<igraph_integer_t>(numVertices), nullptr);
  igraph_add_edges(&graph, &allEdges, nullptr);
  igraph_vector_int_destroy(&allEdges);

  // deserialize vertex attributes
  igraph_strvector_t vnames;
  igraph_strvector_init(&vnames, 1);
  ar(make_nvp("vertex_attr_names", vnames));
  igraph_vector_int_t vtypes;
  igraph_vector_int_init(&vtypes, 1);
  ar(make_nvp("vertex_attr_types", vtypes));

  const igraph_integer_t numVertexAttributes = igraph_vector_int_size(&vtypes);
  for (igraph_integer_t i = 0; i < numVertexAttributes; ++i) {
    std::string attributeName = std::string(igraph_strvector_get(&vnames, i));
    const int attributeType = igraph_vector_int_get(&vtypes, i);
    switch (attributeType) {
      // case IGRAPH_ATTRIBUTE_DEFAULT:
      case IGRAPH_ATTRIBUTE_NUMERIC: {
        igraph_vector_t results;
        igraph_vector_init(&results,
                           static_cast<igraph_integer_t>(numVertices));
        ar(make_nvp("vertex_attr_" + attributeName, results));
        igraph_cattribute_VAN_setv(&graph, attributeName.c_str(), &results);
        igraph_vector_destroy(&results);
      }; break;
      case IGRAPH_ATTRIBUTE_STRING: {
        igraph_strvector_t strresults;
        igraph_strvector_init(&strresults,
                              static_cast<igraph_integer_t>(numVertices));
        ar(make_nvp("vertex_attr_" + attributeName, strresults));
        igraph_cattribute_VAS_setv(&graph, attributeName.c_str(), &strresults);
        igraph_strvector_destroy(&strresults);
      }; break;
      default:
        throw std::runtime_error("This attribute type (" +
                                 std::to_string(attributeType) +
                                 ") is not supported");
    }
  }
  igraph_vector_int_destroy(&vtypes);
  igraph_strvector_destroy(&vnames);

  // and same for edge attributes
  igraph_strvector_t enames;
  igraph_strvector_init(&enames, 1);
  ar(make_nvp("edge_attr_names", enames));
  igraph_vector_int_t etypes;
  igraph_vector_int_init(&etypes, 1);
  ar(make_nvp("edge_attr_types", etypes));

  const igraph_integer_t numEdgeAttributes = igraph_vector_int_size(&etypes);
  for (igraph_integer_t i = 0; i < numEdgeAttributes; ++i) {
    std::string attributeName = std::string(igraph_strvector_get(&enames, i));
    ;
    const int attributeType = igraph_vector_int_get(&etypes, i);
    switch (attributeType) {
      // case IGRAPH_ATTRIBUTE_DEFAULT:
      case IGRAPH_ATTRIBUTE_NUMERIC: {
        igraph_vector_t results;
        igraph_vector_init(&results, 1);
        ar(make_nvp("edge_attr_" + attributeName, results));
        igraph_cattribute_EAN_setv(&graph, attributeName.c_str(), &results);
        igraph_vector_destroy(&results);
      } break;
      case IGRAPH_ATTRIBUTE_STRING: {
        igraph_strvector_t strresults;
        igraph_strvector_init(&strresults, 1);
        ar(make_nvp("edge_attr_" + attributeName, strresults));
        igraph_cattribute_EAS_setv(&graph, attributeName.c_str(), &strresults);
        igraph_strvector_destroy(&strresults);
      } break;
      default:
        throw std::runtime_error("This attribute type (" +
                                 std::to_string(attributeType) +
                                 ") is not supported");
    }
  }
  igraph_vector_int_destroy(&etypes);
  igraph_strvector_destroy(&enames);
}

} // namespace cereal
CEREAL_CLASS_VERSION(igraph_t, 1);

//////////////////////////////////////////////////////////////////////////////////
// MARK: actual utils to serialize fast
namespace pylimer_tools::utils {

template<typename T>
std::string
serializeToString(T obj)
{
  std::stringstream os;
  {
    cereal::BinaryOutputArchive oarchive(os);
    oarchive(obj);
  }
  return os.str();
}

template<typename T>
void
deserializeFromString(T& obj, std::string& in)
{
  std::stringstream is(in);
  cereal::BinaryInputArchive iarchive(is);
  iarchive(obj);
}

template<typename T>
void
serializeToFile(T obj, const std::string file)
{
  std::ofstream os(file);
  if (file.ends_with("json")) {
    cereal::JSONOutputArchive oarchive(os);
    oarchive(obj);
  } else if (file.ends_with("xml")) {
    cereal::XMLOutputArchive oarchive(os);
    oarchive(obj);
  } else {
    cereal::BinaryOutputArchive oarchive(os);
    oarchive(obj);
  }
}

template<typename T>
void
deserializeFromFile(T& obj, const std::string file)
{
  std::ifstream is(file);
  if (file.ends_with("json")) {
    cereal::JSONInputArchive iarchive(is);
    iarchive(obj);
  } else if (file.ends_with("xml")) {
    cereal::XMLInputArchive iarchive(is);
    iarchive(obj);
  } else {
    cereal::BinaryInputArchive iarchive(is);
    iarchive(obj);
  }
}
}

#endif
