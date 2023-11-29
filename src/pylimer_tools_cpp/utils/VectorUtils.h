#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

#include <algorithm>
// #include <iostream>
#include <iterator>
#include <map>
#include <unordered_map>
#include <vector>
extern "C"
{
#include <igraph/igraph.h>
}
#include "utilityMacros.h"
#include <Eigen/Dense>
#include <cassert>

namespace pylimer_tools {
namespace utils {
  template<typename T>
  static inline void removeIfContained(std::vector<T>& vec, const T& value)
  {
    vec.erase(std::remove(vec.begin(), vec.end(), value), vec.end());
  }

  template<typename T>
  static inline bool contains(std::vector<T>& vec, const T value)
  {
    if (std::find(vec.begin(), vec.end(), value) == vec.end()) {
      return false;
    }
    return true;
  }

  template<typename T>
  static inline bool addIfNotContained(std::vector<T>& vec, const T value)
  {
    if (std::find(vec.begin(), vec.end(), value) == vec.end()) {
      vec.push_back(value);
      return true;
    }
    return false;
  }

  template<typename T>
  static inline T max_element(std::vector<T>& vec, const T defaultMax)
  {
    if (vec.size() == 0) {
      return defaultMax;
    }
    if (vec.size() == 1) {
      return vec[0];
    }
    T value = *std::max_element(vec.begin(), vec.end());
    return value;
  }

  typedef Eigen::Array<bool, Eigen::Dynamic, 1> ArrayXb;

  /**
   * @brief Remove a row from an Eigen vector
   *
   * @param vec
   * @param rowToRemove
   */
#define MAKE_REMOVE_ROW(EIGEN_TYPE)                                            \
  static inline void removeRow(                                                \
    EIGEN_TYPE& vec, unsigned int rowToRemove, bool noResize = false)          \
  {                                                                            \
    INVALIDARG_EXP_IFN(vec.size() > rowToRemove,                               \
                       "Cannot remove row " + std::to_string(rowToRemove) +    \
                         " from vector with size " +                           \
                         std::to_string(vec.size()) + "!");                    \
    unsigned int numRows = vec.size() - 1;                                     \
    vec.segment(rowToRemove, numRows - rowToRemove) =                          \
      vec.segment(rowToRemove + 1, numRows - rowToRemove);                     \
    if (!noResize) {                                                           \
      vec.conservativeResize(numRows);                                         \
    }                                                                          \
  }

  MAKE_REMOVE_ROW(Eigen::VectorXd);
  MAKE_REMOVE_ROW(Eigen::VectorXi);
  MAKE_REMOVE_ROW(Eigen::ArrayXi);
  MAKE_REMOVE_ROW(Eigen::ArrayXd);
  MAKE_REMOVE_ROW(ArrayXb);

  /**
   * @brief Remove sequential rows from an Eigen vector
   *
   * @param vec
   * @param rowToRemove
   */
#define MAKE_REMOVE_ROWS(EIGEN_TYPE)                                           \
  static inline void removeRows(EIGEN_TYPE& vec,                               \
                                unsigned int rowToStartRemove,                 \
                                unsigned int nrOfRowsToRemove,                 \
                                bool noResize = false)                         \
  {                                                                            \
    INVALIDARG_EXP_IFN(                                                        \
      vec.size() >= rowToStartRemove + nrOfRowsToRemove,                       \
      "Cannot remove rows " + std::to_string(nrOfRowsToRemove) + " from " +    \
        std::to_string(rowToStartRemove) + " from vector with size " +         \
        std::to_string(vec.size()) + "!");                                     \
    unsigned int numRows = vec.size() - nrOfRowsToRemove;                      \
    vec.segment(rowToStartRemove, numRows - rowToStartRemove) = vec.segment(   \
      rowToStartRemove + nrOfRowsToRemove, numRows - rowToStartRemove);        \
    if (!noResize) {                                                           \
      vec.conservativeResize(numRows);                                         \
    }                                                                          \
  }

  MAKE_REMOVE_ROWS(Eigen::VectorXd);
  MAKE_REMOVE_ROWS(Eigen::VectorXi);
  MAKE_REMOVE_ROWS(Eigen::ArrayXi);
  MAKE_REMOVE_ROWS(Eigen::ArrayXd);
  MAKE_REMOVE_ROWS(ArrayXb);

  template<typename T>
  static inline T last(const std::vector<T>& v)
  {
    return v[v.size() - 1];
  }
  /**
   * @brief Find whether a map contains a value
   *
   * @param map T0<T1, T2>
   * @param value
   * @return true|false
   */
  template<typename T0, typename T1>
  static inline bool set_has_key(const T0 map, const T1 key)
  {
#if __cplusplus >= 202002L
    // C++20 (and later) code
    return map.contains(key);
#else
    return map.count(key) > 0;
#endif
  }

  /**
   * @brief Find whether a map contains a value
   *
   * @param map T0<T1, T2>
   * @param value
   * @return true|false
   */
  template<typename T0, typename T1>
  static inline bool map_has_key(const T0 map, const T1 key)
  {
#if __cplusplus >= 202002L
    // C++20 (and later) code
    return map.contains(key);
#else
    return map.find(key) != map.end();
#endif
  }

  template<typename IN>
  static inline std::vector<IN> interleave(const std::vector<IN>& in1,
                                           const std::vector<IN>& in2)
  {
    size_t size = in1.size();
    assert(size == in2.size());
    std::vector<IN> out;
    out.reserve(2 * size);
    // interleave until at least one container is done
    for (size_t i = 0; i < size; ++i) {
      out.push_back(in1[i]);
      out.push_back(in2[i]);
    }

    return out; // both done
  }

  template<typename IN>
  static inline bool vector_has_duplicates(const std::vector<IN>& vec)
  {
    std::vector<IN> vecSorted;
    vecSorted.reserve(vec.size());
    std::copy(vec.begin(), vec.end(), std::back_inserter(vecSorted));
    std::sort(vecSorted.begin(), vecSorted.end());
    return std::adjacent_find(vecSorted.begin(), vecSorted.end()) !=
           vecSorted.end();
  }

  template<typename IN>
  static inline void eraseIndices(std::vector<IN> from,
                                  std::vector<long int>& indices)
  {
    for (auto index : indices) {
      from.erase(index);
    }
  }

#define MAKE_CONVERSION_FROM_STD_VEC_TO_IGRAPH(IGRAPH_VEC)                     \
  template<typename IN1>                                                       \
  static inline void StdVectorToIgraphVectorT(IN1& vectR, IGRAPH_VEC##_t* v)   \
  {                                                                            \
    size_t n = vectR.size();                                                   \
                                                                               \
    /* Make sure that there is enough space for the items in v */              \
    IGRAPH_VEC##_resize(v, n);                                                 \
                                                                               \
    /* Copy all the items */                                                   \
    for (size_t i = 0; i < n; ++i) {                                           \
      IGRAPH_VEC##_set(v, i, vectR[i]);                                        \
    }                                                                          \
  }

  MAKE_CONVERSION_FROM_STD_VEC_TO_IGRAPH(igraph_vector);
  MAKE_CONVERSION_FROM_STD_VEC_TO_IGRAPH(igraph_vector_int);

  static inline void StdVectorToIgraphVectorT(std::vector<std::string>& vectR,
                                              igraph_strvector_t* v)
  {
    size_t n = vectR.size();
    igraph_strvector_resize(v, n);
    for (size_t i = 0; i < n; ++i) {
      igraph_strvector_set(v, i, vectR[i].c_str());
    }
  }

  // MAKE_CONVERSION_FROM_STD_VEC_TO_IGRAPH(igraph_strvector);

#define MAKE_CONVERSION_FROM_IGRAPH_VEC_TO_STD(IGRAPH_VEC)                     \
  template<typename IN>                                                        \
  static inline void igraphVectorTToStdVector(IGRAPH_VEC##_t* v,               \
                                              std::vector<IN>& vectL)          \
  {                                                                            \
    long n = IGRAPH_VEC##_size(v);                                             \
                                                                               \
    /* Make sure that there is enough space for the items in v */              \
    vectL.clear();                                                             \
    vectL.reserve(n);                                                          \
                                                                               \
    /* Copy all the items */                                                   \
    for (size_t i = 0; i < n; ++i) {                                           \
      vectL.push_back(IGRAPH_VEC##_get(v, i));                                 \
    }                                                                          \
  }

  MAKE_CONVERSION_FROM_IGRAPH_VEC_TO_STD(igraph_vector);
  MAKE_CONVERSION_FROM_IGRAPH_VEC_TO_STD(igraph_vector_int);
  MAKE_CONVERSION_FROM_IGRAPH_VEC_TO_STD(igraph_strvector);

  template<typename IN>
  static inline std::vector<IN> initializeWithValue(size_t n, IN value)
  {
    std::vector<IN> result = std::vector<IN>(n, value);
    // for (size_t i = 0; i < n; ++i) {
    //   result.push_back(value);
    // }
    return result;
  }
} // namespace utils
} // namespace pylimer_tools

#endif
