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

  /**
   * @brief Remove a row from an Eigen vector
   *
   * @param vec
   * @param rowToRemove
   */
  template<typename T>
  static inline void removeRow(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec,
                 unsigned int rowToRemove)
  {
    INVALIDARG_EXP_IFN(vec.size() > rowToRemove,
                       "Cannot remove row " + std::to_string(rowToRemove) +
                         " from vector with size " +
                         std::to_string(vec.size()) + "!");
    unsigned int numRows = vec.size() - 1;
    vec.segment(rowToRemove, numRows - rowToRemove) =
      vec.segment(rowToRemove + 1, numRows - rowToRemove);
    vec.conservativeResize(numRows);
  }

  /**
   * @brief Remove sequential rows from an Eigen vector
   *
   * @param vec
   * @param rowToRemove
   */
  template<typename T>
  static inline void removeRows(Eigen::Matrix<T, Eigen::Dynamic, 1>& vec,
                  unsigned int rowToStartRemove,
                  unsigned int nrOfRowsToRemove)
  {
    INVALIDARG_EXP_IFN(
      vec.size() > rowToStartRemove + nrOfRowsToRemove,
      "Cannot remove rows " + std::to_string(nrOfRowsToRemove) + " from " +
        std::to_string(rowToStartRemove) + " from vector with size " +
        std::to_string(vec.size()) + "!");
    unsigned int numRows = vec.size() - nrOfRowsToRemove;
    vec.segment(rowToStartRemove, numRows - rowToStartRemove) =
      vec.segment(rowToStartRemove + 1, numRows - rowToStartRemove);
    vec.conservativeResize(numRows);
  }

  /**
   * @brief Remove a row from an Eigen vector
   *
   * @param vec
   * @param rowToRemove
   */
  template<typename T>
  static inline void removeRow(Eigen::Array<T, Eigen::Dynamic, 1>& vec,
                 unsigned int rowToRemove)
  {
    INVALIDARG_EXP_IFN(vec.size() > rowToRemove,
                       "Cannot remove row " + std::to_string(rowToRemove) +
                         " from vector with size " +
                         std::to_string(vec.size()) + "!");
    unsigned int numRows = vec.size() - 1;
    vec.segment(rowToRemove, numRows - rowToRemove) =
      vec.segment(rowToRemove + 1, numRows - rowToRemove);
    vec.conservativeResize(numRows);
  }

  /**
   * @brief Remove sequential rows from an Eigen vector
   *
   * @param vec
   * @param rowToRemove
   */
  template<typename T>
  static inline void removeRows(Eigen::Array<T, Eigen::Dynamic, 1>& vec,
                  unsigned int rowToStartRemove,
                  unsigned int nrOfRowsToRemove)
  {
    INVALIDARG_EXP_IFN(
      vec.size() > rowToStartRemove + nrOfRowsToRemove,
      "Cannot remove rows " + std::to_string(nrOfRowsToRemove) + " from " +
        std::to_string(rowToStartRemove) + " from vector with size " +
        std::to_string(vec.size()) + "!");
    unsigned int numRows = vec.size() - nrOfRowsToRemove;
    vec.segment(rowToStartRemove, numRows - rowToStartRemove) =
      vec.segment(rowToStartRemove + 1, numRows - rowToStartRemove);
    vec.conservativeResize(numRows);
  }

  template<typename T>
  static inline T last(std::vector<T>& v)
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
  static inline bool map_has_key(T0 map, T1 key)
  {
#if __cplusplus >= 202002L
    // C++20 (and later) code
    return map.contains(key);
#else
    return map.find(key) != map.end();
#endif
  }

  template<typename IN>
  static inline std::vector<IN> interleave(std::vector<IN> in1,
                                           std::vector<IN> in2)
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
  static inline bool vector_has_duplicates(std::vector<IN> vec)
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
                                  std::vector<long int> indices)
  {
    for (auto index : indices) {
      from.erase(index);
    }
  }

  template<typename IN1>
  static inline void StdVectorToIgraphVectorT(IN1& vectR, igraph_vector_t* v)
  {
    size_t n = vectR.size();

    /* Make sure that there is enough space for the items in v */
    igraph_vector_resize(v, n);

    /* Copy all the items */
    for (size_t i = 0; i < n; ++i) {
      igraph_vector_set(v, i, vectR[i]);
    }
  }

  template<typename IN1>
  static inline void StdVectorToIgraphVectorT(IN1& vectR,
                                              igraph_vector_int_t* v)
  {
    size_t n = vectR.size();

    /* Make sure that there is enough space for the items in v */
    igraph_vector_int_resize(v, n);

    /* Copy all the items */
    for (size_t i = 0; i < n; ++i) {
      igraph_vector_int_set(v, i, vectR[i]);
    }
  }

  template<typename IN>
  static inline void igraphVectorTToStdVector(igraph_vector_t* v,
                                              std::vector<IN>& vectL)
  {
    long n = igraph_vector_size(v);

    vectL.clear();
    vectL.reserve(n);

    for (long i = 0; i < n; ++i) {
      vectL.push_back(igraph_vector_get(v, i));
    }
  }

  template<typename IN>
  static inline void igraphVectorTToStdVector(igraph_vector_int_t* v,
                                              std::vector<IN>& vectL)
  {
    long n = igraph_vector_int_size(v);

    vectL.clear();
    vectL.reserve(n);

    for (long i = 0; i < n; ++i) {
      vectL.push_back(igraph_vector_int_get(v, i));
    }
  }

  template<typename IN>
  static inline std::vector<IN> initializeWithValue(size_t n, IN value)
  {
    std::vector<IN> result;
    result.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      result.push_back(value);
    }
    return result;
  }
} // namespace utils
} // namespace pylimer_tools

#endif
