#pragma once

#include <algorithm>
// #include <iostream>
#include "ExtraEigenTypes.h"
#include "utilityMacros.h"
#include <Eigen/Dense>
#include <cassert>
#include <iterator>
#include <map>
#include <type_traits>
#include <vector>

extern "C"
{
#include <igraph/igraph.h>
}

namespace pylimer_tools::utils {
/**
 * @brief Creates a vector with n times the specified value.
 *
 * @tparam IN the type of the elements in the vector
 * @param n the number of elements to initialize
 * @param value the value to initialize with
 * @return the vector of `n` times the `value`
 */
template<typename IN>
static inline std::vector<IN>
initializeWithValue(size_t n, IN value)
{
  std::vector<IN> result = std::vector<IN>(n, value);
  // for (size_t i = 0; i < n; ++i) {
  //   result.push_back(value);
  // }
  return result;
}

/**
 * @brief Gets the last element of a vector. Throws an exception if the vector
 * is empty.
 *
 * @tparam T the vector's type
 * @param v the vector to get the last element from
 * @return the last element of the vector
 */
template<typename T>
static inline T
last(const std::vector<T>& v)
{
  if (v.empty()) {
    throw std::runtime_error("Cannot get last element from an empty vector.");
  }
  return v[v.size() - 1];
}

/**
 * @brief Removes all occurrences of a specific value from a vector.
 *
 * This function uses the erase-remove idiom to efficiently remove all
 * elements that match the specified value from the vector. The vector is
 * modified in place.
 *
 * @tparam T The type of elements in the vector.
 * @param vec The vector from which to remove elements.
 * @param value The value to be removed from the vector.
 */
template<typename T>
static inline void
removeIfContained(std::vector<T>& vec, const T& value)
{
  vec.erase(std::remove(vec.begin(), vec.end(), value), vec.end());
}

/**
 *  @brief Checks if a vector contains a specific value.
 *
 * This function searches for the specified value in the vector and returns
 * whether it is present or not.
 *
 * @tparam T The type of elements in the vector.
 * @param vec The vector to search in.
 * @param value The value to search for.
 * @return true If the value is found in the vector.
 * @return false If the value is not found in the vector.
 */
template<typename T>
static inline bool
contains(std::vector<T>& vec, const T value)
{
  if (std::find(vec.begin(), vec.end(), value) == vec.end()) {
    return false;
  }
  return true;
}

/**
 * @brief Adds a value to a vector if it is not already present.
 *
 * This function checks if the specified value exists in the vector. If the
 * value is not found, it adds the value to the end of the vector and returns
 * true. If the value is already present, the function does nothing and
 * returns false.
 *
 * @tparam T The type of elements in the vector.
 * @param vec The vector to which the value may be added.
 * @param value The value to add to the vector if not already present.
 * @return true If the value was added to the vector.
 * @return false If the value was already present in the vector.
 */
template<typename T>
static inline bool
addIfNotContained(std::vector<T>& vec, const T value)
{
  if (!pylimer_tools::utils::contains(vec, value)) {
    vec.push_back(value);
    return true;
  }
  return false;
}

/**
 * @brief Finds the maximum element in a vector.
 *
 * This function returns the maximum element in the given vector. If the
 * vector is empty, it returns the provided default value. If the vector
 * contains only one element, it returns that element.
 *
 * @tparam T The type of elements in the vector.
 * @param vec The vector to search for the maximum element.
 * @param defaultMax The default value to return if the vector is empty.
 * @return T The maximum element in the vector, or the default value if the
 * vector is empty.
 */
template<typename T>
static inline T
max_element(std::vector<T>& vec, const T defaultMax)
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

/**
 * @brief Finds the index of a given value in a vector.
 *
 * This function searches for the first occurrence of a specified value in the
 * vector and returns its index. If the value is not found, it throws an
 * exception.
 *
 * @tparam T The type of elements in the vector.
 * @param vec The vector to search in.
 * @param val The value to search for.
 * @return size_t The index of the first occurrence of the value in the
 * vector.
 * @throws std::invalid_argument If the value is not found in the vector.
 */
template<typename T>
static inline size_t
index_of(std::vector<T>& vec, const T val)
{
  auto it = std::find(vec.begin(), vec.end(), val);
  if (it == vec.end()) {
    throw std::invalid_argument("Value must be present, not found in vector");
  }
  return std::distance(vec.begin(), it);
}

/**
 * @brief Compares two vectors for equality.
 *
 * This function checks if two vectors have the same size and contain
 * identical elements in the same order.
 *
 * @tparam T The type of elements in the vectors.
 * @param vec1 The first vector to compare.
 * @param vec2 The second vector to compare.
 * @return true If both vectors have the same size and all elements are equal.
 * @return false If the vectors have different sizes or any elements are not
 * equal.
 */
template<typename T>
static inline bool
equal(const std::vector<T>& vec1, const std::vector<T>& vec2)
{
  if (vec1.size() != vec2.size()) {
    return false;
  }
  for (size_t i = 0; i < vec1.size(); i++) {
    if (vec1[i] != vec2[i]) {
      return false;
    }
  }
  return true;
}

/**
 * @brief Calculates the Euclidean norm of each segment of a vector.
 *
 * This function divides the input vector into segments of a specified size
 * and calculates the Euclidean norm (L2 norm) of each segment. The results
 * are returned as a vector of norms.
 *
 * @param vecs The input vector to be segmented and processed.
 * @param segmentSize The size of each segment. Default is 3, which is common
 *                    for 3D vectors.
 * @return std::vector<double> A vector containing the norm of each segment.
 * @throws std::invalid_argument If segmentSize is not positive or if the
 * input vector's size is not a multiple of segmentSize.
 */
static inline std::vector<double>
segmentwise_norm(const Eigen::VectorXd& vecs, const size_t segmentSize = 3)
{
  INVALIDARG_EXP_IFN(segmentSize > 0, "Segmentwise requires a usable size");
  const size_t vecSize = static_cast<size_t>(vecs.size());
  INVALIDARG_EXP_IFN(vecSize % segmentSize == 0,
                     "The size of the supplied vector, " +
                       std::to_string(vecSize) +
                       " is not a multiple of the segment size, " +
                       std::to_string(segmentSize) + ".");
  std::vector<double> results;
  results.reserve(vecSize / segmentSize);
  for (size_t i = 0; i < vecSize / segmentSize; i++) {
    results.push_back(vecs
                        .segment(static_cast<Eigen::Index>(segmentSize * i),
                                 static_cast<Eigen::Index>(segmentSize))
                        .norm());
  }
  return results;
}

/**
 * @brief Calculates the maximum Euclidean norm among all segments of a
 * vector.
 *
 * This function divides the input vector into segments of a specified size
 * and calculates the Euclidean norm (L2 norm) of each segment. It then
 * returns the maximum norm value found across all segments.
 *
 * @param vecs The input vector to be segmented and processed.
 * @param segmentSize The size of each segment. Default is 3, which is common
 *                    for 3D vectors.
 * @return double The maximum norm value among all segments.
 * @throws std::invalid_argument If segmentSize is not positive or if the
 * input vector's size is not a multiple of segmentSize.
 */
static inline double
segmentwise_norm_max(const Eigen::VectorXd& vecs, const size_t segmentSize = 3)
{
  INVALIDARG_EXP_IFN(segmentSize > 0, "Segmentwise requires a usable size");
  const size_t vecSize = static_cast<size_t>(vecs.size());
  INVALIDARG_EXP_IFN(vecSize % segmentSize == 0,
                     "The size of the supplied vector, " +
                       std::to_string(vecSize) +
                       " is not a multiple of the segment size, " +
                       std::to_string(segmentSize) + ".");
  double result = 0.; //-DBL_MAX;

  for (size_t i = 0; i < vecSize / segmentSize; i++) {
    result = std::max(vecs
                        .segment(static_cast<Eigen::Index>(segmentSize * i),
                                 static_cast<Eigen::Index>(segmentSize))
                        .norm(),
                      result);
  }
  return result;
}

/**
 * @brief Calculates the mean Euclidean norm of all segments of a vector.
 *
 * This function divides the input vector into segments of a specified size
 * and calculates the Euclidean norm (L2 norm) of each segment. It then
 * returns the mean (average) of all segment norms.
 *
 * @param vecs The input vector to be segmented and processed.
 * @param segmentSize The size of each segment. Default is 3, which is common
 *                    for 3D vectors.
 * @return double The mean norm value across all segments.
 * @throws std::invalid_argument If segmentSize is not positive or if the
 * input vector's size is not a multiple of segmentSize.
 */
static inline double
segmentwise_norm_mean(const Eigen::VectorXd& vecs, const size_t segmentSize = 3)
{
  INVALIDARG_EXP_IFN(segmentSize > 0, "Segmentwise requires a usable size");
  const size_t vecSize = static_cast<size_t>(vecs.size());
  INVALIDARG_EXP_IFN(vecSize % segmentSize == 0,
                     "The size of the supplied vector, " +
                       std::to_string(vecSize) +
                       " is not a multiple of the segment size, " +
                       std::to_string(segmentSize) + ".");
  double result = 0.; //-DBL_MAX;
  const double denominator = 1. / static_cast<double>(vecSize / segmentSize);

  for (size_t i = 0; i < vecSize / segmentSize; i++) {
    const double norm = vecs
                          .segment(static_cast<Eigen::Index>(segmentSize * i),
                                   static_cast<Eigen::Index>(segmentSize))
                          .norm();
    result += (norm * denominator);
  }
  return result;
}

/**
 * @brief Checks if all components of a vector are finite numbers.
 *
 * This function iterates through each element of the input vector and checks
 * if it is a finite number (not infinity or NaN) using std::isfinite.
 * It returns false as soon as a non-finite element is encountered.
 *
 * @tparam IN The type of the vector or vector-like container.
 * @param vec The vector or vector-like container to check.
 * @return true If all elements in the vector are finite numbers.
 * @return false If at least one element in the vector is not a finite number.
 */

// For Eigen types - using a simpler SFINAE approach
template<typename IN>
static inline typename std::enable_if<
  std::is_class<IN>::value &&
    !std::is_same<IN, std::vector<typename IN::value_type>>::value,
  bool>::type
all_components_finite(const IN& vec)
{
  for (typename IN::Index i = 0; i < vec.size(); ++i) {
    if (!std::isfinite(vec[i])) {
      return false;
    }
  }
  return true;
}

// For std::vector types
template<typename IN>
static inline typename std::enable_if<
  std::is_same<IN, std::vector<typename IN::value_type>>::value,
  bool>::type
all_components_finite(const IN& vec)
{
  for (size_t i = 0; i < vec.size(); ++i) {
    if (!std::isfinite(vec[i])) {
      return false;
    }
  }
  return true;
}

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
MAKE_REMOVE_ROW(Eigen::ArrayXb);
#undef MAKE_REMOVE_ROW

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
    INVALIDINDEX_EXP_IFN(                                                      \
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
MAKE_REMOVE_ROWS(Eigen::ArrayXb);

#undef MAKE_REMOVE_ROWS

/**
 * @brief Remove multiple sequential rows from an Eigen vector based on
 * indices
 *
 * @param vec the Eigen vector to remove the rows from
 * @param indicesToRemove the indices of the rows to remove
 * @param isSorted whether the indices are already sorted in descending order
 * @param isUnique whether the indices are unique or not
 */
#define MAKE_REMOVE_ROWS(EIGEN_TYPE)                                           \
  static inline void removeRows(EIGEN_TYPE& vec,                               \
                                std::vector<size_t>& indicesToRemove,          \
                                bool isSorted = false,                         \
                                bool isUnique = false)                         \
  {                                                                            \
    if (indicesToRemove.empty())                                               \
      return;                                                                  \
                                                                               \
    /* Sort indices in descending order */                                     \
    if (!isSorted) {                                                           \
      std::ranges::sort(indicesToRemove, std::greater<size_t>());              \
    }                                                                          \
                                                                               \
    /* Remove duplicates */                                                    \
    if (!isUnique) {                                                           \
      indicesToRemove.erase(std::ranges::unique(indicesToRemove).begin(),      \
                            indicesToRemove.end());                            \
    }                                                                          \
    /* Check if the largest index is valid */                                  \
    if (indicesToRemove[0] >= static_cast<size_t>(vec.size())) {               \
      throw std::out_of_range("Index out of range");                           \
    }                                                                          \
                                                                               \
    size_t j = 0;                                                              \
    for (size_t i = 0; i < static_cast<size_t>(vec.size()); ++i) {             \
      if (j < indicesToRemove.size() &&                                        \
          i == indicesToRemove[indicesToRemove.size() - j - 1]) {              \
        ++j;                                                                   \
      } else {                                                                 \
        vec[static_cast<Eigen::Index>(i - j)] =                                \
          vec[static_cast<Eigen::Index>(i)];                                   \
      }                                                                        \
    }                                                                          \
                                                                               \
    vec.conservativeResize(static_cast<Eigen::Index>(                          \
      static_cast<size_t>(vec.size()) - indicesToRemove.size()));              \
  }

MAKE_REMOVE_ROWS(Eigen::VectorXd);
MAKE_REMOVE_ROWS(Eigen::VectorXi);
MAKE_REMOVE_ROWS(Eigen::ArrayXi);
MAKE_REMOVE_ROWS(Eigen::ArrayXd);
MAKE_REMOVE_ROWS(Eigen::ArrayXb);

#undef MAKE_REMOVE_ROWS

/**
 * @brief Remove multiple sequential rows from an Eigen vector based on
 * indices
 *
 * @param vec the Eigen vector to remove the rows from
 * @param indicesToRemove the indices of the rows to remove
 * @param isSorted whether the indices are already sorted in descending order
 * @param isUnique whether the indices are unique or not
 */
template<typename T>
static inline void
removeRows(std::vector<T>& vec,
           std::vector<size_t>& indicesToRemove,
           const bool isSorted = false,
           const bool isUnique = false)
{
  if (indicesToRemove.empty())
    return;

  /* Sort indices in descending order */
  if (!isSorted) {
    std::ranges::sort(indicesToRemove, std::greater<size_t>());
  }

  /* Remove duplicates */
  if (!isUnique) {
    indicesToRemove.erase(std::ranges::unique(indicesToRemove).begin(),
                          indicesToRemove.end());
  }
  /* Check if the largest index is valid */
  if (indicesToRemove[0] >= vec.size()) {
    throw std::out_of_range("Index " + std::to_string(indicesToRemove[0]) +
                            " is out of range for a vector of size " +
                            std::to_string(vec.size()) + ".");
  }

  size_t j = 0;
  for (size_t i = 0; i < vec.size(); ++i) {
    if (j < indicesToRemove.size() &&
        i == indicesToRemove[indicesToRemove.size() - j - 1]) {
      ++j;
    } else {
      vec[i - j] = vec[i];
    }
  }

  vec.resize(vec.size() - indicesToRemove.size());
}

static inline std::vector<long int>
getMappingForRenumbering(const std::vector<size_t>& removedValues,
                         const size_t nRemovableValues)
{
  // make sure things are sorted
  for (long int i = static_cast<long int>(removedValues.size()) - 2; i >= 0;
       --i) {
    INVALIDARG_EXP_IFN(
      removedValues[static_cast<size_t>(i + 1)] <
        removedValues[static_cast<size_t>(i)],
      "Values to remove must be sorted descending and unique, got values " +
        std::to_string(removedValues[static_cast<size_t>(i)]) + "@" +
        std::to_string(i) + " and " +
        std::to_string(removedValues[static_cast<size_t>(i + 1)]) + "@" +
        std::to_string(i + 1) + ".");
  }

  // for performance reasons, first assemble a new mapping
  std::vector<long int> newMapping =
    initializeWithValue<long int>(nRemovableValues, -1);
  long int idxInDeletedStrands =
    static_cast<long int>(removedValues.size()) - 1;
  long int nDeletedSoFar = 0;
  for (size_t i = 0; i < nRemovableValues; ++i) {
    if (idxInDeletedStrands >= 0 &&
        i == removedValues[static_cast<size_t>(idxInDeletedStrands)]) {
      idxInDeletedStrands -= 1;
      nDeletedSoFar += 1;
    } else {
      newMapping[i] = static_cast<long int>(i) - nDeletedSoFar;
    }
  }

  return newMapping;
}

/**
 * In a vector of vectors of indices to another structure,
 * renumber the indices to compensate for indices that have been removed in
 * the other structure.
 *
 * @param v the vector of vectors to renumber
 * @param removedValues the numbers that have been removed
 * @param nRemovableValues the maximum of the numbers
 */
static inline void
renumberWithMapping(std::vector<std::vector<size_t>>& v,
                    const std::vector<long int>& newMapping)
{
  // then, apply this mapping to all strands
  for (size_t linkI = 0; linkI < v.size(); ++linkI) {
    for (size_t& strandIdx : v[linkI]) {
      assert(strandIdx < newMapping.size());
      assert(newMapping[strandIdx] >= 0);
      strandIdx = static_cast<size_t>(newMapping[strandIdx]);
    }
  }
}

// For Eigen-like types - SFINAE to detect Eigen-like types
template<typename VecType>
static inline typename std::enable_if<
  std::is_class<VecType>::value &&
    !std::is_same<VecType, std::vector<typename VecType::value_type>>::value,
  void>::type
renumberWithMapping(VecType& v, const std::vector<long int>& newMapping)
{
  for (size_t i = 0; i < static_cast<size_t>(v.size()); ++i) {
    auto idx = static_cast<typename VecType::Index>(i);
    assert(v[idx] >= 0 && static_cast<size_t>(v[idx]) < newMapping.size());
    // the following assertion fails if the mapping is < 0,
    // which may happen if the value is one that should have been removed
    assert(newMapping[static_cast<size_t>(v[idx])] >= 0);
    v[idx] = static_cast<typename VecType::Scalar>(
      newMapping[static_cast<size_t>(v[idx])]);
  }
}

// For std::vector types
template<typename T>
static inline void
renumberWithMapping(std::vector<T>& v, const std::vector<long int>& newMapping)
{
  for (size_t i = 0; i < v.size(); ++i) {
    assert(v[i] >= 0 && static_cast<size_t>(v[i]) < newMapping.size());
    // the following assertion fails if the mapping is < 0,
    // which may happen if the value is one that should have been removed
    assert(newMapping[static_cast<size_t>(v[i])] >= 0);
    v[i] = static_cast<T>(newMapping[static_cast<size_t>(v[i])]);
  }
}

/**
 * @brief Appends elements from another vector to the end of the current
 * vector.
 *
 * @tparam T The type of elements in the vectors.
 * @param v The vector to which to append elements.
 * @param other The vector from which to append elements.
 * @param startOffset The starting index in the other vector to append.
 * @param endOffset The offset to the end of the other vector to append.
 * @note The order of elements in the current vector remains unchanged.
 */
template<typename T>
static inline void
append(std::vector<T>& v,
       const std::vector<T>& other,
       const int startOffset = 0,
       const int endOffset = 0)
{
  INVALIDARG_EXP_IFN(startOffset >= 0 && endOffset >= 0,
                     "Invalid offset, needs to be non-negative.");
  INVALIDARG_EXP_IFN(startOffset + endOffset <= other.size(),
                     "Invalid offset, expected total to at most skip adding "
                     "the whole vector, got " +
                       std::to_string(startOffset) + " + " +
                       std::to_string(endOffset) + " vs. a total size of " +
                       std::to_string(other.size()) + ".");
  v.insert(v.end(), other.begin() + startOffset, other.end() - endOffset);
}

/**
 * @brief Appends elements from the end another vector to the end of the
 * current vector.
 *
 * @tparam T The type of elements in the vectors.
 * @param v The vector to which to append elements.
 * @param other The vector from which to append elements.
 * @param startOffset The starting index (from the end) in the other vector to
 * append.
 * @param endOffset The offset to the end (from the start) of the other vector
 * to append.
 * @note The order of elements in the current vector remains unchanged.
 */
template<typename T>
static inline void
append_inverse(std::vector<T>& v,
               const std::vector<T>& other,
               const int startOffset = 0,
               const int endOffset = 0)
{
  INVALIDARG_EXP_IFN(startOffset >= 0 && endOffset >= 0,
                     "Invalid offset, needs to be non-negative.");
  INVALIDARG_EXP_IFN(startOffset + endOffset <= other.size(),
                     "Invalid offset, expected total to at most skip adding "
                     "the whole vector, got " +
                       std::to_string(startOffset) + " + " +
                       std::to_string(endOffset) + " vs. a total size of " +
                       std::to_string(other.size()) + ".");
  v.insert(v.end(), other.rbegin() + startOffset, other.rend() - endOffset);
}

/**
 * @brief Prepends elements from another vector to the beginning of the
 * current vector.
 *
 * @tparam T The type of elements in the vectors.
 * @param v The vector to which to prepend elements.
 * @param startOffset The starting index in the vector to prepend from.
 * @param endOffset The offset to the end of the other vector to prepend from.
 * @param other The vector from which to prepend elements.
 */
template<typename T>
static inline void
prepend(std::vector<T>& v,
        const std::vector<T>& other,
        const int startOffset = 0,
        const int endOffset = 0)
{
  INVALIDARG_EXP_IFN(startOffset >= 0 && endOffset >= 0,
                     "Invalid offset, needs to be non-negative.");
  INVALIDARG_EXP_IFN(startOffset + endOffset <= other.size(),
                     "Invalid offset, expected total to at most skip adding "
                     "the whole vector, got " +
                       std::to_string(startOffset) + " + " +
                       std::to_string(endOffset) + " vs. a total size of " +
                       std::to_string(other.size()) + ".");
  v.insert(v.begin(), other.begin() + startOffset, other.end() - endOffset);
}

/**
 * @brief Prepends elements from the end of another vector to the beginning of
 * the current vector.
 *
 * @tparam T The type of elements in the vectors.
 * @param v The vector to which to prepend elements.
 * @param startOffset The starting index (from the end) in the other vector to
 * append.
 * @param endOffset The offset to the end (from the start) of the other vector
 * to append.
 * @param other The vector from which to prepend elements.
 */
template<typename T>
static inline void
prepend_inverse(std::vector<T>& v,
                const std::vector<T>& other,
                const int startOffset = 0,
                const int endOffset = 0)
{
  INVALIDARG_EXP_IFN(startOffset >= 0 && endOffset >= 0,
                     "Invalid offset, needs to be non-negative.");
  INVALIDARG_EXP_IFN(startOffset + endOffset <= other.size(),
                     "Invalid offset, expected total to at most skip adding "
                     "the whole vector, got " +
                       std::to_string(startOffset) + " + " +
                       std::to_string(endOffset) + " vs. a total size of " +
                       std::to_string(other.size()) + ".");
  v.insert(v.begin(), other.rbegin() + startOffset, other.rend() - endOffset);
}

template<typename T>
static inline void
sort_remove_duplicates(std::vector<T>& v)
{
  std::ranges::sort(v);
  v.erase(std::ranges::unique(v).begin(), v.end());
}

/**
 * @brief Finds the index of the first occurrence of a value in a vector,
 * starting from a specified position.
 *
 * If the value is found at the start position, the function searches
 * backwards to find the first index. Otherwise, it searches forward from the
 * start position.
 *
 * @tparam T The type of elements in the vector.
 * @param v The vector to search in.
 * @param value The value to search for.
 * @param start The starting index for the search (default is 0).
 * @return size_t The index of the first occurrence of the value, or the size
 * of the vector if not found.
 */
template<typename T>
static inline size_t
first_occuring_index(const std::vector<T>& v, const T value, size_t start = 0)
{
  if (v[start] == value) {
    for (long int i = static_cast<long int>(start); i >= 0; i--) {
      if (v[static_cast<size_t>(i)] != value) {
        return static_cast<size_t>(i + 1);
      }
    }
    return 0;
  } else {
    for (size_t i = start; i < v.size(); i++) {
      if (v[i] == value) {
        return i;
      }
    }
  }
  return v.size();
}

/**
 * @brief Find whether a map contains a value
 *
 * @param map T0<T1, T2>
 * @param value
 * @return true|false
 */
template<typename T0, typename T1>
static inline bool
set_has_key(const T0 map, const T1 key)
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
static inline bool
map_has_key(const T0 map, const T1 key)
{
#if __cplusplus >= 202002L
  // C++20 (and later) code
  return map.contains(key);
#else
  return map.find(key) != map.end();
#endif
}

template<typename IN>
static inline std::vector<IN>
interleave(const std::vector<IN>& in1, const std::vector<IN>& in2)
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
static inline bool
vector_has_duplicates(const std::vector<IN>& vec)
{
  std::vector<IN> vecSorted;
  vecSorted.reserve(vec.size());
  std::copy(vec.begin(), vec.end(), std::back_inserter(vecSorted));
  std::sort(vecSorted.begin(), vecSorted.end());
  return std::adjacent_find(vecSorted.begin(), vecSorted.end()) !=
         vecSorted.end();
}

template<typename IN>
static inline bool
vector_approx_equal(const IN& v1,
                    const IN& v2,
                    const double absEps = 1e-12,
                    const bool echo = false)
{
  if (v1.size() != v2.size()) {
    return false;
  }
  for (size_t i = 0; i < v1.size(); ++i) {
    if (!APPROX_EQUAL(v1[i], v2[i], absEps)) {
      if (echo) {
        std::cout << "Detected unequality: v1[" << i << "] = " << v1[i]
                  << ", v2[" << i << "] = " << v2[i] << std::endl;
      }
      return false;
    }
  }
  return true;
}

template<typename IN>
static inline bool
vector_approx_rel_equal(const IN& v1,
                        const IN& v2,
                        const double eps = 1e-12,
                        const bool echo = false)
{
  if (v1.size() != v2.size()) {
    return false;
  }
  for (size_t i = 0; i < v1.size(); ++i) {
    if (!APPROX_REL_EQUAL(v1[i], v2[i], eps)) {
      if (echo) {
        std::cout << "Detected unequality: v1[" << i << "] = " << v1[i]
                  << ", v2[" << i << "] = " << v2[i] << std::endl;
      }
      return false;
    }
  }
  return true;
}

template<typename IN>
static inline void
eraseIndices(std::vector<IN> from, std::vector<long int>& indices)
{
  for (const auto index : indices) {
    from.erase(static_cast<size_t>(index));
  }
}

#define MAKE_CONVERSION_FROM_STD_VEC_TO_IGRAPH(IGRAPH_VEC)                     \
  template<typename IN1>                                                       \
  static inline void StdVectorToIgraphVectorT(IN1& vectR, IGRAPH_VEC##_t* v)   \
  {                                                                            \
    igraph_integer_t n = static_cast<igraph_integer_t>(vectR.size());          \
    IGRAPH_VEC##_resize(v, n);                                                 \
    for (igraph_integer_t i = 0; i < n; ++i) {                                 \
      IGRAPH_VEC##_set(v, i, vectR[static_cast<size_t>(i)]);                   \
    }                                                                          \
  }

MAKE_CONVERSION_FROM_STD_VEC_TO_IGRAPH(igraph_vector);
MAKE_CONVERSION_FROM_STD_VEC_TO_IGRAPH(igraph_vector_int);

static inline void
StdVectorToIgraphVectorT(std::vector<std::string>& vectR, igraph_strvector_t* v)
{
  const igraph_integer_t n = static_cast<igraph_integer_t>(vectR.size());
  igraph_strvector_resize(v, n);
  for (igraph_integer_t i = 0; i < n; ++i) {
    igraph_strvector_set(v, i, vectR[static_cast<size_t>(i)].c_str());
  }
}

// MAKE_CONVERSION_FROM_STD_VEC_TO_IGRAPH(igraph_strvector);

#define MAKE_CONVERSION_FROM_IGRAPH_VEC_TO_STD(IGRAPH_VEC)                     \
  template<typename IN>                                                        \
  static inline void igraphVectorTToStdVector(IGRAPH_VEC##_t* v,               \
                                              std::vector<IN>& vectL)          \
  {                                                                            \
    igraph_integer_t n = IGRAPH_VEC##_size(v);                                 \
                                                                               \
    /* Make sure that there is enough space for the items in v */              \
    vectL.clear();                                                             \
    vectL.reserve(static_cast<size_t>(n));                                     \
                                                                               \
    /* Copy all the items */                                                   \
    for (igraph_integer_t i = 0; i < n; ++i) {                                 \
      vectL.push_back(IGRAPH_VEC##_get(v, i));                                 \
    }                                                                          \
  }

MAKE_CONVERSION_FROM_IGRAPH_VEC_TO_STD(igraph_vector);
MAKE_CONVERSION_FROM_IGRAPH_VEC_TO_STD(igraph_vector_int);
MAKE_CONVERSION_FROM_IGRAPH_VEC_TO_STD(igraph_strvector);

/**
 * @brief Add an element to a (assumed) sorted vector where it belongs
 *
 * @tparam IN the vector's template parameter
 * @param vec the vector to add the element to
 * @param value the element to add
 */
template<typename IN>
static inline void
addToSorted(std::vector<IN>& vec, IN value)
{
  auto it = std::upper_bound(vec.begin(), vec.end(), value);
  vec.insert(it, value);
}
}
