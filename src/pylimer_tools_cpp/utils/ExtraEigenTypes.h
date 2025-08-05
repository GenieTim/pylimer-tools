#pragma once

#include <fstream>

#include <Eigen/Dense>
#include <Eigen/SparseCore>

namespace Eigen {

typedef Array<long int, 3, 1> Array3li;
typedef Array<long int, Dynamic, 1> ArrayXli;
typedef Array<size_t, Dynamic, 1> ArrayXst;
typedef Array<bool, Dynamic, 1> ArrayXb;

/**
 * @brief Equality comparison between Eigen vector and std::vector
 *
 * @tparam Derived The derived Eigen type
 * @tparam T The scalar type of the std::vector
 * @param eigenVec The Eigen vector
 * @param stdVec The std::vector
 * @return true if vectors are equal
 * @return false otherwise
 */
template<typename Derived, typename T>
bool
operator==(const Eigen::DenseBase<Derived>& eigenVec,
           const std::vector<T>& stdVec)
{
  if (eigenVec.size() != static_cast<Eigen::Index>(stdVec.size())) {
    return false;
  }

  for (Eigen::Index i = 0; i < eigenVec.size(); ++i) {
    if (eigenVec(i) != stdVec[static_cast<size_t>(i)]) {
      return false;
    }
  }

  return true;
}

/**
 * @brief Equality comparison between std::vector and Eigen vector
 *
 * @tparam T The scalar type of the std::vector
 * @tparam Derived The derived Eigen type
 * @param stdVec The std::vector
 * @param eigenVec The Eigen vector
 * @return true if vectors are equal
 * @return false otherwise
 */
template<typename T, typename Derived>
bool
operator==(const std::vector<T>& stdVec,
           const Eigen::DenseBase<Derived>& eigenVec)
{
  return eigenVec == stdVec;
}

template<typename Derived>
typename Derived::Scalar
median(Eigen::DenseBase<Derived>& d)
{
  auto r{ d.reshaped() };
  std::sort(r.begin(), r.end());
  return r.size() % 2 == 0 ? r.segment((r.size() - 2) / 2, 2).mean()
                           : r(r.size() / 2);
}

template<typename Derived>
typename Derived::Scalar
median(const Eigen::DenseBase<Derived>& d)
{
  typename Derived::PlainObject m{ d.replicate(1, 1) };
  return median(m);
}

static inline bool
isSelfAdjoint(const SparseMatrix<double>& mat, const double tol = 1e-9)
{
  if (mat.rows() != mat.cols())
    return false; // Must be square

  for (Eigen::Index k = 0; k < mat.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
      const Eigen::Index row = it.row();
      const Eigen::Index col = it.col();
      const double value = it.value();

      // Check symmetry within tolerance
      if (std::abs(value - mat.coeff(col, row)) > tol) {
        return false;
      }
    }
  }
  return true;
}

static inline void
saveSparseMatrix(const Eigen::SparseMatrix<double>& mat,
                 const std::string& filename)
{
  std::ofstream file(filename);
  if (!file.is_open())
    return;

  file << "%%MatrixMarket matrix coordinate real general\n";
  file << mat.rows() << " " << mat.cols() << " " << mat.nonZeros() << "\n";

  for (int k = 0; k < mat.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it) {
      file << (it.row() + 1) << " " << (it.col() + 1) << " " << it.value()
           << "\n";
    }
  }

  file.close();
}

static inline void
saveDenseVector(const Eigen::VectorXd& vec, const std::string& filename)
{
  std::ofstream file(filename);
  if (!file.is_open())
    return;

  for (const double i : vec) {
    file << i << "\n";
  }

  file.close();
}

}

#ifdef OPENMP_FOUND
#pragma omp declare reduction(+ : Eigen::VectorXd : omp_out =                  \
                                omp_out + omp_in)                              \
  initializer(omp_priv = Eigen::VectorXd::Zero(omp_orig.size()))
#pragma omp declare reduction(+ : Eigen::Matrix3d : omp_out =                  \
                                omp_out + omp_in)                              \
  initializer(omp_priv = Eigen::Matrix3d::Zero())
#endif
