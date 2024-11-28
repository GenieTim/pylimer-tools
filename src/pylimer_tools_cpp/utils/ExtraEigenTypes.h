#ifndef EXTRA_EIGEN_TYPES_H
#define EXTRA_EIGEN_TYPES_H

#include <Eigen/Dense>

namespace Eigen {

typedef Array<long int, 3, 1> Array3li;
typedef Array<long int, Dynamic, 1> ArrayXli;
typedef Array<size_t, Dynamic, 1> ArrayXst;
typedef Array<bool, Dynamic, 1> ArrayXb;

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
}

#ifdef OPENMP_FOUND
#pragma omp declare reduction(+ : Eigen::VectorXd : omp_out =                  \
omp_out + omp_in)                              \
initializer(omp_priv = Eigen::VectorXd::Zero(omp_orig.size()))
#pragma omp declare reduction(+ : Eigen::Matrix3d : omp_out =                  \
omp_out + omp_in)                              \
initializer(omp_priv = Eigen::Matrix3d::Zero())
#endif

#endif
