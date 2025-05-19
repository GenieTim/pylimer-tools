
#include "MEHPForceEvaluator.h"
#include "MEHPForceRelaxation.h"
#include "MEHPUtilityStructures.h"
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>
#include <math.h> // fabs, log, copysign, fma
#include <nlopt.hpp>

namespace pylimer_tools::sim::mehp {

double
MEHPForceEvaluator::evaluateForceSetGradient(const size_t n,
                                             const Eigen::VectorXd& u,
                                             double* grad,
                                             void* f_data) const
{
  assert(n == this->net.nrOfNodes * 3);
  assert(u.size() == this->net.coordinates.size());
  const Eigen::VectorXd springDistances =
    -1. *
    MEHPForceRelaxation::evaluateSpringDistances(&this->net, u, this->is2D);
  assert(n == this->net.nrOfNodes * 3);
  assert(u.size() == this->net.coordinates.size());

  return evaluateForceSetGradient(n, springDistances, grad);
}

double
SimpleSpringMEHPForceEvaluator::evaluateForceSetGradient(
  const size_t n,
  const Eigen::VectorXd& springDistances,
  double* grad) const
{
  assert(n == this->net.nrOfNodes * 3);
  assert(springDistances.size() == 3 * this->net.nrOfSprings);

  double s2 = 0.0;
  for (size_t i = 0; i < this->net.nrOfSprings; ++i) {
    s2 += springDistances.segment(3 * i, 3).squaredNorm() /
          this->net.springsContourLength[i];
  }
  if (grad != nullptr) {
    const double constantMultiplier = this->kappa; // * 0.5 / s2;
    const int nrOfDim = this->is2D ? 2 : 3;
    for (size_t j = 0; j < n; ++j) {
      grad[j] = 0.0;
    }
    for (size_t j = 0; j < this->net.nrOfSprings; ++j) {
      const int a = this->net.springIndexA[j];
      const int b = this->net.springIndexB[j];
      const double Nterm = 1.0 / this->net.springsContourLength[j];
      for (size_t dir = 0; dir < nrOfDim; ++dir) {
        grad[3 * a + dir] +=
          springDistances[3 * j + dir] * constantMultiplier * Nterm;
        grad[3 * b + dir] -=
          springDistances[3 * j + dir] * constantMultiplier * Nterm;
      }
    }
  }
  // std::cout << "Evaluated force to " << std::setprecision(15)
  //           << 0.5 * kappa * s2 << " with kappa " << this->kappa
  //           << std::endl;
  return 0.5 * this->kappa * s2;
};

#ifdef USE_FMA
#ifdef FP_FAST_FMA
#define MYFMA FP_FAST_FMA(p, t, s)
#else
#define MYFMA fma(p, t, s)
#endif
#else
#define MYFMA(p, t, s) (p * t + s)
#endif

/**
 * @brief Compute inverse Langevin function accurate to almost machine
 * precision
 *
 *  USE_FMA == 0: max. ulp error < 4.27, max. relative error < 4.43e-7
 *  USE_FMA == 1: max. ulp error < 3.64, max. relative error < 3.84e-7
 * @source: https://scicomp.stackexchange.com/a/30251
 */
// double langevin_inv(const double x)
// {
//   double p, r, t;
//   if (std::fabs(x) > 0.99999) {
//     // TODO: do better.
//     // we have two problems: The value must be larger than whatever the
//     // langevin should return, and second, the value should be small
//     enough
//     // to prevent overflow when summing them up.
//     return 1e5 * x * x;
//     //} else if ((std::fabs(x) > 0.890625f) && (std::fabs(x) <= 1.0f)) {
//     // r = copysignf(1.0f / (std::fabs(x) - 1.0f), x);
//   } else {
//     t = MYFMA(x, 0.0 - x, 1.0); // compute 1-x*x accurately
//     t = log(t);
//     p = 2.18808651e-4;               //  0x1.cae000p-13
//     p = MYFMA(p, t, -7.90076610e-3); // -0x1.02e46ep-7
//     p = MYFMA(p, t, -7.12909698e-2); // -0x1.240200p-4
//     p = MYFMA(p, t, -2.40409270e-1); // -0x1.ec5bb2p-3
//     p = MYFMA(p, t, -4.14386481e-1); // -0x1.a854eep-2
//     p = MYFMA(p, t, -4.05752033e-1); // -0x1.9f7d76p-2
//     p = MYFMA(p, t, -2.56382942e-1); // -0x1.068940p-2
//     p = MYFMA(p, t, -1.22061931e-1); // -0x1.f3f736p-4
//     p = MYFMA(p, t, 5.00488468e-2);  //  0x1.9a000ap-5
//     p = MYFMA(p, t, -1.84208602e-1); // -0x1.79425cp-3
//     p = MYFMA(p, t, 3.98338169e-1);  //  0x1.97e5f6p-2
//     p = MYFMA(p, t, -9.00006115e-1); // -0x1.cccd9ap-1
//     p = MYFMA(p, t, 5.00000000e-1);  //  0x1.000000p-1
//     t = x + x;
//     r = MYFMA(p, t, t);
//   }
//   return r;
// }

/**
 * @brief Compute inverse Langevin function approximation
 *
 * @source
 * https://www.sciencedirect.com/science/article/pii/S0377025715001007?via%3Dihub#e0275
 *
 * @param x
 * @return double
 */
double
langevin_inv(const double x)
{
  if (x < 0.0) {
    return langevin_inv(-x);
  }
  if (x > 0.99999) {
    return 1e5 * (1 + x * x);
  }
  const double x2 = x * x;
  return 3. * x / ((1. - x2) * (1. + 0.5 * x2));
}

/**
 * @brief Compute inverse Langevin function approximation
 *
 * @source
 * https://www.sciencedirect.com/science/article/pii/S0377025715001007?via%3Dihub#e0275
 *
 * @param x
 * @return double
 */
// double langevin_inv(const double x)
// {
//   if (x < 0.0) {
//     return langevin_inv(-x);
//   }
//   if (x > 0.99999) {
//     return 1e5 * (1 + x * x);
//   }
//   const double x2 = x * x;
//   const double x4 = x2 * x2;
//   const double x6 = x4 * x2;
//   return ((3 * x - (x / 5.) * (6. * x2 + x4 + 2 * x6)) / (1. - x2));
// }

double
csch(const double x)
{
  return 1. / sinh(x);
};

double
NonGaussianSpringForceEvaluator::evaluateForceSetGradient(
  const size_t n,
  const Eigen::VectorXd& springDistances,
  double* grad) const
{
  // this->springForceEvaluator.setNetwork(this->net);
  // this->springForceEvaluator.setIs2D(this->is2D);

  double force = 0.0; // this->springForceEvaluator.evaluateForceSetGradient(n,
  // springDistances, u, grad);

  if (grad != nullptr) {
    // no reset of grad to 0 needed, will already be done by
    // springForceEvaluator
    for (size_t j = 0; j < n; ++j) {
      grad[j] = 0.0;
    }
    const int nrOfDim = this->is2D ? 2 : 3;
    for (size_t i = 0; i < this->net.nrOfSprings; ++i) {
      const int a = this->net.springIndexA[i];
      const int b = this->net.springIndexB[i];
      const double r =
        sqrt(springDistances[3 * i] * springDistances[3 * i] +
             springDistances[3 * i + 1] * springDistances[3 * i + 1] +
             springDistances[3 * i + 2] * springDistances[3 * i + 2]);
      const double linv =
        langevin_inv(r * this->oneOverl / (this->net.springsContourLength[i]));
      if (std::isnan(linv) || std::isinf(linv)) {
        std::cerr << "Got " << linv << " for spring " << i << " and distance "
                  << r << std::endl;
      }
      // dF/dr
      const double fr = (this->oneOverl) * linv;

      for (size_t dir = 0; dir < nrOfDim; ++dir) {
        const double springDistance = springDistances[3 * i + dir];
        // dr/dui
        double gradTerm = (springDistance / r) * fr;
        if (r == 0.0) {
          gradTerm = 0.0;
        }
        if (std::isnan(gradTerm)) {
          std::cerr << "Got " << gradTerm << " grad term, with " << linv
                    << " for spring " << i << " and distance " << r
                    << std::endl;
        }

        grad[3 * a + dir] += gradTerm;
        grad[3 * b + dir] -= gradTerm;
      }
    }
  }

  for (size_t i = 0; i < this->net.nrOfSprings; ++i) {
    const double r =
      std::sqrt(springDistances[3 * i] * springDistances[3 * i] +
                springDistances[3 * i + 1] * springDistances[3 * i + 1] +
                springDistances[3 * i + 2] * springDistances[3 * i + 2]);
    const double rOverNl =
      r * this->oneOverl / (this->net.springsContourLength[i]);
    const double beta = langevin_inv(rOverNl);
    const double cschTerm =
      csch(beta) *
      beta; // can be inf or 0. In case of zero, std::log() returns -inf.
    if (beta > 0.0 && !std::isinf(cschTerm) && cschTerm > 0.0) {
      force += this->net.springsContourLength[i] *
               (rOverNl * beta + std::log(cschTerm));
    } else if (beta > 0.0 && (std::isinf(cschTerm) || cschTerm == 0.0)) {
      force += this->net.springsContourLength[i] * rOverNl * beta;
    }
  }

  return force;
}
}
