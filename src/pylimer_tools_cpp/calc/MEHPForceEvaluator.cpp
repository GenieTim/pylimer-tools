
#include "MEHPForceEvaluator.h"
#include "MEHPForceRelaxation.h"
#include "MEHPUtilityStructures.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h> // fabs, log, copysign, fma
#include <nlopt.hpp>
#include <string>
#include <tuple>
#include <vector>

namespace pylimer_tools {
namespace calc {
  namespace mehp {

    double MEHPForceEvaluator::evaluateForceSetGradient(
      const size_t n,
      const Eigen::VectorXd& u,
      double* grad,
      void* f_data) const
    {
      assert(n == this->net.nrOfNodes * 3);
      assert(u.size() == this->net.coordinates.size());
      Eigen::VectorXd springDistances =
        MEHPForceRelaxation::evaluateSpringDistances(&this->net, u, this->is2D);
      assert(n == this->net.nrOfNodes * 3);
      assert(u.size() == this->net.coordinates.size());

      return evaluateForceSetGradient(n, springDistances, u, grad);
    }

    double SimpleSpringMEHPForceEvaluator::evaluateForceSetGradient(
      const size_t n,
      const Eigen::VectorXd& springDistances,
      const Eigen::VectorXd& u,
      double* grad) const
    {
      assert(n == this->net.nrOfNodes * 3);
      assert(u.size() == this->net.coordinates.size());
      assert(n == this->net.nrOfNodes * 3);
      assert(u.size() == this->net.coordinates.size());

      double s2 = springDistances.squaredNorm();
      if (grad != nullptr) {
        const double constantMultiplier = this->kappa; // * 0.5 / s2;
        const int nrOfDim = this->is2D ? 2 : 3;
        for (size_t j = 0; j < n; ++j) {
          grad[j] = 0.0;
        }
        for (size_t j = 0; j < this->net.nrOfSprings; ++j) {
          const int a = this->net.springIndexA[j];
          const int b = this->net.springIndexB[j];
          for (size_t dir = 0; dir < nrOfDim; ++dir) {
            grad[3 * a + dir] +=
              springDistances[3 * j + dir] * constantMultiplier;
            grad[3 * b + dir] -=
              springDistances[3 * j + dir] * constantMultiplier;
          }
        }
      }
      // std::cout << "Evaluated force to " << std::setprecision(15)
      //           << 0.5 * kappa * s2 << " with kappa " << this->kappa
      //           << std::endl;
      return 0.5 * this->kappa * s2;
    };

#ifdef USE_FMA
#ifdef FP_FAST_FMAF
#define MYFMAF FP_FAST_FMAF(p, t, s)
#else
#define MYFMAF fmaf(p, t, s)
#endif
#else
#define MYFMAF(p, t, s) (p * t + s)
#endif

    /**
     * @brief Compute inverse Langevin function accurate to almost machine
     * precision
     *
     *  USE_FMA == 0: max. ulp error < 4.27, max. relative error < 4.43e-7
     *  USE_FMA == 1: max. ulp error < 3.64, max. relative error < 3.84e-7
     * @source: https://scicomp.stackexchange.com/a/30251
     */
    float langevin_invf(const float x)
    {
      float p, r, t;
      if (fabsf(x) > 0.99999) {
        // TODO: do better.
        // we have two problems: the value must be larger than whatever the
        // langevin should return, and second, the value should be small enough
        // to prevent overflow when summing them up.
        return 1e5 * x * x;
        //} else if ((fabsf(x) > 0.890625f) && (fabsf(x) <= 1.0f)) {
        // r = copysignf(1.0f / (fabsf(x) - 1.0f), x);
      } else {
        t = MYFMAF(x, 0.0f - x, 1.0f); // compute 1-x*x accurately
        t = logf(t);
        p = 2.18808651e-4f;                //  0x1.cae000p-13
        p = MYFMAF(p, t, -7.90076610e-3f); // -0x1.02e46ep-7
        p = MYFMAF(p, t, -7.12909698e-2f); // -0x1.240200p-4
        p = MYFMAF(p, t, -2.40409270e-1f); // -0x1.ec5bb2p-3
        p = MYFMAF(p, t, -4.14386481e-1f); // -0x1.a854eep-2
        p = MYFMAF(p, t, -4.05752033e-1f); // -0x1.9f7d76p-2
        p = MYFMAF(p, t, -2.56382942e-1f); // -0x1.068940p-2
        p = MYFMAF(p, t, -1.22061931e-1f); // -0x1.f3f736p-4
        p = MYFMAF(p, t, 5.00488468e-2f);  //  0x1.9a000ap-5
        p = MYFMAF(p, t, -1.84208602e-1f); // -0x1.79425cp-3
        p = MYFMAF(p, t, 3.98338169e-1f);  //  0x1.97e5f6p-2
        p = MYFMAF(p, t, -9.00006115e-1f); // -0x1.cccd9ap-1
        p = MYFMAF(p, t, 5.00000000e-1f);  //  0x1.000000p-1
        t = x + x;
        r = MYFMAF(p, t, t);
      }
      return r;
    }

    double csch(double x)
    {
      return 1. / sinh(x);
    };

    double NonGaussianSpringForceEvaluator::evaluateForceSetGradient(
      const size_t n,
      const Eigen::VectorXd& springDistances,
      const Eigen::VectorXd& u,
      double* grad) const
    {
      if (grad != nullptr) {
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
          const double linv = static_cast<double>(
            langevin_invf(static_cast<float>(r * this->oneOverNl)));
          if (std::isnan(linv) || std::isinf(linv)) {
            std::cerr << "Got " << linv << " for spring " << i
                     << " and distance " << r << std::endl;
          }
          for (size_t dir = 0; dir < nrOfDim; ++dir) {
            double springDistance = springDistances[3 * i + dir];
            double cschTerm = csch(linv);
            double gradTerm =
              -this->kappa * this->oneOverl * this->oneOverNl * springDistance;
            gradTerm =
              gradTerm / (r * (-cschTerm * cschTerm + 1. / (linv * linv)));
            // fixes to prevent nans
            if (std::isinf(cschTerm) || springDistance == 0.0) {
              gradTerm = 0.0;
            }
            if (std::isnan(gradTerm) || std::isinf(gradTerm)) {
              std::cerr << "Got " << gradTerm << " for gradTerm for spring " << i
                       << ", dir " << dir << " and distance " << r
                       << ", scsch term " << cschTerm << std::endl;
            }
            grad[3 * a + dir] += gradTerm;
            grad[3 * b + dir] -= gradTerm;
          }
        }
      }

      double force = 0.0;
      for (size_t i = 0; i < this->net.nrOfSprings; ++i) {
        double r =
          sqrt(springDistances[3 * i] * springDistances[3 * i] +
               springDistances[3 * i + 1] * springDistances[3 * i + 1] +
               springDistances[3 * i + 2] * springDistances[3 * i + 2]);
        double linv = static_cast<double>(
          langevin_invf(static_cast<float>(r * this->oneOverNl)));
        force += 0.5 * (this->kappa * this->oneOverl) * linv;
      }

      return force;
    }
  }
}
}
