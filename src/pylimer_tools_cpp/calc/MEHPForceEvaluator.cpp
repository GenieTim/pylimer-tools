
#include "MEHPForceEvaluator.h"
#include "MEHPForceRelaxation.h"
#include "MEHPUtilityStructures.h"
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
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
        double constantMultiplier = this->kappa; // * 0.5 / s2;
        for (size_t j = 0; j < n; ++j) {
          grad[j] = 0.0;
        }
        for (size_t j = 0; j < this->net.nrOfSprings; ++j) {
          int a = this->net.springIndexA[j];
          int b = this->net.springIndexB[j];
          int nrOfDim = this->is2D ? 2 : 3;
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
     */
    float langevin_invf(float x)
    {
      float p, r, t;
      if ((fabsf(x) > 0.890625f) && (fabsf(x) <= 1.0f)) {
        r = copysignf(1.0f / (fabsf(x) - 1.0f), x);
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

    double NonGaussianSpringForceEvaluator::evaluateForceSetGradient(
      const size_t n,
      const Eigen::VectorXd& springDistances,
      const Eigen::VectorXd& u,
      double* grad) const
    {
      if (grad != nullptr) {
        throw std::runtime_error(
          "NonGaussianSpringForceEvaluator::evaluateForceSetGradient does not "
          "support setting the gradient.");
      }

      double s2 = springDistances.squaredNorm();

      return 0.5 * (this->kappa / this->l) * langevin_invf(s2 / (N * l));
    }
  }
}
}
