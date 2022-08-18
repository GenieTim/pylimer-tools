#ifndef MEHP_FORCE_EVAL_H
#define MEHP_FORCE_EVAL_H
#include "MEHPUtilityStructures.h"
#include <iostream>

namespace pylimer_tools {
namespace calc {
  namespace mehp {
    double langevin_inv(double x);

    double csch(double x);

    // abstract class for having different force evaluations
    class MEHPForceEvaluator
    {
    protected:
      Network net;
      bool is2D = false;

    public:
      virtual ~MEHPForceEvaluator() = default;
      void setNetwork(Network& net) { this->net = net; }
      Network getNetwork() const { return this->net; }
      void setIs2D(bool is2D) { this->is2D = is2D; }
      bool getIs2D() { return this->is2D; };
      double evaluateForceSetGradient(const size_t n,
                                      const double* x,
                                      double* grad,
                                      void* f_data) const
      {
        Eigen::Map<const Eigen::VectorXd> u =
          Eigen::Map<const Eigen::VectorXd>(x, n);
        return evaluateForceSetGradient(n, u, grad, f_data);
      }

      double evaluateForceSetGradient(const size_t n,
                                      const Eigen::VectorXd& u,
                                      double* grad,
                                      void* f_data) const;

      virtual void prepareForEvaluations() = 0;
      virtual double evaluateForceSetGradient(
        const size_t n,
        const Eigen::VectorXd& springDistances,
        const Eigen::VectorXd& u,
        double* grad) const = 0;

      virtual double evaluateStressContribution(double springDistances[3],
                                                size_t i,
                                                size_t j) const = 0;
    };

    // example implementation of MEHPForceRelaxation for simple spring (phantom
    // systems)
    class SimpleSpringMEHPForceEvaluator : public MEHPForceEvaluator
    {
    protected:
      double kappa = 1.0;

    public:
      using MEHPForceEvaluator::getIs2D;
      using MEHPForceEvaluator::getNetwork;
      using MEHPForceEvaluator::setIs2D;
      using MEHPForceEvaluator::setNetwork;
      SimpleSpringMEHPForceEvaluator(double kappa = 1.0)
      {
        this->kappa = kappa;
      }

      double evaluateForceSetGradient(const size_t n,
                                      const Eigen::VectorXd& springDistances,
                                      const Eigen::VectorXd& u,
                                      double* grad) const override;
      double evaluateStressContribution(double springDistances[3],
                                        size_t i,
                                        size_t j) const override
      {
        return this->kappa * springDistances[i] * springDistances[j];
      }

      void prepareForEvaluations() override{};
    };

    class NonGaussianSpringForceEvaluator : public MEHPForceEvaluator
    {
    protected:
      double kappa = 1.0;
      double oneOverNl = 1.0;
      double oneOverl = 1.0;
      double N = 1.0;
      SimpleSpringMEHPForceEvaluator springForceEvaluator;

    public:
      using MEHPForceEvaluator::getIs2D;
      using MEHPForceEvaluator::getNetwork;
      using MEHPForceEvaluator::setIs2D;
      using MEHPForceEvaluator::setNetwork;
      NonGaussianSpringForceEvaluator(double kappa = 1.0,
                                      double N = 1.0,
                                      double l = 1.0)
      {
        this->kappa = kappa;
        this->springForceEvaluator = SimpleSpringMEHPForceEvaluator(kappa);
        assert(l > 0);
        assert(l * N > 0);
        this->N = N;
        this->oneOverNl = 1.0 / (N * l);
        this->oneOverl = 1.0 / l;
      }

      double evaluateForceSetGradient(const size_t n,
                                      const Eigen::VectorXd& springDistances,
                                      const Eigen::VectorXd& u,
                                      double* grad) const override;
      double evaluateStressContribution(double springDistances[3],
                                        size_t i,
                                        size_t j) const override
      {
        return this->kappa * springDistances[i] * springDistances[j];
      }

      void prepareForEvaluations() override
      {
        // propagate network and other config to decorated force evaluator
        this->springForceEvaluator.setNetwork(this->net);
        this->springForceEvaluator.setIs2D(this->is2D);
      }
    };
  }
}
}

#endif
