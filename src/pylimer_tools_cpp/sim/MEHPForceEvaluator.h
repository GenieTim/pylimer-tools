#ifndef MEHP_FORCE_EVAL_H
#define MEHP_FORCE_EVAL_H
#include "./MEHPUtilityStructures.h"
// #include <iostream>

namespace pylimer_tools::sim::mehp {
double
langevin_inv(double x);

double
csch(double x);

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
  void setIs2D(const bool is2D) { this->is2D = is2D; }
  bool getIs2D() const { return this->is2D; };
  double evaluateForceSetGradient(const size_t n,
                                  const double* x,
                                  double* grad,
                                  void* f_data) const
  {
    const Eigen::Map<const Eigen::VectorXd> u =
      Eigen::Map<const Eigen::VectorXd>(x, static_cast<Eigen::Index>(n));
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
    double* grad) const = 0;

  virtual double evaluateStressContribution(double springDistances[3],
                                            size_t i,
                                            size_t j,
                                            size_t spring_index) const = 0;
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
  SimpleSpringMEHPForceEvaluator(const double kappa = 1.0)
  {
    this->kappa = kappa;
  }

  double evaluateForceSetGradient(const size_t n,
                                  const Eigen::VectorXd& springDistances,
                                  double* grad) const override;
  double evaluateStressContribution(double springDistances[3],
                                    const size_t i,
                                    const size_t j,
                                    const size_t spring_index) const override
  {
    return this->kappa * springDistances[i] * springDistances[j] /
           this->net
             .springsContourLength[static_cast<Eigen::Index>(spring_index)];
  }

  void prepareForEvaluations() override {};
};

class NonGaussianSpringForceEvaluator : public MEHPForceEvaluator
{
protected:
  double kappa = 1.0;
  double oneOverl = 1.0;
  double l = 1.0;
  SimpleSpringMEHPForceEvaluator springForceEvaluator;

public:
  using MEHPForceEvaluator::getIs2D;
  using MEHPForceEvaluator::getNetwork;
  using MEHPForceEvaluator::setIs2D;
  using MEHPForceEvaluator::setNetwork;
  NonGaussianSpringForceEvaluator(const double kappa = 1.0,
                                  double N = 1.0,
                                  const double l = 1.0)
  {
    this->kappa = kappa;
    this->springForceEvaluator = SimpleSpringMEHPForceEvaluator(kappa);
    assert(l > 0);
    this->l = l;
    this->oneOverl = 1.0 / l;
  }

  double evaluateForceSetGradient(const size_t n,
                                  const Eigen::VectorXd& springDistances,
                                  double* grad) const override;
  double evaluateStressContribution(double springDistances[3],
                                    const size_t i,
                                    const size_t j,
                                    const size_t spring_index) const override
  {
    return this->kappa * springDistances[i] * springDistances[j] /
           this->net
             .springsContourLength[static_cast<Eigen::Index>(spring_index)];
  }

  void prepareForEvaluations() override
  {
    // propagate network and other config to decorated force evaluator
    this->springForceEvaluator.setNetwork(this->net);
    this->springForceEvaluator.setIs2D(this->is2D);
  }
};
}

#endif
