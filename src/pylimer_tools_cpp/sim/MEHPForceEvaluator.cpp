#include "MEHPForceEvaluator.h"
#include "../utils/utilityMacros.h"
#include "MEHPNetwork.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace pylimer_tools::sim::mehp {

namespace {

  void validateContourLength(const double contourLength)
  {
    INVALIDARG_EXP_IFN(std::isfinite(contourLength) && contourLength > 0.0,
                       "Spring contour length must be finite and positive.");
  }

  void validateNetworkAndSpringVectors(const MEHPNetwork& net,
                                       const Eigen::VectorXd& springVectors,
                                       const double* const grad)
  {
    INVALIDARG_EXP_IFN(
      springVectors.size() == static_cast<Eigen::Index>(3 * net.nrOfSprings),
      "springVectors size (" + std::to_string(springVectors.size()) +
        ") does not match 3 * net.nrOfSprings (" +
        std::to_string(3 * net.nrOfSprings) + ").");
    INVALIDARG_EXP_IFN(
      net.springContourLength.size() ==
        static_cast<Eigen::Index>(net.nrOfSprings),
      "net.springContourLength size does not match net.nrOfSprings.");
    INVALIDARG_EXP_IFN(
      net.springIndexA.size() == static_cast<Eigen::Index>(net.nrOfSprings) &&
        net.springIndexB.size() == static_cast<Eigen::Index>(net.nrOfSprings),
      "net.springIndexA / springIndexB sizes do not match net.nrOfSprings.");
    INVALIDARG_EXP_IFN(
      net.springContourLength.allFinite() &&
        (net.springContourLength.array() > 0.0).all(),
      "All spring contour lengths must be finite and positive.");
    if (grad != nullptr) {
      for (Eigen::Index s = 0; s < static_cast<Eigen::Index>(net.nrOfSprings);
           ++s) {
        const int a = net.springIndexA[s];
        const int b = net.springIndexB[s];
        INVALIDARG_EXP_IFN(
          a >= 0 && static_cast<size_t>(a) < net.nrOfLinks && b >= 0 &&
            static_cast<size_t>(b) < net.nrOfLinks,
          "Spring endpoint link index out of range [0, net.nrOfLinks).");
      }
    }
  }

  [[nodiscard]] double computeNonGaussianSpringEnergyScalar(
    const double r,
    const double contourLength,
    const double oneOverl)
  {
    if (r <= 0.0) {
      return 0.0;
    }
    const double rOverNl = (r * oneOverl) / contourLength;
    const double beta = langevin_inv(rOverNl);
    if (beta <= 0.0) {
      return 0.0;
    }
    double logCschTerm = 0.0;
    if (beta < 1e-8) {
      logCschTerm = -(beta * beta) / 6.0;
    } else if (beta <= 500.0) {
      const double cschTerm = csch(beta) * beta;
      logCschTerm = (cschTerm > 0.0 && std::isfinite(cschTerm))
                      ? std::log(cschTerm)
                      : (std::log(2.0 * beta) - beta);
    } else if (std::isfinite(beta)) {
      logCschTerm = std::log(2.0 * beta) - beta;
    }
    return contourLength * ((rOverNl * beta) + logCschTerm);
  }

  [[nodiscard]] double computeNonGaussianForceOverR(const double r,
                                                    const double contourLength,
                                                    const double oneOverl)
  {
    if (r <= 0.0) {
      return 0.0;
    }
    const double x = (r * oneOverl) / contourLength;
    if (x < 1e-8) {
      return (3.0 * oneOverl * oneOverl) / contourLength;
    }
    return (oneOverl * langevin_inv(x)) / r;
  }

  [[nodiscard]] double langevinInvDerivative(const double x)
  {
    if (x > 0.99999) {
      return 2.0e5 * x;
    }
    const double x2 = x * x;
    const double denom = (1.0 - x2) * (1.0 + (0.5 * x2));
    return (3.0 * (1.0 + (0.5 * x2) + (1.5 * x2 * x2))) / (denom * denom);
  }

} // namespace

double
langevin_inv(const double x)
{
  if (x < 0.0) {
    return -langevin_inv(-x);
  }
  if (x > 0.99999) {
    return 1e5 * (1.0 + (x * x));
  }
  const double x2 = x * x;
  return (3.0 * x) / ((1.0 - x2) * (1.0 + (0.5 * x2)));
}

double
csch(const double x)
{
  if (x == 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  return 1.0 / std::sinh(x);
}

Eigen::Matrix3d
MEHPForceEvaluator::evaluateSpringStiffness(const Eigen::Vector3d& springVec,
                                            const double contourLength) const
{
  validateContourLength(contourLength);
  const double h = 1e-6 * std::max(1.0, springVec.norm());
  Eigen::Matrix3d stiffness = Eigen::Matrix3d::Zero();
  for (Eigen::Index col = 0; col < 3; ++col) {
    Eigen::Vector3d drPlus = springVec;
    Eigen::Vector3d drMinus = springVec;
    drPlus[col] += h;
    drMinus[col] -= h;
    const Eigen::Vector3d fPlus =
      this->evaluateSpringForce(drPlus, contourLength);
    const Eigen::Vector3d fMinus =
      this->evaluateSpringForce(drMinus, contourLength);
    stiffness.col(col) = (fPlus - fMinus) / (2.0 * h);
  }
  return stiffness;
}

double
MEHPForceEvaluator::evaluateStressContribution(const Eigen::Vector3d& springVec,
                                               const size_t i,
                                               const size_t j,
                                               const double contourLength) const
{
  INVALIDARG_EXP_IFN(i < 3 && j < 3,
                     "Stress tensor indices must be in [0, 2].");
  validateContourLength(contourLength);
  const Eigen::Vector3d force =
    this->evaluateSpringForce(springVec, contourLength);
  return springVec[static_cast<Eigen::Index>(i)] *
         force[static_cast<Eigen::Index>(j)];
}

double
MEHPForceEvaluator::evaluateNetworkEnergyAndGradient(
  const MEHPNetwork& net,
  const Eigen::VectorXd& springVectors,
  const bool is2D,
  double* const grad) const
{
  validateNetworkAndSpringVectors(net, springVectors, grad);
  if (grad != nullptr) {
    std::fill(grad, grad + (3 * net.nrOfLinks), 0.0);
  }
  const Eigen::Index nrOfDim = is2D ? 2 : 3;
  double totalEnergy = 0.0;
  for (Eigen::Index s = 0; s < static_cast<Eigen::Index>(net.nrOfSprings);
       ++s) {
    const double contourLength = net.springContourLength[s];
    const Eigen::Vector3d springVec(springVectors[3 * s],
                                    springVectors[3 * s + 1],
                                    is2D ? 0.0 : springVectors[3 * s + 2]);
    totalEnergy += this->evaluateSpringEnergy(springVec, contourLength);
    if (grad != nullptr) {
      const Eigen::Vector3d force =
        this->evaluateSpringForce(springVec, contourLength);
      const size_t a = static_cast<size_t>(net.springIndexA[s]);
      const size_t b = static_cast<size_t>(net.springIndexB[s]);
      for (Eigen::Index d = 0; d < nrOfDim; ++d) {
        grad[3 * a + static_cast<size_t>(d)] -= force[d];
        grad[3 * b + static_cast<size_t>(d)] += force[d];
      }
    }
  }
  return totalEnergy;
}

Eigen::Matrix3d
MEHPForceEvaluator::evaluateStressTensor(const MEHPNetwork& net,
                                         const Eigen::VectorXd& springVectors,
                                         const double volume,
                                         const bool is2D) const
{
  INVALIDARG_EXP_IFN(std::isfinite(volume) && volume > 0.0,
                     "Volume must be finite and positive.");
  validateNetworkAndSpringVectors(net, springVectors, nullptr);
  Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();
  for (Eigen::Index s = 0; s < static_cast<Eigen::Index>(net.nrOfSprings);
       ++s) {
    const Eigen::Vector3d springVec(springVectors[3 * s],
                                    springVectors[3 * s + 1],
                                    is2D ? 0.0 : springVectors[3 * s + 2]);
    const Eigen::Vector3d force =
      this->evaluateSpringForce(springVec, net.springContourLength[s]);
    stress += springVec * force.transpose();
  }
  if (is2D) {
    stress.row(2).setZero();
    stress.col(2).setZero();
  }
  stress /= volume;
  return stress;
}

HookeanForceEvaluator::HookeanForceEvaluator(const double kappaVal)
  : kappa(kappaVal)
{
  INVALIDARG_EXP_IFN(std::isfinite(kappaVal) && kappaVal >= 0.0,
                     "Spring constant kappa must be finite and non-negative.");
}

double
HookeanForceEvaluator::evaluateSpringEnergy(const Eigen::Vector3d& springVec,
                                            const double contourLength) const
{
  validateContourLength(contourLength);
  return 0.5 * this->kappa * springVec.squaredNorm() / contourLength;
}

Eigen::Vector3d
HookeanForceEvaluator::evaluateSpringForce(const Eigen::Vector3d& springVec,
                                           const double contourLength) const
{
  validateContourLength(contourLength);
  return (this->kappa / contourLength) * springVec;
}

Eigen::Matrix3d
HookeanForceEvaluator::evaluateSpringStiffness(
  const Eigen::Vector3d& /* springVec */,
  const double contourLength) const
{
  validateContourLength(contourLength);
  return (this->kappa / contourLength) * Eigen::Matrix3d::Identity();
}

double
HookeanForceEvaluator::evaluateStressContribution(
  const Eigen::Vector3d& springVec,
  const size_t i,
  const size_t j,
  const double contourLength) const
{
  INVALIDARG_EXP_IFN(i < 3 && j < 3,
                     "Stress tensor indices must be in [0, 2].");
  validateContourLength(contourLength);
  return (this->kappa / contourLength) *
         springVec[static_cast<Eigen::Index>(i)] *
         springVec[static_cast<Eigen::Index>(j)];
}

double
HookeanForceEvaluator::evaluateNetworkEnergyAndGradient(
  const MEHPNetwork& net,
  const Eigen::VectorXd& springVectors,
  const bool is2D,
  double* const grad) const
{
  validateNetworkAndSpringVectors(net, springVectors, grad);
  const Eigen::Index nSprings = static_cast<Eigen::Index>(net.nrOfSprings);
  if (nSprings == 0) {
    if (grad != nullptr) {
      std::fill(grad, grad + (3 * net.nrOfLinks), 0.0);
    }
    return 0.0;
  }
  if (grad == nullptr) {
    const auto springMat = springVectors.reshaped(3, nSprings);
    if (is2D) {
      return 0.5 * this->kappa *
             (springMat.topRows<2>().colwise().squaredNorm().array() /
              net.springContourLength.transpose().array())
               .sum();
    }
    return 0.5 * this->kappa *
           (springMat.colwise().squaredNorm().array() /
            net.springContourLength.transpose().array())
             .sum();
  }
  std::fill(grad, grad + (3 * net.nrOfLinks), 0.0);
  double totalEnergy = 0.0;
  for (Eigen::Index s = 0; s < nSprings; ++s) {
    const double dx = springVectors[3 * s];
    const double dy = springVectors[3 * s + 1];
    const double dz = is2D ? 0.0 : springVectors[3 * s + 2];
    const double kOverN = this->kappa / net.springContourLength[s];
    totalEnergy += 0.5 * kOverN * ((dx * dx) + (dy * dy) + (dz * dz));
    if (grad != nullptr) {
      const double fx = kOverN * dx;
      const double fy = kOverN * dy;
      const size_t a = static_cast<size_t>(net.springIndexA[s]);
      const size_t b = static_cast<size_t>(net.springIndexB[s]);
      grad[3 * a] -= fx;
      grad[3 * a + 1] -= fy;
      grad[3 * b] += fx;
      grad[3 * b + 1] += fy;
      if (!is2D) {
        const double fz = kOverN * dz;
        grad[3 * a + 2] -= fz;
        grad[3 * b + 2] += fz;
      }
    }
  }
  return totalEnergy;
}

NonGaussianSpringForceEvaluator::NonGaussianSpringForceEvaluator(
  const double kappaVal,
  const double lVal)
  : kappa(kappaVal)
  , l(lVal)
  , oneOverl(lVal > 0.0 ? (1.0 / lVal) : 1.0)
{
  INVALIDARG_EXP_IFN(std::isfinite(kappaVal) && kappaVal >= 0.0,
                     "Spring constant kappa must be finite and non-negative.");
  INVALIDARG_EXP_IFN(std::isfinite(lVal) && lVal > 0.0,
                     "Kuhn segment length l must be finite and positive.");
}

double
NonGaussianSpringForceEvaluator::evaluateSpringEnergy(
  const Eigen::Vector3d& springVec,
  const double contourLength) const
{
  validateContourLength(contourLength);
  return computeNonGaussianSpringEnergyScalar(
    springVec.norm(), contourLength, this->oneOverl);
}

Eigen::Vector3d
NonGaussianSpringForceEvaluator::evaluateSpringForce(
  const Eigen::Vector3d& springVec,
  const double contourLength) const
{
  validateContourLength(contourLength);
  const double r = springVec.norm();
  if (r <= 0.0) {
    return Eigen::Vector3d::Zero();
  }
  return computeNonGaussianForceOverR(r, contourLength, this->oneOverl) *
         springVec;
}

Eigen::Matrix3d
NonGaussianSpringForceEvaluator::evaluateSpringStiffness(
  const Eigen::Vector3d& springVec,
  const double contourLength) const
{
  validateContourLength(contourLength);
  const double r = springVec.norm();
  if (r < 1e-12) {
    return (3.0 / (contourLength * this->l * this->l)) *
           Eigen::Matrix3d::Identity();
  }
  const double x = (r * this->oneOverl) / contourLength;
  const double fr = this->oneOverl * langevin_inv(x);
  const double dlinvDx = langevinInvDerivative(x);
  const double dfDr =
    ((this->oneOverl * this->oneOverl) / contourLength) * dlinvDx;
  const Eigen::Vector3d rHat = springVec / r;
  const Eigen::Matrix3d rHatOuter = rHat * rHat.transpose();
  return ((fr / r) * (Eigen::Matrix3d::Identity() - rHatOuter)) +
         (dfDr * rHatOuter);
}

double
NonGaussianSpringForceEvaluator::evaluateNetworkEnergyAndGradient(
  const MEHPNetwork& net,
  const Eigen::VectorXd& springVectors,
  const bool is2D,
  double* const grad) const
{
  validateNetworkAndSpringVectors(net, springVectors, grad);
  if (grad != nullptr) {
    std::fill(grad, grad + (3 * net.nrOfLinks), 0.0);
  }
  double totalEnergy = 0.0;
  for (Eigen::Index s = 0; s < static_cast<Eigen::Index>(net.nrOfSprings);
       ++s) {
    const double dx = springVectors[3 * s];
    const double dy = springVectors[3 * s + 1];
    const double dz = is2D ? 0.0 : springVectors[3 * s + 2];
    const double r = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
    const double contourLength = net.springContourLength[s];
    totalEnergy +=
      computeNonGaussianSpringEnergyScalar(r, contourLength, this->oneOverl);
    if (grad != nullptr && r > 0.0) {
      const double fOverR =
        computeNonGaussianForceOverR(r, contourLength, this->oneOverl);
      const double fx = fOverR * dx;
      const double fy = fOverR * dy;
      const size_t a = static_cast<size_t>(net.springIndexA[s]);
      const size_t b = static_cast<size_t>(net.springIndexB[s]);
      grad[3 * a] -= fx;
      grad[3 * a + 1] -= fy;
      grad[3 * b] += fx;
      grad[3 * b + 1] += fy;
      if (!is2D) {
        const double fz = fOverR * dz;
        grad[3 * a + 2] -= fz;
        grad[3 * b + 2] += fz;
      }
    }
  }
  return totalEnergy;
}

} // namespace pylimer_tools::sim::mehp
