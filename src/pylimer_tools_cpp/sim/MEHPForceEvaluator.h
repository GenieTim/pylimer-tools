#ifndef MEHP_FORCE_EVAL_H
#define MEHP_FORCE_EVAL_H

#include "./MEHPUtilityStructures.h"
#include "../utils/utilityMacros.h"
#include <Eigen/Dense>
#include <cstddef>

#ifdef CEREALIZABLE
#include "../utils/CerealUtils.h"
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

namespace pylimer_tools::sim::mehp {

struct MEHPNetwork;

/**
 * @brief Compute the Padé approximation of the inverse Langevin function.
 *
 * @param x Fractional extension \f$r / (N l)\f$.
 * @return Approximation to \f$\mathcal{L}^{-1}(x)\f$.
 */
[[nodiscard]] double
langevin_inv(double x);

/**
 * @brief Compute the hyperbolic cosecant \f$\operatorname{csch}(x) = 1 / \sinh(x)\f$.
 *
 * @param x Input argument.
 * @return \f$\operatorname{csch}(x)\f$, or \f$+\infty\f$ when @p x is 0.
 */
[[nodiscard]] double
csch(double x);

/**
 * @brief Stateless constitutive model for evaluating spring energies, forces,
 * stiffnesses, and virial stress tensors in MEHP networks.
 *
 * Geometry and sign convention:
 * - @p springVec is \f$\Delta\mathbf{r}_s = \mathbf{r}_B - \mathbf{r}_A + \mathbf{o}_s\f$.
 * - @ref evaluateSpringForce returns \f$\mathbf{F}_s = +\partial E_s / \partial(\Delta\mathbf{r}_s)\f$.
 * - Consequently, for link displacements \f$\mathbf{u}\f$,
 *   \f$\partial E_s / \partial\mathbf{u}_A = -\mathbf{F}_s\f$ and
 *   \f$\partial E_s / \partial\mathbf{u}_B = +\mathbf{F}_s\f$.
 * - Virial stress uses the positive convention
 *   \f$\boldsymbol{\sigma} = +\frac{1}{V} \sum_s \Delta\mathbf{r}_s \mathbf{F}_s^T\f$.
 */
class MEHPForceEvaluator
{
protected:
  MEHPForceEvaluator(const MEHPForceEvaluator&) = default;
  MEHPForceEvaluator& operator=(const MEHPForceEvaluator&) = default;
  MEHPForceEvaluator(MEHPForceEvaluator&&) noexcept = default;
  MEHPForceEvaluator& operator=(MEHPForceEvaluator&&) noexcept = default;

public:
  MEHPForceEvaluator() = default;
  virtual ~MEHPForceEvaluator() = default;

  /**
   * @brief Indicate whether the constitutive force law is linear in displacement.
   *
   * @return @c true for linear (Hookean) springs, @c false otherwise.
   */
  [[nodiscard]] virtual bool isLinear() const { return false; }

  /**
   * @brief Return the spring stiffness parameter \f$\kappa\f$.
   *
   * @return Stiffness multiplier \f$\kappa\f$.
   */
  [[nodiscard]] virtual double getKappa() const { return 1.0; }

  /**
   * @brief Evaluate the strain energy of a single spring.
   *
   * @param springVec Spring end-to-end vector \f$\Delta\mathbf{r}_s\f$.
   * @param contourLength Effective strand contour length \f$N_s > 0\f$.
   * @return Spring energy \f$E_s\f$.
   */
  [[nodiscard]] virtual double evaluateSpringEnergy(
    const Eigen::Vector3d& springVec,
    double contourLength) const = 0;

  /**
   * @brief Evaluate the spring force \f$\mathbf{F}_s = +\partial E_s / \partial(\Delta\mathbf{r}_s)\f$.
   *
   * @param springVec Spring end-to-end vector \f$\Delta\mathbf{r}_s\f$.
   * @param contourLength Effective strand contour length \f$N_s > 0\f$.
   * @return Force vector \f$\mathbf{F}_s\f$.
   */
  [[nodiscard]] virtual Eigen::Vector3d evaluateSpringForce(
    const Eigen::Vector3d& springVec,
    double contourLength) const = 0;

  /**
   * @brief Evaluate the \f$3 \times 3\f$ tangent stiffness matrix
   * \f$\partial \mathbf{F}_s / \partial(\Delta\mathbf{r}_s)\f$.
   *
   * Default implementation uses symmetric central finite differences of
   * @ref evaluateSpringForce with step \f$h = 10^{-6}\max(1, \|\Delta\mathbf{r}_s\|)\f$.
   *
   * @param springVec Spring end-to-end vector \f$\Delta\mathbf{r}_s\f$.
   * @param contourLength Effective strand contour length \f$N_s > 0\f$.
   * @return \f$3 \times 3\f$ stiffness matrix.
   */
  [[nodiscard]] virtual Eigen::Matrix3d evaluateSpringStiffness(
    const Eigen::Vector3d& springVec,
    double contourLength) const;

  /**
   * @brief Evaluate the \f$(i, j)\f$ virial stress contribution of a single spring.
   *
   * @param springVec Spring end-to-end vector \f$\Delta\mathbf{r}_s\f$.
   * @param i First Cartesian index in @c [0, 2].
   * @param j Second Cartesian index in @c [0, 2].
   * @param contourLength Effective strand contour length \f$N_s > 0\f$.
   * @return Stress contribution \f$\Delta r_{s,i} F_{s,j}\f$.
   * @throws std::invalid_argument if @p i >= 3, @p j >= 3, or @p contourLength <= 0.
   */
  [[nodiscard]] virtual double evaluateStressContribution(
    const Eigen::Vector3d& springVec,
    size_t i,
    size_t j,
    double contourLength) const;

  /**
   * @brief Evaluate total network elastic energy and optional link gradient.
   *
   * @param net Network topology and contour lengths.
   * @param springVectors Flattened spring vectors of size \f$3 \cdot \text{nrOfSprings}\f$.
   * @param is2D If @c true, treat the system as planar (ignore \f$z\f$-components).
   * @param grad Optional output array of size \f$3 \cdot \text{nrOfLinks}\f$; if non-null,
   * zeroed and populated with \f$\partial E / \partial \mathbf{u}\f$.
   * @return Total network elastic energy.
   * @throws std::invalid_argument on mismatched vector sizes, invalid contour lengths,
   * or out-of-range link indices.
   */
  [[nodiscard]] virtual double evaluateNetworkEnergyAndGradient(
    const MEHPNetwork& net,
    const Eigen::VectorXd& springVectors,
    bool is2D,
    double* grad) const;

  /**
   * @brief Evaluate the network virial stress tensor
   * \f$\boldsymbol{\sigma} = \frac{1}{V} \sum_s \Delta\mathbf{r}_s \mathbf{F}_s^T\f$.
   *
   * @param net Network topology and contour lengths.
   * @param springVectors Flattened spring vectors of size \f$3 \cdot \text{nrOfSprings}\f$.
   * @param volume Simulation box volume \f$V > 0\f$.
   * @param is2D If @c true, ignore \f$z\f$-components and zero row/column 2.
   * @return \f$3 \times 3\f$ virial stress tensor.
   * @throws std::invalid_argument if @p volume <= 0, sizes mismatch, or contour lengths are invalid.
   */
  [[nodiscard]] Eigen::Matrix3d evaluateStressTensor(
    const MEHPNetwork& net,
    const Eigen::VectorXd& springVectors,
    double volume,
    bool is2D) const;

#ifdef CEREALIZABLE
  friend class cereal::access;

  template<class Archive>
  void serialize(Archive& /* ar */)
  {
  }
#endif
};

/**
 * @brief Linear Gaussian (Hookean) spring constitutive evaluator:
 * \f$E_s = \frac{\kappa}{2 N_s} \|\Delta\mathbf{r}_s\|^2\f$.
 */
class HookeanForceEvaluator final : public MEHPForceEvaluator
{
private:
  double kappa = 1.0;

public:
  /**
   * @brief Construct a Hookean spring evaluator.
   *
   * @param kappa Spring stiffness prefactor (\f$\ge 0\f$).
   * @throws std::invalid_argument if @p kappa is negative or non-finite.
   */
  explicit HookeanForceEvaluator(double kappa = 1.0);

  [[nodiscard]] bool isLinear() const override { return true; }

  [[nodiscard]] double getKappa() const override { return this->kappa; }

  [[nodiscard]] double evaluateSpringEnergy(
    const Eigen::Vector3d& springVec,
    double contourLength) const override;

  [[nodiscard]] Eigen::Vector3d evaluateSpringForce(
    const Eigen::Vector3d& springVec,
    double contourLength) const override;

  [[nodiscard]] Eigen::Matrix3d evaluateSpringStiffness(
    const Eigen::Vector3d& springVec,
    double contourLength) const override;

  [[nodiscard]] double evaluateStressContribution(
    const Eigen::Vector3d& springVec,
    size_t i,
    size_t j,
    double contourLength) const override;

  [[nodiscard]] double evaluateNetworkEnergyAndGradient(
    const MEHPNetwork& net,
    const Eigen::VectorXd& springVectors,
    bool is2D,
    double* grad) const override;

#ifdef CEREALIZABLE
  friend class cereal::access;

  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(cereal::base_class<MEHPForceEvaluator>(this), this->kappa);
  }
#endif
};

/**
 * @brief Non-Gaussian freely-jointed-chain (inverse Langevin) spring evaluator.
 */
class NonGaussianSpringForceEvaluator final : public MEHPForceEvaluator
{
private:
  double kappa = 1.0;
  double l = 1.0;
  double oneOverl = 1.0;

public:
  /**
   * @brief Construct a non-Gaussian inverse-Langevin spring evaluator.
   *
   * @param kappa Stiffness parameter (\f$\ge 0\f$).
   * @param l Kuhn segment length (\f$> 0\f$).
   * @throws std::invalid_argument if parameters are out of domain or non-finite.
   */
  explicit NonGaussianSpringForceEvaluator(double kappa = 1.0, double l = 1.0);

  [[nodiscard]] double getKappa() const override { return this->kappa; }

  /**
   * @brief Return the Kuhn segment length \f$l\f$.
   *
   * @return Segment length \f$l\f$.
   */
  [[nodiscard]] double getSegmentLength() const { return this->l; }

  [[nodiscard]] double evaluateSpringEnergy(
    const Eigen::Vector3d& springVec,
    double contourLength) const override;

  [[nodiscard]] Eigen::Vector3d evaluateSpringForce(
    const Eigen::Vector3d& springVec,
    double contourLength) const override;

  [[nodiscard]] Eigen::Matrix3d evaluateSpringStiffness(
    const Eigen::Vector3d& springVec,
    double contourLength) const override;

  [[nodiscard]] double evaluateNetworkEnergyAndGradient(
    const MEHPNetwork& net,
    const Eigen::VectorXd& springVectors,
    bool is2D,
    double* grad) const override;

#ifdef CEREALIZABLE
  friend class cereal::access;

  template<class Archive>
  void serialize(Archive& ar)
  {
    ar(cereal::base_class<MEHPForceEvaluator>(this),
       this->kappa,
       this->l,
       this->oneOverl);
  }
#endif
};

} // namespace pylimer_tools::sim::mehp

#ifdef CEREALIZABLE
CEREAL_REGISTER_TYPE(pylimer_tools::sim::mehp::HookeanForceEvaluator);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
  pylimer_tools::sim::mehp::MEHPForceEvaluator,
  pylimer_tools::sim::mehp::HookeanForceEvaluator);

CEREAL_REGISTER_TYPE(pylimer_tools::sim::mehp::NonGaussianSpringForceEvaluator);
CEREAL_REGISTER_POLYMORPHIC_RELATION(
  pylimer_tools::sim::mehp::MEHPForceEvaluator,
  pylimer_tools::sim::mehp::NonGaussianSpringForceEvaluator);
#endif

#endif
