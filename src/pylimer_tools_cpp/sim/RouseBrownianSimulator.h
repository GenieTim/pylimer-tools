#ifndef ROUSE_BROWNIAN_SIMULATOR_H
#define ROUSE_BROWNIAN_SIMULATOR_H

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

namespace pylimer_tools::sim {

/**
 * @brief Available linear solver algorithms for implicit Rouse/Brownian
 * time-stepping.
 */
enum class RouseSolverType
{
  DIRECT_LDLT,
  DIRECT_LLT,
  CONJUGATE_GRADIENT
};

/**
 * @brief Output container for Rouse dynamics simulation runs.
 */
struct RouseSimulationResult
{
  Eigen::VectorXd times{};
  Eigen::VectorXd msd{};
  Eigen::VectorXd msdCenterOfMass{};
  Eigen::VectorXd endToEndCorrelation{};
  std::vector<Eigen::MatrixXd> trajectory{};
};

/**
 * @brief Implicit time-stepping Rouse and Brownian dynamics simulator for
 * polymer chains and general spring networks.
 *
 * Discretization:
 * Backward Euler: (I + dt * (k / gamma) * L) r_{t+dt} = r_t + sqrt(2 * D * dt)
 * * xi where L is the scalar graph Laplacian, D = (kB * T) / gamma, and
 * coordinates x, y, z advance simultaneously against the 3 RHS column vectors
 * using a single factorization of the shifted scalar matrix A = I + dt * (k /
 * gamma) * L.
 */
class RouseBrownianSimulator
{
public:
  /**
   * @brief Construct from edge lists (springFrom, springTo), friction gamma,
   *        spring constant k, temperature T, time step dt, and optional
   * numNodes.
   */
  RouseBrownianSimulator(const std::vector<size_t>& springFrom,
                         const std::vector<size_t>& springTo,
                         double gamma = 1.0,
                         double k = 1.0,
                         double T = 1.0,
                         double dt = 0.01,
                         size_t numNodes = 0);

  /**
   * @brief Construct from pre-assembled sparse scalar graph Laplacian matrix L.
   */
  explicit RouseBrownianSimulator(const Eigen::SparseMatrix<double>& L,
                                  double gamma = 1.0,
                                  double k = 1.0,
                                  double T = 1.0,
                                  double dt = 0.01);

  // Factory methods for canonical topologies
  static RouseBrownianSimulator createLinearChain(size_t numNodes,
                                                  double gamma = 1.0,
                                                  double k = 1.0,
                                                  double T = 1.0,
                                                  double dt = 0.01);

  static RouseBrownianSimulator createRing(size_t numNodes,
                                           double gamma = 1.0,
                                           double k = 1.0,
                                           double T = 1.0,
                                           double dt = 0.01);

  static RouseBrownianSimulator createStar(size_t numArms,
                                           size_t armLength = 1,
                                           double gamma = 1.0,
                                           double k = 1.0,
                                           double T = 1.0,
                                           double dt = 0.01);

  // Parameter getters
  size_t getNrOfNodes() const { return numNodes_; }
  double getGamma() const { return gamma_; }
  double getK() const { return k_; }
  double getTemperature() const { return T_; }
  double getDt() const { return dt_; }
  double getKB() const { return kB_; }
  double getDiffusionCoefficient() const { return (kB_ * T_) / gamma_; }
  RouseSolverType getSolverType() const { return solverType_; }

  // Parameter setters
  void setParameters(double gamma, double k, double T, double dt);
  void setGamma(double gamma);
  void setK(double k);
  void setTemperature(double T);
  void setDt(double dt);
  void setKB(double kB);
  void setSolverType(RouseSolverType type);
  void setSeed(uint64_t seed);

  // Position management (N x 3 matrix or 3N vector)
  void setPositions(const Eigen::MatrixXd& positions);
  void setPositionsFlattened(const Eigen::VectorXd& r, bool interleaved = true);
  const Eigen::MatrixXd& getPositions() const { return positions_; }
  Eigen::VectorXd getPositionsFlattened(bool interleaved = true) const;

  void initializeLinearPositions(double bondLength = 1.0, int axis = 0);
  void initializeRandomWalk(double bondLength = 1.0, uint64_t seed = 42);

  // Access system matrices
  const Eigen::SparseMatrix<double>& getLaplacian() const { return L_; }
  const Eigen::SparseMatrix<double>& getShiftedMatrix() const
  {
    ensureFactorized();
    return A_;
  }

  // Time-stepping methods
  void step(size_t numSteps = 1);
  void stepWithExternalForces(const Eigen::MatrixXd& forces,
                              size_t numSteps = 1);
  void stepDeterministic(size_t numSteps = 1);
  void stepDeterministicWithExternalForces(const Eigen::MatrixXd& forces,
                                           size_t numSteps = 1);

  // Physical observables at current state
  Eigen::Vector3d getCenterOfMass() const;
  Eigen::Vector3d getEndToEndVector(
    size_t nodeA = 0,
    size_t nodeB = static_cast<size_t>(-1)) const;
  double getEndToEndDistance(size_t nodeA = 0,
                             size_t nodeB = static_cast<size_t>(-1)) const;
  double getRadiusOfGyration() const;
  double getPotentialEnergy() const;

  // Run simulation and track dynamics
  RouseSimulationResult run(size_t totalSteps,
                            size_t sampleInterval = 1,
                            bool saveTrajectory = false,
                            size_t endToEndNodeA = 0,
                            size_t endToEndNodeB = static_cast<size_t>(-1));

  // Trajectory post-processing
  static Eigen::VectorXd computeMSD(
    const std::vector<Eigen::MatrixXd>& trajectory,
    bool timeAveraged = false);

  static Eigen::VectorXd computeCenterOfMassMSD(
    const std::vector<Eigen::MatrixXd>& trajectory,
    bool timeAveraged = false);

  static Eigen::VectorXd computeEndToEndAutocorrelation(
    const std::vector<Eigen::MatrixXd>& trajectory,
    size_t nodeA = 0,
    size_t nodeB = static_cast<size_t>(-1),
    bool timeAveraged = false);

private:
  size_t numNodes_{ 0 };
  double gamma_{ 1.0 };
  double k_{ 1.0 };
  double T_{ 1.0 };
  double dt_{ 0.01 };
  double kB_{ 1.0 };
  RouseSolverType solverType_{ RouseSolverType::DIRECT_LDLT };

  Eigen::SparseMatrix<double> L_{};
  Eigen::SparseMatrix<double> A_{};
  Eigen::MatrixXd positions_{};

  mutable Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldltSolver_{};
  mutable Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> lltSolver_{};
  mutable Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                                   Eigen::Lower | Eigen::Upper>
    cgSolver_{};
  mutable bool solverFactorized_{ false };

  std::mt19937_64 rng_{ 42 };
  std::normal_distribution<double> normalDist_{ 0.0, 1.0 };

  void initializeFromEdges(const std::vector<size_t>& springFrom,
                           const std::vector<size_t>& springTo,
                           size_t numNodes);
  void updateShiftedMatrixAndFactorization();
  void ensureFactorized() const;
  void solveShiftedRHS(const Eigen::MatrixXd& rhs, Eigen::MatrixXd& result);
};

} // namespace pylimer_tools::sim

#endif // ROUSE_BROWNIAN_SIMULATOR_H
