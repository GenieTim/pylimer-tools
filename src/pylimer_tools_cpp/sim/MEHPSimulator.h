#ifndef MEHP_SIMULATOR_H
#define MEHP_SIMULATOR_H

#include "../entities/Box.h"
#include "../entities/Universe.h"
#include "../topo/EntanglementDetector.h"
#include "../utils/utilityMacros.h"
#include "./MEHPForceEvaluator.h"
#include "./MEHPLinearSystemSolver.h"
#include "./MEHPNetwork.h"
#include "./MEHPNonlinearOptimizationSolver.h"
#include "./MEHPSlipLinkRelaxationSolver.h"
#include "./MEHPSolver.h"
#include "./MEHPUtilityStructures.h"
#include "./OutputSupportingSimulation.h"
#include <Eigen/Dense>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef CEREALIZABLE
#include "../utils/CerealUtils.h"
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/memory.hpp>
#endif

namespace pylimer_tools::sim::mehp {

/**
 * @brief Unified Mechanical Equilibrium of Heterogeneous Polymers (MEHP)
 * simulator coordinating a canonical @ref MEHPNetwork, a stateless
 * @ref MEHPForceEvaluator, and an interchangeable @ref MEHPSolver strategy.
 */
class MEHPSimulator final : public pylimer_tools::sim::OutputSupportingSimulation
{
private:
#ifdef CEREALIZABLE
  MEHPSimulator() = default;
  friend class cereal::access;
#endif

  MEHPNetwork network;
  Eigen::VectorXd currentDisplacements;
  std::shared_ptr<MEHPForceEvaluator> forceEvaluator =
    std::make_shared<HookeanForceEvaluator>(1.0);
  bool is2D = false;
  int crossLinkerType = 2;
  int entanglementType = 3;
  MEHPSolveResult lastSolveResult;
  bool simulationHasRun = false;

  void initializeAfterNetworkConstruction();

public:
  // ===========================================================================
  // Constructors & Factories
  // ===========================================================================

  explicit MEHPSimulator(
    const MEHPNetwork& net,
    std::shared_ptr<MEHPForceEvaluator> evaluator = nullptr,
    bool is2D = false);

  explicit MEHPSimulator(
    const pylimer_tools::entities::Universe& universe,
    const pylimer_tools::topo::entanglement_detection::AtomPairEntanglements&
      entanglements = {},
    int crossLinkerType = 2,
    EntanglementMode entanglementMode = EntanglementMode::MERGED_LINKS,
    std::shared_ptr<MEHPForceEvaluator> evaluator = nullptr,
    bool is2D = false,
    bool keep2FunctionalCrosslinkers = false,
    bool removePrimaryLoops = false,
    bool assumeBoxLargeEnough = false);

  MEHPSimulator(
    const pylimer_tools::entities::Box& box,
    const CrosslinkerUniverse& cu,
    const pylimer_tools::topo::entanglement_detection::StrandPairEntanglements&
      entanglements = {},
    int crossLinkerType = 2,
    EntanglementMode entanglementMode = EntanglementMode::MERGED_LINKS,
    std::shared_ptr<MEHPForceEvaluator> evaluator = nullptr,
    bool is2D = false,
    bool removeDanglingChains = false,
    bool includeFreeChains = true,
    bool removePrimaryLoops = false,
    bool assumeBoxLargeEnough = false);

  [[nodiscard]] static MEHPSimulator constructWithRandomEntanglements(
    const pylimer_tools::entities::Universe& universe,
    size_t nrOfEntanglementsToSample,
    double upperSamplingCutoff,
    double lowerSamplingCutoff = 0.0,
    size_t minNrOfEntanglements = 0,
    double sameStrandCutoff = 0.0,
    const std::string& seed = "",
    int crossLinkerType = 2,
    EntanglementMode entanglementMode = EntanglementMode::MERGED_LINKS,
    bool withCrosslinkers = true,
    bool filterDanglingAndSoluble = false,
    std::shared_ptr<MEHPForceEvaluator> evaluator = nullptr,
    bool is2D = false,
    bool keep2FunctionalCrosslinkers = false,
    bool removePrimaryLoops = false,
    bool assumeBoxLargeEnough = false);

  // ===========================================================================
  // Configuration & Mutators
  // ===========================================================================

  void setForceEvaluator(std::shared_ptr<MEHPForceEvaluator> evaluator);

  [[nodiscard]] std::shared_ptr<MEHPForceEvaluator> getForceEvaluator() const
  {
    return this->forceEvaluator;
  }

  void configMeanBondLength(double meanBondLength);
  void configMeanSquaredBondLength(double r02);
  void configAssumeBoxLargeEnough(bool assume);
  void configAssumeNetworkIsComplete(bool assume);
  void configRemovePrimaryLoops(bool remove);
  void setIs2D(bool newIs2D);

  [[nodiscard]] bool getIs2D() const { return this->is2D; }

  void deformTo(const pylimer_tools::entities::Box& targetBox);

  size_t addSlipLinks(
    const pylimer_tools::topo::entanglement_detection::StrandPairEntanglements&
      entanglements,
    bool clampAlpha = false);

  size_t addSlipLinks(const std::vector<size_t>& strandIdx1,
                      const std::vector<size_t>& strandIdx2,
                      const Eigen::VectorXd& x,
                      const Eigen::VectorXd& y,
                      const Eigen::VectorXd& z,
                      const std::vector<double>& alpha1,
                      const std::vector<double>& alpha2,
                      bool clampAlpha = false);

  size_t addSlipLinks(const std::vector<size_t>& strandIdx1,
                      const std::vector<size_t>& strandIdx2,
                      const Eigen::VectorXd& x,
                      const Eigen::VectorXd& y,
                      const Eigen::VectorXd& z,
                      bool clampAlpha = false);

  size_t randomlyAddSlipLinks(size_t nrOfEntanglementsToSample,
                              double cutoff = 2.0,
                              size_t minNrOfEntanglements = 0,
                              double sameStrandCutoff = 2.0,
                              const std::string& seed = "",
                              bool clampAlpha = false);

  size_t breakTooLongStrands(double breakingLength);
  size_t removeInactiveLinks(double tolerance = 1e-6);
  size_t unlinkBifunctionalLinks(bool entanglementsOnly = false);
  size_t simplifyNetwork(StructureSimplificationMode mode,
                         double tolerance = 1e-6,
                         double breakingLength = -1.0);

  // ===========================================================================
  // Unified Solver Execution
  // ===========================================================================

  MEHPSolveResult run(const MEHPSolver& solver,
                      const MEHPSolveCallbacks& callbacks = {});

  MEHPSolveResult run(const MEHPLinearSolverConfig& config,
                      const MEHPSolveCallbacks& callbacks = {});

  MEHPSolveResult run(const MEHPSlipLinkSolverConfig& config,
                      const MEHPSolveCallbacks& callbacks = {});

  MEHPSolveResult run(const MEHPNonlinearSolverConfig& config,
                      const MEHPSolveCallbacks& callbacks = {});

  MEHPSolveResult run();

  // ===========================================================================
  // Network & Diagnostic Getters
  // ===========================================================================

  [[nodiscard]] const MEHPNetwork& getNetwork() const { return this->network; }

  [[nodiscard]] const Eigen::VectorXd& getCurrentDisplacements() const
  {
    return this->currentDisplacements;
  }

  void setCurrentDisplacements(const Eigen::VectorXd& newDisplacements);

  [[nodiscard]] Eigen::VectorXd getInitialCoordinates() const
  {
    return this->network.coordinates;
  }

  [[nodiscard]] Eigen::VectorXd getSpringVectors() const
  {
    return this->network.evaluateSpringVectors(this->currentDisplacements,
                                               this->is2D);
  }

  [[nodiscard]] Eigen::VectorXd getSpringLengths() const
  {
    return this->network.evaluateSpringLengths(this->currentDisplacements,
                                               this->is2D);
  }

  [[nodiscard]] double getAverageSpringLength() const;

  [[nodiscard]] Eigen::VectorXd getSpringContourLength() const
  {
    return this->network.springContourLength;
  }

  [[nodiscard]] Eigen::VectorXd getSpringPartitions() const
  {
    return this->network.springPartitions;
  }

  [[nodiscard]] Eigen::VectorXd getStrandContourLength() const
  {
    return this->network.strandContourLength;
  }

  [[nodiscard]] double getAverageContourLength() const
  {
    return this->network.meanSpringContourLength;
  }

  [[nodiscard]] size_t getNrOfLinks() const { return this->network.nrOfLinks; }
  [[nodiscard]] size_t getNrOfNodes() const { return this->network.nrOfNodes; }
  [[nodiscard]] size_t getNrOfStrands() const
  {
    return this->network.nrOfStrands;
  }
  [[nodiscard]] size_t getNrOfSprings() const
  {
    return this->network.nrOfSprings;
  }

  [[nodiscard]] size_t getNrOfActiveNodes(double tolerance = 1e-3) const;
  [[nodiscard]] size_t getNrOfActiveSprings(double tolerance = 1e-3) const;
  [[nodiscard]] size_t getNrOfActiveStrands(double tolerance = 1e-3) const;

  [[nodiscard]] pylimer_tools::entities::Universe getCrosslinkerVerse() const
  {
    return this->network.getCrosslinkerVerse(this->currentDisplacements);
  }

  [[nodiscard]] double getDefaultR0Square() const
  {
    return this->network.defaultR02;
  }

  [[nodiscard]] double getDefaultBondLength() const
  {
    return this->network.defaultBondLength;
  }

  void validateNetwork() const
  {
    this->network.validate(this->currentDisplacements);
  }

  [[nodiscard]] const MEHPSolveResult& getLastSolveResult() const
  {
    return this->lastSolveResult;
  }

  [[nodiscard]] ExitReason getExitReason() const
  {
    return this->lastSolveResult.exitReason;
  }

  [[nodiscard]] size_t getNrOfIterations() const
  {
    return this->lastSolveResult.nrOfIterations;
  }

  [[nodiscard]] bool getSimulationHasRun() const
  {
    return this->simulationHasRun;
  }

  // ===========================================================================
  // Physical & Topological Observables
  // ===========================================================================

  [[nodiscard]] double getEnergy() const;
  [[nodiscard]] Eigen::VectorXd getResiduals() const;
  [[nodiscard]] double getResidualNorm() const;
  [[nodiscard]] double getPressure() const;

  [[nodiscard]] Eigen::VectorXd getGammaFactors(double b02 = -1.0) const
  {
    return this->network.getGammaFactors(
      this->currentDisplacements, b02, this->is2D);
  }

  [[nodiscard]] Eigen::VectorXd getGammaFactorsInDir(int dir,
                                                     double b02 = -1.0) const
  {
    return this->network.getGammaFactorsInDir(
      this->currentDisplacements, b02, dir, this->is2D);
  }

  [[nodiscard]] double getGammaFactor(double b02 = -1.0,
                                      int nrOfChains = -1) const
  {
    return this->network.getGammaFactor(
      this->currentDisplacements, b02, nrOfChains, this->is2D);
  }

  [[nodiscard]] double getSolubleWeightFraction(double tolerance = 1e-3) const
  {
    return this->network.computeSolubleWeightFraction(
      this->currentDisplacements, tolerance, this->is2D);
  }

  [[nodiscard]] double getDanglingWeightFraction(double tolerance = 1e-3) const
  {
    return this->network.computeDanglingWeightFraction(
      this->currentDisplacements, tolerance, this->is2D);
  }

  [[nodiscard]] double getActiveWeightFraction(double tolerance = 1e-3) const
  {
    return this->network.computeActiveWeightFraction(
      this->currentDisplacements, tolerance, this->is2D);
  }

  [[nodiscard]] std::unordered_map<long int, int>
  getEffectiveFunctionalityOfAtoms(double tolerance = 1e-3) const
  {
    return this->network.getEffectiveFunctionalityOfAtoms(
      this->currentDisplacements, tolerance, this->is2D);
  }

  [[nodiscard]] double countActiveClusteredAtoms(double tolerance = 1e-3) const
  {
    return this->network.countActiveClusteredAtoms(
      this->currentDisplacements, tolerance, this->is2D);
  }

  [[nodiscard]] std::vector<int> getIndicesOfActiveNodes(
    double tolerance = 1e-3) const
  {
    return this->network.getIndicesOfActiveNodes(
      this->currentDisplacements, tolerance, this->is2D);
  }

  [[nodiscard]] size_t getNumIntraChainSlipLinks() const
  {
    return this->network.getNumIntraChainSlipLinks();
  }

  // ===========================================================================
  // OutputSupportingSimulation Overrides
  // ===========================================================================

  double getCurrentTime(double currentStep) override;
  double getTemperature() override;
  double getGamma() override { return this->getGammaFactor(1.0, -1); }
  double getTimestep() override { return 1.0; }
  double getVolume() override { return this->network.vol; }
  double getResidual() override { return this->getResidualNorm(); }
  Eigen::Matrix3d getStressTensor() override;
  Eigen::VectorXd getBondLengths() override;
  Eigen::VectorXd getCoordinates() override;
  int getNumRelocations() override
  {
    return this->lastSolveResult.nrOfRelocations;
  }
  int getNumShifts() override { return this->lastSolveResult.nrOfShifts; }
  long int getNumBondsToForm() override { return 0; }
  size_t getNumAtoms() override { return this->network.nrOfNodes; }
  size_t getNumBonds() override { return this->network.nrOfStrands; }
  size_t getNumExtraAtoms() override
  {
    return this->network.nrOfLinks - this->network.nrOfNodes;
  }
  size_t getNumExtraBonds() override
  {
    return static_cast<size_t>(this->network.springIsEntanglement.count());
  }
  size_t getNumParticles() override { return this->network.nrOfNodes; }

#ifdef CEREALIZABLE
  template<class Archive>
  void serialize(Archive& ar, std::uint32_t const /* version */)
  {
    ar(cereal::virtual_base_class<OutputSupportingSimulation>(this));
    ar(network,
       currentDisplacements,
       forceEvaluator,
       is2D,
       crossLinkerType,
       entanglementType,
       lastSolveResult,
       simulationHasRun);
  }

  void writeRestartFile(std::string& filename) override
  {
    pylimer_tools::utils::serializeToFile<MEHPSimulator>(*this, filename);
  }

  [[nodiscard]] static MEHPSimulator readRestartFile(std::string& filename)
  {
    MEHPSimulator res;
    pylimer_tools::utils::deserializeFromFile(res, filename);
    return res;
  }

  [[nodiscard]] std::string serializeToString() const
  {
    return pylimer_tools::utils::serializeToString<MEHPSimulator>(*this);
  }

  [[nodiscard]] static MEHPSimulator constructFromString(const std::string& s)
  {
    MEHPSimulator res;
    std::string copy = s;
    pylimer_tools::utils::deserializeFromString(res, copy);
    return res;
  }
#endif
};

} // namespace pylimer_tools::sim::mehp

#ifdef CEREALIZABLE
CEREAL_CLASS_VERSION(pylimer_tools::sim::mehp::MEHPSimulator, 1);
#endif

#endif // MEHP_SIMULATOR_H
