#ifndef PYBIND_SIM_H
#define PYBIND_SIM_H

#include <functional>

#include "../entities/Universe.h"
#include "../sim/DPDSimulator.h"
#include "../sim/MEHPForceEvaluator.h"
#include "../sim/MEHPLinearSystemSolver.h"
#include "../sim/MEHPNetwork.h"
#include "../sim/MEHPNonlinearOptimizationSolver.h"
#include "../sim/MEHPSimulator.h"
#include "../sim/MEHPSlipLinkRelaxationSolver.h"
#include "../sim/MEHPSolver.h"
#include "../sim/RouseBrownianSimulator.h"

#include <algorithm>
#include <cassert>
#include <memory>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/native_enum.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;

using namespace pylimer_tools::sim;

namespace pylimer_tools::sim::mehp {
class PyMEHPForceEvaluator
  : public MEHPForceEvaluator
  , public py::trampoline_self_life_support
{
public:
  using MEHPForceEvaluator::MEHPForceEvaluator;

  ~PyMEHPForceEvaluator() override = default;

  [[nodiscard]] bool isLinear() const override
  {
    PYBIND11_OVERRIDE_NAME(bool, MEHPForceEvaluator, "is_linear", isLinear, );
  }

  [[nodiscard]] double getKappa() const override
  {
    PYBIND11_OVERRIDE_NAME(double, MEHPForceEvaluator, "get_kappa", getKappa, );
  }

  [[nodiscard]] double evaluateSpringEnergy(
    const Eigen::Vector3d& springVec,
    const double contourLength) const override
  {
    PYBIND11_OVERRIDE_PURE_NAME(double,
                                MEHPForceEvaluator,
                                "evaluate_spring_energy",
                                evaluateSpringEnergy,
                                springVec,
                                contourLength);
  }

  [[nodiscard]] Eigen::Vector3d evaluateSpringForce(
    const Eigen::Vector3d& springVec,
    const double contourLength) const override
  {
    PYBIND11_OVERRIDE_PURE_NAME(Eigen::Vector3d,
                                MEHPForceEvaluator,
                                "evaluate_spring_force",
                                evaluateSpringForce,
                                springVec,
                                contourLength);
  }

  [[nodiscard]] Eigen::Matrix3d evaluateSpringStiffness(
    const Eigen::Vector3d& springVec,
    const double contourLength) const override
  {
    PYBIND11_OVERRIDE_NAME(Eigen::Matrix3d,
                           MEHPForceEvaluator,
                           "evaluate_spring_stiffness",
                           evaluateSpringStiffness,
                           springVec,
                           contourLength);
  }

  [[nodiscard]] double evaluateStressContribution(
    const Eigen::Vector3d& springVec,
    const size_t i,
    const size_t j,
    const double contourLength) const override
  {
    PYBIND11_OVERRIDE_NAME(double,
                           MEHPForceEvaluator,
                           "evaluate_stress_contribution",
                           evaluateStressContribution,
                           springVec,
                           i,
                           j,
                           contourLength);
  }

  [[nodiscard]] double evaluateNetworkEnergyAndGradient(
    const MEHPNetwork& net,
    const Eigen::VectorXd& springVectors,
    const bool is2D,
    double* const grad) const override
  {
    py::gil_scoped_acquire gil;
    const py::function override =
      py::get_override(static_cast<const MEHPForceEvaluator*>(this),
                       "evaluate_network_energy_and_gradient");
    if (override) {
      const py::object result =
        override(std::cref(net), springVectors, is2D, grad != nullptr);
      if (!py::isinstance<py::tuple>(result)) {
        throw std::runtime_error(
          "evaluate_network_energy_and_gradient must return a tuple of "
          "(energy, gradient)");
      }
      const py::tuple t = result.cast<py::tuple>();
      if (t.size() != 2) {
        throw std::runtime_error(
          "evaluate_network_energy_and_gradient must return a tuple of "
          "(energy, gradient)");
      }
      const double energy = t[0].cast<double>();
      if (grad != nullptr) {
        if (t[1].is_none()) {
          throw std::runtime_error(
            "evaluate_network_energy_and_gradient returned None for gradient "
            "when compute_gradient=True");
        }
        const Eigen::VectorXd gradVec = t[1].cast<Eigen::VectorXd>();
        const size_t expectedSize = 3 * net.nrOfLinks;
        if (static_cast<size_t>(gradVec.size()) != expectedSize) {
          throw std::runtime_error("Gradient size mismatch: expected " +
                                   std::to_string(expectedSize) + ", got " +
                                   std::to_string(gradVec.size()));
        }
        for (size_t i = 0; i < expectedSize; ++i) {
          grad[i] = gradVec[static_cast<Eigen::Index>(i)];
        }
      }
      return energy;
    }
    return MEHPForceEvaluator::evaluateNetworkEnergyAndGradient(
      net, springVectors, is2D, grad);
  }
};
} // namespace pylimer_tools::sim::mehp

namespace {
[[nodiscard]] pylimer_tools::sim::mehp::MEHPSolveCallbacks
makePythonInterruptCallbacks()
{
  pylimer_tools::sim::mehp::MEHPSolveCallbacks callbacks;
  callbacks.shouldInterrupt = []() {
    py::gil_scoped_acquire gil;
    return PyErr_CheckSignals() != 0;
  };
  callbacks.cleanupInterrupt = []() {
    py::gil_scoped_acquire gil;
    throw py::error_already_set();
  };
  return callbacks;
}

void
writeBackNumpyVectorIfCompatible(const py::object& uObj,
                                 const Eigen::VectorXd& uVec)
{
  if (py::isinstance<py::array_t<double, py::array::c_style>>(uObj)) {
    py::array_t<double, py::array::c_style> arr =
      py::reinterpret_borrow<py::array_t<double, py::array::c_style>>(uObj);
    if (arr.writeable() && arr.ndim() == 1 && arr.size() == uVec.size()) {
      std::copy(uVec.data(), uVec.data() + uVec.size(), arr.mutable_data());
    }
  }
}

[[nodiscard]] pylimer_tools::sim::mehp::MEHPSolveResult
solveWithPythonDisplacements(
  const pylimer_tools::sim::mehp::MEHPSolver& solver,
  pylimer_tools::sim::mehp::MEHPNetwork& net,
  const py::object& uObj,
  const pylimer_tools::sim::mehp::MEHPForceEvaluator& forceEvaluator,
  const bool is2D)
{
  Eigen::VectorXd uVec = uObj.cast<Eigen::VectorXd>();
  const pylimer_tools::sim::mehp::MEHPSolveResult result = solver.solve(
    net, uVec, forceEvaluator, is2D, makePythonInterruptCallbacks());
  writeBackNumpyVectorIfCompatible(uObj, uVec);
  return result;
}

template<typename ClassType>
void
bindSlipLinkInsertionMethods(ClassType& cls,
                             const char* addName,
                             const char* randomAddName)
{
  cls
    .def(addName,
         py::overload_cast<const pylimer_tools::topo::entanglement_detection::
                             StrandPairEntanglements&,
                           bool>(&mehp::MEHPSimulator::addSlipLinks),
         py::arg("entanglements"),
         py::arg("clamp_alpha") = false)
    .def(addName,
         py::overload_cast<const std::vector<size_t>&,
                           const std::vector<size_t>&,
                           const Eigen::VectorXd&,
                           const Eigen::VectorXd&,
                           const Eigen::VectorXd&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           bool>(&mehp::MEHPSimulator::addSlipLinks),
         py::arg("strand_idx_1"),
         py::arg("strand_idx_2"),
         py::arg("x"),
         py::arg("y"),
         py::arg("z"),
         py::arg("alpha_1"),
         py::arg("alpha_2"),
         py::arg("clamp_alpha") = false)
    .def(addName,
         py::overload_cast<const std::vector<size_t>&,
                           const std::vector<size_t>&,
                           const Eigen::VectorXd&,
                           const Eigen::VectorXd&,
                           const Eigen::VectorXd&,
                           bool>(&mehp::MEHPSimulator::addSlipLinks),
         py::arg("strand_idx_1"),
         py::arg("strand_idx_2"),
         py::arg("x"),
         py::arg("y"),
         py::arg("z"),
         py::arg("clamp_alpha") = false)
    .def(randomAddName,
         &mehp::MEHPSimulator::randomlyAddSlipLinks,
         py::arg("nr_of_entanglements_to_sample"),
         py::arg("cutoff") = 2.0,
         py::arg("min_nr_of_entanglements") = 0,
         py::arg("same_strand_cutoff") = 2.0,
         py::arg("seed") = "",
         py::arg("clamp_alpha") = false);
}
} // namespace

void
init_pylimer_bound_sim(py::module_& m)
{
  ////////////////////////////////////////////////////////////////
  // MARK: Output Quantities

  py::native_enum<ComputedIntValues>(
    m, "ComputedIntValues", "enum.IntEnum", "Integer output quantities")
#define X(e, n)                                                                \
  .value(#e, ComputedIntValues::e, "Results in the output column " #n ".")
    COMPUTED_INT_VALUES
#undef X
      .finalize();

  py::native_enum<ComputedDoubleValues>(m,
                                        "ComputedDoubleValues",
                                        "enum.IntEnum",
                                        "Floating point output quantities")
#define X(e, n)                                                                \
  .value(#e, ComputedDoubleValues::e, "Results in the output column " #n ".")
    COMPUTED_DOUBLE_VALUES
#undef X
      .finalize();

  py::class_<OutputConfiguration, py::smart_holder>(m,
                                                    "OutputConfiguration",
                                                    py::module_local(),
                                                    R"pbdoc(
     A configuration object to configure the output values and frequency
     for simulation classes in this package.

     This class specifies which quantities to output and how often to write them
     during simulations.
    )pbdoc")
    .def(py::init<>(),
         R"pbdoc(
          Create a new OutputConfiguration instance.

          :return: A new OutputConfiguration object with default settings
         )pbdoc")
    .def_readwrite("int_values",
                   &OutputConfiguration::intValues,
                   R"pbdoc(
                    List of integer-valued quantities to output.
                    
                    Use :class:`~pylimer_tools_cpp.ComputedIntValues` enum to specify which integer quantities
                    should be computed and written to output.
                   )pbdoc")
    .def_readwrite("double_values",
                   &OutputConfiguration::doubleValues,
                   R"pbdoc(
                    List of double-valued quantities to output.
                    
                    Use :class:`~pylimer_tools_cpp.ComputedDoubleValues` enum to specify which floating-point quantities
                    should be computed and written to output.
                   )pbdoc")
    .def_readwrite("use_every",
                   &OutputConfiguration::useEvery,
                   R"pbdoc(
     For autocorrelation and averaging, how often to include values.

     Use a value of 1 to take average of or autocorrelate, respectively,
     all values encountered during the simulation or optimization procedure.
     )pbdoc")
    .def_readwrite("append",
                   &OutputConfiguration::append,
                   R"pbdoc(
     Whether to append to the file or truncate it
     )pbdoc")
    .def_readwrite("filename",
                   &OutputConfiguration::filename,
                   R"pbdoc(
      The path and name of the file to write to.
      An empty string ("") means standard output (console).
     )pbdoc")
    .def_readwrite("output_every",
                   &OutputConfiguration::outputEvery,
                   R"pbdoc(
     How often to write the values to the output.
     For averages, this value also says how many values will be averaged.
     )pbdoc");

  /**
   * ////////////////////////////////////////////////////////////////
   * MEHP
   * ////////////////////////////////////////////////////////////////
   */

  py::native_enum<mehp::ExitReason>(m,
                                    "ExitReason",
                                    "enum.IntEnum",
                                    R"pbdoc(
An enum representing the reason for exiting
the simulation or optimization procedure.)pbdoc")
#define X(e, n) .value(#e, mehp::ExitReason::e, "Exit reason: " #n ".")
    EXIT_REASONS
#undef X
      .finalize();

  m.def("inverse_langevin",
        &mehp::langevin_inv,
        R"pbdoc(
     A somewhat accurate (for :math:`x \in (-1, 1)`) implementation of the inverse Langevin function.

     Source: https://scicomp.stackexchange.com/a/30251

     :param x: Input value in the range (-1, 1)
     :return: Inverse Langevin function value
  )pbdoc",
        py::arg("x"));

  ////////////////////////////////////////////////////////////////
  // MARK: Configuration Enums
  py::native_enum<mehp::StructureSimplificationMode>(
    m,
    "StructureSimplificationMode",
    "enum.Enum",
    "How the structure shall be simplified during the optimization in order to "
    "remove non-trapped entanglement links.")
#define X(e, name)                                                             \
  .value(STRINGINFY(e), mehp::StructureSimplificationMode::e, name)
    STRUCTURE_SIMPLIFICATION_MODES
#undef X
      .finalize();

  py::native_enum<mehp::LinkSwappingMode>(
    m,
    "LinkSwappingMode",
    "enum.Enum",
    "How slip-links may act when they reach each-other or even a crosslink.")
#define X(e, name) .value(STRINGINFY(e), mehp::LinkSwappingMode::e, name)
    LINK_SWAPPING_MODES
#undef X
      .finalize();

  py::native_enum<mehp::SLESolver>(
    m, "SLESolver", "enum.Enum", "Solver for sparse linear equation systems")
#define X(e, name) .value(STRINGINFY(e), mehp::SLESolver::e, name)
    SLE_SOLVERS
#undef X
      .finalize();

  py::native_enum<mehp::EntanglementMode>(
    m,
    "EntanglementMode",
    "enum.Enum",
    "Representation of topological entanglements in an MEHP network.")
    .value("NONE",
           mehp::EntanglementMode::NONE,
           "Phantom network (crosslinks only).")
    .value(
      "MERGED_LINKS",
      mehp::EntanglementMode::MERGED_LINKS,
      "Entanglement links lie directly on strand contour (no extra spring).")
    .value("SPRING",
           mehp::EntanglementMode::SPRING,
           "Entangled beads coupled by compliant Hookean spring.")
    .finalize();

  py::native_enum<mehp::SlipMode>(
    m,
    "SlipMode",
    "enum.Enum",
    "Mobility of entanglement links along polymer strands.")
    .value("FIXED",
           mehp::SlipMode::FIXED,
           "Entanglement links remain fixed at initial strand positions.")
    .value("SLIDING",
           mehp::SlipMode::SLIDING,
           "Entanglement links can slide along polymer strands.")
    .finalize();

  py::native_enum<mehp::NLoptAlgorithm>(
    m,
    "NLoptAlgorithm",
    "enum.Enum",
    "Supported NLopt optimization algorithms for MEHP force relaxation.")
    .value("LD_MMA", mehp::NLoptAlgorithm::LD_MMA, "NLopt LD_MMA algorithm.")
    .value(
      "LD_LBFGS", mehp::NLoptAlgorithm::LD_LBFGS, "NLopt LD_LBFGS algorithm.")
    .value(
      "LD_SLSQP", mehp::NLoptAlgorithm::LD_SLSQP, "NLopt LD_SLSQP algorithm.")
    .value("LD_TNEWTON",
           mehp::NLoptAlgorithm::LD_TNEWTON,
           "NLopt LD_TNEWTON algorithm.")
    .value("LD_TNEWTON_RESTART",
           mehp::NLoptAlgorithm::LD_TNEWTON_RESTART,
           "NLopt LD_TNEWTON_RESTART algorithm.")
    .value("LD_TNEWTON_PRECOND",
           mehp::NLoptAlgorithm::LD_TNEWTON_PRECOND,
           "NLopt LD_TNEWTON_PRECOND algorithm.")
    .value("LD_TNEWTON_PRECOND_RESTART",
           mehp::NLoptAlgorithm::LD_TNEWTON_PRECOND_RESTART,
           "NLopt LD_TNEWTON_PRECOND_RESTART algorithm.")
    .value("LD_VAR1", mehp::NLoptAlgorithm::LD_VAR1, "NLopt LD_VAR1 algorithm.")
    .value("LD_VAR2", mehp::NLoptAlgorithm::LD_VAR2, "NLopt LD_VAR2 algorithm.")
    .value(
      "LN_SBPLX", mehp::NLoptAlgorithm::LN_SBPLX, "NLopt LN_SBPLX algorithm.")
    .value("LN_BOBYQA",
           mehp::NLoptAlgorithm::LN_BOBYQA,
           "NLopt LN_BOBYQA algorithm.")
    .value("LN_NELDERMEAD",
           mehp::NLoptAlgorithm::LN_NELDERMEAD,
           "NLopt LN_NELDERMEAD algorithm.")
    .value("LN_COBYLA",
           mehp::NLoptAlgorithm::LN_COBYLA,
           "NLopt LN_COBYLA algorithm.")
    .value("LN_NEWUOA",
           mehp::NLoptAlgorithm::LN_NEWUOA,
           "NLopt LN_NEWUOA algorithm.")
    .value("LN_NEWUOA_BOUND",
           mehp::NLoptAlgorithm::LN_NEWUOA_BOUND,
           "NLopt LN_NEWUOA_BOUND algorithm.")
    .value("LN_PRAXIS",
           mehp::NLoptAlgorithm::LN_PRAXIS,
           "NLopt LN_PRAXIS algorithm.")
    .value("GN_DIRECT",
           mehp::NLoptAlgorithm::GN_DIRECT,
           "NLopt GN_DIRECT algorithm.")
    .value("GN_CRS2_LM",
           mehp::NLoptAlgorithm::GN_CRS2_LM,
           "NLopt GN_CRS2_LM algorithm.")
    .value(
      "GD_STOGO", mehp::NLoptAlgorithm::GD_STOGO, "NLopt GD_STOGO algorithm.")
    .finalize();

  ////////////////////////////////////////////////////////////////
  // MARK: Canonical MEHP Network
  py::class_<mehp::MEHPNetwork, py::smart_holder>(m,
                                                  "MEHPNetwork",
                                                  R"pbdoc(
      Canonical polymer network representation shared across MEHP linear-system,
      slip-link relaxation, and nonlinear optimization solvers.
    )pbdoc")
    .def(py::init<>(), "Construct an empty MEHPNetwork.")
    .def_static(
      "from_universe",
      &mehp::MEHPNetwork::fromUniverse,
      R"pbdoc(
        Build a canonical MEHPNetwork from an atomistic Universe using DFS strand contraction.

        :param universe: Atomistic polymer Universe.
        :param entanglements: Optional pre-sampled atom-pair entanglements.
        :param crosslinker_type: Atom type identifying crosslinker beads.
        :param entanglement_mode: Entanglement representation mode.
        :param keep_2functional_crosslinkers: Keep degree-2 crosslinkers as explicit network nodes.
        :param remove_primary_loops: Remove polymer primary loops upon construction.
        :param assume_box_large_enough: Whether springs can be evaluated with minimum-image PBC without baked offsets.
        :return: Constructed MEHPNetwork.
      )pbdoc",
      py::arg("universe"),
      py::arg_v(
        "entanglements",
        pylimer_tools::topo::entanglement_detection::AtomPairEntanglements(),
        "AtomPairEntanglements()"),
      py::arg("crosslinker_type") = 2,
      py::arg_v("entanglement_mode",
                mehp::EntanglementMode::MERGED_LINKS,
                "EntanglementMode.MERGED_LINKS"),
      py::arg("keep_2functional_crosslinkers") = false,
      py::arg("remove_primary_loops") = false,
      py::arg("assume_box_large_enough") = false)
    .def_static("from_universe_with_random_entanglements",
                &mehp::MEHPNetwork::fromUniverseWithRandomEntanglements,
                R"pbdoc(
        Sample random entanglements from an atomistic Universe and build a canonical MEHPNetwork.

        :param universe: Atomistic polymer Universe.
        :param nr_of_entanglements_to_sample: Target number of entanglements to sample.
        :param upper_sampling_cutoff: Upper distance cutoff for pair sampling.
        :param lower_sampling_cutoff: Lower distance cutoff for pair sampling.
        :param min_nr_of_entanglements: Minimum entanglements required if possible.
        :param same_strand_cutoff: Minimum chemical distance for same-strand pairs.
        :param seed: PRNG seed string (empty for random seed).
        :param crosslinker_type: Atom type identifying crosslinker beads.
        :param entanglement_mode: Entanglement representation mode.
        :param with_crosslinkers: Exclude crosslinker beads as entanglement targets when True.
        :param filter_dangling_and_soluble: Omit entanglements touching dangling/soluble strands during sampling.
        :param keep_2functional_crosslinkers: Keep degree-2 crosslinkers as explicit nodes.
        :param remove_primary_loops: Remove polymer primary loops upon construction.
        :param assume_box_large_enough: Assume box large enough for minimum-image PBC.
        :return: Constructed MEHPNetwork.
      )pbdoc",
                py::arg("universe"),
                py::arg("nr_of_entanglements_to_sample"),
                py::arg("upper_sampling_cutoff"),
                py::arg("lower_sampling_cutoff") = 0.0,
                py::arg("min_nr_of_entanglements") = 0,
                py::arg("same_strand_cutoff") = 0.0,
                py::arg("seed") = "",
                py::arg("crosslinker_type") = 2,
                py::arg_v("entanglement_mode",
                          mehp::EntanglementMode::MERGED_LINKS,
                          "EntanglementMode.MERGED_LINKS"),
                py::arg("with_crosslinkers") = true,
                py::arg("filter_dangling_and_soluble") = false,
                py::arg("keep_2functional_crosslinkers") = false,
                py::arg("remove_primary_loops") = false,
                py::arg("assume_box_large_enough") = false)
    .def_static(
      "from_crosslinker_universe",
      &mehp::MEHPNetwork::fromCrosslinkerUniverse,
      R"pbdoc(
        Build a canonical MEHPNetwork from a coarse-grained CrosslinkerUniverse and optional StrandPairEntanglements.

        :param box: Simulation box.
        :param crosslinker_universe: Crosslinker & strand representation.
        :param entanglements: Optional pre-sampled strand-pair entanglements.
        :param crosslinker_type: Atom type identifying crosslinkers.
        :param entanglement_mode: Entanglement representation mode.
        :param remove_dangling_chains: Omit dangling chains when True.
        :param include_free_chains: Include free/soluble chains when True.
        :param remove_primary_loops: Remove polymer primary loops upon construction.
        :param assume_box_large_enough: Assume box large enough for minimum-image PBC.
        :return: Constructed MEHPNetwork.
      )pbdoc",
      py::arg("box"),
      py::arg("crosslinker_universe"),
      py::arg_v(
        "entanglements",
        pylimer_tools::topo::entanglement_detection::StrandPairEntanglements(),
        "StrandPairEntanglements()"),
      py::arg("crosslinker_type") = 2,
      py::arg_v("entanglement_mode",
                mehp::EntanglementMode::MERGED_LINKS,
                "EntanglementMode.MERGED_LINKS"),
      py::arg("remove_dangling_chains") = false,
      py::arg("include_free_chains") = true,
      py::arg("remove_primary_loops") = false,
      py::arg("assume_box_large_enough") = false)
    .def_readwrite("nr_of_links", &mehp::MEHPNetwork::nrOfLinks)
    .def_readwrite("nr_of_nodes", &mehp::MEHPNetwork::nrOfNodes)
    .def_readwrite("nr_of_strands", &mehp::MEHPNetwork::nrOfStrands)
    .def_readwrite("nr_of_springs", &mehp::MEHPNetwork::nrOfSprings)
    .def_readwrite("L", &mehp::MEHPNetwork::L)
    .def_readwrite("vol", &mehp::MEHPNetwork::vol)
    .def_readonly("box", &mehp::MEHPNetwork::box)
    .def_readwrite("coordinates", &mehp::MEHPNetwork::coordinates)
    .def_readwrite("spring_contour_length",
                   &mehp::MEHPNetwork::springContourLength)
    .def_readwrite("strand_contour_length",
                   &mehp::MEHPNetwork::strandContourLength)
    .def_readwrite("spring_partitions", &mehp::MEHPNetwork::springPartitions)
    .def_readwrite("spring_index_a", &mehp::MEHPNetwork::springIndexA)
    .def_readwrite("spring_index_b", &mehp::MEHPNetwork::springIndexB)
    .def_readwrite("spring_box_offset", &mehp::MEHPNetwork::springBoxOffset)
    .def_readwrite("strand_index_of_spring",
                   &mehp::MEHPNetwork::strandIndexOfSpring)
    .def_readwrite("link_is_entanglement",
                   &mehp::MEHPNetwork::linkIsEntanglement)
    .def_readwrite("spring_is_entanglement",
                   &mehp::MEHPNetwork::springIsEntanglement)
    .def_readwrite("old_atom_ids", &mehp::MEHPNetwork::oldAtomIds)
    .def_readwrite("old_atom_types", &mehp::MEHPNetwork::oldAtomTypes)
    .def_readwrite("assume_box_large_enough",
                   &mehp::MEHPNetwork::assumeBoxLargeEnough)
    .def_readwrite("assume_complete", &mehp::MEHPNetwork::assumeComplete)
    .def_readwrite("remove_primary_loops",
                   &mehp::MEHPNetwork::removePrimaryLoops)
    .def_readwrite("mean_spring_contour_length",
                   &mehp::MEHPNetwork::meanSpringContourLength)
    .def_readwrite("mean_strand_contour_length",
                   &mehp::MEHPNetwork::meanStrandContourLength)
    .def_readwrite("default_bond_length", &mehp::MEHPNetwork::defaultBondLength)
    .def_readwrite("default_r02", &mehp::MEHPNetwork::defaultR02)
    .def_readwrite("nr_of_crosslink_swaps_endured",
                   &mehp::MEHPNetwork::nrOfCrosslinkSwapsEndured)
    .def_readwrite("link_indices_of_strand",
                   &mehp::MEHPNetwork::linkIndicesOfStrand)
    .def_readwrite("spring_indices_of_strand",
                   &mehp::MEHPNetwork::springIndicesOfStrand)
    .def_readwrite("strand_indices_of_link",
                   &mehp::MEHPNetwork::strandIndicesOfLink)
    .def("validate",
         &mehp::MEHPNetwork::validate,
         "Validate network dimensions, index mappings, and contour invariants.",
         py::arg("u") = Eigen::VectorXd())
    .def("deform_to",
         py::overload_cast<const pe::Box&>(&mehp::MEHPNetwork::deformTo),
         "Affinely deform network coordinates and baked box offsets to "
         "target_box.",
         py::arg("target_box"))
    .def(
      "deform_to",
      [](mehp::MEHPNetwork& self,
         const pe::Box& targetBox,
         const py::object& uObj) {
        Eigen::VectorXd uVec = uObj.cast<Eigen::VectorXd>();
        self.deformTo(targetBox, uVec);
        writeBackNumpyVectorIfCompatible(uObj, uVec);
      },
      "Affinely deform network coordinates, displacements u, and baked box "
      "offsets to target_box.",
      py::arg("target_box"),
      py::arg("u"))
    .def("evaluate_spring_vectors",
         &mehp::MEHPNetwork::evaluateSpringVectors,
         "Evaluate all flattened spring end-to-end vectors (size 3 * "
         "nr_of_springs).",
         py::arg("u"),
         py::arg("is_2d") = false)
    .def("evaluate_spring_lengths",
         &mehp::MEHPNetwork::evaluateSpringLengths,
         "Evaluate Euclidean lengths of all springs.",
         py::arg("u"),
         py::arg("is_2d") = false)
    .def("find_active_springs",
         &mehp::MEHPNetwork::findActiveSprings,
         "Return boolean mask of active springs.",
         py::arg("u"),
         py::arg("tolerance") = 1e-3,
         py::arg("is_2d") = false)
    .def("find_active_strands",
         &mehp::MEHPNetwork::findActiveStrands,
         "Return boolean mask of active strands.",
         py::arg("u"),
         py::arg("tolerance") = 1e-3,
         py::arg("is_2d") = false)
    .def("compute_active_weight_fraction",
         &mehp::MEHPNetwork::computeActiveWeightFraction,
         "Compute the weight fraction of active strands.",
         py::arg("u"),
         py::arg("tolerance") = 1e-3,
         py::arg("is_2d") = false)
    .def("compute_soluble_weight_fraction",
         &mehp::MEHPNetwork::computeSolubleWeightFraction,
         "Compute the weight fraction of soluble strands.",
         py::arg("u"),
         py::arg("tolerance") = 1e-3,
         py::arg("is_2d") = false)
    .def("compute_dangling_weight_fraction",
         &mehp::MEHPNetwork::computeDanglingWeightFraction,
         "Compute the weight fraction of dangling strands.",
         py::arg("u"),
         py::arg("tolerance") = 1e-3,
         py::arg("is_2d") = false)
    .def("get_gamma_factor",
         &mehp::MEHPNetwork::getGammaFactor,
         R"pbdoc(
Compute the average topological Gamma factor ``sum(gamma_factors) / denom``.

:param u: Displacement vector.
:param b02: Mean-squared bond length :math:`b_0^2` (defaults to network bond length squared).
:param nr_of_chains: Denominator for averaging. When ``< 1`` (default ``-1``), divides by
                     ``nr_of_strands`` — note that in ``EntanglementMode.SPRING``,
                     ``nr_of_strands`` includes entanglement coupling springs. Pass ``1`` (or
                     sum ``get_gamma_factors``) for the unnormalized sum ``sum(gamma_factors)``.
:param is_2d: Whether to evaluate in 2D.
)pbdoc",
         py::arg("u"),
         py::arg("b02") = -1.0,
         py::arg("nr_of_chains") = -1,
         py::arg("is_2d") = false)
    .def("get_crosslinker_verse",
         &mehp::MEHPNetwork::getCrosslinkerVerse,
         "Reconstruct a coarse-grained Universe of crosslinker nodes bonded by "
         "polymer strands.",
         py::arg("u") = Eigen::VectorXd())
    .def(
      "simplify",
      [](mehp::MEHPNetwork& self,
         const py::object& uObj,
         const mehp::StructureSimplificationMode mode,
         const double tolerance,
         const bool is2D,
         const double breakingLength) -> size_t {
        Eigen::VectorXd uVec = uObj.cast<Eigen::VectorXd>();
        const size_t removed =
          self.simplify(uVec, mode, tolerance, is2D, breakingLength);
        writeBackNumpyVectorIfCompatible(uObj, uVec);
        return removed;
      },
      "Run iterative topology simplification until convergence.",
      py::arg("u"),
      py::arg("mode"),
      py::arg("tolerance") = 1e-9,
      py::arg("is_2d") = false,
      py::arg("breaking_length") = -1.0)
#ifdef CEREALIZABLE
    .def(py::pickle(
      [](const mehp::MEHPNetwork& net) {
        return py::make_tuple(
          py::bytes(pylimer_tools::utils::serializeToString(net)));
      },
      [](const py::tuple& t) {
        if (t.size() != 1) {
          throw std::runtime_error("Invalid pickle state for MEHPNetwork.");
        }
        std::string in = t[0].cast<std::string>();
        mehp::MEHPNetwork net;
        pylimer_tools::utils::deserializeFromString(net, in);
        return net;
      }))
#endif
    ;

  ////////////////////////////////////////////////////////////////
  // MARK: Force Evaluators
  py::class_<mehp::MEHPForceEvaluator,
             mehp::PyMEHPForceEvaluator,
             py::smart_holder>(m,
                               "MEHPForceEvaluator",
                               R"pbdoc(
      Stateless constitutive model interface for evaluating spring energies, forces,
      stiffnesses, and virial stress tensors in MEHP networks.

      To create a custom force evaluator in Python, subclass this class and implement:

      - :meth:`evaluate_spring_energy`: Compute strain energy for one spring
      - :meth:`evaluate_spring_force`: Compute force vector for one spring
    )pbdoc")
    .def(py::init<>())
    .def(
      "is_linear",
      &mehp::MEHPForceEvaluator::isLinear,
      "Indicate whether the constitutive force law is linear in displacement.")
    .def("get_kappa",
         &mehp::MEHPForceEvaluator::getKappa,
         "Return the spring stiffness parameter kappa.")
    .def("evaluate_spring_energy",
         &mehp::MEHPForceEvaluator::evaluateSpringEnergy,
         "Evaluate strain energy of a single spring.",
         py::arg("spring_vec"),
         py::arg("contour_length"))
    .def("evaluate_spring_force",
         &mehp::MEHPForceEvaluator::evaluateSpringForce,
         "Evaluate spring force vector dE/d(spring_vec) for a single spring.",
         py::arg("spring_vec"),
         py::arg("contour_length"))
    .def("evaluate_spring_stiffness",
         &mehp::MEHPForceEvaluator::evaluateSpringStiffness,
         "Evaluate 3x3 tangent stiffness matrix dF/d(spring_vec) for a single "
         "spring.",
         py::arg("spring_vec"),
         py::arg("contour_length"))
    .def("evaluate_stress_contribution",
         &mehp::MEHPForceEvaluator::evaluateStressContribution,
         "Evaluate the (i, j) virial stress contribution of a single spring.",
         py::arg("spring_vec"),
         py::arg("i"),
         py::arg("j"),
         py::arg("contour_length"))
    .def(
      "evaluate_network_energy_and_gradient",
      [](const mehp::MEHPForceEvaluator& self,
         const mehp::MEHPNetwork& net,
         const Eigen::VectorXd& springVectors,
         const bool is2D,
         const bool computeGradient) -> py::tuple {
        if (computeGradient) {
          Eigen::VectorXd gradVec =
            Eigen::VectorXd::Zero(static_cast<Eigen::Index>(3 * net.nrOfLinks));
          const double energy = self.evaluateNetworkEnergyAndGradient(
            net, springVectors, is2D, gradVec.data());
          return py::make_tuple(energy, gradVec);
        }
        const double energy = self.evaluateNetworkEnergyAndGradient(
          net, springVectors, is2D, nullptr);
        return py::make_tuple(energy, py::none());
      },
      R"pbdoc(
        Evaluate total network elastic energy and optional link gradient.

        :param net: Network topology and contour lengths.
        :param spring_vectors: Flattened spring vectors (size 3 * nr_of_springs).
        :param is_2d: If True, treat the system as planar (ignore z components).
        :param compute_gradient: If True, compute and return link gradient ndarray.
        :return: Tuple (energy, gradient_or_none).
      )pbdoc",
      py::arg("net"),
      py::arg("spring_vectors"),
      py::arg("is_2d") = false,
      py::arg("compute_gradient") = true)
    .def("evaluate_stress_tensor",
         &mehp::MEHPForceEvaluator::evaluateStressTensor,
         "Evaluate the 3x3 virial stress tensor for the network.",
         py::arg("net"),
         py::arg("spring_vectors"),
         py::arg("volume"),
         py::arg("is_2d") = false);

  py::class_<mehp::HookeanForceEvaluator,
             mehp::MEHPForceEvaluator,
             py::smart_holder>(m,
                               "HookeanForceEvaluator",
                               R"pbdoc(
      Linear Gaussian (Hookean) spring constitutive evaluator:
      :math:`E_s = \frac{\kappa}{2 N_s} \|\Delta\mathbf{r}_s\|^2`.

      :param kappa: Spring stiffness prefactor (>= 0).
    )pbdoc")
    .def(py::init<double>(), py::arg("kappa") = 1.0);

  py::class_<mehp::NonGaussianSpringForceEvaluator,
             mehp::MEHPForceEvaluator,
             py::smart_holder>(m,
                               "NonGaussianSpringForceEvaluator",
                               R"pbdoc(
      Non-Gaussian freely-jointed-chain (inverse Langevin) spring evaluator.

      :param kappa: Stiffness parameter (>= 0).
      :param l: Kuhn segment length (> 0).
    )pbdoc")
    .def(py::init<double, double>(), py::arg("kappa") = 1.0, py::arg("l") = 1.0)
    .def_property_readonly(
      "l",
      &mehp::NonGaussianSpringForceEvaluator::getSegmentLength,
      "Kuhn segment length l.")
    .def("get_segment_length",
         &mehp::NonGaussianSpringForceEvaluator::getSegmentLength,
         "Return the Kuhn segment length l.");

  ////////////////////////////////////////////////////////////////
  // MARK: Solver Configs, Result, and Solvers
  py::class_<mehp::MEHPLinearSolverConfig, py::smart_holder>(
    m,
    "MEHPLinearSolverConfig",
    "Configuration for the sparse linear-system MEHP equilibrium solver.")
    .def(py::init([](const mehp::SLESolver solver,
                     const mehp::StructureSimplificationMode simplificationMode,
                     const double inactiveRemovalCutoff,
                     const double residualReduction,
                     const int maxIterations,
                     const double springBreakingLength,
                     const double inPlaceRebuildThreshold,
                     const int cgFallbackMaxIterations,
                     const double gradientDescentLearningRate,
                     const mehp::SlipMode slipMode) {
           return mehp::MEHPLinearSolverConfig{ solver,
                                                simplificationMode,
                                                inactiveRemovalCutoff,
                                                residualReduction,
                                                maxIterations,
                                                springBreakingLength,
                                                inPlaceRebuildThreshold,
                                                cgFallbackMaxIterations,
                                                gradientDescentLearningRate,
                                                slipMode };
         }),
         py::arg_v("solver", mehp::SLESolver::DEFAULT, "SLESolver.DEFAULT"),
         py::arg_v("simplification_mode",
                   mehp::StructureSimplificationMode::NO_SIMPLIFICATION,
                   "StructureSimplificationMode.NO_SIMPLIFICATION"),
         py::arg("inactive_removal_cutoff") = 1e-6,
         py::arg("residual_reduction") = 1e-6,
         py::arg("max_iterations") = 50000,
         py::arg("spring_breaking_length") = -1.0,
         py::arg("in_place_rebuild_threshold") = 0.15,
         py::arg("cg_fallback_max_iterations") = 20000,
         py::arg("gradient_descent_learning_rate") = 0.01,
         py::arg_v("slip_mode", mehp::SlipMode::FIXED, "SlipMode.FIXED"))
    .def_readwrite("solver", &mehp::MEHPLinearSolverConfig::solver)
    .def_readwrite("simplification_mode",
                   &mehp::MEHPLinearSolverConfig::simplificationMode)
    .def_readwrite("inactive_removal_cutoff",
                   &mehp::MEHPLinearSolverConfig::inactiveRemovalCutoff)
    .def_readwrite("residual_reduction",
                   &mehp::MEHPLinearSolverConfig::residualReduction)
    .def_readwrite("max_iterations",
                   &mehp::MEHPLinearSolverConfig::maxIterations)
    .def_readwrite("spring_breaking_length",
                   &mehp::MEHPLinearSolverConfig::springBreakingLength)
    .def_readwrite("in_place_rebuild_threshold",
                   &mehp::MEHPLinearSolverConfig::inPlaceRebuildThreshold)
    .def_readwrite("cg_fallback_max_iterations",
                   &mehp::MEHPLinearSolverConfig::cgFallbackMaxIterations)
    .def_readwrite("gradient_descent_learning_rate",
                   &mehp::MEHPLinearSolverConfig::gradientDescentLearningRate)
    .def_readwrite("slip_mode", &mehp::MEHPLinearSolverConfig::slipMode);

  py::class_<mehp::MEHPSlipLinkSolverConfig, py::smart_holder>(
    m,
    "MEHPSlipLinkSolverConfig",
    "Configuration for alternating coordinate/partition slip-link relaxation.")
    .def(py::init([](const long int maxNrOfSteps,
                     const double xtol,
                     const double initialResidualToUse,
                     const mehp::StructureSimplificationMode simplificationMode,
                     const double inactiveRemovalCutoff,
                     const bool doInnerIterations,
                     const mehp::LinkSwappingMode linkSwappingMode,
                     const int swappingFrequency,
                     const double oneOverSpringPartitionUpperLimit,
                     const int nrOfCrosslinkSwapsAllowedPerSlipLink,
                     const mehp::SlipMode slipMode,
                     const double springBreakingLength) {
           return mehp::MEHPSlipLinkSolverConfig{
             maxNrOfSteps,
             xtol,
             initialResidualToUse,
             simplificationMode,
             inactiveRemovalCutoff,
             doInnerIterations,
             linkSwappingMode,
             swappingFrequency,
             oneOverSpringPartitionUpperLimit,
             nrOfCrosslinkSwapsAllowedPerSlipLink,
             slipMode,
             springBreakingLength
           };
         }),
         py::arg("max_nr_of_steps") = 50000,
         py::arg("xtol") = 1e-5,
         py::arg("initial_residual_to_use") = -1.0,
         py::arg_v("simplification_mode",
                   mehp::StructureSimplificationMode::NO_SIMPLIFICATION,
                   "StructureSimplificationMode.NO_SIMPLIFICATION"),
         py::arg("inactive_removal_cutoff") = 1e-3,
         py::arg("do_inner_iterations") = false,
         py::arg_v("link_swapping_mode",
                   mehp::LinkSwappingMode::NO_SWAPPING,
                   "LinkSwappingMode.NO_SWAPPING"),
         py::arg("swapping_frequency") = 10,
         py::arg("one_over_spring_partition_upper_limit") = 1.0,
         py::arg("nr_of_crosslink_swaps_allowed_per_sliplink") = -1,
         py::arg_v("slip_mode", mehp::SlipMode::SLIDING, "SlipMode.SLIDING"),
         py::arg("spring_breaking_length") = -1.0)
    .def_readwrite("max_nr_of_steps",
                   &mehp::MEHPSlipLinkSolverConfig::maxNrOfSteps)
    .def_readwrite("xtol", &mehp::MEHPSlipLinkSolverConfig::xtol)
    .def_readwrite("initial_residual_to_use",
                   &mehp::MEHPSlipLinkSolverConfig::initialResidualToUse)
    .def_readwrite("simplification_mode",
                   &mehp::MEHPSlipLinkSolverConfig::simplificationMode)
    .def_readwrite("inactive_removal_cutoff",
                   &mehp::MEHPSlipLinkSolverConfig::inactiveRemovalCutoff)
    .def_readwrite("do_inner_iterations",
                   &mehp::MEHPSlipLinkSolverConfig::doInnerIterations)
    .def_readwrite("link_swapping_mode",
                   &mehp::MEHPSlipLinkSolverConfig::linkSwappingMode)
    .def_readwrite("swapping_frequency",
                   &mehp::MEHPSlipLinkSolverConfig::swappingFrequency)
    .def_readwrite(
      "one_over_spring_partition_upper_limit",
      &mehp::MEHPSlipLinkSolverConfig::oneOverSpringPartitionUpperLimit)
    .def_readwrite(
      "nr_of_crosslink_swaps_allowed_per_sliplink",
      &mehp::MEHPSlipLinkSolverConfig::nrOfCrosslinkSwapsAllowedPerSlipLink)
    .def_readwrite("slip_mode", &mehp::MEHPSlipLinkSolverConfig::slipMode)
    .def_readwrite("spring_breaking_length",
                   &mehp::MEHPSlipLinkSolverConfig::springBreakingLength);

  py::class_<mehp::MEHPNonlinearSolverConfig, py::smart_holder>(
    m,
    "MEHPNonlinearSolverConfig",
    "Configuration for general nonlinear optimization MEHP equilibrium solving "
    "via NLopt.")
    .def(py::init([](const mehp::NLoptAlgorithm algorithm,
                     const long int maxNrOfSteps,
                     const double xtol,
                     const double ftol,
                     const mehp::StructureSimplificationMode simplificationMode,
                     const double inactiveRemovalCutoff,
                     const double springBreakingLength,
                     const mehp::SlipMode slipMode) {
           return mehp::MEHPNonlinearSolverConfig{ algorithm,
                                                   maxNrOfSteps,
                                                   xtol,
                                                   ftol,
                                                   simplificationMode,
                                                   inactiveRemovalCutoff,
                                                   springBreakingLength,
                                                   slipMode };
         }),
         py::arg_v(
           "algorithm", mehp::NLoptAlgorithm::LD_MMA, "NLoptAlgorithm.LD_MMA"),
         py::arg("max_nr_of_steps") = 250000,
         py::arg("xtol") = 1e-12,
         py::arg("ftol") = 1e-9,
         py::arg_v("simplification_mode",
                   mehp::StructureSimplificationMode::NO_SIMPLIFICATION,
                   "StructureSimplificationMode.NO_SIMPLIFICATION"),
         py::arg("inactive_removal_cutoff") = 1e-6,
         py::arg("spring_breaking_length") = -1.0,
         py::arg_v("slip_mode", mehp::SlipMode::FIXED, "SlipMode.FIXED"))
    .def_readwrite("algorithm", &mehp::MEHPNonlinearSolverConfig::algorithm)
    .def_readwrite("max_nr_of_steps",
                   &mehp::MEHPNonlinearSolverConfig::maxNrOfSteps)
    .def_readwrite("xtol", &mehp::MEHPNonlinearSolverConfig::xtol)
    .def_readwrite("ftol", &mehp::MEHPNonlinearSolverConfig::ftol)
    .def_readwrite("simplification_mode",
                   &mehp::MEHPNonlinearSolverConfig::simplificationMode)
    .def_readwrite("inactive_removal_cutoff",
                   &mehp::MEHPNonlinearSolverConfig::inactiveRemovalCutoff)
    .def_readwrite("spring_breaking_length",
                   &mehp::MEHPNonlinearSolverConfig::springBreakingLength)
    .def_readwrite("slip_mode", &mehp::MEHPNonlinearSolverConfig::slipMode);

  py::class_<mehp::MEHPSolveResult, py::smart_holder>(
    m,
    "MEHPSolveResult",
    "Diagnostic summary returned by an MEHPSolver execution.")
    .def(py::init<>())
    .def_readonly("exit_reason", &mehp::MEHPSolveResult::exitReason)
    .def_readonly("nr_of_iterations", &mehp::MEHPSolveResult::nrOfIterations)
    .def_readonly("residual", &mehp::MEHPSolveResult::residual)
    .def_readonly("energy", &mehp::MEHPSolveResult::energy)
    .def_readonly("cg_fallback_triggered",
                  &mehp::MEHPSolveResult::cgFallbackTriggered)
    .def_readonly("effective_sle_solver",
                  &mehp::MEHPSolveResult::effectiveSLESolver)
    .def_readonly("nr_of_relocations", &mehp::MEHPSolveResult::nrOfRelocations)
    .def_readonly("nr_of_shifts", &mehp::MEHPSolveResult::nrOfShifts)
    .def_readonly("one_over_spring_partition_upper_limit_used",
                  &mehp::MEHPSolveResult::oneOverSpringPartitionUpperLimitUsed);

  py::class_<mehp::MEHPSolver, py::smart_holder>(
    m,
    "MEHPSolver",
    "Abstract strategy interface for solving mechanical equilibrium on an "
    "MEHPNetwork.")
    .def("solve",
         &solveWithPythonDisplacements,
         "Solve for equilibrium link displacements on net starting from u.",
         py::arg("net"),
         py::arg("u"),
         py::arg("force_evaluator"),
         py::arg("is_2d") = false);

  py::class_<mehp::LinearSystemMEHPSolver, mehp::MEHPSolver, py::smart_holder>(
    m,
    "LinearSystemMEHPSolver",
    "Equilibrium solver for linear (Hookean) MEHP networks via sparse linear "
    "system assembly.")
    .def(py::init<const mehp::MEHPLinearSolverConfig&>(),
         py::arg_v("config",
                   mehp::MEHPLinearSolverConfig{},
                   "MEHPLinearSolverConfig()"))
    .def_property("config",
                  &mehp::LinearSystemMEHPSolver::getConfig,
                  &mehp::LinearSystemMEHPSolver::setConfig)
    .def("get_config", &mehp::LinearSystemMEHPSolver::getConfig)
    .def("set_config",
         &mehp::LinearSystemMEHPSolver::setConfig,
         py::arg("config"));

  py::
    class_<mehp::SlipLinkRelaxationSolver, mehp::MEHPSolver, py::smart_holder>(
      m,
      "SlipLinkRelaxationSolver",
      "Iterative alternating coordinate and contour-partition relaxation "
      "solver for slip-link networks.")
      .def(py::init<const mehp::MEHPSlipLinkSolverConfig&>(),
           py::arg_v("config",
                     mehp::MEHPSlipLinkSolverConfig{},
                     "MEHPSlipLinkSolverConfig()"))
      .def_property("config",
                    &mehp::SlipLinkRelaxationSolver::getConfig,
                    &mehp::SlipLinkRelaxationSolver::setConfig)
      .def("get_config", &mehp::SlipLinkRelaxationSolver::getConfig)
      .def("set_config",
           &mehp::SlipLinkRelaxationSolver::setConfig,
           py::arg("config"));

  py::class_<mehp::NonlinearOptimizationSolver,
             mehp::MEHPSolver,
             py::smart_holder>(m,
                               "NonlinearOptimizationSolver",
                               "General nonlinear energy minimization solver "
                               "for MEHP networks using NLopt.")
    .def(py::init<const mehp::MEHPNonlinearSolverConfig&>(),
         py::arg_v("config",
                   mehp::MEHPNonlinearSolverConfig{},
                   "MEHPNonlinearSolverConfig()"))
    .def_property("config",
                  &mehp::NonlinearOptimizationSolver::getConfig,
                  &mehp::NonlinearOptimizationSolver::setConfig)
    .def("get_config", &mehp::NonlinearOptimizationSolver::getConfig)
    .def("set_config",
         &mehp::NonlinearOptimizationSolver::setConfig,
         py::arg("config"));

  ////////////////////////////////////////////////////////////////
  // MARK: Unified MEHP Simulator
  py::class_<mehp::MEHPSimulator, py::smart_holder> mehpSimulatorClass(
    m,
    "MEHPSimulator",
    R"pbdoc(
      Unified Mechanical Equilibrium of Heterogeneous Polymers (MEHP) simulator
      coordinating a canonical MEHPNetwork, a stateless MEHPForceEvaluator,
      and an interchangeable MEHPSolver strategy.
    )pbdoc");

  mehpSimulatorClass
    .def(py::init<const mehp::MEHPNetwork&,
                  std::shared_ptr<mehp::MEHPForceEvaluator>,
                  bool>(),
         "Instantiate MEHPSimulator from an existing MEHPNetwork.",
         py::arg("network"),
         py::arg("force_evaluator") = nullptr,
         py::arg("is_2d") = false)
    .def(py::init<const pe::Universe&,
                  const pylimer_tools::topo::entanglement_detection::
                    AtomPairEntanglements&,
                  int,
                  mehp::EntanglementMode,
                  std::shared_ptr<mehp::MEHPForceEvaluator>,
                  bool,
                  bool,
                  bool,
                  bool>(),
         "Instantiate MEHPSimulator from an atomistic Universe.",
         py::arg("universe"),
         py::arg_v(
           "entanglements",
           pylimer_tools::topo::entanglement_detection::AtomPairEntanglements(),
           "AtomPairEntanglements()"),
         py::arg("crosslinker_type") = 2,
         py::arg_v("entanglement_mode",
                   mehp::EntanglementMode::MERGED_LINKS,
                   "EntanglementMode.MERGED_LINKS"),
         py::arg("force_evaluator") = nullptr,
         py::arg("is_2d") = false,
         py::arg("keep_2functional_crosslinkers") = false,
         py::arg("remove_primary_loops") = false,
         py::arg("assume_box_large_enough") = false)
    .def(
      py::init<const pe::Box&,
               const mehp::CrosslinkerUniverse&,
               const pylimer_tools::topo::entanglement_detection::
                 StrandPairEntanglements&,
               int,
               mehp::EntanglementMode,
               std::shared_ptr<mehp::MEHPForceEvaluator>,
               bool,
               bool,
               bool,
               bool,
               bool>(),
      "Instantiate MEHPSimulator from a Box and a CrosslinkerUniverse.",
      py::arg("box"),
      py::arg("crosslinker_universe"),
      py::arg_v(
        "entanglements",
        pylimer_tools::topo::entanglement_detection::StrandPairEntanglements(),
        "StrandPairEntanglements()"),
      py::arg("crosslinker_type") = 2,
      py::arg_v("entanglement_mode",
                mehp::EntanglementMode::MERGED_LINKS,
                "EntanglementMode.MERGED_LINKS"),
      py::arg("force_evaluator") = nullptr,
      py::arg("is_2d") = false,
      py::arg("remove_dangling_chains") = false,
      py::arg("include_free_chains") = true,
      py::arg("remove_primary_loops") = false,
      py::arg("assume_box_large_enough") = false)
    .def_static("construct_with_random_entanglements",
                &mehp::MEHPSimulator::constructWithRandomEntanglements,
                "Construct MEHPSimulator from a Universe with randomly sampled "
                "entanglements.",
                py::arg("universe"),
                py::arg("nr_of_entanglements_to_sample"),
                py::arg("upper_sampling_cutoff"),
                py::arg("lower_sampling_cutoff") = 0.0,
                py::arg("min_nr_of_entanglements") = 0,
                py::arg("same_strand_cutoff") = 0.0,
                py::arg("seed") = "",
                py::arg("crosslinker_type") = 2,
                py::arg_v("entanglement_mode",
                          mehp::EntanglementMode::MERGED_LINKS,
                          "EntanglementMode.MERGED_LINKS"),
                py::arg("with_crosslinkers") = true,
                py::arg("filter_dangling_and_soluble") = false,
                py::arg("force_evaluator") = nullptr,
                py::arg("is_2d") = false,
                py::arg("keep_2functional_crosslinkers") = false,
                py::arg("remove_primary_loops") = false,
                py::arg("assume_box_large_enough") = false)
    .def(
      "__copy__",
      [](const mehp::MEHPSimulator& self) { return mehp::MEHPSimulator(self); })
    // Configuration & Mutators
    .def_property("force_evaluator",
                  &mehp::MEHPSimulator::getForceEvaluator,
                  &mehp::MEHPSimulator::setForceEvaluator)
    .def("set_force_evaluator",
         &mehp::MEHPSimulator::setForceEvaluator,
         py::arg("force_evaluator"))
    .def("get_force_evaluator", &mehp::MEHPSimulator::getForceEvaluator)
    .def("config_mean_bond_length",
         &mehp::MEHPSimulator::configMeanBondLength,
         py::arg("mean_bond_length"))
    .def("config_mean_squared_bond_length",
         &mehp::MEHPSimulator::configMeanSquaredBondLength,
         py::arg("r02"))
    .def("config_assume_box_large_enough",
         &mehp::MEHPSimulator::configAssumeBoxLargeEnough,
         py::arg("assume") = true)
    .def("config_assume_network_is_complete",
         &mehp::MEHPSimulator::configAssumeNetworkIsComplete,
         py::arg("assume") = false)
    .def("config_remove_primary_loops",
         &mehp::MEHPSimulator::configRemovePrimaryLoops,
         py::arg("remove") = false)
    .def_property(
      "is_2d", &mehp::MEHPSimulator::getIs2D, &mehp::MEHPSimulator::setIs2D)
    .def("set_is_2d", &mehp::MEHPSimulator::setIs2D, py::arg("is_2d"))
    .def("get_is_2d", &mehp::MEHPSimulator::getIs2D)
    .def("deform_to", &mehp::MEHPSimulator::deformTo, py::arg("target_box"))
    // Simplification
    .def("break_too_long_strands",
         &mehp::MEHPSimulator::breakTooLongStrands,
         py::arg("breaking_length"))
    .def("remove_inactive_links",
         &mehp::MEHPSimulator::removeInactiveLinks,
         py::arg("tolerance") = 1e-6)
    .def("unlink_bifunctional_links",
         &mehp::MEHPSimulator::unlinkBifunctionalLinks,
         py::arg("entanglements_only") = false)
    .def("simplify_network",
         &mehp::MEHPSimulator::simplifyNetwork,
         py::arg("mode"),
         py::arg("tolerance") = 1e-6,
         py::arg("breaking_length") = -1.0)
    // Unified Solver Execution
    .def("run",
         [](mehp::MEHPSimulator& self) {
           const std::shared_ptr<mehp::MEHPForceEvaluator> eval =
             self.getForceEvaluator();
           if (eval == nullptr) {
             return self.run();
           }
           const mehp::MEHPSolveCallbacks cb = makePythonInterruptCallbacks();
           if (eval->isLinear()) {
             return self.run(mehp::MEHPLinearSolverConfig{}, cb);
           }
           return self.run(mehp::MEHPNonlinearSolverConfig{}, cb);
         })
    .def(
      "run",
      [](mehp::MEHPSimulator& self, const mehp::MEHPSolver& solver) {
        return self.run(solver, makePythonInterruptCallbacks());
      },
      py::arg("solver"))
    .def(
      "run",
      [](mehp::MEHPSimulator& self,
         const mehp::MEHPLinearSolverConfig& config) {
        return self.run(config, makePythonInterruptCallbacks());
      },
      py::arg("config"))
    .def(
      "run",
      [](mehp::MEHPSimulator& self,
         const mehp::MEHPSlipLinkSolverConfig& config) {
        return self.run(config, makePythonInterruptCallbacks());
      },
      py::arg("config"))
    .def(
      "run",
      [](mehp::MEHPSimulator& self,
         const mehp::MEHPNonlinearSolverConfig& config) {
        return self.run(config, makePythonInterruptCallbacks());
      },
      py::arg("config"))
    // Network & Diagnostic Getters
    .def_property_readonly("network", &mehp::MEHPSimulator::getNetwork)
    .def("get_network",
         &mehp::MEHPSimulator::getNetwork,
         py::return_value_policy::reference_internal)
    .def("get_current_displacements",
         &mehp::MEHPSimulator::getCurrentDisplacements)
    .def("set_current_displacements",
         &mehp::MEHPSimulator::setCurrentDisplacements,
         py::arg("new_displacements"))
    .def("get_initial_coordinates", &mehp::MEHPSimulator::getInitialCoordinates)
    .def("get_spring_vectors", &mehp::MEHPSimulator::getSpringVectors)
    .def("get_spring_lengths", &mehp::MEHPSimulator::getSpringLengths)
    .def("get_average_spring_length",
         &mehp::MEHPSimulator::getAverageSpringLength)
    .def("get_spring_contour_length",
         &mehp::MEHPSimulator::getSpringContourLength)
    .def("get_spring_partitions", &mehp::MEHPSimulator::getSpringPartitions)
    .def("get_strand_contour_length",
         &mehp::MEHPSimulator::getStrandContourLength)
    .def("get_average_contour_length",
         &mehp::MEHPSimulator::getAverageContourLength)
    .def("get_nr_of_links", &mehp::MEHPSimulator::getNrOfLinks)
    .def("get_nr_of_nodes", &mehp::MEHPSimulator::getNrOfNodes)
    .def("get_nr_of_strands", &mehp::MEHPSimulator::getNrOfStrands)
    .def("get_nr_of_springs", &mehp::MEHPSimulator::getNrOfSprings)
    .def("get_nr_of_active_nodes",
         &mehp::MEHPSimulator::getNrOfActiveNodes,
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_springs",
         &mehp::MEHPSimulator::getNrOfActiveSprings,
         py::arg("tolerance") = 1e-3)
    .def("get_nr_of_active_strands",
         &mehp::MEHPSimulator::getNrOfActiveStrands,
         py::arg("tolerance") = 1e-3)
    .def("get_crosslinker_universe", &mehp::MEHPSimulator::getCrosslinkerVerse)
    .def("get_crosslinker_verse", &mehp::MEHPSimulator::getCrosslinkerVerse)
    .def("get_default_r0_square", &mehp::MEHPSimulator::getDefaultR0Square)
    .def("get_default_bond_length", &mehp::MEHPSimulator::getDefaultBondLength)
    .def("validate_network", &mehp::MEHPSimulator::validateNetwork)
    .def_property_readonly("last_solve_result",
                           &mehp::MEHPSimulator::getLastSolveResult)
    .def("get_last_solve_result", &mehp::MEHPSimulator::getLastSolveResult)
    .def("get_exit_reason", &mehp::MEHPSimulator::getExitReason)
    .def("get_nr_of_iterations", &mehp::MEHPSimulator::getNrOfIterations)
    .def("get_simulation_has_run", &mehp::MEHPSimulator::getSimulationHasRun)
    // Physical & Topological Observables
    .def("get_energy", &mehp::MEHPSimulator::getEnergy)
    .def("get_residuals", &mehp::MEHPSimulator::getResiduals)
    .def("get_residual_norm", &mehp::MEHPSimulator::getResidualNorm)
    .def("get_pressure", &mehp::MEHPSimulator::getPressure)
    .def("get_gamma_factors",
         &mehp::MEHPSimulator::getGammaFactors,
         R"pbdoc(
Compute the per-spring topological Gamma factors :math:`\Gamma_i = \mathbf{r}_i^2 / (N_i b_0^2)`.

Summing this vector (``np.sum(sim.get_gamma_factors(b02))``) gives the unnormalized total
needed to compute the shear modulus :math:`G = (k_B T / V) \sum_i \Gamma_i`.

:param b02: Mean-squared bond length :math:`b_0^2` (defaults to network bond length squared).
)pbdoc",
         py::arg("b02") = -1.0)
    .def("get_gamma_factors_in_dir",
         &mehp::MEHPSimulator::getGammaFactorsInDir,
         R"pbdoc(
Compute the per-spring directional Gamma factors :math:`\Gamma_{i,\alpha} = r_{i,\alpha}^2 / (N_i b_0^2)`.

:param dir: Cartesian direction index (``0``, ``1``, or ``2``).
:param b02: Mean-squared bond length :math:`b_0^2` (defaults to network bond length squared).
)pbdoc",
         py::arg("dir"),
         py::arg("b02") = -1.0)
    .def("get_gamma_factor",
         &mehp::MEHPSimulator::getGammaFactor,
         R"pbdoc(
Compute the average topological Gamma factor ``sum(gamma_factors) / denom``.

:param b02: Mean-squared bond length :math:`b_0^2` (defaults to network bond length squared).
:param nr_of_chains: Denominator for averaging. When ``< 1`` (default ``-1``), divides by
                     ``network.nr_of_strands`` — note that in ``EntanglementMode.SPRING``,
                     ``nr_of_strands`` includes entanglement coupling springs, so comparing
                     ``MERGED_LINKS`` vs. ``SPRING`` with the default ``-1`` is not apples-to-apples.
                     Pass ``1`` or use ``np.sum(get_gamma_factors(b02))`` when computing shear modulus.
)pbdoc",
         py::arg("b02") = -1.0,
         py::arg("nr_of_chains") = -1)
    .def("get_soluble_weight_fraction",
         &mehp::MEHPSimulator::getSolubleWeightFraction,
         py::arg("tolerance") = 1e-3)
    .def("get_dangling_weight_fraction",
         &mehp::MEHPSimulator::getDanglingWeightFraction,
         py::arg("tolerance") = 1e-3)
    .def("get_active_weight_fraction",
         &mehp::MEHPSimulator::getActiveWeightFraction,
         py::arg("tolerance") = 1e-3)
    .def("get_effective_functionality_of_atoms",
         &mehp::MEHPSimulator::getEffectiveFunctionalityOfAtoms,
         py::arg("tolerance") = 1e-3)
    .def("count_active_clustered_atoms",
         &mehp::MEHPSimulator::countActiveClusteredAtoms,
         py::arg("tolerance") = 1e-3)
    .def("get_indices_of_active_nodes",
         &mehp::MEHPSimulator::getIndicesOfActiveNodes,
         py::arg("tolerance") = 1e-3)
    .def("get_num_intra_chain_slip_links",
         &mehp::MEHPSimulator::getNumIntraChainSlipLinks)
    .def("get_num_intra_chain_sliplinks",
         &mehp::MEHPSimulator::getNumIntraChainSlipLinks)
    .def("get_nr_of_intra_chain_sliplinks",
         &mehp::MEHPSimulator::getNumIntraChainSlipLinks)
    // OutputSupportingSimulation methods & overrides
    .def("config_step_output",
         &mehp::MEHPSimulator::configStepOutput,
         py::arg("output_configuration"))
    .def("config_average_output",
         &mehp::MEHPSimulator::configAverageOutput,
         py::arg("values"))
    .def("config_auto_correlator_output",
         &mehp::MEHPSimulator::configAutoCorrelatorOutput,
         py::arg("values"),
         py::arg("num_corr_in") = 32,
         py::arg("p") = 16,
         py::arg("m") = 2)
    .def("config_restart_output",
         &mehp::MEHPSimulator::configRestartOutput,
         py::arg("file"),
         py::arg("output_every") = 50000)
    .def("has_outputs", &mehp::MEHPSimulator::hasOutputs)
    .def("get_current_time",
         &mehp::MEHPSimulator::getCurrentTime,
         py::arg("current_step") = 0.0)
    .def("get_temperature", &mehp::MEHPSimulator::getTemperature)
    .def("get_gamma", &mehp::MEHPSimulator::getGamma)
    .def("get_timestep", &mehp::MEHPSimulator::getTimestep)
    .def("get_volume", &mehp::MEHPSimulator::getVolume)
    .def("get_residual", &mehp::MEHPSimulator::getResidual)
    .def("get_stress_tensor", &mehp::MEHPSimulator::getStressTensor)
    .def("get_bond_lengths", &mehp::MEHPSimulator::getBondLengths)
    .def("get_coordinates", &mehp::MEHPSimulator::getCoordinates)
    .def("get_num_relocations", &mehp::MEHPSimulator::getNumRelocations)
    .def("get_num_shifts", &mehp::MEHPSimulator::getNumShifts)
    .def("get_num_bonds_to_form", &mehp::MEHPSimulator::getNumBondsToForm)
    .def("get_nr_of_atoms", &mehp::MEHPSimulator::getNumAtoms)
    .def("get_nr_of_extra_atoms", &mehp::MEHPSimulator::getNumExtraAtoms)
    .def("get_nr_of_bonds", &mehp::MEHPSimulator::getNumBonds)
    .def("get_nr_of_extra_bonds", &mehp::MEHPSimulator::getNumExtraBonds)
    .def("get_nr_of_particles", &mehp::MEHPSimulator::getNumParticles)
#ifdef CEREALIZABLE
    .def(
      "write_restart_file",
      [](mehp::MEHPSimulator& self, const std::string& filename) {
        std::string filenameCopy = filename;
        self.writeRestartFile(filenameCopy);
      },
      py::arg("file"))
    .def_static(
      "read_restart_file",
      [](const std::string& filename) {
        std::string filenameCopy = filename;
        return mehp::MEHPSimulator::readRestartFile(filenameCopy);
      },
      py::arg("file"))
    .def("serialize_to_string",
         [](const mehp::MEHPSimulator& sim) {
           return py::bytes(sim.serializeToString());
         })
    .def_static("construct_from_string",
                &mehp::MEHPSimulator::constructFromString,
                py::arg("s"))
    .def(py::pickle(
      [](const mehp::MEHPSimulator& sim) {
        return py::make_tuple(py::bytes(sim.serializeToString()));
      },
      [](const py::tuple& t) {
        if (t.size() != 1) {
          throw std::runtime_error("Invalid pickle state for MEHPSimulator.");
        }
        const std::string in = t[0].cast<std::string>();
        return mehp::MEHPSimulator::constructFromString(in);
      }))
#endif
    ;

  // Bind slip-link insertion with both underscore and legacy naming
  bindSlipLinkInsertionMethods(
    mehpSimulatorClass, "add_slip_links", "randomly_add_slip_links");
  bindSlipLinkInsertionMethods(
    mehpSimulatorClass, "add_sliplinks", "randomly_add_sliplinks");

  /**
   * DPD Simulations
   */

  ////////////////////////////////////////////////////////////////
  // MARK: DPD Simulator
  py::class_<dpd::DPDSimulator, py::smart_holder>(m,
                                                  "DPDSimulator",
                                                  R"pbdoc(
     A quick-and-dirty implementation of the dissipative particle dynamics (DPD) simulation
     with slip-springs as presented by :cite:t:`langeloth_recovering_2013` 
     and :cite:t:`schneider_simulation_2021`.
     )pbdoc")
    .def(py::init<const pe::Universe&,
                  const int,
                  const int,
                  const bool,
                  const std::string&,
                  const bool>(),
         "Get an instance of this class",
         py::arg("universe"),
         py::arg("crosslinker_type") = 2,
         py::arg("slipspring_bond_type") = 9,
         py::arg("is_2d") = false,
         py::arg("seed") = "",
         py::arg("assume_box_large_enough") = false)
    .def(
      "run_simulation",
      [](dpd::DPDSimulator& sim,
         const int nSteps,
         const double dt,
         const bool withMC) {
        sim.configTimeStep(dt);
        return sim.runSimulation(
          nSteps,
          withMC,
          []() { return PyErr_CheckSignals() != 0; },
          []() { throw py::error_already_set(); });
      },
      py::arg("n_steps"),
      py::arg("dt") = 0.06,
      py::arg("with_MC") = false)
    .def("assume_box_large_enough",
         &dpd::DPDSimulator::configAssumeBoxLargeEnough,
         R"pbdoc(
          Configure whether to assume minimum-image PBC is sufficient for all bonds.

          If your bonds could exceed half the box length, keep this ``False`` (default) so per-bond periodic box offsets are tracked.
          Set to ``True`` to apply minimum-image PBC on every bond evaluation.
         )pbdoc",
         py::arg("assume") = true)
    .def("config_assume_box_large_enough",
         &dpd::DPDSimulator::configAssumeBoxLargeEnough,
         R"pbdoc(
          Configure whether to assume minimum-image PBC is sufficient for all bonds.

          If your bonds could exceed half the box length, keep this ``False`` (default) so per-bond periodic box offsets are tracked.
          Set to ``True`` to apply minimum-image PBC on every bond evaluation.
         )pbdoc",
         py::arg("assume") = true)
    .def("assumes_box_large_enough",
         &dpd::DPDSimulator::assumesBoxLargeEnough,
         R"pbdoc(
          Return ``True`` if minimum-image PBC is assumed sufficient for bonds, ``False`` if per-bond periodic box offsets are tracked.
         )pbdoc")
    .def("create_slip_springs",
         &dpd::DPDSimulator::createSlipSprings,
         R"pbdoc(
          Randomly add up to the specified number of slip-springs to eligible neighbours within the configured cut-offs.

          :param num: Maximum number of slip-springs to create
          :param bond_type: Bond type for created slip-springs (``0`` inherits ``slipspring_bond_type``)
          :return: Number of slip-springs actually created
     )pbdoc",
         py::arg("num"),
         py::arg("bond_type") = 0)
    .def("shift_slip_springs",
         &dpd::DPDSimulator::shiftSlipSprings,
         R"pbdoc(
          Perform one sweep of Monte Carlo slip-spring shift moves along polymer chains.

          :param kbt: Thermal energy scale :math:`k_B T` for Metropolis acceptance (must be ``>= 0``)
          :return: Number of accepted shift moves
     )pbdoc",
         py::arg("kbt") = 1.0)
    .def("relocate_slip_springs",
         &dpd::DPDSimulator::relocateSlipSprings,
         R"pbdoc(
          Perform one sweep of Monte Carlo slip-spring chain-end relocation moves.

          :param kbt: Thermal energy scale :math:`k_B T` for Metropolis acceptance (must be ``>= 0``)
          :return: Number of accepted relocation moves
     )pbdoc",
         py::arg("kbt") = 1.0)
    .def("config_a",
         &dpd::DPDSimulator::configA,
         R"pbdoc(
          Configure the force-field (pair-style) parameter `A`.
     )pbdoc",
         py::arg("A") = 25.)
    .def("config_sigma",
         &dpd::DPDSimulator::configSigma,
         R"pbdoc(
          Configure the force-field (pair-style) parameter `\sigma` (and set `\gamma = 0.5 \sigma^2`).
     )pbdoc",
         py::arg("sigma") = 3.)
    .def("config_gamma",
         &dpd::DPDSimulator::configGamma,
         R"pbdoc(
          Configure the DPD dissipative friction parameter `\gamma` independently of `\sigma`.
     )pbdoc",
         py::arg("gamma") = 4.5)
    .def("config_spring_constant",
         &dpd::DPDSimulator::configSpringConstant,
         R"pbdoc(
          Configure the force-field (bond-style) parameter `k`, the spring constant.
     )pbdoc",
         py::arg("k") = 2.)
    .def("config_lambda",
         &dpd::DPDSimulator::configLambda,
         R"pbdoc(
          Configure the modified velocity verlet integration parameter `\lambda`.
     )pbdoc",
         py::arg("l") = 0.65)
    .def("config_slipspring_high_cutoff",
         &dpd::DPDSimulator::configSlipspringHighCutoff,
         R"pbdoc(
          Configure the higher cut-off of how far a pair may be distanced for a slip-spring to be created.
     )pbdoc",
         py::arg("cutoff") = 2.)
    .def("config_slipspring_low_cutoff",
         &dpd::DPDSimulator::configSlipspringLowCutoff,
         R"pbdoc(
          Configure the lower cut-off of how far a pair may be distanced for a slip-spring to be created.
     )pbdoc",
         py::arg("cutoff") = 0.5)
    .def("config_box_deformation",
         &dpd::DPDSimulator::configBoxDeformation,
         R"pbdoc(
          Configure where to (incrementally) deform the box to during the next simulation run.
     )pbdoc",
         py::arg("target_box"))
#ifdef CEREALIZABLE
    .def_static(
      "read_restart_file",
      [](const std::string& filename) {
        return dpd::DPDSimulator::readRestartFile(filename);
      },
      R"pbdoc(
          Read a restart file in order to continue a simulation.

          :param file: The file path to the restart file to read
     )pbdoc",
      py::arg("file"))
    .def("config_restart_output",
         &dpd::DPDSimulator::configRestartOutput,
         R"pbdoc(
          Set when to output a restart file.

          Note:
               The filename determines the type of serialization:
               .json, .xml are supported; other file endings will lead to binary serialization (fastest!).

          Caution:
               This method may not be backwards- nor forward-compatible.
               Use the same version of pylimer-tools if you want to be sure that things work.

          :param file: The file path to the restart file to write
          :param output_every: How often to write the restart file (default: 50000)
     )pbdoc",
         py::arg("file"),
         py::arg("output_every") = 50000)
    .def(
      "write_restart_file",
      [](dpd::DPDSimulator& self, const std::string& filename) {
        std::string filenameCopy = filename;
        self.writeRestartFile(filenameCopy);
      },
      R"pbdoc(
          Explicitly force the writing of a restart file, now!

          :param file: The file path and name of the restart file to be written.
                       Can end in .xml, .json or anything else (-> binary)
     )pbdoc",
      py::arg("file"))
#endif
    .def("config_average_output",
         &dpd::DPDSimulator::configAverageOutput,
         R"pbdoc(
          Set which values to compute averages for.

          :param values: A list of OutputConfiguration structs specifying what to average
     )pbdoc",
         py::arg("values"))
    .def("config_auto_correlator_output",
         &dpd::DPDSimulator::configAutoCorrelatorOutput,
         R"pbdoc(
          Set which values to compute multiple-tau autocorrelation for.
          If you use this, you should cite :cite:t:`ramirez_efficient_2010`.

          :param values: a list of OutputConfiguration structs
          :param num_corr_in: Number of correlations in
          :param p: Parameter p for the autocorrelator
          :param m: Parameter m for the autocorrelator
     )pbdoc",
         py::arg("values"),
         py::arg("num_corr_in") = 32,
         py::arg("p") = 16,
         py::arg("m") = 2)
    .def("config_step_output",
         &dpd::DPDSimulator::configStepOutput,
         R"pbdoc(
          Set which values to log.

          :param values: a list of OutputConfiguration structs
     )pbdoc",
         py::arg("values"))
    .def("config_shift_possibility_empty",
         &dpd::DPDSimulator::configShiftPossibilityEmpty,
         R"pbdoc(
          Configure the possibility of shifting to empty positions.

          This setting affects Monte Carlo moves in the simulation.
         )pbdoc",
         py::arg("shift_possibility_empty") = true)
    .def("config_shift_one_at_a_time",
         &dpd::DPDSimulator::configShiftOneAtATime,
         R"pbdoc(
          Configure whether to shift atoms one at a time.

          This setting affects Monte Carlo move behavior in the simulation.
         )pbdoc",
         py::arg("shift_one_at_a_time") = false)
    .def("config_num_steps_mc",
         &dpd::DPDSimulator::configNumStepsMC,
         R"pbdoc(
          Configure the number of steps to do in one MC sequence.

          :param num_steps: Number of Monte Carlo steps per sequence
     )pbdoc",
         py::arg("num_steps") = 500)
    .def("config_num_steps_dpd",
         &dpd::DPDSimulator::configNumStepsDPD,
         R"pbdoc(
          Configure the number of steps to do in one DPD sequence.

          :param num_steps: Number of DPD steps per sequence
     )pbdoc",
         py::arg("num_steps") = 500)
    .def("config_bond_formation",
         &dpd::DPDSimulator::configBondFormation,
         R"pbdoc(
          Configure how to do bond formation during the run.

          :param num_bonds_to_form (int): The nr of bonds to form in total. Use 0 to stop bond formation.
          :param num_bonds_per_atom_type (dict): The nr of bonds each atom type may have at most (e.g., 2 for strand atoms, 4 for a tertiary crosslinkers)
          :param bond_formation_dist (float): The maximum distance allowed to form bonds
          :param attempt_bond_formation_every (int): attempt to form bonds every this many steps during the simulation run
          :param atom_type_form_from (int): The atom type to start forming bonds from.
          :param atom_type_form_to (int): The atom type to start forming bonds to.
         )pbdoc",
         py::arg("num_bonds_to_form"),
         py::arg("max_bonds_per_atom_type"),
         py::arg("bond_formation_dist") = 1.0,
         py::arg("attempt_bond_formation_every") = 50,
         py::arg("atom_type_form_from") = 2,
         py::arg("atom_type_form_to") = 1)
    .def("get_nr_of_bonds_to_form",
         &dpd::DPDSimulator::getNrOfBondsToForm,
         R"pbdoc(
          Get the number of bonds that are configured to have to be formed.
     )pbdoc")
    .def("config_allow_relocation_in_network",
         &dpd::DPDSimulator::configAllowRelocationInNetwork,
         R"pbdoc(
          Configure whether a relocation step may happen when a slip-spring has ended at a crosslink.

          Side-effect: if true, the relocations may also happen *to* a slip-spring next to a crosslink.

          :param allow_relocation_in_network (bool): Whether to allow relocation in the network or not.
         )pbdoc",
         py::arg("allow_relocation_in_network") = false)
    .def("start_measuring_msd_for_atoms",
         &dpd::DPDSimulator::startMeasuringMSDForAtoms,
         R"pbdoc(
          Set a new origin for measuing the mean square displacement for a specified set of atoms
         )pbdoc",
         py::arg("atom_ids"))
    .def("get_universe",
         &dpd::DPDSimulator::getUniverse,
         R"pbdoc(
     Get a universe instance from the current coordinates (and connectivity).

     :param with_slipsprings: Whether to include slip-springs in the returned universe (default: True)
    )pbdoc",
         py::arg("with_slipsprings") = true)
    .def("refresh_current_state",
         &dpd::DPDSimulator::refreshCurrentState,
         R"pbdoc(
          After re-configuring the force-field parameters,
          this method should be called to update the current stress tensor etc.
     )pbdoc")
    .def("get_timestep",
         &dpd::DPDSimulator::getTimestep,
         R"pbdoc(
          Get the timestep used in the simulation.

          :return: The simulation timestep value
         )pbdoc")
    .def("get_current_timestep",
         &dpd::DPDSimulator::getCurrentTimestep,
         R"pbdoc(
          Get the current timestep number.

          :return: The current timestep index
         )pbdoc")
    .def("get_temperature",
         &dpd::DPDSimulator::getTemperature,
         R"pbdoc(
          Get the current system temperature.

          :return: The current temperature of the system
         )pbdoc")
    .def("get_bond_lengths",
         &dpd::DPDSimulator::getBondLengths,
         R"pbdoc(
          Get the lengths of all bonds in the system.

          :return: Vector containing the length of each bond
         )pbdoc")
    .def("get_coordinates",
         &dpd::DPDSimulator::getCoordinates,
         R"pbdoc(
          Get the current particle coordinates.

          :return: Vector of particle coordinates (x1,y1,z1,x2,y2,z2,...)
         )pbdoc")
    .def("get_spring_constant",
         &dpd::DPDSimulator::getSpringConstant,
         R"pbdoc(
          Get the current spring constant value.

          :return: The current spring constant for bond interactions
         )pbdoc")
    .def("get_shift_one_at_a_time",
         &dpd::DPDSimulator::getShiftOneAtATime,
         R"pbdoc(
          Get whether slip-springs are shifted one at a time.

          :return: True if shifting one at a time, False otherwise
         )pbdoc")
    .def("get_nr_of_slip_springs",
         &dpd::DPDSimulator::getNumSlipSprings,
         R"pbdoc(
          Get the current number of slip-springs in the system.

          :return: Number of slip-springs
         )pbdoc")
    .def("get_nr_of_atoms",
         &dpd::DPDSimulator::getNumAtoms,
         R"pbdoc(
          Get the total number of atoms in the system.

          :return: Total number of atoms
         )pbdoc")
    .def("get_nr_of_extra_atoms",
         &dpd::DPDSimulator::getNumExtraAtoms,
         R"pbdoc(
          Get the number of extra atoms (always 0 for DPD simulations).

          :return: Number of extra atoms
         )pbdoc")
    .def("get_nr_of_bonds",
         &dpd::DPDSimulator::getNumBonds,
         R"pbdoc(
          Get the number of regular bonds (excluding slip-springs).

          :return: Number of regular bonds
         )pbdoc")
    .def("get_nr_of_extra_bonds",
         &dpd::DPDSimulator::getNumExtraBonds,
         R"pbdoc(
          Get the number of extra bonds (slip-springs).

          :return: Number of slip-springs
         )pbdoc")
    .def("get_stress_tensor",
         &dpd::DPDSimulator::getStressTensor,
         R"pbdoc(
          Get the current stress tensor.

          :return: 3x3 stress tensor matrix
         )pbdoc")
    .def("get_nr_of_steps_dpd",
         &dpd::DPDSimulator::getNumStepsDPD,
         R"pbdoc(
          Get the configured number of DPD steps per sequence.

          :return: Number of DPD steps
         )pbdoc")
    .def("get_nr_of_steps_mc",
         &dpd::DPDSimulator::getNumStepsMC,
         R"pbdoc(
          Get the configured number of Monte Carlo steps per sequence.

          :return: Number of MC steps
         )pbdoc")
    .def("get_volume",
         &dpd::DPDSimulator::getVolume,
         R"pbdoc(
          Get the current system volume.

          :return: Current simulation box volume
         )pbdoc")
    .def("get_slip_spring_bond_type",
         &dpd::DPDSimulator::getSlipSpringBondType,
         R"pbdoc(
          Get the bond type identifier used for slip-springs.

          :return: Bond type for slip-springs
         )pbdoc")
    .def("get_velocities",
         &dpd::DPDSimulator::getVelocities,
         R"pbdoc(
          Get the current particle velocities.

          :return: Vector of particle velocities (vx1,vy1,vz1,vx2,vy2,vz2,...)
         )pbdoc")
    .def("get_target_kbt",
         &dpd::DPDSimulator::getTargetKbT,
         R"pbdoc(
          Get the target thermal energy :math:`k_B T = \sigma^2 / (2 \gamma)` configured on the DPD thermostat.

          :return: Target thermal energy :math:`k_B T`
         )pbdoc")
    .def("get_sigma",
         &dpd::DPDSimulator::getSigma,
         R"pbdoc(
          Get the configured DPD random force noise amplitude `\sigma`.
         )pbdoc")
    .def("get_dpd_gamma",
         &dpd::DPDSimulator::getDpdGamma,
         R"pbdoc(
          Get the configured DPD dissipative friction coefficient `\gamma`.
         )pbdoc")
    .def("get_num_shifts",
         &dpd::DPDSimulator::getNumShifts,
         R"pbdoc(
          Get the number of accepted slip-spring shift moves in the most recent MC sequence.
         )pbdoc")
    .def("get_num_relocations",
         &dpd::DPDSimulator::getNumRelocations,
         R"pbdoc(
          Get the number of accepted slip-spring relocation moves in the most recent MC sequence.
         )pbdoc")
    .def("get_shift_possibility_empty",
         &dpd::DPDSimulator::getShiftPossibilityEmpty,
         R"pbdoc(
          Get whether shifting to empty positions is allowed.

          :return: True if shifting to empty positions is allowed
         )pbdoc")
    .def("validate_neighbour_list",
         &dpd::DPDSimulator::validateNeighbourlist,
         R"pbdoc(
          Validate the neighbor list consistency for debugging purposes.

          :param cutoff: Cutoff distance for validation
         )pbdoc",
         py::arg("cutoff"))
    .def("validate_state",
         &dpd::DPDSimulator::validateState,
         R"pbdoc(
          Validate the current simulation state for debugging purposes.
          
          Checks internal data structure consistency and throws exceptions if issues are found.
         )pbdoc");

  py::native_enum<RouseSolverType>(
    m,
    "RouseSolverType",
    "enum.IntEnum",
    "Linear solver algorithm for implicit Rouse/Brownian time-stepping.")
    .value("DIRECT_LDLT",
           RouseSolverType::DIRECT_LDLT,
           "Direct sparse SimplicialLDLT solver")
    .value("DIRECT_LLT",
           RouseSolverType::DIRECT_LLT,
           "Direct sparse SimplicialLLT solver")
    .value("CONJUGATE_GRADIENT",
           RouseSolverType::CONJUGATE_GRADIENT,
           "Iterative Conjugate Gradient solver")
    .finalize();

  py::class_<RouseSimulationResult, py::smart_holder>(m,
                                                      "RouseSimulationResult",
                                                      R"pbdoc(
      Output container for Rouse dynamics simulation runs.
    )pbdoc")
    .def_readonly(
      "times", &RouseSimulationResult::times, "Sampled simulation times.")
    .def_readonly("msd",
                  &RouseSimulationResult::msd,
                  "Bead mean-squared displacement at sampled times.")
    .def_readonly("msd_center_of_mass",
                  &RouseSimulationResult::msdCenterOfMass,
                  "Center-of-mass mean-squared displacement at sampled times.")
    .def_readonly(
      "end_to_end_correlation",
      &RouseSimulationResult::endToEndCorrelation,
      "Normalized end-to-end vector autocorrelation at sampled times.")
    .def_readonly(
      "trajectory",
      &RouseSimulationResult::trajectory,
      "List of (N, 3) position matrices if save_trajectory was enabled.");

  py::class_<RouseBrownianSimulator, py::smart_holder>(m,
                                                       "RouseBrownianSimulator",
                                                       R"pbdoc(
      Implicit time-stepping Rouse and Brownian dynamics simulator for polymer chains
      and general spring networks.
    )pbdoc")
    .def(py::init<const std::vector<size_t>&,
                  const std::vector<size_t>&,
                  double,
                  double,
                  double,
                  double,
                  size_t,
                  const std::vector<double>&>(),
         R"pbdoc(
           Construct from edge lists (spring_from, spring_to), friction gamma,
           spring constant k, temperature T, time step dt, optional num_nodes,
           and optional edge weights.

           :param spring_from: Source node indices of springs
           :param spring_to: Target node indices of springs
           :param gamma: Bead friction coefficient
           :param k: Hookean spring constant
           :param T: Temperature
           :param dt: Integration time step
           :param num_nodes: Total number of nodes (deduced if 0)
           :param weights: Optional positive edge weights (defaults to 1.0)
         )pbdoc",
         py::arg("spring_from"),
         py::arg("spring_to"),
         py::arg("gamma") = 1.0,
         py::arg("k") = 1.0,
         py::arg("T") = 1.0,
         py::arg("dt") = 0.01,
         py::arg("num_nodes") = 0,
         py::arg("weights") = std::vector<double>{})
    .def(py::init<const Eigen::SparseMatrix<double>&,
                  double,
                  double,
                  double,
                  double>(),
         R"pbdoc(
           Construct from a pre-assembled sparse scalar graph Laplacian matrix L.

           :param L: Symmetric positive semidefinite graph Laplacian matrix
           :param gamma: Bead friction coefficient
           :param k: Hookean spring constant
           :param T: Temperature
           :param dt: Integration time step
         )pbdoc",
         py::arg("L"),
         py::arg("gamma") = 1.0,
         py::arg("k") = 1.0,
         py::arg("T") = 1.0,
         py::arg("dt") = 0.01)
    .def(py::init<const pe::Universe&, double, double, double, double>(),
         R"pbdoc(
           Construct from a Universe's atom-bond graph and initial coordinates.

           :param universe: Universe instance
           :param gamma: Bead friction coefficient
           :param k: Hookean spring constant
           :param T: Temperature
           :param dt: Integration time step
         )pbdoc",
         py::arg("universe"),
         py::arg("gamma") = 1.0,
         py::arg("k") = 1.0,
         py::arg("T") = 1.0,
         py::arg("dt") = 0.01)
    .def_static(
      "from_matrix",
      [](const Eigen::SparseMatrix<double>& L,
         double gamma,
         double k,
         double T,
         double dt) {
        return std::unique_ptr<RouseBrownianSimulator>(
          new RouseBrownianSimulator(
            RouseBrownianSimulator::fromMatrix(L, gamma, k, T, dt)));
      },
      R"pbdoc(
        Create a RouseBrownianSimulator from a sparse graph Laplacian matrix L.

        :param L: Symmetric positive semidefinite graph Laplacian matrix
        :param gamma: Bead friction coefficient
        :param k: Hookean spring constant
        :param T: Temperature
        :param dt: Integration time step
        :return: RouseBrownianSimulator instance
      )pbdoc",
      py::arg("L"),
      py::arg("gamma") = 1.0,
      py::arg("k") = 1.0,
      py::arg("T") = 1.0,
      py::arg("dt") = 0.01)
    .def_static(
      "from_universe",
      [](const pe::Universe& universe,
         double gamma,
         double k,
         double T,
         double dt) {
        return std::unique_ptr<RouseBrownianSimulator>(
          new RouseBrownianSimulator(
            RouseBrownianSimulator::fromUniverse(universe, gamma, k, T, dt)));
      },
      R"pbdoc(
        Create a RouseBrownianSimulator from a Universe's atom-bond graph and coordinates.

        :param universe: Universe instance
        :param gamma: Bead friction coefficient
        :param k: Hookean spring constant
        :param T: Temperature
        :param dt: Integration time step
        :return: RouseBrownianSimulator instance
      )pbdoc",
      py::arg("universe"),
      py::arg("gamma") = 1.0,
      py::arg("k") = 1.0,
      py::arg("T") = 1.0,
      py::arg("dt") = 0.01)
    .def_static(
      "create_linear_chain",
      [](size_t num_nodes, double gamma, double k, double T, double dt) {
        return std::unique_ptr<RouseBrownianSimulator>(
          new RouseBrownianSimulator(RouseBrownianSimulator::createLinearChain(
            num_nodes, gamma, k, T, dt)));
      },
      R"pbdoc(
        Create a RouseBrownianSimulator for a linear chain of num_nodes beads.
      )pbdoc",
      py::arg("num_nodes"),
      py::arg("gamma") = 1.0,
      py::arg("k") = 1.0,
      py::arg("T") = 1.0,
      py::arg("dt") = 0.01)
    .def_static(
      "create_ring",
      [](size_t num_nodes, double gamma, double k, double T, double dt) {
        return std::unique_ptr<RouseBrownianSimulator>(
          new RouseBrownianSimulator(
            RouseBrownianSimulator::createRing(num_nodes, gamma, k, T, dt)));
      },
      R"pbdoc(
        Create a RouseBrownianSimulator for a ring polymer of num_nodes beads.
      )pbdoc",
      py::arg("num_nodes"),
      py::arg("gamma") = 1.0,
      py::arg("k") = 1.0,
      py::arg("T") = 1.0,
      py::arg("dt") = 0.01)
    .def_static(
      "create_complete",
      [](size_t num_nodes, double gamma, double k, double T, double dt) {
        return std::unique_ptr<RouseBrownianSimulator>(
          new RouseBrownianSimulator(RouseBrownianSimulator::createComplete(
            num_nodes, gamma, k, T, dt)));
      },
      R"pbdoc(
        Create a RouseBrownianSimulator for a complete graph K_N of num_nodes beads.

        :param num_nodes: Number of nodes in the complete graph (at least 2)
        :param gamma: Bead friction coefficient
        :param k: Hookean spring constant
        :param T: Temperature
        :param dt: Integration time step
        :return: RouseBrownianSimulator instance
      )pbdoc",
      py::arg("num_nodes"),
      py::arg("gamma") = 1.0,
      py::arg("k") = 1.0,
      py::arg("T") = 1.0,
      py::arg("dt") = 0.01)
    .def_static(
      "create_star",
      [](size_t num_arms,
         size_t arm_length,
         double gamma,
         double k,
         double T,
         double dt) {
        return std::unique_ptr<RouseBrownianSimulator>(
          new RouseBrownianSimulator(RouseBrownianSimulator::createStar(
            num_arms, arm_length, gamma, k, T, dt)));
      },
      R"pbdoc(
        Create a RouseBrownianSimulator for a star polymer with num_arms arms of arm_length bonds each.
      )pbdoc",
      py::arg("num_arms"),
      py::arg("arm_length") = 1,
      py::arg("gamma") = 1.0,
      py::arg("k") = 1.0,
      py::arg("T") = 1.0,
      py::arg("dt") = 0.01)
    .def("get_nr_of_nodes",
         &RouseBrownianSimulator::getNrOfNodes,
         "Get the number of nodes (beads) in the simulator.")
    .def("get_gamma",
         &RouseBrownianSimulator::getGamma,
         "Get the friction coefficient gamma.")
    .def("get_k",
         &RouseBrownianSimulator::getK,
         "Get the Hookean spring constant k.")
    .def("get_temperature",
         &RouseBrownianSimulator::getTemperature,
         "Get the temperature T.")
    .def("get_dt",
         &RouseBrownianSimulator::getDt,
         "Get the integration time step dt.")
    .def("get_kb",
         &RouseBrownianSimulator::getKB,
         "Get the Boltzmann constant kB.")
    .def("get_diffusion_coefficient",
         &RouseBrownianSimulator::getDiffusionCoefficient,
         "Get the bead diffusion coefficient D = kB * T / gamma.")
    .def("get_solver_type",
         &RouseBrownianSimulator::getSolverType,
         "Get the active linear solver type.")
    .def("set_parameters",
         &RouseBrownianSimulator::setParameters,
         "Set gamma, k, T, and dt simultaneously.",
         py::arg("gamma"),
         py::arg("k"),
         py::arg("T"),
         py::arg("dt"))
    .def("set_gamma",
         &RouseBrownianSimulator::setGamma,
         "Set the friction coefficient gamma.",
         py::arg("gamma"))
    .def("set_k",
         &RouseBrownianSimulator::setK,
         "Set the Hookean spring constant k.",
         py::arg("k"))
    .def("set_temperature",
         &RouseBrownianSimulator::setTemperature,
         "Set the temperature T.",
         py::arg("T"))
    .def("set_dt",
         &RouseBrownianSimulator::setDt,
         "Set the integration time step dt.",
         py::arg("dt"))
    .def("set_kb",
         &RouseBrownianSimulator::setKB,
         "Set the Boltzmann constant kB.",
         py::arg("kb"))
    .def("set_solver_type",
         &RouseBrownianSimulator::setSolverType,
         "Set the linear solver algorithm used for implicit time-stepping.",
         py::arg("solver_type"))
    .def("set_seed",
         &RouseBrownianSimulator::setSeed,
         "Set the random number generator seed for thermal noise.",
         py::arg("seed"))
    .def("set_positions",
         &RouseBrownianSimulator::setPositions,
         "Set the (N, 3) bead positions matrix.",
         py::arg("positions"))
    .def("set_positions_flattened",
         &RouseBrownianSimulator::setPositionsFlattened,
         "Set bead positions from a flattened 3N vector.",
         py::arg("positions"),
         py::arg("interleaved") = true)
    .def("get_positions",
         &RouseBrownianSimulator::getPositions,
         "Get the (N, 3) bead positions matrix.")
    .def("get_positions_flattened",
         &RouseBrownianSimulator::getPositionsFlattened,
         "Get bead positions as a flattened 3N vector.",
         py::arg("interleaved") = true)
    .def("initialize_linear_positions",
         &RouseBrownianSimulator::initializeLinearPositions,
         "Initialize bead positions spaced along a coordinate axis.",
         py::arg("bond_length") = 1.0,
         py::arg("axis") = 0)
    .def("initialize_random_walk",
         &RouseBrownianSimulator::initializeRandomWalk,
         "Initialize bead positions as a 3D random walk.",
         py::arg("bond_length") = 1.0,
         py::arg("seed") = 42)
    .def("get_laplacian",
         &RouseBrownianSimulator::getLaplacian,
         "Get the sparse graph Laplacian matrix L.")
    .def(
      "get_shifted_matrix",
      &RouseBrownianSimulator::getShiftedMatrix,
      "Get the shifted implicit stepping matrix A = I + dt * (k / gamma) * L.")
    .def("step",
         &RouseBrownianSimulator::step,
         "Advance the simulation by num_steps implicit Brownian steps.",
         py::arg("num_steps") = 1)
    .def("step_with_external_forces",
         &RouseBrownianSimulator::stepWithExternalForces,
         "Advance the simulation by num_steps under external forces and "
         "thermal noise.",
         py::arg("forces"),
         py::arg("num_steps") = 1)
    .def("step_deterministic",
         &RouseBrownianSimulator::stepDeterministic,
         "Advance the simulation by num_steps deterministic relaxation steps "
         "(without thermal noise).",
         py::arg("num_steps") = 1)
    .def("step_deterministic_with_external_forces",
         &RouseBrownianSimulator::stepDeterministicWithExternalForces,
         "Advance the simulation by num_steps deterministic steps under "
         "external forces.",
         py::arg("forces"),
         py::arg("num_steps") = 1)
    .def("get_center_of_mass",
         &RouseBrownianSimulator::getCenterOfMass,
         "Compute the 3D center of mass of all nodes.")
    .def(
      "get_end_to_end_vector",
      &RouseBrownianSimulator::getEndToEndVector,
      "Compute the 3D vector from node_a to node_b (None defaults to N - 1).",
      py::arg("node_a") = 0,
      py::arg("node_b") = std::nullopt)
    .def("get_end_to_end_distance",
         &RouseBrownianSimulator::getEndToEndDistance,
         "Compute the Euclidean distance between node_a and node_b (None "
         "defaults to N - 1).",
         py::arg("node_a") = 0,
         py::arg("node_b") = std::nullopt)
    .def("get_radius_of_gyration",
         &RouseBrownianSimulator::getRadiusOfGyration,
         "Compute the radius of gyration of the current bead configuration.")
    .def("get_potential_energy",
         &RouseBrownianSimulator::getPotentialEnergy,
         "Compute the harmonic spring potential energy of the current "
         "configuration.")
    .def(
      "run",
      &RouseBrownianSimulator::run,
      "Run the simulation for total_steps and record dynamical observables "
      "every sample_interval steps (end_to_end_node_b=None defaults to N - 1).",
      py::arg("total_steps"),
      py::arg("sample_interval") = 1,
      py::arg("save_trajectory") = false,
      py::arg("end_to_end_node_a") = 0,
      py::arg("end_to_end_node_b") = std::nullopt)
    .def_static("compute_msd",
                &RouseBrownianSimulator::computeMSD,
                "Compute bead mean-squared displacement from a trajectory of "
                "(N, 3) position matrices.",
                py::arg("trajectory"),
                py::arg("time_averaged") = false)
    .def_static("compute_center_of_mass_msd",
                &RouseBrownianSimulator::computeCenterOfMassMSD,
                "Compute center-of-mass mean-squared displacement from a "
                "trajectory of (N, 3) position matrices.",
                py::arg("trajectory"),
                py::arg("time_averaged") = false)
    .def_static("compute_end_to_end_autocorrelation",
                &RouseBrownianSimulator::computeEndToEndAutocorrelation,
                "Compute normalized end-to-end vector autocorrelation from a "
                "trajectory of (N, 3) position matrices (node_b=None defaults "
                "to N - 1).",
                py::arg("trajectory"),
                py::arg("node_a") = 0,
                py::arg("node_b") = std::nullopt,
                py::arg("time_averaged") = false);
}

#endif /* PYBIND_CALC_H */
