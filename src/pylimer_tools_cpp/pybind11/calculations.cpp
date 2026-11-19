#ifndef PYBIND_CALC_H
#define PYBIND_CALC_H

#include "../calc/Correlator.h"
#include "../calc/MMTanalysis.h"
#include "../calc/NormalModeAnalyzer.h"
#include "../entities/Universe.h"

#include <pybind11/eigen.h>
#include <pybind11/native_enum.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>
#include <variant>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;
namespace calc = pylimer_tools::calc;

using namespace pylimer_tools::calc;

void
init_pylimer_bound_calc(py::module_& m)
{
  py::native_enum<calc::NormalModeSolver>(
    m,
    "NormalModeSolver",
    "enum.Enum",
    "Exact eigenvalue solver algorithms for normal mode analysis")
    .value("AUTO", calc::NormalModeSolver::AUTO)
    .value("DENSE_LAPACK", calc::NormalModeSolver::DENSE_LAPACK)
    .value("DENSE_EIGEN", calc::NormalModeSolver::DENSE_EIGEN)
    .value("SPARSE_SLICING", calc::NormalModeSolver::SPARSE_SLICING)
    .value("METAL", calc::NormalModeSolver::METAL)
    .value("CUDA", calc::NormalModeSolver::CUDA)
    .finalize();

  py::native_enum<calc::ComputeBackend>(
    m,
    "ComputeBackend",
    "enum.Enum",
    "Hardware compute backend for approximate spectral methods")
    .value("AUTO", calc::ComputeBackend::AUTO)
    .value("CPU", calc::ComputeBackend::CPU)
    .value("METAL", calc::ComputeBackend::METAL)
    .value("CUDA", calc::ComputeBackend::CUDA)
    .finalize();

  auto backendRepr = [](const ComputeBackend b) -> const char* {
    switch (b) {
      case ComputeBackend::AUTO:
        return "ComputeBackend.AUTO";
      case ComputeBackend::CPU:
        return "ComputeBackend.CPU";
      case ComputeBackend::METAL:
        return "ComputeBackend.METAL";
      case ComputeBackend::CUDA:
        return "ComputeBackend.CUDA";
      default:
        return "ComputeBackend.AUTO";
    }
  };

  py::class_<SLQConfig>(
    m,
    "SLQConfig",
    py::module_local(),
    "Configuration for Stochastic Lanczos Quadrature (SLQ).")
    .def(
      py::init([](const size_t numSteps,
                  const size_t numRandomVectors,
                  const ComputeBackend backend) {
        return SLQConfig{ numSteps, numRandomVectors, backend };
      }),
      py::arg("num_steps") = 200,
      py::arg("num_random_vectors") = 20,
      py::arg_v("backend", calc::ComputeBackend::AUTO, "ComputeBackend.AUTO"))
    .def_readwrite("num_steps", &SLQConfig::numSteps)
    .def_readwrite("num_random_vectors", &SLQConfig::numRandomVectors)
    .def_readwrite("backend", &SLQConfig::backend)
    .def("__repr__", [backendRepr](const SLQConfig& c) {
      std::ostringstream oss;
      oss << "SLQConfig(num_steps=" << c.numSteps
          << ", num_random_vectors=" << c.numRandomVectors
          << ", backend=" << backendRepr(c.backend) << ")";
      return oss.str();
    });

  py::class_<HybridConfig>(
    m,
    "HybridConfig",
    py::module_local(),
    "Configuration for Hybrid exact low-modes + deflated SLQ.")
    .def(
      py::init([](const size_t numSteps,
                  const size_t numRandomVectors,
                  const size_t numExactLowModes,
                  const ComputeBackend backend) {
        return HybridConfig{
          numSteps, numRandomVectors, numExactLowModes, backend
        };
      }),
      py::arg("num_steps") = 200,
      py::arg("num_random_vectors") = 20,
      py::arg("num_exact_low_modes") = 100,
      py::arg_v("backend", calc::ComputeBackend::AUTO, "ComputeBackend.AUTO"))
    .def_readwrite("num_steps", &HybridConfig::numSteps)
    .def_readwrite("num_random_vectors", &HybridConfig::numRandomVectors)
    .def_readwrite("num_exact_low_modes", &HybridConfig::numExactLowModes)
    .def_readwrite("backend", &HybridConfig::backend)
    .def("__repr__", [backendRepr](const HybridConfig& c) {
      std::ostringstream oss;
      oss << "HybridConfig(num_steps=" << c.numSteps
          << ", num_random_vectors=" << c.numRandomVectors
          << ", num_exact_low_modes=" << c.numExactLowModes
          << ", backend=" << backendRepr(c.backend) << ")";
      return oss.str();
    });

  py::class_<KPMConfig>(
    m,
    "KPMConfig",
    py::module_local(),
    "Configuration for Kernel Polynomial Method (KPM) spectral approximation.")
    .def(
      py::init([](const size_t numMoments,
                  const size_t numRandomVectors,
                  const size_t numQuadraturePoints,
                  const size_t batchSize,
                  const ComputeBackend backend) {
        return KPMConfig{
          numMoments, numRandomVectors, numQuadraturePoints, batchSize, backend
        };
      }),
      py::arg("num_moments") = 1024,
      py::arg("num_random_vectors") = 20,
      py::arg("num_quadrature_points") = 2048,
      py::arg("batch_size") = 4,
      py::arg_v("backend", calc::ComputeBackend::AUTO, "ComputeBackend.AUTO"))
    .def_readwrite("num_moments", &KPMConfig::numMoments)
    .def_readwrite("num_random_vectors", &KPMConfig::numRandomVectors)
    .def_readwrite("num_quadrature_points", &KPMConfig::numQuadraturePoints)
    .def_readwrite("batch_size", &KPMConfig::batchSize)
    .def_readwrite("backend", &KPMConfig::backend)
    .def("__repr__", [backendRepr](const KPMConfig& c) {
      std::ostringstream oss;
      oss << "KPMConfig(num_moments=" << c.numMoments
          << ", num_random_vectors=" << c.numRandomVectors
          << ", num_quadrature_points=" << c.numQuadraturePoints
          << ", batch_size=" << c.batchSize
          << ", backend=" << backendRepr(c.backend) << ")";
      return oss.str();
    });

  py::class_<NormalModeSpectrum>(m,
                                 "NormalModeSpectrum",
                                 py::module_local(),
                                 R"pbdoc(
    Immutable normal mode spectrum and viscoelastic observable evaluator.

    Represents any discrete or quadrature-weighted spectrum :math:`\{(\lambda_i, w_i)\}`.
    )pbdoc")
    .def(py::init<Eigen::VectorXd,
                  std::optional<Eigen::VectorXd>,
                  std::optional<Eigen::MatrixXd>,
                  bool>(),
         R"pbdoc(
         Construct a NormalModeSpectrum from eigenvalues and optional weights/eigenvectors.

         :param eigenvalues: 1D array of eigenvalues
         :param weights: Optional 1D array of quadrature weights (defaults to all 1.0)
         :param eigenvectors: Optional 2D array of eigenvectors (columns correspond to modes)
         :param is_complete: Whether this spectrum covers all modes of the system
         )pbdoc",
         py::arg("eigenvalues"),
         py::arg("weights") = std::nullopt,
         py::arg("eigenvectors") = std::nullopt,
         py::arg("is_complete") = true)
    .def_property_readonly("eigenvalues",
                           &NormalModeSpectrum::getEigenvalues,
                           "Eigenvalues of the spectrum.")
    .def_property_readonly("weights",
                           &NormalModeSpectrum::getWeights,
                           "Quadrature weights of the spectrum.")
    .def_property_readonly(
      "eigenvectors",
      &NormalModeSpectrum::getEigenvectors,
      "Eigenvectors matrix (columns are modes), or None if not computed.")
    .def_property_readonly("is_complete",
                           &NormalModeSpectrum::isComplete,
                           "Whether this spectrum is complete.")
    .def_property_readonly("num_modes",
                           &NormalModeSpectrum::getNumModes,
                           "Number of modes/quadrature nodes in the spectrum.")
    .def_property_readonly("num_zero_modes",
                           &NormalModeSpectrum::getNumZeroModes,
                           "Number of zero modes in the spectrum.")
    .def("__len__", &NormalModeSpectrum::getNumModes)
    .def("stress_autocorrelation",
         &NormalModeSpectrum::stressAutocorrelation,
         R"pbdoc(
         Evaluate stress autocorrelation :math:`G(t)`.

         :param t: Time points at which to evaluate
         :return: Stress autocorrelation values
         )pbdoc",
         py::arg("t"))
    .def("storage_modulus",
         &NormalModeSpectrum::storageModulus,
         R"pbdoc(
         Evaluate the storage modulus :math:`G'(\omega)` (without prefactor).

         :param omega: Angular frequencies
         :return: Storage modulus values
         )pbdoc",
         py::arg("omega"))
    .def("loss_modulus",
         &NormalModeSpectrum::lossModulus,
         R"pbdoc(
         Evaluate the loss modulus :math:`G''(\omega)` (without prefactor).

         :param omega: Angular frequencies
         :return: Loss modulus values
         )pbdoc",
         py::arg("omega"))
    .def("moduli",
         &NormalModeSpectrum::moduli,
         R"pbdoc(
         Evaluate both storage and loss moduli :math:`(G'(\omega), G''(\omega))`
         in a single pass.

         :param omega: Angular frequencies
         :return: Tuple of (storage_modulus, loss_modulus)
         )pbdoc",
         py::arg("omega"))
#ifdef CEREALIZABLE
    .def(py::pickle(
      [](const NormalModeSpectrum& s) {
        return py::make_tuple(
          py::bytes(pylimer_tools::utils::serializeToString(s)));
      },
      [](py::tuple t) {
        const std::string in = t[0].cast<std::string>();
        return NormalModeSpectrum::fromString(in);
      }))
#endif
    ;

  py::class_<NormalModeAnalyzer>(m,
                                 "NormalModeAnalyzer",
                                 py::module_local(),
                                 R"pbdoc(
    Stateless graph Laplacian builder and normal mode solver.

    Please cite :cite:t:`gusev_molecular_2024` if you use this method in your work.
    )pbdoc")
    .def(py::init<const std::vector<size_t>&,
                  const std::vector<size_t>&,
                  size_t,
                  const std::optional<std::vector<double>>&>(),
         R"pbdoc(
         Initialize the NormalModeAnalyzer from bond/spring edge lists.

         :param spring_from: Starting node indices for springs/bonds
         :param spring_to: Ending node indices for springs/bonds
         :param nr_of_nodes: Optional total number of nodes (to include isolated vertices)
         :param weights: Optional per-spring weights (defaults to unit weight 1.0)
         )pbdoc",
         py::arg("spring_from"),
         py::arg("spring_to"),
         py::arg("nr_of_nodes") = 0,
         py::arg("weights") = std::nullopt)
    .def_static("from_matrix",
                &NormalModeAnalyzer::fromMatrix,
                R"pbdoc(
                Construct a NormalModeAnalyzer directly from a sparse Laplacian/connectivity matrix.

                :param matrix: Square sparse connectivity matrix
                :return: NormalModeAnalyzer instance
                )pbdoc",
                py::arg("matrix"))
    .def_static("from_universe",
                &NormalModeAnalyzer::fromUniverse,
                R"pbdoc(
                Construct a NormalModeAnalyzer from a Universe's atom-bond graph.

                :param universe: Universe instance
                :return: NormalModeAnalyzer instance
                )pbdoc",
                py::arg("universe"))
    .def_property_readonly("size",
                           &NormalModeAnalyzer::getSize,
                           "Dimension of the connectivity matrix.")
    .def_property_readonly("matrix",
                           &NormalModeAnalyzer::getMatrix,
                           "Assembled sparse connectivity matrix.")
    .def_property_readonly(
      "num_soluble_clusters",
      &NormalModeAnalyzer::getNumSolubleClusters,
      "Number of connected components (soluble clusters / zero modes).")
    .def("solve_exact",
         &NormalModeAnalyzer::solveExact,
         R"pbdoc(
         Compute exact eigenvalues (and optionally eigenvectors).

         :param solver: Solver algorithm (NormalModeSolver enum; default: NormalModeSolver.AUTO). Ignored when ``num_modes < size``.
         :param compute_eigenvectors: Whether to also compute eigenvectors (default: False)
         :param decompose_components: Whether to decompose into connected components (default: True). Ignored when ``num_modes < size``.
         :param num_modes: Optional number of smallest modes to compute via Spectra sparse solver
         :return: NormalModeSpectrum instance
         )pbdoc",
         py::arg_v(
           "solver", calc::NormalModeSolver::AUTO, "NormalModeSolver.AUTO"),
         py::arg("compute_eigenvectors") = false,
         py::arg("decompose_components") = true,
         py::arg("num_modes") = std::nullopt)
    .def(
      "solve_approx",
      [](const NormalModeAnalyzer& self,
         const std::optional<std::variant<SLQConfig, HybridConfig, KPMConfig>>&
           config) {
        if (!config.has_value()) {
          return self.solveApprox();
        }
        return std::visit(
          [&self](const auto& cfg) { return self.solveApprox(cfg); }, *config);
      },
      R"pbdoc(
      Compute an approximate / stochastic normal mode spectrum.

      :param config: SLQConfig, HybridConfig, KPMConfig, or None for AUTO
      :return: NormalModeSpectrum instance
      )pbdoc",
      py::arg("config") = std::nullopt)
    .def("chebyshev_moments",
         &NormalModeAnalyzer::chebyshevMoments,
         R"pbdoc(
         Compute Jackson-damped Chebyshev moments using KPM.

         :param config: KPMConfig instance
         :return: Vector of damped Chebyshev moments
         )pbdoc",
         py::arg_v("config", KPMConfig{}, "KPMConfig()"))
    .def("density_of_states",
         &NormalModeAnalyzer::densityOfStates,
         R"pbdoc(
         Compute continuous density of states (DOS) on lambda_grid using KPM.

         :param lambda_grid: Eigenvalue points at which to evaluate DOS
         :param config: KPMConfig instance
         :return: Density of states values on the grid
         )pbdoc",
         py::arg("lambda_grid"),
         py::arg_v("config", KPMConfig{}, "KPMConfig()"))
#ifdef CEREALIZABLE
    .def(py::pickle(
      [](const NormalModeAnalyzer& u) {
        return py::make_tuple(
          py::bytes(pylimer_tools::utils::serializeToString(u)));
      },
      [](py::tuple t) {
        const std::string in = t[0].cast<std::string>();
        return NormalModeAnalyzer::fromString(in);
      }));
#else
    ;
#endif

  py::class_<Correlator, py::smart_holder>(m,
                                           "Correlator",
                                           py::module_local(),
                                           R"pbdoc(
    Multi-tau scalar autocorrelator computing :math:`f(\tau) = \langle A(t) A(t + \tau) \rangle`.

    Uses a hierarchy of coarse-grained shift registers (Likhtman/Ramirez multi-tau algorithm)
    to evaluate time correlations across exponentially spanning lag times in O(1) amortized time per sample.
    )pbdoc")
    .def(py::init<unsigned int, unsigned int, unsigned int>(),
         R"pbdoc(
         Construct a multi-tau scalar Correlator.

         :param num_correlators: Number of hierarchical correlator levels (default: 32)
         :param p: Number of lag points per correlator level (default: 16)
         :param m: Averaging factor between consecutive levels; recommended ``p % m == 0`` (default: 2)
         )pbdoc",
         py::arg("num_correlators") = 32,
         py::arg("p") = 16,
         py::arg("m") = 2)
    .def("set_size",
         &Correlator::setsize,
         R"pbdoc(
         Reconfigure correlator dimensions and reset accumulated state.

         :param num_correlators: Number of hierarchical correlator levels (default: 32)
         :param p: Number of lag points per correlator level (default: 16)
         :param m: Averaging factor between consecutive levels (default: 2)
         )pbdoc",
         py::arg("num_correlators") = 32,
         py::arg("p") = 16,
         py::arg("m") = 2)
    .def("reset",
         &Correlator::reset,
         "Reset accumulated state while keeping the current correlator "
         "dimensions.")
    .def("add",
         &Correlator::add,
         R"pbdoc(
         Add a scalar observation to correlator level ``level``.

         :param value: Finite scalar observation
         :param level: Correlator level index (default: 0)
         )pbdoc",
         py::arg("value"),
         py::arg("level") = 0)
    .def("add_batch",
         &Correlator::addBatch,
         R"pbdoc(
         Add a 1D array of scalar observations at level 0.

         :param values: 1D array of finite scalar observations
         )pbdoc",
         py::arg("values"))
    .def("evaluate",
         &Correlator::evaluate,
         R"pbdoc(
         Evaluate the autocorrelation function at all populated lag times.

         :param norm: If True, subtract the squared sample mean :math:`\langle A \rangle^2` (default: False)
         )pbdoc",
         py::arg("norm") = false)
    .def_property_readonly("times",
                           &Correlator::getTimes,
                           "Evaluated lag times (call evaluate() first).")
    .def_property_readonly(
      "correlation",
      &Correlator::getCorrelation,
      "Evaluated correlation values (call evaluate() first).")
    .def_property_readonly("num_correlators",
                           &Correlator::getNumCorrelators,
                           "Number of correlator levels configured.")
    .def_property_readonly("num_samples",
                           &Correlator::getNumSamples,
                           "Total number of samples added at level 0.")
    .def_property_readonly(
      "mean", &Correlator::getMean, "Running sample mean at level 0.")
    .def_readonly("npcorr",
                  &Correlator::npcorr,
                  "Number of populated correlation points after evaluate().")
    .def_readonly(
      "p", &Correlator::p, "Number of lag points per correlator level.")
    .def_readonly(
      "m", &Correlator::m, "Averaging factor between consecutive levels.")
    .def_readonly(
      "accval", &Correlator::accval, "Accumulated sum of level-0 samples.")
    .def_readonly(
      "raw_t", &Correlator::t, "Full pre-allocated lag time buffer.")
    .def_readonly(
      "raw_f", &Correlator::f, "Full pre-allocated correlation buffer.")
#ifdef CEREALIZABLE
    .def(py::pickle(
      [](const Correlator& c) {
        return py::make_tuple(
          py::bytes(pylimer_tools::utils::serializeToString(c)));
      },
      [](py::tuple t) {
        if (t.size() != 1) {
          throw std::runtime_error("Invalid Correlator pickle state");
        }
        const std::string in = t[0].cast<std::string>();
        return Correlator::fromString(in);
      }));
#else
    ;
#endif
}

#endif /* PYBIND_CALC_H */
