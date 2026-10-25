#ifndef PYBIND_CALC_H
#define PYBIND_CALC_H

#include "../calc/MMTanalysis.h"
#include "../calc/NormalModeAnalyzer.h"

#include <pybind11/eigen.h>
#include <pybind11/native_enum.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <sstream>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;
namespace calc = pylimer_tools::calc;

using namespace pylimer_tools::calc;

void
init_pylimer_bound_calc(py::module_& m)
{
  //   m.def("predict_gelation_point",
  //         &mmt::predictGelationPoint,
  //         "Predict the gelation point of a Universe");
  // m.def("computeExtentOfReaction", &mmt::computeExtentOfReaction, "Compute
  // extent of reaction");
  //   m.def("compute_stoichiometric_imbalance",
  //         &mmt::computeStoichiometricInbalance,
  //         "Compute stoichiometric imbalance");

  py::native_enum<calc::NormalModeSolver>(
    m,
    "NormalModeSolver",
    "enum.Enum",
    "Solver algorithms for normal mode analysis")
    .value("AUTO", calc::NormalModeSolver::AUTO)
    .value("DENSE_LAPACK", calc::NormalModeSolver::DENSE_LAPACK)
    .value("DENSE_EIGEN", calc::NormalModeSolver::DENSE_EIGEN)
    .value("SPARSE_SLICING", calc::NormalModeSolver::SPARSE_SLICING)
    .finalize();

  py::class_<NormalModeAnalyzer>(m,
                                 "NormalModeAnalyzer",
                                 py::module_local(),
                                 R"pbdoc(
    Compute the normal modes and predict the loss/storage moduli.

    Please cite :cite:t:`gusev_molecular_2024` if you use this method in your work.
    )pbdoc")
    .def(py::init<const std::vector<size_t>, const std::vector<size_t>>(),
         R"pbdoc(
         Initialize the NormalModeAnalyzer with the bonds (edges).

         Constructs the connectivity matrix from the given edges.
         
         :param spring_from: Vector of starting node indices for springs/bonds
         :param spring_to: Vector of ending node indices for springs/bonds
         )pbdoc",
         py::arg("spring_from"),
         py::arg("spring_to"))
    .def("get_matrix_size",
         &NormalModeAnalyzer::getMatrixSize,
         R"pbdoc(
         Get the size of the matrix (the maximum number of eigenvalues that could be queried).
         
         :return: Size of the connectivity matrix
         )pbdoc")
    .def("get_matrix",
         &NormalModeAnalyzer::getAssembledConnectivityMatrix,
         R"pbdoc(
         Get the assembled connectivity matrix.
         
         :return: The connectivity matrix
         )pbdoc")
    .def("find_sparse_eigenvalues",
         &NormalModeAnalyzer::findSparseEigenvalues,
         R"pbdoc(
         Find the k smallest eigenvalues using a sparse solver.
         
         :param nr_of_eigenvalues: Number of smallest eigenvalues to find
         :param compute_eigenvectors: Whether to also compute eigenvectors (default: False)
         :return: True if computation was successful
         )pbdoc",
         py::arg("nr_of_eigenvalues"),
         py::arg("compute_eigenvectors") = false)
    .def(
      "find_all_eigenvalues",
      [](NormalModeAnalyzer& self,
         const py::object& solver_obj,
         bool compute_eigenvectors,
         const bool decompose_components) {
        calc::NormalModeSolver solver = calc::NormalModeSolver::AUTO;
        if (py::isinstance<py::bool_>(solver_obj)) {
          // Backward compatibility: legacy signature was
          // find_all_eigenvalues(compute_eigenvectors=False)
          compute_eigenvectors = solver_obj.cast<bool>();
          solver = calc::NormalModeSolver::AUTO;
        } else if (py::isinstance<py::str>(solver_obj)) {
          std::string s = solver_obj.cast<std::string>();
          std::transform(s.begin(), s.end(), s.begin(), ::tolower);
          if (s == "auto") {
            solver = calc::NormalModeSolver::AUTO;
          } else if (s == "dense" || s == "dense_lapack") {
            solver = calc::NormalModeSolver::DENSE_LAPACK;
          } else if (s == "dense_eigen") {
            solver = calc::NormalModeSolver::DENSE_EIGEN;
          } else if (s == "sparse" || s == "sparse_slicing" || s == "slicing") {
            solver = calc::NormalModeSolver::SPARSE_SLICING;
          } else {
            throw std::invalid_argument("Unknown NormalModeSolver string: " +
                                        s);
          }
        } else if (py::isinstance<calc::NormalModeSolver>(solver_obj)) {
          solver = solver_obj.cast<calc::NormalModeSolver>();
        } else {
          throw py::type_error(
            "solver must be a NormalModeSolver enum or string ('auto', "
            "'dense_lapack', 'dense_eigen', 'sparse_slicing')");
        }
        self.computeAllEigenvalues(
          solver, compute_eigenvectors, decompose_components);
      },
      R"pbdoc(
         Find all eigenvalues using the specified solver.

         :param solver: Solver method to use (NormalModeSolver enum or string: 'auto', 'dense_lapack', 'dense_eigen', 'sparse_slicing'; default: NormalModeSolver.AUTO)
         :param compute_eigenvectors: Whether to also compute eigenvectors (default: False)
         :param decompose_components: Whether to decompose the graph into connected components (default: True)
         )pbdoc",
      py::arg("solver") = "auto",
      py::arg("compute_eigenvectors") = false,
      py::arg("decompose_components") = true)
#ifdef USE_CUDA
    .def("find_all_eigenvalues_gpu",
         &NormalModeAnalyzer::computeAllEigenvaluesGPU,
         R"pbdoc(
         Find all eigenvalues using a dense GPU solver (cuSOLVER on CUDA).
         
         :param compute_eigenvectors: Whether to also compute eigenvectors (default: False)
         :return: True if computation was successful
         )pbdoc",
         py::arg("compute_eigenvectors") = false)
#endif
    .def("get_eigenvalues",
         &NormalModeAnalyzer::getEigenvalues,
         R"pbdoc(
         Get the eigenvalues.
         
         :return: Vector of eigenvalues
         )pbdoc")
    .def("set_eigenvalues",
         &NormalModeAnalyzer::setEigenvalues,
         R"pbdoc(
         Set the eigenvalues, e.g. if you use an external solver.
         
         :param eigenvalues: Vector of eigenvalues to set
         )pbdoc",
         py::arg("eigenvalues"))
    .def("get_eigenvectors",
         &NormalModeAnalyzer::getEigenvectors,
         R"pbdoc(
         Get eigenvectors.
         
         :return: Matrix of eigenvectors
         )pbdoc")
    .def("set_eigenvectors",
         &NormalModeAnalyzer::setEigenvectors,
         R"pbdoc(
         Set eigenvectors, e.g. if you use an external solver.
         
         :param eigenvectors: Matrix of eigenvectors to set
         )pbdoc",
         py::arg("eigenvectors"))
    .def("evaluate_stress_autocorrelation",
         &NormalModeAnalyzer::evaluateStressAutocorrelation,
         R"pbdoc(
         Evaluate stress autocorrelation :math:`C(t)`.
         
         :param t: The time at which to evaluate the stress autocorrelation
         :return: Stress autocorrelation values
         )pbdoc",
         py::arg("t"))
    .def("evaluate_storage_modulus",
         &NormalModeAnalyzer::evaluateStorageModulus,
         R"pbdoc(
         Evaluate the storage modulus :math:`G'(\omega)`. Yet misses the conversion factor.
         
         :param omega: Angular frequencies
         :return: Storage modulus values
         )pbdoc",
         py::arg("omega"))
    .def("evaluate_loss_modulus",
         &NormalModeAnalyzer::evaluateLossModulus,
         R"pbdoc(
         Evaluate the loss modulus :math:`G''(\omega)`. Yet misses the conversion factor.
         
         :param omega: Angular frequencies
         :return: Loss modulus values
         )pbdoc",
         py::arg("omega"))
    .def("get_nr_of_soluble_clusters",
         &NormalModeAnalyzer::getNrOfSolubleClusters,
         R"pbdoc(
         Get the number of soluble clusters (eigenvalues = 0).
         
         :return: Number of soluble clusters
         )pbdoc")
    .def("compute_chebyshev_moments_kpm",
         &NormalModeAnalyzer::computeChebyshevMomentsKPM,
         R"pbdoc(
         Compute Chebyshev moments using the Kernel Polynomial Method (KPM).

         :param num_moments: Number of moments to compute (default: 1024)
         :param num_random_vectors: Number of random Rademacher vectors for trace estimation (default: 20)
         :return: Vector of damped Chebyshev moments
         )pbdoc",
         py::arg("num_moments") = 1024,
         py::arg("num_random_vectors") = 20)
    .def("compute_density_of_states_kpm",
         &NormalModeAnalyzer::computeDensityOfStatesKPM,
         R"pbdoc(
         Compute continuous density of states (DOS) using KPM.

         :param lambda_grid: Eigenvalue points at which to evaluate DOS
         :param num_moments: Number of Chebyshev moments (default: 1024)
         :param num_random_vectors: Number of random vectors (default: 20)
         :return: Density of states values on the grid
         )pbdoc",
         py::arg("lambda_grid"),
         py::arg("num_moments") = 1024,
         py::arg("num_random_vectors") = 20)
    .def("evaluate_stress_autocorrelation_kpm",
         &NormalModeAnalyzer::evaluateStressAutocorrelationKPM,
         R"pbdoc(
         Evaluate stress autocorrelation directly using KPM without computing eigenvalues.

         :param t: Time points to evaluate
         :param num_moments: Number of Chebyshev moments (default: 1024)
         :param num_random_vectors: Number of random vectors (default: 20)
         :param num_quadrature_points: Number of Gauss-Chebyshev quadrature points (default: 2048)
         :return: Stress autocorrelation values
         )pbdoc",
         py::arg("t"),
         py::arg("num_moments") = 1024,
         py::arg("num_random_vectors") = 20,
         py::arg("num_quadrature_points") = 2048)
    .def("evaluate_storage_modulus_kpm",
         &NormalModeAnalyzer::evaluateStorageModulusKPM,
         R"pbdoc(
         Evaluate storage modulus G' directly using KPM without computing eigenvalues.

         :param omega: Angular frequencies to evaluate
         :param num_moments: Number of Chebyshev moments (default: 1024)
         :param num_random_vectors: Number of random vectors (default: 20)
         :param num_quadrature_points: Number of Gauss-Chebyshev quadrature points (default: 2048)
         :return: Storage modulus values
         )pbdoc",
         py::arg("omega"),
         py::arg("num_moments") = 1024,
         py::arg("num_random_vectors") = 20,
         py::arg("num_quadrature_points") = 2048)
    .def("evaluate_loss_modulus_kpm",
         &NormalModeAnalyzer::evaluateLossModulusKPM,
         R"pbdoc(
         Evaluate loss modulus G'' directly using KPM without computing eigenvalues.

         :param omega: Angular frequencies to evaluate
         :param num_moments: Number of Chebyshev moments (default: 1024)
         :param num_random_vectors: Number of random vectors (default: 20)
         :param num_quadrature_points: Number of Gauss-Chebyshev quadrature points (default: 2048)
         :return: Loss modulus values
         )pbdoc",
         py::arg("omega"),
         py::arg("num_moments") = 1024,
         py::arg("num_random_vectors") = 20,
         py::arg("num_quadrature_points") = 2048)
#ifdef CEREALIZABLE
    .def(py::pickle(
      [](const NormalModeAnalyzer& u) {
        return py::make_tuple(
          py::bytes(pylimer_tools::utils::serializeToString(u)));
      },
      [](py::tuple t) {
        std::string in = t[0].cast<std::string>();
        return NormalModeAnalyzer::fromString(in);
      }));
#else
    ;
#endif
}

#endif /* PYBIND_CALC_H */
