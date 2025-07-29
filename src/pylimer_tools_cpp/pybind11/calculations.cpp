#ifndef PYBIND_CALC_H
#define PYBIND_CALC_H

#include "../calc/MMTanalysis.h"
#include "../calc/NormalModeAnalyzer.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;

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
    .def("find_all_eigenvalues",
         &NormalModeAnalyzer::computeAllEigenvalues,
         R"pbdoc(
         Find all eigenvalues using a dense solver.
         
         :param compute_eigenvectors: Whether to also compute eigenvectors (default: False)
         :return: True if computation was successful
         )pbdoc",
         py::arg("compute_eigenvectors") = false)
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
#ifdef CEREALIZABLE
    .def(py::pickle(
      [](const NormalModeAnalyzer& u) {
        return py::make_tuple(pylimer_tools::utils::serializeToString(u));
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
