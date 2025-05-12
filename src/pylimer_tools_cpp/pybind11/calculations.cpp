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
  m.def("predict_gelation_point",
        &mmt::predictGelationPoint,
        "Predict the gelation point of a Universe");
  // m.def("computeExtentOfReaction", &mmt::computeExtentOfReaction, "Compute
  // extent of reaction");
  m.def("compute_stoichiometric_imbalance",
        &mmt::computeStoichiometricInbalance,
        "Compute stoichiometric imbalance");

  py::class_<NormalModeAnalyzer>(
    m,
    "NormalModeAnalyzer",
    py::module_local(),
    "Compute the normal modes and loss/storage moduli.")
    .def(py::init<const std::vector<size_t>, const std::vector<size_t>>(),
         "Initialize NormalModeAnalyzer",
         py::arg("spring_from"),
         py::arg("spring_to"))
    .def("get_matrix_size",
         &NormalModeAnalyzer::getMatrixSize,
         "Get the size of the matrix (the maximum of eigenvalues that could be "
         "queried)")
    .def("get_matrix",
         &NormalModeAnalyzer::getAssembledConnectivityMatrix,
         "Get the assembled connectivity matrix.")
    .def("find_sparse_eigenvalues",
         &NormalModeAnalyzer::findSparseEigenvalues,
         "Find the `k` smallest eigenvalues using a sparse solver",
         py::arg("nr_of_eigenvalues"),
         py::arg("compute_eigenvectors") = false)
    .def("find_all_eigenvalues",
         &NormalModeAnalyzer::computeAllEigenvalues,
         "Find all eigenvalues using a dense solver",
         py::arg("compute_eigenvectors") = false)
    .def("get_eigenvalues",
         &NormalModeAnalyzer::getEigenvalues,
         "Get the eigenvalues")
    .def("set_eigenvalues",
         &NormalModeAnalyzer::setEigenvalues,
         "Set the eigenvalues, e.g. if you use an external solver",
         py::arg("eigenvalues"))
    .def("get_eigenvectors",
         &NormalModeAnalyzer::getEigenvectors,
         "Get eigenvectors")
    .def("set_eigenvectors",
         &NormalModeAnalyzer::setEigenvectors,
         "Set eigenvectors, e.g. if you use an external solver",
         py::arg("eigenvectors"))
    .def("evaluate_stress_autocorrelation",
         &NormalModeAnalyzer::evaluateStressAutocorrelation,
         "Evaluate stress autocorrelation",
         py::arg("t"))
    .def("evaluate_storage_modulus",
         &NormalModeAnalyzer::evaluateStorageModulus,
         "Evaluate the storage modulus :math:`G'`. Yet misses the conversion "
         "factor.",
         py::arg("omega"))
    .def("evaluate_loss_modulus",
         &NormalModeAnalyzer::evaluateLossModulus,
         "Evaluate the loss modulus :math:`G''`. Yet misses the conversion "
         "factor.",
         py::arg("omega"))
    .def("get_nr_of_soluble_clusters",
         &NormalModeAnalyzer::getNrOfSolubleClusters,
         "Get the number of soluble clusters (Eigenvalues = 0)")
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
