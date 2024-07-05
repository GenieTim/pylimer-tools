#ifndef PYBIND_CALC_H
#define PYBIND_CALC_H

#include "../calc/MEHPanalysis.h"
#include "../calc/MMTanalysis.h"
#include "../calc/NormalModeAnalyzer.h"
#include "../entities/Universe.h"

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

  py::class_<NormalModeAnalyzer>(m, "NormalModeAnalyzer")
    .def(py::init<const std::vector<size_t>, const std::vector<size_t>>(),
         "Initialize NormalModeAnalyzer",
         py::arg("spring_from"),
         py::arg("spring_to"))
    .def("find_sparse_eigenvalues",
         &NormalModeAnalyzer::findSparseEigenvalues,
         "Find the `k` smallest eigenvalues using a sparse solver",
         py::arg("nr_of_eigenvalues"))
    .def("find_all_eigenvalues",
         &NormalModeAnalyzer::computeAllEigenvalues,
         "Find all eigenvalues using a dense solver",
         py::arg("compute_eigenvectors") = false)
    .def(
      "get_eigenvalues", &NormalModeAnalyzer::getEigenvalues, "Get eigenvalues")
    .def(
      "set_eigenvalues", &NormalModeAnalyzer::setEigenvalues, "Set eigenvalues")
    .def("get_eigenvectors",
         &NormalModeAnalyzer::getEigenvectors,
         "Get eigenvectors")
    .def("set_eigenvectors",
         &NormalModeAnalyzer::setEigenvectors,
         "Set eigenvectors")
    .def("evaluate_stress_autocorrelation",
         &NormalModeAnalyzer::evaluateStressAutocorrelation,
         "Evaluate stress autocorrelation",
         py::arg("t"))
    .def("evaluate_storage_modulus",
         &NormalModeAnalyzer::evaluateStorageModulus,
         "Evaluate storage modulus",
         py::arg("omega"))
    .def("evaluate_loss_modulus",
         &NormalModeAnalyzer::evaluateLossModulus,
         "Evaluate loss modulus",
         py::arg("omega"));
}

#endif /* PYBIND_CALC_H */
