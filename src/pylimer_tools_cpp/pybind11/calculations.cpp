#ifndef PYBIND_CALC_H
#define PYBIND_CALC_H

#include "../calc/MEHPForceRelaxation.h"
#include "../calc/MEHPanalysis.h"
#include "../calc/MMTanalysis.h"
#include "../entities/Universe.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pylimer_tools::calc;
namespace pe = pylimer_tools::entities;

void init_pylimer_bound_calc(py::module_ &m) {
  m.def("predictGelationPoint", &mmt::predictGelationPoint,
        "Predict the gelation point of a Universe");
  // m.def("computeExtentOfReaction", &mmt::computeExtentOfReaction, "Compute
  // extent of reaction");
  m.def("computeStoichiometricInbalance", &mmt::computeStoichiometricInbalance,
        "Compute stoichiometric inbalance");

  py::class_<mehp::MEHPForceRelaxation>(m, "MEHPForceRelaxation", R"pbdoc(
    A small simulation tool for quickly minimizing the force between the cross-linker beads.
  )pbdoc")
      .def(py::init<pe::Universe, int>())
      .def("runForceRelaxation", &mehp::MEHPForceRelaxation::runForceRelaxation,
           "Run the simulation")
      .def("getVolume", &mehp::MEHPForceRelaxation::getVolume)
      .def("getNrOfNodes", &mehp::MEHPForceRelaxation::getNrOfNodes)
      .def("getNrOfSprings", &mehp::MEHPForceRelaxation::getNrOfSprings)
      .def("getNrOfActiveNodes", &mehp::MEHPForceRelaxation::getNrOfActiveNodes)
      .def("getNrOfActiveSprings",
           &mehp::MEHPForceRelaxation::getNrOfActiveSprings)
      .def("getAverageSpringLength",
           &mehp::MEHPForceRelaxation::getAverageSpringLength)
      .def("getFinalNrOfLoops", &mehp::MEHPForceRelaxation::getFinalNrOfLoops)
      .def("getInitialShearModulus",
           &mehp::MEHPForceRelaxation::getInitialShearModulus)
      .def("getFinalShearModulus",
           &mehp::MEHPForceRelaxation::getFinalShearModulus)
      .def("getInitialSquareDistance",
           &mehp::MEHPForceRelaxation::getInitialSquareDistance)
      .def("getFinalSquareDistance",
           &mehp::MEHPForceRelaxation::getFinalSquareDistance)
      .def("getInitialSquareRelativeDistance",
           &mehp::MEHPForceRelaxation::getInitialSquareRelativeDistance)
      .def("getFinalSquareRelativeDistance",
           &mehp::MEHPForceRelaxation::getFinalSquareRelativeDistance);
}

#endif /* PYBIND_CALC_H */
