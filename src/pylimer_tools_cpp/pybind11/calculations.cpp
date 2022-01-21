#ifndef PYBIND_CALC_H
#define PYBIND_CALC_H

#include "../calc/MEHPanalysis.h"
#include "../calc/MEHPForceRelaxation.h"
#include "../calc/MMTanalysis.h"
#include "../entities/Universe.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pylimer_tools::calc;
namespace pe = pylimer_tools::entities;

void init_pylimer_bound_calc(py::module_ &m)
{
  m.def("predictGelationPoint", &mmt::predictGelationPoint, "Predict the gelation point of a Universe");
  // m.def("computeExtentOfReaction", &mmt::computeExtentOfReaction, "Compute extent of reaction");
  m.def("computeStoichiometricInbalance", &mmt::computeStoichiometricInbalance, "Compute stoichiometric inbalance");

  py::class_<mehp::MEHPForceRelaxation>(m, "MEHPForceRelaxation", R"pbdoc(
    A small simulation tool for quickly minimizing the force between the cross-linker beads.
  )pbdoc").def(py::init<pe::Universe>())
  .def("runForceRelaxation", &mehp::MEHPForceRelaxation::runForceRelaxation, "Run the simulation");
}

#endif /* PYBIND_CALC_H */
