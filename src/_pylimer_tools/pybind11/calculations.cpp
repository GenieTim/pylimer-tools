#ifndef PYBIND_CALC_H
#define PYBIND_CALC_H

#include "../calc/MEHPanalysis.h"
#include "../calc/MMTanalysis.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pylimer_tools::calc;

void init_pylimer_bound_calc(py::module_ &m)
{
  m.doc() = "Calculations ";

  m.def("predictGelationPoint", &mmt::predictGelationPoint, "Predict the gelation point of a Universe");
  // m.def("computeExtentOfReaction", &mmt::computeExtentOfReaction, "Compute extent of reaction");
  m.def("computeStoichiometricInbalance", &mmt::computeStoichiometricInbalance, "Compute stoichiometric inbalance");
}

#endif /* PYBIND_CALC_H */
