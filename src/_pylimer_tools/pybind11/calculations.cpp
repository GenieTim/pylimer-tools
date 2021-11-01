#ifndef PYBIND_CALC_H
#define PYBIND_CALC_H

#include "../calc/MEHPanalysis.h"
#include "../calc/MMTanalysis.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pylimer_tools::calc;

PYBIND11_MODULE(pylimer_calc, m)
{
  m.def("mmt.predictGelationPoint", &mmt::predictGelationPoint, "Predict the gelation point of a Universe");
  m.def("mmt.computeExtentOfReaction", &mmt::computeExtentOfReaction, "Compute extent of reaction");
  m.def("mmt.computeStoichiometricInbalance", &mmt::computeStoichiometricInbalance, "Compute stoichiometric inbalance");
}

#endif /* PYBIND_CALC_H */
