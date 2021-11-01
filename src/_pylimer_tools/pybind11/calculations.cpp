#ifndef PYBIND_CALC_H
#define PYBIND_CALC_H

#include "../calc/MEHPanalysis.h"
#include "../calc/MMTanalysis.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pylimer_tools::entities;

PYBIND11_MODULE(pylimer_tools, m)
{
  m.def("mmt.predictGelationPoint", &predictGelationPoint);
  m.def("mmt.computeExtentOfReaction", &computeExtentOfReaction);
  m.def("mmt.computeStoichiometricInbalance", &computeStoichiometricInbalance);
}

#endif /* PYBIND_CALC_H */
