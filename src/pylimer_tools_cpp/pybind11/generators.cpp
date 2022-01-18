#ifndef PYBIND_GENERATORS_H
#define PYBIND_GENERATORS_H

#include "../utils/MCUniverseGenerator.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;
using namespace pylimer_tools::utils;

void init_pylimer_bound_generators(py::module_ &m) {

  py::class_<MCUniverseGenerator>(m, "MCUniverseGenerator")
      .def(py::init<const double, const double, const double>())
      .def("setSeed", &MCUniverseGenerator::setSeed, "Set the seed for the random generator.")
      .def("setBeadDistance", &MCUniverseGenerator::setBeadDistance, "Set the optimal distance between beads.")
      .def("addCrosslinkers", &MCUniverseGenerator::addCrosslinkers, "Add the cross-linkers.")
      .def("addSolventChains", &MCUniverseGenerator::addSolventChains, "Randomly distribute additional, free chains.")
      .def("addAndLinkStrands", py::overload_cast<int, std::vector<int>,
                             double, int, int>(&MCUniverseGenerator::addAndLinkStrands), "Actually add strands, link them to the previously added cross-linkers.")
      .def("getUniverse", &MCUniverseGenerator::getUniverse, "Fetch the current (or final) state of the universe.");
}

#endif
