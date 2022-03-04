#ifndef PYBIND_GENERATORS_H
#define PYBIND_GENERATORS_H

#include "../utils/MCUniverseGenerator.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;
using namespace pylimer_tools::utils;

void init_pylimer_bound_generators(py::module_ &m) {

  py::class_<MCUniverseGenerator>(m, "MCUniverseGenerator", R"pbdoc(
       A :obj:`pylimer_tools_cpp.pylimer_tools_cpp.Universe` generator using a Monte-Carlo procedure.
  )pbdoc")
      .def(py::init<const double, const double, const double>(), py::arg("Lx"),
           py::arg("Ly"), py::arg("Lz"))
      .def("setSeed", &MCUniverseGenerator::setSeed,
           "Set the seed for the random generator.", py::arg("seed"))
      .def("setBeadDistance", &MCUniverseGenerator::setBeadDistance,
           "Set the optimal distance between beads.", py::arg("distance"))
      .def("addCrosslinkers", &MCUniverseGenerator::addCrosslinkers, R"pbdoc(
            Add the cross-linkers.
            )pbdoc",
           py::arg("nrOfCrosslinkers"), py::arg("crosslinkerAtomType") = 2)
      .def("addSolventChains", &MCUniverseGenerator::addSolventChains, R"pbdoc(
            Randomly distribute additional, free chains.
            )pbdoc",
           py::arg("nrOfSolventChains"), py::arg("solventChainLength"),
           py::arg("solventAtomType") = 3)
      .def("addAndLinkStrands",
           py::overload_cast<int, std::vector<int>, double, int, int>(
               &MCUniverseGenerator::addAndLinkStrands),
           R"pbdoc(
            Actually add strands, link them to the previously added cross-linkers.
            )pbdoc",
           py::arg("nrOfStrands"), py::arg("strandLengths"),
           py::arg("crosslinkerConversion"),
           py::arg("crosslinkerFunctionality"), py::arg("strandAtomType") = 1)
      .def("getUniverse", &MCUniverseGenerator::getUniverse, R"pbdoc(
            Fetch the current (or final) state of the universe.
            )pbdoc");
}

#endif
