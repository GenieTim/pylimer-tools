
#ifndef PYBIND_WRITERS_H
#define PYBIND_WRITERS_H

#include "../utils/DataFileWriter.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;
using namespace pylimer_tools::utils;

void init_pylimer_bound_writers(py::module_ &m) {
  py::class_<DataFileWriter>(m, "DataFileWriter")
      .def(py::init<pe::Universe>(), py::arg("universe"))
      .def("setUniverseToWrite", &DataFileWriter::setUniverseToWrite,
           py::arg("universe"))
      .def("configIncludeAngles", &DataFileWriter::configIncludeAngles,
           py::arg("includeAngles"))
      .def("configReindexAtoms", &DataFileWriter::configReindexAtoms,
           py::arg("reindexAtoms"))
      .def("configCrosslinkerType", &DataFileWriter::configCrosslinkerType,
           py::arg("crosslinkerType"))
      .def("configMoleculeIdxForSwap",
           &DataFileWriter::configMoleculeIdxForSwap,
           py::arg("enableSwappability"), R"pbdoc(
                Swappable chains implies that their `moleculeIdx` in the LAMMPS data file is not 
                identical per chain, but identical per position in the chain.
                That's how you can have bond swapping with constant chain length distribution.
           )pbdoc")
      .def("writeToFile", &DataFileWriter::writeToFile);
}

#endif
