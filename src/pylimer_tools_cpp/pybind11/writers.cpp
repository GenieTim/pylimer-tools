
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
      .def(py::init<pe::Universe>(), py::arg("universe"), R"pbdoc(
           Initialize the writer with the universe to write.
      )pbdoc")
      .def("setUniverseToWrite", &DataFileWriter::setUniverseToWrite, R"pbdoc(
           Re-set the universe to write.
      )pbdoc",
           py::arg("universe"))
      .def("configIncludeAngles", &DataFileWriter::configIncludeAngles, R"pbdoc(
           Set whether to include the angles from the universe in the file or not.
      )pbdoc",
           py::arg("includeAngles") = true)
      .def("configReindexAtoms", &DataFileWriter::configReindexAtoms, R"pbdoc(
           Set whether to reindex the atoms or not. 
           Re-indexing leads to atom ids being in the range of 1 to the number of atoms.
      )pbdoc",
           py::arg("reindexAtoms") = true)
      .def("configCrosslinkerType", &DataFileWriter::configCrosslinkerType,
           R"pbdoc(
           Set which atom type represents cross-linkers. 
           Needed in case the moleculeIdx in the output file should have any meaning.
           (e.g. with :func:`~pylimer_tools_cpp.pylimer_tools_cpp.DataFileWriter.configMoleculeIdxForSwap`).
      )pbdoc",
           py::arg("crosslinkerType") = 2)
      .def("configMoleculeIdxForSwap",
           &DataFileWriter::configMoleculeIdxForSwap,
           py::arg("enableSwappability") = false, R"pbdoc(
                Swappable chains implies that their `moleculeIdx` in the LAMMPS data file is not 
                identical per chain, but identical per position in the chain.
                That's how you can have bond swapping with constant chain length distribution.
           )pbdoc")
      .def("writeToFile", &DataFileWriter::writeToFile);
}

#endif
