
#ifndef PYBIND_READERS_H
#define PYBIND_READERS_H

#include "../utils/DataFileWriter.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;
using namespace pylimer_tools::utils;

void init_pylimer_bound_writers(py::module_ &m) {
  py::class_<DataFileWriter>(m, "DataFileWriter").def(py::init<>())
  .def("setUniverseToWrite", &DataFileWriter::setUniverseToWrite)
  // .def("configIncludeAngles", &DataFileWriter::configIncludeAngles)
  // .def("configMoleculeIdxForSwap", &DataFileWriter::configMoleculeIdxForSwap)
  .def("writeToFile", &DataFileWriter::writeToFile);
}
