#ifndef PYBIND_READERS_H
#define PYBIND_READERS_H

#include "../utils/DataFileParser.h"
#include "../utils/DumpFileParser.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;
using namespace pylimer_tools::utils;

// struct LazyDumpFileIterator {
//     LazyDumpFileIterator(const DumpFileParser &fileParser, py::object ref) :
//     fileParser(fileParser), ref(ref) {}

//     pe::Universe next()
//     {
//         if (fileParser.isFinishedReading())
//         {
//             throw py::stop_iteration();
//         }
//         return molecule[index++];
//     }

//     const DumpFileParser &fileParser;
//     py::object ref;
//     size_t index;
// };

void init_pylimer_bound_readers(py::module_ &m) {

  py::class_<DumpFileParser>(m, "DumpFileReader", R"pbdoc(
       A reader for LAMMPS's `dump` files.
  )pbdoc")
      .def(py::init<const std::string>(), py::arg("pathOfFileToRead"))
      .def("read", &DumpFileParser::read, "Read the whole file")
      .def("getLength", &DumpFileParser::getLength,
           "Get the number of sections (time-steps) in the file")
      .def("getStringValuesForAt", &DumpFileParser::getStringValuesForAt,
           "Get the values for the section `index`, the main header "
           "`headerKey` and the column (in the header) `column`.",
           py::arg("rowIndex"), py::arg("headerKey"), py::arg("columnIndex"))
      .def("getNumericValuesForAt", &DumpFileParser::getNumericValuesForAt,
           "Get the values for the section `index`, the main header "
           "`headerKey` and the column (in the header) `column`.")
      .def("hasKey", &DumpFileParser::hasKey,
           "Check whether the first section has the header specified",
           py::arg("headerKey"))
      .def("keyHasColumn", &DumpFileParser::keyHasColumn,
           "Check whether the header of the first section has the specified "
           "column",
           py::arg("headerKey"), py::arg("columnName"))
      .def("keyHasDirectionalColumn", &DumpFileParser::keyHasDirectionalColumn,
           "Check whether the header of the first section has all the three "
           "columns `{dirPraefix}{x|y|z}{dirSuffix}`.",
           py::arg("headerKey"), py::arg("dirPraefix") = "",
           py::arg("dirSuffix") = "");

  py::class_<DataFileParser>(m, "DataFileReader", R"pbdoc(
       A reader for LAMMPS's `write_data` files.
  )pbdoc")
      .def(py::init<>())
      .def("read", &DataFileParser::read, py::arg("pathOfFileToRead"))
      .def("getNrOfAtoms", &DataFileParser::getNrOfAtoms)
      .def("getNrOfAtomTypes", &DataFileParser::getNrOfAtomTypes)
      .def("getAtomIds", &DataFileParser::getAtomIds)
      .def("getMoleculeIds", &DataFileParser::getMoleculeIds)
      .def("getAtomTypes", &DataFileParser::getAtomTypes)
      .def("getAtomX", &DataFileParser::getAtomX)
      .def("getAtomY", &DataFileParser::getAtomY)
      .def("getAtomZ", &DataFileParser::getAtomZ)
      .def("getAtomNx", &DataFileParser::getAtomNx)
      .def("getAtomNy", &DataFileParser::getAtomNy)
      .def("getAtomNz", &DataFileParser::getAtomNz)
      .def("getMasses", &DataFileParser::getMasses)
      .def("getNrOfBonds", &DataFileParser::getNrOfBonds)
      .def("getNrOfBondTypes", &DataFileParser::getNrOfBondTypes)
      .def("getBondTypes", &DataFileParser::getBondTypes)
      .def("getBondFrom", &DataFileParser::getBondFrom)
      .def("getBondTo", &DataFileParser::getBondTo)
      .def("getLx", &DataFileParser::getLx)
      .def("getLowX", &DataFileParser::getLowX)
      .def("getHighX", &DataFileParser::getHighX)
      .def("getLy", &DataFileParser::getLy)
      .def("getLowY", &DataFileParser::getLowY)
      .def("getHighY", &DataFileParser::getHighY)
      .def("getLz", &DataFileParser::getLz)
      .def("getLowZ", &DataFileParser::getLowZ)
      .def("getHighZ", &DataFileParser::getHighZ);
}

#endif
