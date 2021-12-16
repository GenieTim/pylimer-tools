#ifndef PYBIND_READERS_H
#define PYBIND_READERS_H

#include "../utils/DumpFileParser.h"
#include "../utils/DataFileParser.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;
using namespace pylimer_tools::utils;

// struct LazyDumpFileIterator {
//     LazyDumpFileIterator(const DumpFileParser &fileParser, py::object ref) : fileParser(fileParser), ref(ref) {}

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

void init_pylimer_bound_readers(py::module_ &m)
{

    py::class_<DumpFileParser>(m, "DumpFileReader")
        .def(py::init<>())
        .def("read", &DumpFileParser::read, "Read a file")
        .def("getLength", &DumpFileParser::getLength, "Get the number of sections in the file")
        .def("getStringValuesForAt", &DumpFileParser::getStringValuesForAt, "Get the values for the section `index`, the main header `headerKey` and the column (in the header) `column`.")
        .def("getNumericValuesForAt", &DumpFileParser::getNumericValuesForAt, "Get the values for the section `index`, the main header `headerKey` and the column (in the header) `column`.")
        .def("hasKey", &DumpFileParser::hasKey, "Check whether the first section has the header specified")
        .def("keyHasColumn", &DumpFileParser::keyHasColumn, "Check whether the header of the first section has the specified column")
        .def("keyHasDirectionalColumn", &DumpFileParser::keyHasDirectionalColumn, "Check whether the header of the first section has all the three columns `{dirPraefix}{x|y|z}{dirSuffix}`.");

    py::class_<DataFileParser>(m, "DataFileReader")
        .def(py::init<>())
        .def("read", &DataFileParser::read)
        .def("getNrOfAtoms", &DataFileParser::getNrOfAtoms)
        .def("getNrOfAtomTypes", &DataFileParser::getNrOfAtomTypes)
        .def("getAtomIds", &DataFileParser::getAtomIds)
        .def("getMoleculeIds", &DataFileParser::getMoleculeIds)
        .def("getAtomTypes", &DataFileParser::getAtomTypes)
        .def("getAtomX", &DataFileParser::getAtomX)
        .def("getAtomY", &DataFileParser::getAtomY)
        .def("getAtomNx", &DataFileParser::getAtomNx)
        .def("getAtomNy", &DataFileParser::getAtomNy)
        .def("getMasses", &DataFileParser::getMasses)
        .def("getNrOfBonds", &DataFileParser::getNrOfBonds)
        .def("getNrOfBondTypes", &DataFileParser::getNrOfBondTypes)
        .def("getBondTypes", &DataFileParser::getBondTypes)
        .def("getBondFrom", &DataFileParser::getBondFrom)
        .def("getBondTo", &DataFileParser::getBondTo)
        .def("getLx", &DataFileParser::getLx)
        .def("getLy", &DataFileParser::getLy)
        .def("getLz", &DataFileParser::getLz);
}

#endif
