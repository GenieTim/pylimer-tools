#ifndef PYBIND_READERS_H
#define PYBIND_READERS_H

#include "../io/DataFileParser.h"
#include "../io/DumpFileParser.h"
#include "../io/CSVSplitter.h"

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

void
init_pylimer_bound_readers(py::module_& m)
{

  py::enum_<AtomStyle>(m, "AtomStyle")
    .value("NONE", AtomStyle::NONE)
    .value("ANGLE", AtomStyle::ANGLE)
    .value("ATOMIC", AtomStyle::ATOMIC)
    .value("BODY", AtomStyle::BODY)
    .value("BOND", AtomStyle::BOND)
    .value("CHARGE", AtomStyle::CHARGE)
    .value("DIELECTRIC", AtomStyle::DIELECTRIC)
    .value("DIPOLE", AtomStyle::DIPOLE)
    .value("DPD", AtomStyle::DPD)
    .value("EDPD", AtomStyle::EDPD)
    .value("ELECTRON", AtomStyle::ELECTRON)
    .value("ELLIPSOID", AtomStyle::ELLIPSOID)
    .value("FULL", AtomStyle::FULL)
    .value("LINE", AtomStyle::LINE)
    .value("MDPD", AtomStyle::MDPD)
    .value("MOLECULAR", AtomStyle::MOLECULAR)
    .value("PERI", AtomStyle::PERI)
    .value("SMD", AtomStyle::SMD)
    .value("SPH", AtomStyle::SPH)
    .value("SPHERE", AtomStyle::SPHERE)
    .value("BPM_SPHERE", AtomStyle::BPM_SPHERE)
    .value("SPIN", AtomStyle::SPIN)
    .value("TDPD", AtomStyle::TDPD)
    .value("TEMPLATE", AtomStyle::TEMPLATE)
    .value("TRI", AtomStyle::TRI)
    .value("WAVEPACKET", AtomStyle::WAVEPACKET)
    .value("HYBRID", AtomStyle::HYBRID);

  py::class_<DumpFileParser>(m, "DumpFileReader", R"pbdoc(
       A reader for LAMMPS's `dump` files.
  )pbdoc")
    .def(py::init<const std::string>(), py::arg("pathOfFileToRead"))
    .def("read", &DumpFileParser::read, "Read the whole file")
    .def("getLength",
         &DumpFileParser::getLength,
         "Get the number of sections (time-steps) in the file")
    .def("getStringValuesForAt",
         &DumpFileParser::getStringValuesForAt,
         "Get the values for the section `index`, the main header "
         "`headerKey` and the column (in the header) `column`.",
         py::arg("rowIndex"),
         py::arg("headerKey"),
         py::arg("columnIndex"))
    .def("getNumericValuesForAt",
         &DumpFileParser::getNumericValuesForAt,
         "Get the values for the section `index`, the main header "
         "`headerKey` and the column (in the header) `column`.")
    .def("hasKey",
         &DumpFileParser::hasKey,
         "Check whether the first section has the header specified",
         py::arg("headerKey"))
    .def("keyHasColumn",
         &DumpFileParser::keyHasColumn,
         "Check whether the header of the first section has the specified "
         "column",
         py::arg("headerKey"),
         py::arg("columnName"))
    .def("keyHasDirectionalColumn",
         &DumpFileParser::keyHasDirectionalColumn,
         "Check whether the header of the first section has all the three "
         "columns `{dirPraefix}{x|y|z}{dirSuffix}`.",
         py::arg("headerKey"),
         py::arg("dirPraefix") = "",
         py::arg("dirSuffix") = "");

  py::class_<DataFileParser>(m, "DataFileReader", R"pbdoc(
       A reader for LAMMPS's `write_data` files.
  )pbdoc")
    .def(py::init<>())
    .def("read",
         &DataFileParser::read,
         R"pbdoc(
       Actually read a LAMMPS's `write_data` file.

       Arguments:
          - `path_of_file_to_read`: The path to the file to read
          - `atom_style`: The format of the "Atoms" section, see https://docs.lammps.org/read_data.html
          - `atom_style2`: The format of the "Atoms" section if the previous parameter is equal to AtomStyle::HYBRID
          - `atom_style3`: The format of the "Atoms" section if the second to last parameter is equal to AtomStyle::HYBRID
  )pbdoc",
         py::arg("path_of_file_to_read"),
         py::arg("atom_style") = AtomStyle::ANGLE,
         py::arg("atom_style2") = AtomStyle::NONE,
         py::arg("atom_style_3") = AtomStyle::NONE)
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

     m.def("splitCSV", &splitCSV, "Read a file containing a number of CSVs. Returns them split up.");
}

#endif
