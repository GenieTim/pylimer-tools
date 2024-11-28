#ifndef PYBIND_READERS_H
#define PYBIND_READERS_H

#include "../io/AveFileReader.h"
#include "../io/CSVSplitter.h"
#include "../io/DataFileParser.h"
#include "../io/DumpFileParser.h"

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

    py::class_<AveFileReader>(m, "AveFileReader", R"pbdoc(
          Alternative implementation of the data file reader implemented in 
          :func:`pylimer_tools.read_lammps_output_file.read_averages_file`.

          This implementation is better for certain use cases, worse for others.
          In the end, only performance and memory usage are different.
          For moderately sized and small files, we recommend to use the Python interface instead.
     )pbdoc")
    .def(py::init<const std::string>(), py::arg("file_path"))
    .def("get_column_names", &AveFileReader::getColumnNames)
    .def("get_nr_of_rows", &AveFileReader::getNrOfRows)
    .def("get_nr_of_columns", &AveFileReader::getNrOfColumns)
    .def("get_data", &AveFileReader::getData)
    .def("autocorrelate_column",
         &AveFileReader::autocorrelateColumn,
         R"pbdoc(
          Do autocorrelation on one particular column for a specified set of delta indices.

          Assumes the data is equally spaced.
     )pbdoc",
         py::arg("column_index"),
         py::arg("delta_indices"))
    .def("autocorrelate_column_difference",
         &AveFileReader::autocorrelateColumnDifference,
         R"pbdoc(
          Do autocorrelation on the difference between two particular columns for a specified set of delta indices.

          Assumes the data is equally spaced.
     )pbdoc",
         py::arg("column_index1"),
         py::arg("column_index2"),
         py::arg("delta_indices"));

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
    .def(py::init<const std::string>(), py::arg("path_of_file_to_read"))
    .def("read", &DumpFileParser::read, "Read the whole file")
    .def("get_length",
         &DumpFileParser::getLength,
         "Get the number of sections (time-steps) in the file")
    .def("get_string_values_for_at",
         &DumpFileParser::getStringValuesForAt,
         "Get the values for the section `index`, the main header "
         "`headerKey` and the column (in the header) `column`.",
         py::arg("rowIndex"),
         py::arg("headerKey"),
         py::arg("columnIndex"))
    .def("get_numeric_values_for_at",
         &DumpFileParser::getNumericValuesForAt,
         "Get the values for the section `index`, the main header "
         "`headerKey` and the column (in the header) `column`.")
    .def("has_key",
         &DumpFileParser::hasKey,
         "Check whether the first section has the header specified",
         py::arg("headerKey"))
    .def("key_has_column",
         &DumpFileParser::keyHasColumn,
         "Check whether the header of the first section has the specified "
         "column",
         py::arg("headerKey"),
         py::arg("columnName"))
    .def("key_has_directional_column",
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
    .def("get_nr_of_atoms", &DataFileParser::getNrOfAtoms)
    .def("get_nr_of_atom_types", &DataFileParser::getNrOfAtomTypes)
    .def("get_atom_ids", &DataFileParser::getAtomIds)
    .def("get_molecule_ids", &DataFileParser::getMoleculeIds)
    .def("get_atom_types", &DataFileParser::getAtomTypes)
    .def("get_atom_x", &DataFileParser::getAtomX)
    .def("get_atom_y", &DataFileParser::getAtomY)
    .def("get_atom_z", &DataFileParser::getAtomZ)
    .def("get_atom_nx", &DataFileParser::getAtomNx)
    .def("get_atom_ny", &DataFileParser::getAtomNy)
    .def("get_atom_nz", &DataFileParser::getAtomNz)
    .def("get_masses", &DataFileParser::getMasses)
    .def("get_nr_of_bonds", &DataFileParser::getNrOfBonds)
    .def("get_nr_of_bond_types", &DataFileParser::getNrOfBondTypes)
    .def("get_bond_types", &DataFileParser::getBondTypes)
    .def("get_bond_from", &DataFileParser::getBondFrom)
    .def("get_bond_to", &DataFileParser::getBondTo)
    .def("get_lx", &DataFileParser::getLx)
    .def("get_low_x", &DataFileParser::getLowX)
    .def("get_high_x", &DataFileParser::getHighX)
    .def("get_ly", &DataFileParser::getLy)
    .def("get_low_y", &DataFileParser::getLowY)
    .def("get_high_y", &DataFileParser::getHighY)
    .def("get_lz", &DataFileParser::getLz)
    .def("get_low_z", &DataFileParser::getLowZ)
    .def("get_high_z", &DataFileParser::getHighZ);

    m.def("split_csv",
          &splitCSV,
          "Read a file containing a number of CSVs. Returns them split up.");
}

#endif
