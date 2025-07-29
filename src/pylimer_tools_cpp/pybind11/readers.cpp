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

  py::enum_<AtomStyle>(
    m, "AtomStyle", "An enumeration of the LAMMPS atom styles.")
#define X(e, n) .value(#e, AtomStyle::e, "LAMMPS atom style '" n "'")
    LAMMPS_ATOM_STYLES
#undef X
    ;

  py::class_<AveFileReader>(m, "AveFileReader", R"pbdoc(
          Alternative implementation of the data file reader implemented in 
          :func:`pylimer_tools.read_lammps_output_file.read_averages_file`.

          This implementation is better for certain use cases, worse for others.
          In the end, only performance and memory usage are different.
          For moderately sized and small files, we recommend to use the Python interface instead.
     )pbdoc")
    .def(py::init<const std::string>(),
         R"pbdoc(
         Initialize the AveFileReader with a file path.
         
         :param file_path: Path to the averages file to read
         )pbdoc",
         py::arg("file_path"))
    .def("get_column_names", &AveFileReader::getColumnNames, R"pbdoc(
         Get the names of all columns in the file.
         
         :return: List of column names
         )pbdoc")
    .def("get_nr_of_rows", &AveFileReader::getNrOfRows, R"pbdoc(
         Get the number of data rows in the file.
         
         :return: Number of rows
         )pbdoc")
    .def("get_nr_of_columns", &AveFileReader::getNrOfColumns, R"pbdoc(
         Get the number of columns in the file.
         
         :return: Number of columns
         )pbdoc")
    .def("get_data", &AveFileReader::getData, R"pbdoc(
         Get all data from the file.
         
         :return: 2D array containing all the numerical data
         )pbdoc")
    .def("autocorrelate_column",
         &AveFileReader::autocorrelateColumn,
         R"pbdoc(
          Do autocorrelation on one particular column for a specified set of delta indices.

          Assumes the data is equally spaced.
          
          :param column_index: Index of the column to autocorrelate
          :param delta_indices: List of delta indices for the autocorrelation
          :return: Autocorrelation values for the specified deltas
     )pbdoc",
         py::arg("column_index"),
         py::arg("delta_indices"))
    .def("autocorrelate_column_difference",
         &AveFileReader::autocorrelateColumnDifference,
         R"pbdoc(
          Do autocorrelation on the difference between two particular columns for a specified set of delta indices.

          Assumes the data is equally spaced.
          
          :param column_index1: Index of the first column
          :param column_index2: Index of the second column  
          :param delta_indices: List of delta indices for the autocorrelation
          :return: Autocorrelation values for the column differences at specified deltas
     )pbdoc",
         py::arg("column_index1"),
         py::arg("column_index2"),
         py::arg("delta_indices"));

  py::class_<DumpFileParser>(m, "DumpFileReader", R"pbdoc(
       A reader for LAMMPS's `dump` files.
  )pbdoc")
    .def(py::init<const std::string>(),
         R"pbdoc(
         Initialize the dump file reader.
         
         :param path_of_file_to_read: Path to the dump file to read
         )pbdoc",
         py::arg("path_of_file_to_read"))
    .def("read", &DumpFileParser::read, R"pbdoc(
         Read the whole file.
                  )pbdoc")
    .def("get_length",
         &DumpFileParser::getLength,
         R"pbdoc(
         Get the number of sections (time-steps) in the file.
         
         :return: Number of timesteps/sections in the dump file
         )pbdoc")
    .def("get_string_values_for_at",
         &DumpFileParser::getStringValuesForAt,
         R"pbdoc(
         Get string values for the section at index, the main header headerKey and the column columnIndex.
         
         :param rowIndex: Index of the section/timestep to query
         :param headerKey: Name of the header section to query
         :param columnIndex: Index of the column within the header
         :return: Vector of string values for the specified location
         )pbdoc",
         py::arg("rowIndex"),
         py::arg("headerKey"),
         py::arg("columnIndex"))
    .def("get_numeric_values_for_at",
         &DumpFileParser::getNumericValuesForAt,
         R"pbdoc(
         Get numeric values for the section at index, the main header headerKey and the column columnIndex.
         
         :param rowIndex: Index of the section/timestep to query
         :param headerKey: Name of the header section to query  
         :param columnIndex: Index of the column within the header
         :return: Vector of numeric values for the specified location
         )pbdoc",
         py::arg("rowIndex"),
         py::arg("headerKey"),
         py::arg("columnIndex"))
    .def("has_key",
         &DumpFileParser::hasKey,
         R"pbdoc(
         Check whether the first section has the header specified.
         
         :param headerKey: Name of the header to check for
         :return: True if the header exists in the first section
         )pbdoc",
         py::arg("headerKey"))
    .def("key_has_column",
         &DumpFileParser::keyHasColumn,
         R"pbdoc(
         Check whether the header of the first section has the specified column.
         
         :param headerKey: Name of the header section to check
         :param columnName: Name of the column to look for
         :return: True if the column exists in the specified header
         )pbdoc",
         py::arg("headerKey"),
         py::arg("columnName"))
    .def("key_has_directional_column",
         &DumpFileParser::keyHasDirectionalColumn,
         R"pbdoc(
         Check whether the header of the first section has all the three columns {dir_prefix}{x|y|z}{dir_suffix}.
         
         :param header_key: Name of the header section to check
         :param dir_prefix: Prefix for the directional columns (default: "")
         :param dir_suffix: Suffix for the directional columns (default: "")
         :return: True if all three directional columns (x, y, z) exist
         )pbdoc",
         py::arg("header_key"),
         py::arg("dir_prefix") = "",
         py::arg("dir_suffix") = "");

  py::class_<DataFileParser>(m, "DataFileReader", R"pbdoc(
       A reader for LAMMPS's `write_data` files.
  )pbdoc")
    .def(py::init<>(), R"pbdoc(
         Initialize a new data file parser.
         )pbdoc")
    .def("read",
         &DataFileParser::read,
         R"pbdoc(
       Actually read a LAMMPS's `write_data` file.

       :param path_of_file_to_read: The path to the file to read
       :param atom_style: The format of the "Atoms" section, see https://docs.lammps.org/read_data.html
       :param atom_style2: The format of the "Atoms" section if the previous parameter is equal to AtomStyle::HYBRID
       :param atom_style_3: The format of the "Atoms" section if the second to last parameter is equal to AtomStyle::HYBRID
  )pbdoc",
         py::arg("path_of_file_to_read"),
         py::arg("atom_style") = AtomStyle::ANGLE,
         py::arg("atom_style2") = AtomStyle::NONE,
         py::arg("atom_style_3") = AtomStyle::NONE)
    .def("get_nr_of_atoms", &DataFileParser::getNrOfAtoms, R"pbdoc(
         Get the number of atoms in the data file.
         
         :return: Number of atoms
         )pbdoc")
    .def("get_nr_of_atom_types", &DataFileParser::getNrOfAtomTypes, R"pbdoc(
         Get the number of atom types in the data file.
         
         :return: Number of atom types
         )pbdoc")
    .def("get_atom_ids", &DataFileParser::getAtomIds, R"pbdoc(
         Get all atom IDs from the data file.
         
         :return: Vector of atom IDs
         )pbdoc")
    .def("get_molecule_ids", &DataFileParser::getMoleculeIds, R"pbdoc(
         Get all molecule IDs from the data file.
         
         :return: Vector of molecule IDs
         )pbdoc")
    .def("get_atom_types", &DataFileParser::getAtomTypes, R"pbdoc(
         Get all atom types from the data file.
         
         :return: Vector of atom types
         )pbdoc")
    .def("get_atom_x", &DataFileParser::getAtomX, R"pbdoc(
         Get x coordinates of all atoms.
         
         :return: Vector of x coordinates
         )pbdoc")
    .def("get_atom_y", &DataFileParser::getAtomY, R"pbdoc(
         Get y coordinates of all atoms.
         
         :return: Vector of y coordinates
         )pbdoc")
    .def("get_atom_z", &DataFileParser::getAtomZ, R"pbdoc(
         Get z coordinates of all atoms.
         
         :return: Vector of z coordinates
         )pbdoc")
    .def("get_atom_nx", &DataFileParser::getAtomNx, R"pbdoc(
         Get periodic image flags in x direction for all atoms.
         
         :return: Vector of nx values (image flags)
         )pbdoc")
    .def("get_atom_ny", &DataFileParser::getAtomNy, R"pbdoc(
         Get periodic image flags in y direction for all atoms.
         
         :return: Vector of ny values (image flags)
         )pbdoc")
    .def("get_atom_nz", &DataFileParser::getAtomNz, R"pbdoc(
         Get periodic image flags in z direction for all atoms.
         
         :return: Vector of nz values (image flags)
         )pbdoc")
    .def("get_masses", &DataFileParser::getMasses, R"pbdoc(
         Get the mass values for each atom type.
         
         :return: Map of atom types to their masses
         )pbdoc")
    .def("get_nr_of_bonds", &DataFileParser::getNrOfBonds, R"pbdoc(
         Get the number of bonds in the data file.
         
         :return: Number of bonds
         )pbdoc")
    .def("get_nr_of_bond_types", &DataFileParser::getNrOfBondTypes, R"pbdoc(
         Get the number of bond types in the data file.
         
         :return: Number of bond types
         )pbdoc")
    .def("get_bond_types", &DataFileParser::getBondTypes, R"pbdoc(
         Get all bond types from the data file.
         
         :return: Vector of bond types
         )pbdoc")
    .def("get_bond_from", &DataFileParser::getBondFrom, R"pbdoc(
         Get starting atom IDs for all bonds.
         
         :return: Vector of starting atom IDs for bonds
         )pbdoc")
    .def("get_bond_to", &DataFileParser::getBondTo, R"pbdoc(
         Get ending atom IDs for all bonds.
         
         :return: Vector of ending atom IDs for bonds
         )pbdoc")
    .def("get_lx", &DataFileParser::getLx, R"pbdoc(
         Get the box length in x direction.
         
         :return: Box length in x direction
         )pbdoc")
    .def("get_low_x", &DataFileParser::getLowX, R"pbdoc(
         Get the lower bound of the box in x direction.
         
         :return: Lower x boundary
         )pbdoc")
    .def("get_high_x", &DataFileParser::getHighX, R"pbdoc(
         Get the upper bound of the box in x direction.
         
         :return: Upper x boundary
         )pbdoc")
    .def("get_ly", &DataFileParser::getLy, R"pbdoc(
         Get the box length in y direction.
         
         :return: Box length in y direction
         )pbdoc")
    .def("get_low_y", &DataFileParser::getLowY, R"pbdoc(
         Get the lower bound of the box in y direction.
         
         :return: Lower y boundary
         )pbdoc")
    .def("get_high_y", &DataFileParser::getHighY, R"pbdoc(
         Get the upper bound of the box in y direction.
         
         :return: Upper y boundary
         )pbdoc")
    .def("get_lz", &DataFileParser::getLz, R"pbdoc(
         Get the box length in z direction.
         
         :return: Box length in z direction
         )pbdoc")
    .def("get_low_z", &DataFileParser::getLowZ, R"pbdoc(
         Get the lower bound of the box in z direction.
         
         :return: Lower z boundary
         )pbdoc")
    .def("get_high_z", &DataFileParser::getHighZ, R"pbdoc(
         Get the upper bound of the box in z direction.
         
         :return: Upper z boundary
         )pbdoc");

  m.def("split_csv",
        &splitCSV,
        R"pbdoc(
        Read a file containing a number of CSVs. Returns them split up.
        
        :param file_path: Path to the file containing CSV data
        :param delimiter: Delimiter used in the CSV file (default: ',')
        :return: Split CSV data structures
        )pbdoc",
        py::arg("file_path"),
        py::arg("delimiter") = ',');
}

#endif
