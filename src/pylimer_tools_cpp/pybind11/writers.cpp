#ifndef PYBIND_WRITERS_H
#define PYBIND_WRITERS_H

#include "../io/DataFileWriter.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pe = pylimer_tools::entities;
using namespace pylimer_tools::utils;

void
init_pylimer_bound_writers(py::module_& m)
{
  py::class_<DataFileWriter>(
    m,
    "DataFileWriter",
    py::module_local(),
    R"pbdoc(A class to write a LAMMPS data file from a universe.)pbdoc")
    .def(py::init<pe::Universe>(), py::arg("universe"), R"pbdoc(
        Initialize the writer with the universe to write.
    )pbdoc")
    .def("set_universe_to_write",
         &DataFileWriter::setUniverseToWrite,
         py::arg("universe"),
         R"pbdoc(
        Re-set the universe to write.
    )pbdoc")
    .def("config_include_angles",
         &DataFileWriter::configIncludeAngles,
         py::arg("include_angles") = true,
         R"pbdoc(
        Set whether to include the angles from the universe in the file or not.

        Default: true.
    )pbdoc")
    .def("config_include_dihedral_angles",
         &DataFileWriter::configIncludeDihedralAngles,
         py::arg("include_dihedral_angles") = true,
         R"pbdoc(
        Set whether to include the dihedral angles from the universe in the file or not.

        Default: true.
    )pbdoc")
    .def("config_include_velocities",
         &DataFileWriter::configIncludeVelocities,
         py::arg("include_velocities") = true,
         R"pbdoc(
        Set whether to include the velocities from the universe (if any) in the file or not.

        Default: true.
    )pbdoc")
    .def("config_reindex_atoms",
         &DataFileWriter::configReindexAtoms,
         py::arg("reindex_atoms") = false,
         R"pbdoc(
        Set whether to reindex the atoms or not.
        Re-indexing leads to atom ids being in the range of 1 to the number of atoms.

        Default: false.
    )pbdoc")
    .def("config_move_into_box",
         &DataFileWriter::configMoveIntoBox,
         py::arg("move_into_box") = false,
         R"pbdoc(
        Set whether to change the output coordinates to lie in the box or not.

        Default: false (used to be true).
    )pbdoc")
    .def("config_attempt_image_reset",
         &DataFileWriter::configAttemptImageReset,
         py::arg("attempt_image_reset") = false,
         R"pbdoc(
        Set whether to attempt to reset image flags so that output coordinates lie in the box.

        Default: false.
    )pbdoc")
    .def("config_atom_style",
         &DataFileWriter::configAtomStyle,
         py::arg("atom_style") = pylimer_tools::utils::AtomStyle::ANGLE,
         R"pbdoc(
        Set the (LAMMPS) atom style to use for writing the atoms.

        Default: AtomStyle.ANGLE.
    )pbdoc")
    .def("set_custom_atom_format",
         &DataFileWriter::setCustomAtomFormat,
         py::arg("atom_format") =
           "\t$atomId\t$moleculeId\t$atomType\t$x\t$y\t$z\t$nx\t$ny\t$nz",
         R"pbdoc(
        Specify a custom format for the atom section.

        Placeholder options are:

          - $atomId
          - $moleculeId
          - $atomType
          - $x
          - $y
          - $z
          - $nx
          - $ny
          - $nz

        Additionally, you can use the keys used in
        :func:`~pylimer_tools_cpp.Universe.setPropertyValue`
        as placeholders (as long as they are alphanumeric only; prefix in the format with '$' as well).
        Specifically useful if you need a different (or hybrid) atom style in LAMMPS.

        Be sure to still call :func:`~pylimer_tools_cpp.DataFileWriter.configAtomStyle`,
        so that the file can be read correctly again.

        Default: tab-separated LAMMPS angle atom style.
    )pbdoc")
    .def("config_crosslinker_type",
         &DataFileWriter::configCrosslinkerType,
         py::arg("crosslinker_type") = 2,
         R"pbdoc(
        Set which atom type represents crosslinkers.
        Needed in case the moleculeIdx in the output file should have any meaning.
        (e.g. with :func:`~pylimer_tools_cpp.DataFileWriter.configMoleculeIdxForSwap`).

        Default: 2.
    )pbdoc")
    .def("config_molecule_idx_for_swap",
         &DataFileWriter::configMoleculeIdxForSwap,
         py::arg("enable_swappability") = false,
         R"pbdoc(
        Swappable chains implies that their `moleculeIdx` in the LAMMPS data file is not
        identical per chain, but identical per position in the chain.
        That's how you can have bond swapping with constant chain length distribution.

        Default: false.
    )pbdoc")
    .def("write_to_file",
         &DataFileWriter::writeToFile,
         py::arg("file"),
         R"pbdoc(
        Actually do the writing to the disk.

        Arguments:
            file (str): The path and file name to write to.
    )pbdoc");
}

#endif
