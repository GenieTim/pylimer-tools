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
  py::class_<DataFileWriter>(m,
                             "DataFileWriter",
                             py::module_local(),
                             R"pbdoc(
    A class to write a LAMMPS data file from a universe.
    
    .. attention::
        The resulting file is not guaranteed to be a completely valid LAMMPS file.
        In particular, this writer does not force you to set masses for all atom types,
        and it does not have limits on the image flags.
        You can use :meth:`~pylimer_tools_cpp.Universe.set_masses` to set the masses for all atom types,
        and either :meth:`~pylimer_tools_cpp.Universe.set_vertex_property` or 
        :meth:`~pylimer_tools_cpp.DataFileWriter.config_attempt_image_reset` and
        :meth:`~pylimer_tools_cpp.DataFileWriter.config_move_into_box` to ensure that the image flags are correct.
    )pbdoc")
    .def(py::init<pe::Universe>(),
         R"pbdoc(
        Initialize the writer with the universe to write.
        
        :param universe: The universe to write to the data file
        )pbdoc",
         py::arg("universe"))
    .def("set_universe_to_write",
         &DataFileWriter::setUniverseToWrite,
         R"pbdoc(
        Re-set the universe to write.
        
        :param universe: The new universe to write to the data file
        )pbdoc",
         py::arg("universe"))
    .def("config_include_angles",
         &DataFileWriter::configIncludeAngles,
         R"pbdoc(
        Set whether to include the angles from the universe in the file or not.

        :param include_angles: Whether to include angles (default: True)
        )pbdoc",
         py::arg("include_angles") = true)
    .def("config_include_dihedral_angles",
         &DataFileWriter::configIncludeDihedralAngles,
         R"pbdoc(
        Set whether to include the dihedral angles from the universe in the file or not.

        :param include_dihedral_angles: Whether to include dihedral angles (default: True)
        )pbdoc",
         py::arg("include_dihedral_angles") = true)
    .def("config_include_velocities",
         &DataFileWriter::configIncludeVelocities,
         R"pbdoc(
        Set whether to include the velocities from the universe (if any) in the file or not.

        :param include_velocities: Whether to include velocities (default: True)
        )pbdoc",
         py::arg("include_velocities") = true)
    .def("config_reindex_atoms",
         &DataFileWriter::configReindexAtoms,
         R"pbdoc(
        Set whether to reindex the atoms or not.
        Re-indexing leads to atom IDs being in the range of 1 to the number of atoms.

        :param reindex_atoms: Whether to reindex atoms (default: False)
        )pbdoc",
         py::arg("reindex_atoms") = false)
    .def("config_move_into_box",
         &DataFileWriter::configMoveIntoBox,
         R"pbdoc(
        Set whether to change the output coordinates to lie in the box or not.

        :param move_into_box: Whether to move coordinates into box (default: False)
        )pbdoc",
         py::arg("move_into_box") = false)
    .def("config_attempt_image_reset",
         &DataFileWriter::configAttemptImageReset,
         R"pbdoc(
        Set whether to attempt to reset image flags so that output coordinates lie in the box.

        :param attempt_image_reset: Whether to attempt image flag reset (default: False)
        )pbdoc",
         py::arg("attempt_image_reset") = false)
    .def("config_atom_style",
         &DataFileWriter::configAtomStyle,
         R"pbdoc(
        Set the (LAMMPS) atom style to use for writing the atoms.

        :param atom_style: The LAMMPS atom style to use (default: AtomStyle.ANGLE)
        )pbdoc",
         py::arg("atom_style") = pylimer_tools::utils::AtomStyle::ANGLE)
    .def("set_custom_atom_format",
         &DataFileWriter::setCustomAtomFormat,
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
        :meth:`~pylimer_tools_cpp.Universe.set_property_value`
        as placeholders (as long as they are alphanumeric only; prefix in the format with '$' as well).
        Other placeholders are available if the universe was read from a LAMMPS data file with an
        atom style with additional data.

        This method is specifically useful if you need a different (or hybrid) atom style in LAMMPS.

        Be sure to still call :meth:`~pylimer_tools_cpp.DataFileWriter.config_atom_style`,
        so that the file can be read correctly again.

        :param atom_format: Custom format string for atoms (default: tab-separated standard format)
        )pbdoc",
         py::arg("atom_format") =
           "\t$atomId\t$moleculeId\t$atomType\t$x\t$y\t$z\t$nx\t$ny\t$nz")
    .def("config_crosslinker_type",
         &DataFileWriter::configCrosslinkerType,
         R"pbdoc(
        Set which atom type represents crosslinkers.
        Needed in case the moleculeIdx in the output file should have any meaning.
        (e.g. with :meth:`~pylimer_tools_cpp.DataFileWriter.config_molecule_idx_for_swap`).

        :param crosslinker_type: The atom type representing crosslinkers (default: 2)
        )pbdoc",
         py::arg("crosslinker_type") = 2)
    .def("config_molecule_idx_for_swap",
         &DataFileWriter::configMoleculeIdxForSwap,
         R"pbdoc(
        Swappable chains implies that their `moleculeIdx` in the LAMMPS data file is not
        identical per chain, but identical per position in the chain.
        That's how you can have bond swapping with constant chain length distribution.

        :param enable_swappability: Whether to enable molecule index swappability (default: False)
        )pbdoc",
         py::arg("enable_swappability") = false)
    .def("write_to_file",
         &DataFileWriter::writeToFile,
         R"pbdoc(
        Actually do the writing to the disk.

        :param file: The path and file name to write to
        )pbdoc",
         py::arg("file"));
}

#endif
