#ifndef PYBIND_ENTITIES_H
#define PYBIND_ENTITIES_H

#include "../entities/Box.h"
#include "../entities/Atom.h"
#include "../entities/Molecule.h"
#include "../entities/Universe.h"
#include "../entities/UniverseSequence.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pylimer_tools::entities;

void init_pylimer_bound_entities(py::module_ &m)
{
    py::class_<Box>(m, "Box", R"pbdoc(
        The box that the simulation is run in.

        *NOTE*: currently, only rectangular boxes are supported.
        )pbdoc")
        .def(py::init<const double, const double, const double>())
        .def("getVolume", &Box::getVolume, R"pbdoc(
            Compute the volume of the box.

            :math:`V = L_x \cdot L_y \cdot L_z`
            )pbdoc")
        .def("getLx", &Box::getLx, "Get the lenght of the box in x direction.")
        .def("getLy", &Box::getLy, "Get the lenght of the box in y direction.")
        .def("getLz", &Box::getLz, "Get the lenght of the box in z direction.")
        .def(py::pickle(
            [](const Box &b) { // __getstate__
                               /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(b.getLx(), b.getLy(), b.getLz());
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 3)
                {
                    throw std::runtime_error("Invalid state!.");
                }

                /* Create a new C++ instance */
                Box b = Box(t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>());

                return b;
            }));

    py::class_<Atom>(m, "Atom")
        .def(py::init<const long int, const int, const double, const double, const double, const int, const int, const int>())
        .def("computeVectorTo", &Atom::computeVectorTo, "Compute the vector to another atom.")
        .def("distanceTo", &Atom::distanceTo, "Compute the distance to another atom.")
        .def("vectorToUnwrapped", &Atom::vectorToUnwrapped, "Compute the vector to another atom respecting the periodic image flag.")
        .def("distanceToUnwrapped", &Atom::distanceToUnwrapped, "Compute the distance to another atom respecting the periodic image flag.")
        .def("getId", &Atom::getId, "Get the id of the atom.")
        .def("getType", &Atom::getType, "Get the type of the atom.")
        .def("getX", &Atom::getX, "Get the x coordinate of the atom.")
        .def("getY", &Atom::getY, "Get the y coordinate of the atom.")
        .def("getZ", &Atom::getZ, "Get the z coordinate of the atom.")
        .def("getNX", &Atom::getNX, "Get the box image that the atom is in in x direction (also known as `ix` or `nx`).")
        .def("getNY", &Atom::getNY, "Get the box image that the atom is in in y direction (also known as `iy` or `ny`).")
        .def("getNZ", &Atom::getNZ, "Get the box image that the atom is in in z direction (also known as `iz` or `nz`).")
        .def(py::pickle(
            [](const Atom &b) { // __getstate__
                                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(b.getId(), b.getType(), b.getX(), b.getY(), b.getZ(), b.getNX(), b.getNY(), b.getNZ());
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 8)
                {
                    throw std::runtime_error("Invalid state!.");
                }

                /* Create a new C++ instance */
                Atom a = Atom(t[0].cast<long int>(), t[1].cast<int>(),
                              t[2].cast<double>(), t[3].cast<double>(), t[4].cast<double>(),
                              t[5].cast<int>(), t[6].cast<int>(), t[7].cast<int>());

                return a;
            }));

    py::enum_<MoleculeType>(m, "MoleculeType")
        .value("UNDEFINED", MoleculeType::UNDEFINED, "This value indicates that either the property was not set or not discovered.")
        .value("NETWORK_STRAND", MoleculeType::NETWORK_STRAND)
        .value("PRIMARY_LOOP", MoleculeType::PRIMARY_LOOP)
        .value("DANGLING_CHAIN", MoleculeType::DANGLING_CHAIN)
        .value("FREE_CHAIN", MoleculeType::FREE_CHAIN)
        .export_values();

    py::class_<Molecule>(m, "Molecule")
        .def(py::init<Box *, igraph_t *, MoleculeType>())
        .def("getLength", &Molecule::getLength, "Counts and returns the number of atoms associated with this molecule.")
        .def("getType", &Molecule::getType, "Get the type of this molecule (see 'MoleculeType' enum).")
        .def("getAtomsWithType", &Molecule::getAtomsWithType, "Get the atoms with the specified type.")
        .def("getAtomTypes", &Molecule::getAtomTypes, "Query all types (each one for each atom) ordered by atom vertex id.")
        .def("getAtomForVertexId", &Molecule::getAtomForVertexId, "Get an atom for a specific vertex.")
        .def("getAtoms", &Molecule::getAtoms, "Returns all atom objects enclosed in this molecule.")
        .def("getNrOfBonds", &Molecule::getNrOfBonds, "Counts and returns the number of bonds associated with this molecule.")
        .def("getNrOfAtoms", &Molecule::getNrOfAtoms, "Counts and returns the number of atoms associated with this molecule.")
        .def("computeBondLengths", &Molecule::computeBondLengths, "Computes the length :math:`b` of each bond in the molecule.")
        .def("computeRadiusOfGyration", &Molecule::computeRadiusOfGyration, R"pbdoc(
            Computes the radius of gyration, :math:`R_g^2` of this molecule.

            **NOTE**: the mass of the atoms is not yet taken into account.
            )pbdoc")
        .def("computeEndToEndDistance", &Molecule::computeEndToEndDistance, R"pbdoc(
            Compute the end-to-end distance (:math:`R_{ee}`) of this molecule. 

            **Caution**:
            Returns 0.0 if the molecule does not have two or more atoms.
            Returns -1.0 if not exactly 2 ends were found.
            )pbdoc");

    py::class_<Universe>(m, "Universe", "Represents a full Polymer Network structure, a collection of molecules.")
        .def(py::init<const double, const double, const double>(), "Instantiate this Universe (Collection of Molecules) providing the box lengths.")
        .def("addAtoms", &Universe::addAtoms, "Add atoms to the Universe, vertices to the underlying graph.")
        .def("addBonds", &Universe::addBonds, "Add bonds to the underlying atoms, edges to the underlying graph. If the connected atoms are not found, the bonds are silently skipped.")
        .def("setMasses", &Universe::setMasses, "Set the mass per type of atom.")
        .def("setTimestep", &Universe::setTimestep, "Set the timestep when this Universe was captured.")
        .def("setBoxLengths", &Universe::setBoxLengths, "Set the box side lengths.")
        .def("setBox", &Universe::setBox, "Override the currently assigned box with the one specified.")
        .def("getMolecules", &Universe::getMolecules, R"pbdoc(
            Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms.
            
            Reduces the Universe to a list of molecules. 
            Specify the crosslinkerType to an existing type id, 
            then those atoms will be omitted, and this function returns chains instead.)pbdoc")
        .def("getChainsWithCrosslinker", &Universe::getChainsWithCrosslinker, R"pbdoc(
            Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms, without omitting the crosslinkers.
            In turn, e.g. for a tetrafunctional crosslinker, it will be 4 times in the resulting molecules.
            
            **NOTE**: Crosslinkers without bonds to non-crosslinkers are not returned.)pbdoc")
        .def("getAtomTypes", &Universe::getAtomTypes, "Get all types (each one for each atom) ordered by atom vertex id.")
        .def("getAtom", &Universe::getAtom, "Find an atom by its ID.")
        .def("getAtomsWithType", &Universe::getAtomsWithType, "Find many atom by their type.")
        .def("getAtomByVertexIdx", &Universe::getAtomByIdx, "Find an atom by the ID of the vertex of the underlying graph.")
        .def("getBox", &Universe::getBox, "Get the underlying bounding box object.")
        .def("getMasses", &Universe::getMasses, "Get the mass of one atom per type")
        .def("getVolume", &Universe::getVolume, "Query the volume of the underlying bounding box.")
        .def("getNrOfAtoms", &Universe::getNrOfAtoms, "Query the number of atoms in this universe.")
        .def("getNrOfBonds", &Universe::getNrOfBonds, "Query the number of bonds associated with this universe.")
        .def("getTimestep", &Universe::getTimestep, "Query the timestep when this universe was captured.")
        .def("getNrOfBondsOfAtom", &Universe::getNrOfBondsOfAtom, "Count the number of immediate neighbours of an atom, specified by its id.")
        .def("getNrOfBondsOfVertex", &Universe::getNrOfBondsOfVertex, "Count the number of immediate neighbours of an atom, specified by its vertex id.")
        .def("determineFunctionalityPerType", &Universe::determineFunctionalityPerType, "Find the maximum functionality of each atom type in the network.");

    py::class_<UniverseSequence>(m, "UniverseSequence", "This class represents a sequence of Universes, with the Universe's data files only being read on request. Dump files are read at once in order to know how many timesteps/universes are available in total.")
        .def(py::init<>())
        .def("initializeFromDumpFile", &UniverseSequence::initializeFromDumpFile, "Reset and initialize the Universes from a Lammps :code:`dump` output.")
        .def("initializeFromDataSequence", &UniverseSequence::initializeFromDataSequence, "Reset and initialize the Universes from an ordered list of Lammps data (:code:`write_data`) files.")
        .def("next", &UniverseSequence::next, "Get the Universe that's next in the sequence.")
        .def("atIndex", &UniverseSequence::atIndex, "Get the Universe at the given index (as of in the sequence given by the dump file).")
        .def("resetIterator", &UniverseSequence::resetIterator, "Reset the internal iterator, such that a subsequent call to :code:`next()` returns the first one again.")
        .def("getLength", &UniverseSequence::getLength, "Get the number of universes in this sequence.")
        .def("getAll", &UniverseSequence::getAll, "Get all universes initialized back in a list");
}

#endif
