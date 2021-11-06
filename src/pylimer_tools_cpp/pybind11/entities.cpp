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
                    throw std::runtime_error("Invalid state!");
                }

                /* Create a new C++ instance */
                Box b = Box(t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>());

                return b;
            }));

    py::class_<Atom>(m, "Atom")
        .def(py::init<const int, const int, const double, const double, const double, const int, const int, const int>())
        .def("vectorTo", &Atom::vectorTo, "Compute the vector to another atom")
        .def("distanceTo", &Atom::distanceTo, "Compute the distance to another atom")
        .def("getId", &Atom::getId, "Get the id of the atom")
        .def("getType", &Atom::getType, "Get the type of the atom")
        .def("getX", &Atom::getX, "Get the x coordinate of the atom")
        .def("getY", &Atom::getY, "Get the y coordinate of the atom")
        .def("getZ", &Atom::getZ, "Get the z coordinate of the atom")
        .def("getNX", &Atom::getNX, "Get the box image that the atom is in in x direction (also known as `ix` or `nx`)")
        .def("getNY", &Atom::getNY, "Get the box image that the atom is in in y direction (also known as `iy` or `ny`)")
        .def("getNZ", &Atom::getNZ, "Get the box image that the atom is in in z direction (also known as `iz` or `nz`)")
        .def(py::pickle(
            [](const Atom &b) { // __getstate__
                                /* Return a tuple that fully encodes the state of the object */
                return py::make_tuple(b.getId(), b.getType(), b.getX(), b.getY(), b.getZ(), b.getNX(), b.getNY(), b.getNZ());
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 8)
                {
                    throw std::runtime_error("Invalid state!");
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
        .def("getAtomForVertexId", &Molecule::getAtomForVertexId, "Get an atom for a specific vertex.")
        .def("getAtoms", &Molecule::getAtoms, "Returns all atom objects enclosed in this molecule.")
        .def("getNrOfBonds", &Molecule::getNrOfBonds, "Counts and returns the number of bonds associated with this molecule.")
        .def("getNrOfAtoms", &Molecule::getNrOfAtoms, "Counts and returns the number of atoms associated with this molecule.")
        .def("computeBondLengths", &Molecule::computeBondLengths, "Computes the length of each bond in the molecule.")
        .def("computeRadiusOfGyration", &Molecule::computeRadiusOfGyration, "Computes the radius of gyration, :math:`R_g^2` of this molecule")
        .def("computeEndToEndDistance", &Molecule::computeEndToEndDistance, R"pbdoc(
        Compute the end-to-end distance of this molecule. 
        
        **Caution**:
        Returns 0.0 if the molecule does not have two or more atoms.
        Returns -1.0 if not exactly 2 ends were found.
        )pbdoc");

    py::class_<Universe>(m, "Universe")
        .def(py::init<const double, const double, const double>())
        .def("addAtoms", &Universe::addAtoms)
        .def("addBonds", &Universe::addBonds)
        .def("setBoxLengths", &Universe::setBoxLengths)
        .def("getMolecules", &Universe::getMolecules, R"pbdoc(
        Reduces the Universe to a list of molecules. 
        Specify the crosslinkerType to an existing type id, 
        then those atoms will be omitted, and this function returns chains instead.)pbdoc")
        .def("getChainsWithCrosslinker", &Universe::getChainsWithCrosslinker)
        .def("getAtom", &Universe::getAtom)
        .def("getAtomsWithType", &Universe::getAtomsWithType)
        .def("getAtomByVertexIdx", &Universe::getAtomByIdx)
        .def("getBox", &Universe::getBox)
        .def("getVolume", &Universe::getVolume)
        .def("getNrOfAtoms", &Universe::getNrOfAtoms)
        .def("determineFunctionalityPerType", &Universe::determineFunctionalityPerType);

    py::class_<UniverseSequence>(m, "UniverseSequence")
        .def(py::init<>())
        .def("initializeFromDumpFile", &UniverseSequence::initializeFromDumpFile)
        .def("initializeFromDataSequence", &UniverseSequence::initializeFromDataSequence)
        .def("next", &UniverseSequence::next)
        .def("atIndex", &UniverseSequence::atIndex)
        .def("resetIterator", &UniverseSequence::resetIterator)
        .def("getLength", &UniverseSequence::getLength);
}

#endif
