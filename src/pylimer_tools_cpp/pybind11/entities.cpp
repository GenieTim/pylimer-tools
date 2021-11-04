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
    py::class_<Box>(m, "Box")
        .def(py::init<const double, const double, const double>())
        .def("getVolume", &Box::getVolume)
        .def("getLx", &Box::getLx)
        .def("getLy", &Box::getLy)
        .def("getLz", &Box::getLz);

    py::class_<Atom>(m, "Atom")
        .def(py::init<const int, const int, const double, const double, const double, const int, const int, const int>())
        .def("vectorTo", &Atom::vectorTo)
        .def("distanceTo", &Atom::distanceTo)
        .def("getX", &Atom::getX)
        .def("getY", &Atom::getY)
        .def("getZ", &Atom::getZ)
        .def("getNX", &Atom::getNX)
        .def("getNY", &Atom::getNY)
        .def("getNZ", &Atom::getNZ);

    py::enum_<MoleculeType>(m, "MoleculeType")
        .value("UNDEFINED", MoleculeType::UNDEFINED)
        .value("NETWORK_STRAND", MoleculeType::NETWORK_STRAND)
        .value("PRIMARY_LOOP", MoleculeType::PRIMARY_LOOP)
        .value("DANGLING_CHAIN", MoleculeType::DANGLING_CHAIN)
        .value("FREE_CHAIN", MoleculeType::FREE_CHAIN)
        .export_values();

    py::class_<Molecule>(m, "Molecule")
        .def(py::init<Box *, igraph_t *, MoleculeType>())
        .def("computeEndToEndDistance", &Molecule::computeEndToEndDistance)
        .def("getLength", &Molecule::getLength)
        .def("getType", &Molecule::getType)
        .def("getAtomForVertexId", &Molecule::getAtomForVertexId)
        .def("getNrOfBonds", &Molecule::getNrOfBonds)
        .def("getNrOfAtoms", &Molecule::getNrOfAtoms)
        .def("computeBondLengths", &Molecule::computeBondLengths);

    py::class_<Universe>(m, "Universe")
        .def(py::init<const double, const double, const double>())
        .def("addAtoms", &Universe::addAtoms)
        .def("addBonds", &Universe::addBonds)
        .def("setBoxLengths", &Universe::setBoxLengths)
        .def("getMolecules", &Universe::getMolecules)
        .def("getChainsWithCrosslinker", &Universe::getChainsWithCrosslinker)
        .def("determineFunctionalityPerType", &Universe::determineFunctionalityPerType)
        .def("getAtom", &Universe::getAtom)
        .def("getAtomsWithType", &Universe::getAtomsWithType)
        .def("getAtomByIdx", &Universe::getAtomByIdx)
        .def("getBox", &Universe::getBox)
        .def("getVolume", &Universe::getVolume)
        .def("getNrOfAtoms", &Universe::getNrOfAtoms);

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
