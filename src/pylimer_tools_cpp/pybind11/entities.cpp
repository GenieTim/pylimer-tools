#ifndef PYBIND_ENTITIES_H
#define PYBIND_ENTITIES_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Molecule.h"
#include "../entities/Universe.h"
#include "../entities/UniverseSequence.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pylimer_tools::entities;

struct MoleculeIterator {
  MoleculeIterator(const Molecule &molecule, py::object ref)
      : molecule(molecule), ref(ref) {}

  Atom next() {
    if (index == molecule.getLength()) {
      throw py::stop_iteration();
    }
    return molecule[index++];
  }

  const Molecule &molecule;
  py::object ref;   // keep a reference
  size_t index = 0; // the index to access
};

void init_pylimer_bound_entities(py::module_ &m) {
  py::class_<Box>(m, "Box", R"pbdoc(
        The box that the simulation is run in.

        NOTE: 
          currently, only rectangular boxes are supported.
        )pbdoc")
      .def(py::init<const double, const double, const double>())
      .def("getVolume", &Box::getVolume, R"pbdoc(
            Compute the volume of the box.

            :math:`V = L_x \cdot L_y \cdot L_z`
            )pbdoc")
      .def("getLx", &Box::getLx, R"pbdoc(
            Get the length of the box in x direction.
            )pbdoc")
      .def("getLowX", &Box::getLowX)
      .def("getHighX", &Box::getHighX)
      .def("getLy", &Box::getLy, R"pbdoc(
            Get the length of the box in y direction.
            )pbdoc")
      .def("getLowY", &Box::getLowY)
      .def("getHighY", &Box::getHighY)
      .def("getLz", &Box::getLz, R"pbdoc(
            Get the length of the box in z direction.
            )pbdoc")
      .def("getLowZ", &Box::getLowZ)
      .def("getHighZ", &Box::getHighZ)
      .def(py::pickle(
               [](const Box &b) { // __getstate__
                 /* Return a tuple that fully encodes the state of the object */
                 return py::make_tuple(b.getLowX(), b.getLowY(), b.getLowZ(),
                                       b.getHighX(), b.getHighY(),
                                       b.getHighZ());
               },
               [](py::tuple t) { // __setstate__
                 if (t.size() != 6) {
                   throw std::runtime_error("Invalid state!.");
                 }

                 /* Create a new C++ instance */
                 Box b = Box(t[0].cast<double>(), t[1].cast<double>(),
                             t[2].cast<double>(), t[3].cast<double>(),
                             t[4].cast<double>(), t[5].cast<double>());

                 return b;
               }),
           "Provides support for pickling.");

  py::class_<Atom>(m, "Atom", R"pbdoc(
       A single bead or atom
  )pbdoc")
      .def(py::init<const long int, const int, const double, const double,
                    const double, const int, const int, const int>(),
           "Construct this atom", py::arg("id"), py::arg("type"), py::arg("x"),
           py::arg("y"), py::arg("z"), py::arg("nx"), py::arg("ny"),
           py::arg("nz"))
      .def("computeVectorTo", &Atom::computeVectorTo, R"pbdoc(
            Compute the vector to another atom.
            )pbdoc")
      .def("distanceTo", &Atom::distanceTo, R"pbdoc(
            Compute the distance to another atom.
            )pbdoc")
      .def("vectorToUnwrapped", &Atom::vectorToUnwrapped,
           "Compute the vector to another atom respecting the periodic image "
           "flag.")
      .def("distanceToUnwrapped", &Atom::distanceToUnwrapped,
           "Compute the distance to another atom respecting the periodic image "
           "flag.")
      .def("getId", &Atom::getId, R"pbdoc(
            Get the id of the atom.
            )pbdoc")
      .def("getType", &Atom::getType, R"pbdoc(
            Get the type of the atom.
            )pbdoc")
      .def("getX", &Atom::getX, R"pbdoc(
            Get the x coordinate of the atom.
            )pbdoc")
      .def("getY", &Atom::getY, R"pbdoc(
            Get the y coordinate of the atom.
            )pbdoc")
      .def("getZ", &Atom::getZ, R"pbdoc(
            Get the z coordinate of the atom.
            )pbdoc")
      .def("getNX", &Atom::getNX,
           "Get the box image that the atom is in in x direction (also known "
           "as `ix` or `nx`).")
      .def("getNY", &Atom::getNY,
           "Get the box image that the atom is in in y direction (also known "
           "as `iy` or `ny`).")
      .def("getNZ", &Atom::getNZ,
           "Get the box image that the atom is in in z direction (also known "
           "as `iz` or `nz`).")
      .def(py::pickle(
               [](const Atom &b) { // __getstate__
                 /* Return a tuple that fully encodes the state of the object */
                 return py::make_tuple(b.getId(), b.getType(), b.getX(),
                                       b.getY(), b.getZ(), b.getNX(), b.getNY(),
                                       b.getNZ());
               },
               [](py::tuple t) { // __setstate__
                 if (t.size() != 8) {
                   throw std::runtime_error("Invalid state!.");
                 }

                 /* Create a new C++ instance */
                 Atom a = Atom(t[0].cast<long int>(), t[1].cast<int>(),
                               t[2].cast<double>(), t[3].cast<double>(),
                               t[4].cast<double>(), t[5].cast<int>(),
                               t[6].cast<int>(), t[7].cast<int>());

                 return a;
               }),
           "Provides support for pickling");

  py::class_<MoleculeIterator>(m, "MoleculeIterator", R"pbdoc(
       An iterator to iterate throught the atoms in :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule`.
  )pbdoc")
      .def("__iter__",
           [](MoleculeIterator &it) -> MoleculeIterator & { return it; })
      .def("__next__", &MoleculeIterator::next);

  py::enum_<MoleculeType>(m, "MoleculeType")
      .value("UNDEFINED", MoleculeType::UNDEFINED,
             "This value indicates that either the property was not set or not "
             "discovered.")
      .value("NETWORK_STRAND", MoleculeType::NETWORK_STRAND, R"pbdoc(
           A network strand is a strand in a network.
      )pbdoc")
      .value("PRIMARY_LOOP", MoleculeType::PRIMARY_LOOP, R"pbdoc(
           A primary loop is a network strand looping from and to the same cross-linker.
      )pbdoc")
      .value("DANGLING_CHAIN", MoleculeType::DANGLING_CHAIN, R"pbdoc(
           A dangling chain is a network strand where only one end is attached to a cross-linker.
      )pbdoc")
      .value("FREE_CHAIN", MoleculeType::FREE_CHAIN, R"pbdoc(
           A free chain is a strand not connected to any cross-linker.
      )pbdoc")
      .export_values();

  py::class_<Molecule>(m, "Molecule", R"pbdoc(
       An (ideally) connected series of atoms/beads.
  )pbdoc")
      .def(py::init<Box *, igraph_t *, MoleculeType, std::map<int, double>>())
      // getters
      .def("getLength", &Molecule::getLength, R"pbdoc(
           Counts and returns the number of atoms associated with this 
           molecule.
      )pbdoc")
      .def("getType", &Molecule::getType, R"pbdoc(
           Get the type of this molecule (see :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType` enum).
      )pbdoc")
      .def("getAtomsWithType", &Molecule::getAtomsWithType, R"pbdoc(
            Get the atoms with the specified type.
            )pbdoc")
      .def("getAtomsWithDegree", &Molecule::getAtomsOfDegree, R"pbdoc(
            Get the atoms that have the specified number of bonds.
            )pbdoc")
      .def("getAtomsConnectedTo", &Molecule::getAtomsConnectedTo, R"pbdoc(
            Get the atoms connected to a specified vertex id.
            )pbdoc")
      .def("getEdges", &Molecule::getEdges, R"pbdoc(
            Get all bonds. Returns a dict with three properties: 'edge_from', 'edge_to' and 'edge_type'.
            The order is not necessarily related to any structural property.
            
            NOTE:
               The integer values returned refer to the vertex ids, not the atom ids.
               Use :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule.getAtomIdByIdx` to translate them to atom ids, or 
               :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule.getBonds` to have that done for you.
            )pbdoc")
      .def("getBonds", &Molecule::getBonds, R"pbdoc(
            Get all bonds. Returns a dict with three properties: 'bond_from', 'bond_to' and 'bond_type'.
            )pbdoc")
      .def(
          "getAtomTypes", &Molecule::getAtomTypes,
          "Query all types (each one for each atom) ordered by atom vertex id.")
      .def("getAtomForVertexId", &Molecule::getAtomByVertexIdx, R"pbdoc(
            Get an atom for a specific vertex.
            )pbdoc")
      .def("getAtoms", &Molecule::getAtoms, R"pbdoc(
            Returns all atom objects enclosed in this molecule.
            )pbdoc")
      .def("getNrOfBonds", &Molecule::getNrOfBonds,
           "Counts and returns the number of bonds associated with this "
           "molecule.")
      .def("getNrOfAtoms", &Molecule::getNrOfAtoms,
           "Counts and returns the number of atoms associated with this "
           "molecule.")
      .def("getAtomIdByIdx", &Molecule::getAtomIdByIdx, R"pbdoc(
            Get the id of the atom by the vertex id of the underlying graph.
            )pbdoc")
      .def("getIdxByAtomId", &Molecule::getIdxByAtomId,
           "Get the vertex id of the underlying graph for an atom with a "
           "specified id.")
      .def("getKey", &Molecule::getKey, R"pbdoc(
            Get a unique identifier for this molecule.
            )pbdoc")
      // computations
      .def("computeWeight", &Molecule::computeWeight, R"pbdoc(
            Computes the total weight of this molecule.
            )pbdoc")
      .def("computeBondLengths", &Molecule::computeBondLengths,
           "Computes the length :math:`b` of each bond in the molecule, "
           "respecting periodic boundaries.")
      .def("computeRadiusOfGyration", &Molecule::computeRadiusOfGyration,
           R"pbdoc(
            Computes the radius of gyration, :math:`R_g^2` of this molecule.

            NOTE: 
               The mass of the atoms is not yet taken into account.
            )pbdoc")
      .def("computeEndToEndDistance", &Molecule::computeEndToEndDistance,
           R"pbdoc(
            Compute the end-to-end distance (:math:`R_{ee}`) of this molecule. 

            CAUTION:
               Returns 0.0 if the molecule does not have two or more atoms.
               Returns -1.0 if not exactly 2 ends were found.
            )pbdoc")
      // operators
      .def(
          "__getitem__",
          [](const Molecule &molecule, size_t index) {
            if (index > molecule.getLength()) {
              throw py::index_error();
            }
            return molecule[index];
          },
          R"pbdoc(
       Access an atom by its vertex index.
  )pbdoc")
      .def("__len__", &Molecule::getLength, R"pbdoc(
       Get the number of atoms in this molecule.
  )pbdoc")
      .def(
          "__iter__",
          [](py::object mol) {
            return MoleculeIterator(mol.cast<const Molecule &>(), mol);
          },
          R"pbdoc(
       Iterate through the atoms in this molecule.
       No specific order is guaranteed.
  )pbdoc")
      // .def(py::pickle(
      //      [](const Molecule &molecule) { // __getstate__
      //        /* Return a tuple that fully encodes the state of the object */
      //           return py::make_tuple(molecule.)
      //      }
      // ))
      ;

  py::class_<Universe>(
      m, "Universe",
      "Represents a full Polymer Network structure, a collection of molecules.")
      .def(py::init<const double, const double, const double>(),
           "Instantiate this Universe (Collection of Molecules) providing the "
           "box lengths.",
           py::arg("Lx"), py::arg("Ly"), py::arg("Lz"))
      // setters
      .def("addAtoms",
           py::overload_cast<std::vector<long int>, std::vector<int>,
                             std::vector<double>, std::vector<double>,
                             std::vector<double>, std::vector<int>,
                             std::vector<int>, std::vector<int>>(
               &Universe::addAtoms),
           "Add atoms to the Universe, vertices to the underlying graph.",
           py::arg("ids"), py::arg("types"), py::arg("x"), py::arg("y"),
           py::arg("z"), py::arg("nx"), py::arg("ny"), py::arg("nz"))
      .def("addBonds",
           py::overload_cast<std::vector<long int>, std::vector<long int>>(
               &Universe::addBonds),
           "Add bonds to the underlying atoms, edges to the underlying graph. "
           "If the connected atoms are not found, the bonds are silently "
           "skipped.",
           py::arg("from"), py::arg("to"))
      .def("addAngles", &Universe::addAngles,
           "Add angles to the Universe. No relation to the underlying graph, "
           "just a method to preserve read & write capabilities",
           py::arg("from"), py::arg("via"), py::arg("to"))
      .def("setMasses", &Universe::setMasses, "Set the mass per type of atom.",
           py::arg("massPerType"))
      .def("setTimestep", &Universe::setTimestep,
           "Set the timestep when this Universe was captured.",
           py::arg("timestep"))
      .def("setBoxLengths", &Universe::setBoxLengths,
           "Set the box side lengths.", py::arg("Lx"), py::arg("Ly"),
           py::arg("Lz"))
      .def("setBox", &Universe::setBox, R"pbdoc(
            Override the currently assigned box with the one specified.
            )pbdoc",
           py::arg("box"))
      // getters
      .def("getClusters", &Universe::getClusters, R"pbdoc(
            Get the components of the universe that are not connected to each other.
            Returns a list of :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule`s.
            Unconnected, free atoms/beads become their own :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule`.
            )pbdoc")
      .def("getMolecules", &Universe::getMolecules, R"pbdoc(
            Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms.
            
            Reduces the Universe to a list of molecules. 
            Specify the crosslinkerType to an existing type id, 
            then those atoms will be omitted, and this function returns chains instead.)pbdoc",
           py::arg("atomTypeToOmit"))
      .def("findLoops", &Universe::findLoops,
           R"pbdoc(
            Decompose the Universe into loops.
            The primary index specifies the degree of the loop.

            CAUTION:
               There are exponentially many paths between two cross-linkers of a network,
               and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
               You can use the maxLength parameter to restrict the algorithm to only search for loops up to a certain length.
               Use a negative value to find all loops and paths.
            )pbdoc",
           py::arg("crosslinkerType"), py::arg("maxLength") = -1,
           py::arg("skipSelfLoops") = false)
      .def("getChainsWithCrosslinker", &Universe::getChainsWithCrosslinker,
           R"pbdoc(
            Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms, without omitting the cross-linkers.
            In turn, e.g. for a tetrafunctional cross-linker, it will be 4 times in the resulting molecules.
            
            NOTE:
               Cross-linkers without bonds to non-cross-linkers are not returned 
               (i.e., cross-linker-cross-linker bonds, or single cross-linkers, are not counted as strands).
           )pbdoc",
           py::arg("crosslinkerType"))
      .def("getNetworkOfCrosslinker", &Universe::getNetworkOfCrosslinker,
           R"pbdoc(
            Reduce the network to contain only cross-linkers, replacing all the strands with a single bond.
            Useful e.g. to reduce the memory useage and runtime of 
            :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.findLoops()` or 
            :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.hasInfiniteStrand()`.
            
            Further use :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.simplify()` to remove primary loops.
          )pbdoc",
           py::arg("crossLinkerType"))
      .def(
          "getAtomTypes", &Universe::getAtomTypes,
          R"pbdoc(Get all types (each one for each atom) ordered by atom vertex id.)pbdoc")
      .def("getAtom", &Universe::getAtom, "Find an atom by its ID.",
           py::arg("atomId"))
      .def("getAtomForVertexId", &Molecule::getAtomByVertexIdx,
           R"pbdoc(Get an atom for a specific vertex.)pbdoc",
           py::arg("vertexId"))
      .def("getAtoms", &Universe::getAtoms, R"pbdoc(
            Get all atoms.
            )pbdoc")
      .def("getAtomsWithType", &Universe::getAtomsWithType, R"pbdoc(
            Find many atom by their type.
            )pbdoc")
      .def(
          "getAtomByVertexIdx", &Universe::getAtomByVertexIdx,
          R"pbdoc(Find an atom by the ID of the vertex of the underlying graph.)pbdoc",
          py::arg("vertexId"))
      .def("getAtomIdByIdx", &Universe::getAtomIdByIdx,
           "Get the id of the atom by the vertex id of the underlying graph.",
           py::arg("atomId"))
      .def("getIdxByAtomId", &Universe::getIdxByAtomId,
           "Get the vertex id of the underlying graph for an atom with a "
           "specified id.",
           py::arg("atomId"))
      .def("getEdges", &Universe::getEdges, R"pbdoc(
            Get all edges. Returns a dict with three properties: 'edge_from', 'edge_to' and 'edge_type'.
            The order is not necessarily related to any structural characteristic.
            
            NOTE:
               The integer values returned refer to the vertex ids, not the atom ids.
               Use :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.getAtomIdByIdx` to translate them to atom ids, or
               :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.getBonds` to have that done for you.
            )pbdoc")
      .def("getBonds", &Universe::getBonds, R"pbdoc(
            Get all bonds. Returns a dict with three properties: 'bond_from', 'bond_to' and 'bond_type'.
            The order is not necessarily related to any structural characteristic.
            )pbdoc")
      .def("getAngles", &Universe::getAngles, R"pbdoc(
           Get all angles added to this network.

           Returns a dict with three properties: 'angle_from', 'angle_via' and 'angle_to'.

            NOTE:
               The integer values returned refer to the the atom ids, not the vertex ids.
               Use :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.getIdxByAtomId` to translate them to vertex ids.
           )pbdoc")
      .def("getBox", &Universe::getBox, R"pbdoc(
            Get the underlying bounding box object.
            )pbdoc")
      .def("getMasses", &Universe::getMasses, R"pbdoc(
            Get the mass of one atom per type
            )pbdoc")
      .def("getVolume", &Universe::getVolume, R"pbdoc(
            Query the volume of the underlying bounding box.
            )pbdoc")
      .def("getNrOfAtoms", &Universe::getNrOfAtoms, R"pbdoc(
            Query the number of atoms in this universe.
            )pbdoc")
      .def("getNrOfBonds", &Universe::getNrOfBonds, R"pbdoc(
            Query the number of bonds associated with this universe.
            )pbdoc")
      .def("getNrOfAngles", &Universe::getNrOfAngles, R"pbdoc(
            Query the number of angles that have been added to this universe.
            )pbdoc")
      .def("getTimestep", &Universe::getTimestep, R"pbdoc(
            Query the timestep when this universe was captured.
            )pbdoc")
      .def(
          "getNrOfBondsOfAtom", &Universe::computeFunctionalityForAtom,
          R"pbdoc(Count the number of immediate neighbours of an atom, specified by its id.)pbdoc")
      .def(
          "getNrOfBondsOfVertex", &Universe::computeFunctionalityForVertex,
          R"pbdoc(Count the number of immediate neighbours of an atom, specified by its vertex id.)pbdoc")
      // computations
      .def("computeBondLengths", &Universe::computeBondLengths,
           "Computes the length :math:`b` of each bond in the molecule, "
           "respecting periodic boundaries.")
      .def("detectAngles", &Universe::detectAngles,
           "Returns just as "
           ":func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.getAngles`, "
           "but "
           "all angles that are found in the "
           "network.")
      .def("hasInfiniteStrand", &Universe::hasInfiniteStrand,
           R"pbdoc(
           Checks whether there is a strand (with cross-linker) in the universe that loops through periodic images without coming back.
           
            CAUTION:
               There are exponentially many paths between two cross-linkers of a network,
               and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
           )pbdoc")
      .def("determineFunctionalityPerType",
           &Universe::determineFunctionalityPerType, R"pbdoc(
            Find the maximum functionality of each atom type in the network.
            )pbdoc")
      .def("determineEffectiveFunctionalityPerType",
           &Universe::determineEffectiveFunctionalityPerType, R"pbdoc(
            Find the average functionality of each atom type in the network.
            )pbdoc")
      .def("computeMeanStrandLength", &Universe::getMeanStrandLength,
           "Compute the mean strand length.", py::arg("crosslinkerType"))
      .def("computeWeightFractions", &Universe::computeWeightFractions, R"pbdoc(
            Compute the weight fractions of each atom type in the network.
            )pbdoc")
      .def("computeDxs", &Universe::computeDxs,
           "Compute the dx distance for certain bonds (length in x direction).")
      .def("computeDys", &Universe::computeDys,
           "Compute the dy distance for certain bonds (length in y direction).")
      .def("computeDzs", &Universe::computeDzs,
           "Compute the dz distance for certain bonds (length in z direction).")
      .def("simplify", &Universe::simplify,
           "Remove self links and double bonds. This function is called "
           "automatically after adding bonds.");

  struct LazyUniverseSequenceIterator {
    LazyUniverseSequenceIterator(UniverseSequence &us, py::object ref)
        : us(us), ref(ref) {}

    Universe next() {
      if (index == us.getLength()) {
        throw py::stop_iteration();
      }
      Universe toReturn = us.atIndex(index);
      us.forgetAtIndex(index);
      index += 1;
      return toReturn;
    }

    UniverseSequence &us;
    py::object ref;   // keep a reference
    size_t index = 0; // the index to access
  };

  py::class_<UniverseSequence>(m, "UniverseSequence", R"pbdoc(
     This class represents a sequence of Universes, with the Universe's data
     only being read on request. Dump files are read at once in order
     to know how many timesteps/universes are available in total 
     (but the universes' data is not read on first look through the file).
     This, while it can lead to two (or more) reads of the whole file, 
     is a measure in order to enable low memory useage if needed (i.e. for large dump files).
     Use Python's iterator to have this UniverseSequence only ever retain one universe in memory.
     Alternatively, use :func:`~pylimer_tools_cpp.pylimer_tools_cpp.UniverseSequence.forgetAtIndex`
     to have the UniverseSequence forget about already read universes.
     )pbdoc")
      .def(py::init<>())
      .def("initializeFromDumpFile", &UniverseSequence::initializeFromDumpFile,
           R"pbdoc(
          Reset and initialize the Universes from a Lammps :code:`dump` output. 
        
          NOTE:
               If you have not output the id of the atoms in the dump file, they will be assigned sequentially. 
               If you have not output the type of the atoms in the dump file, they will be set to -1 if they cannot be infered from the data file.
        )pbdoc")
      .def("initializeFromDataSequence",
           &UniverseSequence::initializeFromDataSequence,
           "Reset and initialize the Universes from an ordered list of Lammps "
           "data (:code:`write_data`) files.")
      .def("next", &UniverseSequence::next,
           R"pbdoc(Get the Universe that's next in the sequence.)pbdoc")
      .def("atIndex", &UniverseSequence::atIndex,
           "Get the Universe at the given index (as of in the sequence given "
           "by the dump file).")
      .def(
          "forgetAtIndex", &UniverseSequence::forgetAtIndex,
          R"pbdoc(Clear the memory of the Universe at the given index (as of in the 
           sequence given by the dump file).)pbdoc")
      .def("resetIterator", &UniverseSequence::resetIterator,
           R"pbdoc(
          Reset the internal iterator, such that a subsequent call to 
          :func:`~pylimer_tools_cpp.pylimer_tools_cpp.UniverseSequence.next` returns the first one again.
          )pbdoc")
      .def("getLength", &UniverseSequence::getLength, R"pbdoc(
            Get the number of universes in this sequence.
            )pbdoc")
      .def("getAll", &UniverseSequence::getAll, R"pbdoc(
            Get all universes initialized back in a list.
            For big dump files or lots of data files, this might lead to memory issues.
            Use :func:`~pylimer_tools_cpp.pylimer_tools_cpp.UniverseSequence.__iter__`
            to have
            or :func:`~pylimer_tools_cpp.pylimer_tools_cpp.UniverseSequence.atIndex`
            and :func:`~pylimer_tools_cpp.pylimer_tools_cpp.UniverseSequence.forgetAtIndex`
            to craft a more memory-efficient retrieval mechanism.
            )pbdoc")
      // operators
      .def(
          "__getitem__",
          [](UniverseSequence &us, size_t index) {
            if (index > us.getLength()) {
              throw py::index_error();
            }
            return us.atIndex(index);
          },
          "Get a universe by its index.")
      .def("__len__", &UniverseSequence::getLength,
           "Get the number of universes")
      .def(
          "__iter__",
          [](py::object us) {
            return LazyUniverseSequenceIterator(us.cast<UniverseSequence &>(),
                                                us);
          },
          R"pbdoc(
           Lazily (memory-efficiently) iterate through all the universes in this sequence.
           This is the standard Python iteration way. 
           
           Example:

           .. code::
           
               for (universe in universeSequence):
                    # do something with the universe
                    pass
           

           Note: 
               this iterator is supposed to be memory-efficient. Therefore, no cache is kept;
               iterating twice will lead to the file(s) being read twice 
               (plus, for dump files, a third time initially to determine the number of universes in the file).
      )pbdoc");
}

#endif
