#ifndef PYBIND_ENTITIES_H
#define PYBIND_ENTITIES_H

#include "../entities/Atom.h"
#include "../entities/Box.h"
#include "../entities/Molecule.h"
#include "../entities/NeighbourList.h"
#include "../entities/Universe.h"
#include "../entities/UniverseSequence.h"

#include "../utils/CerealUtils.h"

#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pylimer_tools::entities;

struct MoleculeIterator
{
  MoleculeIterator(const Molecule& molecule, py::object ref)
    : molecule(molecule)
    , ref(ref)
  {
  }

  Atom next()
  {
    if (index == molecule.getLength()) {
      throw py::stop_iteration();
    }
    return molecule[index++];
  }

  const Molecule& molecule;
  py::object ref;   // keep a reference
  size_t index = 0; // the index to access
};

void
init_pylimer_bound_entities(py::module_& m)
{
  py::class_<Box>(m, "Box", R"pbdoc(
        The box that the simulation is run in.

        NOTE: 
          currently, only rectangular boxes are supported.
        )pbdoc")
    .def(py::init<const double, const double, const double>())
    .def(py::init<const double,
                  const double,
                  const double,
                  const double,
                  const double,
                  const double>())
    .def("apply_simple_shear",
         &Box::applySimpleShear,
         R"pbdoc(
          Apply a simple shear to the box.

          CAUTION:
            currently, this is not supported for all operations.

          For shear magnitude, you specify the angle :math:`\gamma`.

          For the shearDirection parameter, you can specify `0` for x, `1` for y and `2` for z, respectively.
          Specify another integer to disable the shear.
         )pbdoc",
         py::arg("shear_magnitude"),
         py::arg("shear_direction") = 0)
    .def("get_volume", &Box::getVolume, R"pbdoc(
            Compute the volume of the box.

            :math:`V = L_x \cdot L_y \cdot L_z`
            )pbdoc")
    .def("get_lx", &Box::getLx, R"pbdoc(
            Get the length of the box in x direction.
            )pbdoc")
    .def("get_low_x", &Box::getLowX)
    .def("get_high_x", &Box::getHighX)
    .def("get_ly", &Box::getLy, R"pbdoc(
            Get the length of the box in y direction.
            )pbdoc")
    .def("get_low_y", &Box::getLowY)
    .def("get_high_y", &Box::getHighY)
    .def("get_lz", &Box::getLz, R"pbdoc(
            Get the length of the box in z direction.
            )pbdoc")
    .def("get_low_z", &Box::getLowZ)
    .def("get_high_z", &Box::getHighZ)
    .def("get_offset",
         &Box::getOffset,
         R"pbdoc(
     Compute the offset required to compensate for periodic boundary conditions.

     Useful e.g. if you are using absolute coordinates for distances, but 
     still need an infinite network, 
     e.g., if the bonds need to be able to get longer than half the box.
    )pbdoc",
         py::arg("distances"))
    .def(
      "apply_pbc",
      [](const Box& box, const Eigen::VectorXd& distances) {
        Eigen::VectorXd dist = distances;
        box.handlePBC(dist);
        return dist;
      },
      R"pbdoc(
      Apply periodic boundary conditions (PBC): adjust the specified distances to fit into this box.
      )pbdoc",
      py::arg("distances"))
    .def(py::pickle(
           [](const Box& b) { // __getstate__
             /* Return a tuple that fully encodes the state of the object */
             return py::make_tuple(b.getLowX(),
                                   b.getLowY(),
                                   b.getLowZ(),
                                   b.getHighX(),
                                   b.getHighY(),
                                   b.getHighZ(),
                                   b.getShearMagnitude(),
                                   b.getShearDirection());
           },
           [](py::tuple t) { // __setstate__
             if (t.size() != 6) {
               throw std::runtime_error("Invalid state!.");
             }

             /* Create a new C++ instance */
             Box b = Box(t[0].cast<double>(),
                         t[1].cast<double>(),
                         t[2].cast<double>(),
                         t[3].cast<double>(),
                         t[4].cast<double>(),
                         t[5].cast<double>());
             if (t.size() > 6) {
               b.applySimpleShear(t[6].cast<double>(), t[7].cast<int>());
             }

             return b;
           }),
         "Provides support for pickling.");

  py::class_<Atom>(m, "Atom", R"pbdoc(
       A single bead or atom
  )pbdoc")
    .def(py::init<const long int,
                  const int,
                  const double,
                  const double,
                  const double,
                  const int,
                  const int,
                  const int>(),
         "Construct this atom",
         py::arg("id"),
         py::arg("type"),
         py::arg("x"),
         py::arg("y"),
         py::arg("z"),
         py::arg("nx"),
         py::arg("ny"),
         py::arg("nz"))
    .def("compute_vector_to",
         &Atom::vectorTo,
         R"pbdoc(
            Compute the vector to another atom.
            )pbdoc",
         py::arg("to_atom"),
         py::arg("pbc_box"))
    .def("distance_to",
         &Atom::distanceTo,
         R"pbdoc(
            Compute the distance to another atom.
            )pbdoc",
         py::arg("to_atom"),
         py::arg("pbc_box"))
    .def("vector_to_unwrapped",
         &Atom::vectorToUnwrapped,
         "Compute the vector to another atom respecting the periodic image "
         "flags.")
    .def("distance_to_unwrapped",
         &Atom::distanceToUnwrapped,
         "Compute the distance to another atom respecting the periodic image "
         "flag.")
    .def("get_id", &Atom::getId, R"pbdoc(
            Get the id of the atom.
            )pbdoc")
    .def("get_type", &Atom::getType, R"pbdoc(
            Get the type of the atom.
            )pbdoc")
    .def("get_x", &Atom::getX, R"pbdoc(
            Get the x coordinate of the atom.
            )pbdoc")
    .def("get_y", &Atom::getY, R"pbdoc(
            Get the y coordinate of the atom.
            )pbdoc")
    .def("get_z", &Atom::getZ, R"pbdoc(
            Get the z coordinate of the atom.
            )pbdoc")
    .def("get_nx",
         &Atom::getNX,
         "Get the box image that the atom is in in x direction (also known "
         "as `ix` or `nx`).")
    .def("get_ny",
         &Atom::getNY,
         "Get the box image that the atom is in in y direction (also known "
         "as `iy` or `ny`).")
    .def("get_nz",
         &Atom::getNZ,
         "Get the box image that the atom is in in z direction (also known "
         "as `iz` or `nz`).")
    .def("get_coordinates", [](const Atom& a) { return a.getCoordinates(); })
    .def("get_unwrapped_coordinates",
         [](const Atom& a, const Box& box) {
           return a.getUnwrappedCoordinates(box);
         })
    .def(pybind11::self == pybind11::self)
    //     .def(pybind11::self != pybind11::self)
    .def(py::pickle(
           [](const Atom& b) { // __getstate__
                               // TODO: support extra data
             /* Return a tuple that fully encodes the state of the object */
             return py::make_tuple(b.getId(),
                                   b.getType(),
                                   b.getX(),
                                   b.getY(),
                                   b.getZ(),
                                   b.getNX(),
                                   b.getNY(),
                                   b.getNZ());
           },
           [](py::tuple t) { // __setstate__
             if (t.size() != 8) {
               throw std::runtime_error("Invalid state!.");
             }

             /* Create a new C++ instance */
             Atom a = Atom(t[0].cast<long int>(),
                           t[1].cast<int>(),
                           t[2].cast<double>(),
                           t[3].cast<double>(),
                           t[4].cast<double>(),
                           t[5].cast<int>(),
                           t[6].cast<int>(),
                           t[7].cast<int>());

             return a;
           }),
         "Provides support for pickling");

  py::class_<MoleculeIterator>(m, "MoleculeIterator", R"pbdoc(
       An iterator to iterate throught the atoms in :obj:`~pylimer_tools_cpp.Molecule`.
  )pbdoc")
    .def("__iter__",
         [](MoleculeIterator& it) -> MoleculeIterator& { return it; })
    .def("__next__", &MoleculeIterator::next);

  py::enum_<MoleculeType>(m, "MoleculeType")
    .value("UNDEFINED",
           MoleculeType::UNDEFINED,
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
    .def(py::init<Box&, igraph_t*, MoleculeType, std::map<int, double>>())
    // getters
    .def("get_nr_of_bonds",
         &Molecule::getNrOfBonds,
         "Counts and returns the number of bonds associated with this "
         "molecule.")
    .def("get_nr_of_atoms",
         &Molecule::getNrOfAtoms,
         "Counts and returns the number of atoms associated with this "
         "molecule.")
    .def("get_strand_type", &Molecule::getType, R"pbdoc(
           Get the type of this molecule (see :obj:`~pylimer_tools_cpp.MoleculeType` enum).

           Note that this type might be unset; currently, only 
           :func:`~pylimer_tools_cpp.Universe.get_chains_with_crosslinker` assigns them automatically.
      )pbdoc")
    .def("get_strand_ends",
         &Molecule::getChainEnds,
         R"pbdoc(
          Get the ends of the given strand (= molecule).
          In case of a primary loop, the cross-link is returned, if there is one.
          Use the argument `close_loop` to decide, whether this should be returned once or twice.

          NOTE: 
               Currently only works for linear strands.
     )pbdoc",
         py::arg("crosslinker_type") = 2,
         py::arg("close_loop") = false)
    .def("get_atoms", &Molecule::getAtoms, R"pbdoc(
            Returns all atom objects enclosed in this molecule, ordered by vertex id.
            )pbdoc")
    .def("get_atoms_lined_up",
         &Molecule::getAtomsLinedUp,
         R"pbdoc(
            Returns all atom objects enclosed in this molecule based on the connectivity.

            This method works only for lone chains, atoms and loops, 
            as it throws an error if the molecule does not allow such a "line-up", 
            for example because of crosslinkers.

            Use the `crosslinker_type` parameter to force the atoms in a primary loop 
            to start with the cross-link.
            )pbdoc",
         py::arg("crosslinker_type") = 2,
         py::arg("assumed_coordinates") = false,
         py::arg("close_loop") = false)
    .def("get_atoms_by_type",
         &Molecule::getAtomsOfType,
         R"pbdoc(
            Get the atoms with the specified type.
            )pbdoc",
         py::arg("type"))
    .def("get_atoms_by_degree",
         &Molecule::getAtomsOfDegree,
         R"pbdoc(
            Get the atoms that have the specified number of bonds.
            )pbdoc",
         py::arg("degree"))
    .def("get_atoms_connected_to_vertex",
         &Molecule::getAtomsConnectedTo,
         R"pbdoc(
            Get the atoms connected to a specified vertex id.
            )pbdoc",
         py::arg("vertex_idx"))
    .def("get_atoms_connected_to",
         &Molecule::getConnectedAtoms,
         R"pbdoc(
            Get the atoms connected to a specified atom.

            Internally uses :func:`~pylimer_tools_cpp.Molecule.getAtomsConnectedTo`
            )pbdoc",
         py::arg("atom"))
    .def("get_edges", &Molecule::getEdges, R"pbdoc(
            Get all bonds. Returns a dict with three properties: 'edge_from', 'edge_to' and 'edge_type'.
            The order is not necessarily related to any structural property.
            
            NOTE:
               The integer values returned refer to the vertex ids, not the atom ids.
               Use :func:`~pylimer_tools_cpp.Molecule.getAtomIdByIdx` to translate them to atom ids, or 
               :func:`~pylimer_tools_cpp.Molecule.getBonds` to have that done for you.
            )pbdoc")
    .def("get_bonds", &Molecule::getBonds, R"pbdoc(
            Get all bonds. Returns a dict with three properties: 'bond_from', 'bond_to' and 'bond_type'.
            )pbdoc")
    .def("get_atom_types",
         &Molecule::getAtomTypes,
         "Query all types (each one for each atom) ordered by atom vertex id.")
    .def("get_atom_by_vertex_idx",
         &Molecule::getAtomByVertexIdx,
         R"pbdoc(
            Get an atom for a specific vertex.
            )pbdoc",
         py::arg("vertex_idx"))
    .def(
      "get_atom_by_id",
      [](const Molecule& molecule, size_t id) {
        return molecule.getAtomByVertexIdx(molecule.getIdxByAtomId(id));
      },
      R"pbdoc(
            Get an atom by its id.
            )pbdoc",
      py::arg("atom_id"))
    .def("get_atom_id_by_vertex_idx",
         &Molecule::getAtomIdByIdx,
         R"pbdoc(
            Get the id of the atom by the vertex id of the underlying graph.
            )pbdoc",
         py::arg("vertex_idx"))
    .def("get_vertex_idx_by_atom_id",
         &Molecule::getIdxByAtomId,
         "Get the vertex index of the underlying graph for an atom with a "
         "specified id.",
         py::arg("atom_id"))
    .def("get_key", &Molecule::getKey, R"pbdoc(
            Get a unique identifier for this molecule.
            )pbdoc")
    // computations
    .def("compute_total_mass", &Molecule::computeTotalMass, R"pbdoc(
            Computes the total mass of this molecule.
            )pbdoc")
    .def("compute_bond_lengths",
         &Molecule::computeBondLengths,
         "Computes the length :math:`b` of each bond in the molecule, "
         "respecting periodic boundaries.")
    .def("compute_radius_of_gyration",
         &Molecule::computeRadiusOfGyration,
         R"pbdoc(
            Computes the radius of gyration, :math:`R_g^2` of this molecule.
            
            :math:`{R_g}^2 = \frac{1}{M} \sum_i m_i (r_i - r_{cm})^2`,
            where :math:`M` is the total mass of the molecule, :math:`r_{cm}`
            are the coordinates of the center of mass of the molecule and the
            sum is over all contained atoms.
            )pbdoc")
    .def("compute_radius_of_gyration_with_derived_image_flags",
         &Molecule::computeRadiusOfGyrationWithDerivedImageFlags,
         R"pbdoc(
            Computes the radius of gyration, :math:`R_g^2` of this molecule,
            but ignoring the image flags attached to the atoms.
            This only works for Molecules that can be lined up with 
            :func:`~pylimer_tools_cpp.Molecule.getAtomsLinedUp()`,
            as it needs the atoms sorted such that the periodic box can still be respected somewhat.
            In other words, this function computes the radius of gyration 
            assuming the distance between two lined-up beads 
            is smaller than half the periodic box in each direction.
            
            See also: :func:`~pylimer_tools_cpp.Molecule.computeRadiusOfGyration()`.
            )pbdoc")
    .def("compute_end_to_end_distance",
         &Molecule::computeEndToEndDistance,
         R"pbdoc(
            Compute the end-to-end distance (:math:`R_{ee}`) of this molecule. 

            CAUTION:
               Returns 0.0 if the molecule does not have two or more atoms.
               Returns -1.0 if not exactly 2 ends were found.
               Computes the distance between 2 atoms with functionality 1, 
               ignoring whether they are cross-linkers or not.
            )pbdoc")
    .def("compute_end_to_end_distance_with_derived_image_flags",
         &Molecule::computeEndToEndDistanceWithDerivedImageFlags,
         R"pbdoc(
            Compute the end-to-end distance (:math:`R_{ee}`) of this molecule,
            but ignoring the image flags attached to the atoms. 
            This only works for Molecules that can be lined up with 
            :func:`~pylimer_tools_cpp.Molecule.getAtomsLinedUp()`,
            as it needs the atoms sorted such that the periodic box can still be respected somewhat.

            CAUTION:
               Returns 0.0 if the molecule does not have two or more atoms.
               Requires bonds to be shorter than half the box length.
               Computes the distance between 2 atoms with functionality 1, 
               ignoring whether they are cross-linkers or not.
            )pbdoc")
    .def("compute_total_length", &Molecule::computeTotalLength, R"pbdoc(
     Computes the sum of the lengths of all bonds.
     In most cases, this is equal to the contour length.
    )pbdoc")
    // operators
    .def(pybind11::self == pybind11::self)
    .def(
      "__getitem__",
      [](const Molecule& molecule, size_t index) {
        if (index > molecule.getLength()) {
          throw py::index_error();
        }
        return molecule[index];
      },
      R"pbdoc(
       Access an atom by its vertex index.
  )pbdoc")
    .def("__contains__", &Molecule::containsAtom, R"pbdoc(
          Check whether a particular atom is contained in this molecule.
     )pbdoc")
    .def("__len__", &Molecule::getLength, R"pbdoc(
       Get the number of atoms in this molecule.
  )pbdoc")
    .def(
      "__iter__",
      [](py::object mol) {
        return MoleculeIterator(mol.cast<const Molecule&>(), mol);
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
    .def("__copy__",
         [](const Molecule& molecule) { return Molecule(molecule); });

  py::class_<NeighbourList>(
    m,
    "NeighbourList",
    "Gives access to somewhat fast queries on the neighbourhood of atoms")
    .def(py::init<const std::vector<pylimer_tools::entities::Atom>&,
                  const pylimer_tools::entities::Box&,
                  double>(),
         "Instantiates a new neighbour list",
         py::arg("atoms"),
         py::arg("box"),
         py::arg("cutoff"))
    .def("get_atoms_close_to",
         py::overload_cast<const Atom&, double, double, bool>(
           &NeighbourList::getAtomsCloseTo),
         R"pbdoc(
          List all atoms that are close to a given one. 

          It is possible to request it within a new cutoff, 
          though the underlying neighbour list will not be regenerated.
          For performance reasons, it is recommended to initialize a 
          new NeighbourList if you require a different cutoff, depending on your use case.

          You can use a negative value for the newCutoff to use the cutoff used for 
          filling the neighbour list buckets.


         )pbdoc",
         py::arg("atom"),
         py::arg("upper_cutoff") = 1.0,
         py::arg("lower_cutoff") = 0.0,
         py::arg("unwrapped") = true)
    .def("remove_atom",
         &NeighbourList::removeAtom,
         R"pbdoc(
          Remove an atom from this neighbour list.
          It will not show up when querying for neighbours, 
          but its neighbours cannot be queried either.
         )pbdoc",
         py::arg("atom"),
         py::arg("debug_hint") = "");

  py::class_<Universe>(
    m,
    "Universe",
    "Represents a full Polymer Network structure, a collection of molecules.")
    .def(py::init<const double, const double, const double>(),
         "Instantiate this Universe (Collection of Molecules) providing the "
         "box lengths.",
         py::arg("Lx"),
         py::arg("Ly"),
         py::arg("Lz"))
    // setters
    .def("add_atoms",
         py::overload_cast<const std::vector<long int>&,
                           const std::vector<int>&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           const std::vector<double>&,
                           const std::vector<int>&,
                           const std::vector<int>&,
                           const std::vector<int>&>(&Universe::addAtoms),
         "Add atoms to the Universe, vertices to the underlying graph.",
         py::arg("ids"),
         py::arg("types"),
         py::arg("x"),
         py::arg("y"),
         py::arg("z"),
         py::arg("nx"),
         py::arg("ny"),
         py::arg("nz"))
    .def("remove_atoms",
         &Universe::removeAtoms,
         R"pbdoc(
          Remove atoms and all associated bonds by their atom ids. 
          )pbdoc",
         py::arg("atom_ids"))
    .def("replace_atom",
         &Universe::replaceAtom,
         R"pbdoc(
          Replace the properties of an atom with the properties of another given atom.
          )pbdoc",
         py::arg("atom_id"),
         py::arg("replacement_atom"))
    .def("resample_velocities",
         &Universe::resampleVelocities,
         R"pbdoc()pbdoc",
         py::arg("mean"),
         py::arg("variance"),
         py::arg("seed") = "",
         py::arg("is_2d") = false)
    .def("add_bonds",
         py::overload_cast<const std::vector<long int>&,
                           const std::vector<long int>&>(&Universe::addBonds),
         "Add bonds to the underlying atoms, edges to the underlying graph. "
         "If the connected atoms are not found, the bonds are silently "
         "skipped.",
         py::arg("bonds_from"),
         py::arg("bonds_to"))
    .def("add_bonds",
         py::overload_cast<const size_t,
                           const std::vector<long int>&,
                           const std::vector<long int>&,
                           const std::vector<int>&,
                           const bool,
                           const bool>(&Universe::addBonds),
         "Add bonds to the underlying atoms, edges to the underlying graph. ",
         py::arg("nr_of_bonds"),
         py::arg("bonds_from"),
         py::arg("bonds_to"),
         py::arg("bond_types"),
         py::arg("ignore_non_existent_atoms") = false,
         py::arg("simplify_universe") = true)
    .def("add_bonds_with_types",
         py::overload_cast<const std::vector<long int>&,
                           const std::vector<long int>&,
                           const std::vector<int>&>(&Universe::addBonds),
         "Add bonds to the underlying atoms, edges to the underlying graph. "
         "If the connected atoms are not found, the bonds are silently "
         "skipped.",
         py::arg("bonds_from"),
         py::arg("bonds_to"),
         py::arg("bond_types"))
    .def("remove_bonds",
         &Universe::removeBonds,
         R"pbdoc(
          Remove bonds by their connected atom ids. 
          )pbdoc",
         py::arg("bonds_from"),
         py::arg("bonds_to"))
    .def("remove_bonds_by_type",
         &Universe::removeBondsOfType,
         R"pbdoc(
          Remove bonds with a specific type. 
          )pbdoc",
         py::arg("bond_type"))
    .def("add_angles",
         &Universe::addAngles,
         "Add angles to the Universe. No relation to the underlying graph, "
         "just a method to preserve read & write capabilities",
         py::arg("angles_from"),
         py::arg("angles_via"),
         py::arg("angles_to"),
         py::arg("angle_types"))
    .def("add_dihedral_angles",
         &Universe::addDihedralAngles,
         "Add dihedral angles to the Universe. No relation to the underlying "
         "graph, "
         "just a method to preserve read & write capabilities",
         py::arg("angles_from"),
         py::arg("angles_via1"),
         py::arg("angles_via2"),
         py::arg("angles_to"),
         py::arg("angle_types"))
    .def("remove_all_angles", &Universe::removeAllAngles, R"pbdoc(

     )pbdoc")
    .def(
      "remove_all_dihedral_angles", &Universe::removeAllDihedralAngles, R"pbdoc(
          
     )pbdoc")
    .def("hash_angle_type",
         &Universe::hashAngleType,
         R"pbdoc(
          Convert the three integers to one long number/hash.
          Used internally for duplicate detection.
     )pbdoc",
         py::arg("angle_from"),
         py::arg("angle_via"),
         py::arg("angle_to"))
    .def("hash_dihedral_angle_type",
         &Universe::hashDihedralAngleType,
         R"pbdoc(
          Convert the four integers to one long number/hash.
          Used internally for duplicate detection.
     )pbdoc",
         py::arg("angle_from"),
         py::arg("angle_via1"),
         py::arg("angle_via2"),
         py::arg("angle_to"))
    .def("set_masses",
         &Universe::setMasses,
         "Set the mass per type of atom.",
         py::arg("massPerType"))
    .def("set_timestep",
         &Universe::setTimestep,
         "Set the time-step when this Universe was captured.",
         py::arg("timestep"))
    .def("set_box_lengths",
         &Universe::setBoxLengths,
         R"pbdoc(
          Override the currently assigned box with one with the side lengths specified.
          )pbdoc",
         py::arg("lx"),
         py::arg("ly"),
         py::arg("lz"),
         py::arg("rescale_atoms") = false)
    .def("set_box",
         &Universe::setBox,
         R"pbdoc(
          Override the currently assigned box with the one specified.
          )pbdoc",
         py::arg("box"),
         py::arg("rescale_atoms") = false)
    .def("set_vertex_property",
         &Universe::setPropertyValue<double>,
         R"pbdoc(
          Set a specific property for a specific vertex.
          )pbdoc",
         py::arg("vertex_id"),
         py::arg("property_name"),
         py::arg("value"))
    // getters
    .def("get_clusters", &Universe::getClusters, R"pbdoc(
          Get the components of the universe that are not connected to each other.
          Returns a list of :obj:`~pylimer_tools_cpp.Universe`s.
          Unconnected, free atoms/beads become their own :obj:`~pylimer_tools_cpp.Universe`.
          )pbdoc")
    .def("get_molecules",
         &Universe::getMolecules,
         R"pbdoc(
          Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms.
          
          Reduces the Universe to a list of molecules. 
          Specify the crosslinker_type to an existing type id, 
          then those atoms will be omitted, and this function returns chains instead.
          )pbdoc",
         py::arg("atom_type_to_omit"))
    .def("get_atoms_connected_to_vertex",
         &Universe::getAtomsConnectedTo,
         R"pbdoc(
            Get the atoms connected to a specified vertex id.
            )pbdoc",
         py::arg("vertex_idx"))
    .def("get_atoms_connected_to",
         &Universe::getConnectedAtoms,
         R"pbdoc(
            Get the atoms connected to a specified atom.

            Internally uses :func:`~pylimer_tools_cpp.Universe.getAtomsConnectedTo`
            )pbdoc",
         py::arg("atom"))
    .def("get_atoms_by_degree",
         &Universe::getAtomsOfDegree,
         R"pbdoc(
            Get the atoms that have the specified number of bonds.
            )pbdoc",
         py::arg("functionality"))
    .def("find_loops",
         &Universe::findLoopsOfAtoms,
         R"pbdoc(
            Decompose the Universe into loops.
            The primary index specifies the degree of the loop.

            CAUTION:
               There are exponentially many paths between two cross-linkers of a network,
               and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
               You can use the maxLength parameter to restrict the algorithm to only search for loops up to a certain length.
               Use a negative value to find all loops and paths.
            )pbdoc",
         py::arg("crosslinker_type"),
         py::arg("max_length") = -1,
         py::arg("skip_self_loops") = false)
    .def("find_minimal_order_loop_from",
         &Universe::findMinimalOrderLoopFrom,
         R"pbdoc(
            Decompose the Universe into loops.
            The primary index specifies the degree of the loop.

            CAUTION:
               There are exponentially many paths between two cross-linkers of a network,
               and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
               You can use the maxLength parameter to restrict the algorithm to only search for loops up to a certain length.
               Use a negative value to find all loops and paths.
            )pbdoc",
         py::arg("loop_start"),
         py::arg("loop_step1"),
         py::arg("max_length") = -1,
         py::arg("skip_self_loops") = false)
    .def("count_atom_types",
         &Universe::countAtomTypes,
         R"pbdoc(
          Count how often each atom type is present.
     )pbdoc")
    .def("count_atoms_in_skin_distance",
         &Universe::countAtomsInSkinDistance,
         R"pbdoc(
          This is a function that may help you to compute the radial distribution function.
          It loops the 

          Parameters:
               - distances: the edges of the bins
               - unwrapped: whether to measure the distance in unwrapped coordinates or as PBC-corrected distance
     )pbdoc",
         py::arg("distances"),
         py::arg("unwrapped") = false)
    .def("get_chains_with_crosslinker",
         &Universe::getChainsWithCrosslinker,
         R"pbdoc(
            Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms, without omitting the cross-linkers.
            In turn, e.g. for a tetrafunctional cross-linker, it will be 4 times in the resulting molecules.
            
            NOTE:
               Cross-linkers without bonds to non-cross-linkers are not returned 
               (i.e., cross-linker-cross-linker bonds, or single cross-linkers, are not counted as strands).
           )pbdoc",
         py::arg("crosslinker_type"))
    .def("get_network_of_crosslinker",
         &Universe::getNetworkOfCrosslinker,
         R"pbdoc(
            Reduce the network to contain only cross-linkers, replacing all the strands with a single bond.
            Useful e.g. to reduce the memory useage and runtime of 
            :func:`~pylimer_tools_cpp.Universe.findLoops()` or 
            :func:`~pylimer_tools_cpp.Universe.hasInfiniteStrand()`.
            
            Further use :func:`~pylimer_tools_cpp.Universe.simplify()` to remove primary loops.
          )pbdoc",
         py::arg("crosslinker_type"))
    .def("get_atom_types",
         &Universe::getAtomTypes,
         R"pbdoc(
          Get all types (each one for each atom) ordered by atom vertex id.
          )pbdoc")
    .def("get_atom",
         &Universe::getAtom,
         "Find an atom by its ID.",
         py::arg("atom_id"))
    .def(
      "get_atom_by_vertex_id",
      &Molecule::getAtomByVertexIdx,
      R"pbdoc(Find an atom by the ID of the vertex of the underlying graph.)pbdoc",
      py::arg("vertex_id"))
    .def("get_atoms", &Universe::getAtoms, R"pbdoc(
            Get all atoms.
            )pbdoc")
    .def("get_atoms_by_type",
         &Universe::getAtomsOfType,
         R"pbdoc(
            Query all atoms by their type.
            )pbdoc",
         py::arg("atom_type"))
    .def("get_atom_id_by_vertex_idx",
         &Universe::getAtomIdByIdx,
         "Get the id of the atom by the vertex id of the underlying graph.",
         py::arg("vertex_id"))
    .def("get_vertex_idx_by_atom_id",
         &Universe::getIdxByAtomId,
         "Get the vertex id of the underlying graph for an atom with a "
         "specified id.",
         py::arg("atom_id"))
    .def("get_edges", &Universe::getEdges, R"pbdoc(
            Get all edges. Returns a dict with three properties: 'edge_from', 'edge_to' and 'edge_type'.
            The order is not necessarily related to any structural characteristic.
            
            NOTE:
               The integer values returned refer to the vertex ids, not the atom ids.
               Use :func:`~pylimer_tools_cpp.Universe.getAtomIdByIdx` to translate them to atom ids, or
               :func:`~pylimer_tools_cpp.Universe.getBonds` to have that done for you.
            )pbdoc")
    .def("get_bonds", &Universe::getBonds, R"pbdoc(
            Get all bonds. Returns a dict with three properties: 'bond_from', 'bond_to' and 'bond_type'.
            The order is not necessarily related to any structural characteristic.
            )pbdoc")
    .def("get_angles", &Universe::getAngles, R"pbdoc(
           Get all angles added to this network.

           Returns a dict with three properties: 'angle_from', 'angle_via' and 'angle_to'.

           NOTE:
               The integer values returned refer to the the atom ids, not the vertex ids.
               Use :func:`~pylimer_tools_cpp.Universe.getIdxByAtomId` to translate them to vertex ids.
           )pbdoc")
    .def("get_box", &Universe::getBox, R"pbdoc(
            Get the underlying bounding box object.
            )pbdoc")
    .def("get_masses", &Universe::getMasses, R"pbdoc(
            Get the mass of one atom per type
            )pbdoc")
    .def("get_volume", &Universe::getVolume, R"pbdoc(
            Query the volume of the underlying bounding box.
            )pbdoc")
    .def("get_nr_of_atoms", &Universe::getNrOfAtoms, R"pbdoc(
            Query the number of atoms in this universe.
            )pbdoc")
    .def("get_nr_of_bonds", &Universe::getNrOfBonds, R"pbdoc(
            Query the number of bonds associated with this universe.
            )pbdoc")
    .def("get_nr_of_angles", &Universe::getNrOfAngles, R"pbdoc(
            Query the number of angles that have been added to this universe.
            )pbdoc")
    .def("get_nr_of_dihedral_angles", &Universe::getNrOfDihedralAngles, R"pbdoc(
            Query the number of dihedralangles that have been added to this universe.
            )pbdoc")
    .def("get_timestep", &Universe::getTimestep, R"pbdoc(
            Query the timestep when this universe was captured.
            )pbdoc")
    .def(
      "get_nr_of_bonds_of_atom",
      &Universe::computeFunctionalityForAtom,
      R"pbdoc(Count the number of immediate neighbours of an atom, specified by its id.)pbdoc")
    .def(
      "get_nr_of_bonds_of_vertex",
      &Universe::computeFunctionalityForVertex,
      R"pbdoc(Count the number of immediate neighbours of an atom, specified by its vertex id.)pbdoc")
    // computations
    .def("compute_bond_lengths",
         &Universe::computeBondLengths,
         "Computes the length :math:`b` of each bond in the molecule, "
         "respecting periodic boundaries.")
    .def("detect_angles",
         &Universe::detectAngles,
         R"pbdoc(Returns just as 
          :func:`~pylimer_tools_cpp.Universe.getAngles`, 
          but all angles that are detected in the network, rather than the one already set.
          Note that the angle types are determined by 
          :func:`~pylimer_tools_cpp.Universe.hashAngleType`,
          which serves angle types that should be mapped by you back to smaller numbers, 
          before serving them to :func:`~pylimer_tools_cpp.Universe.addAngles`.
         )pbdoc")
    .def("detect_dihedral_angles",
         &Universe::detectDihedralAngles,
         R"pbdoc(Returns just as 
          :func:`~pylimer_tools_cpp.Universe.getDihedralAngles`, 
          but all dihedral angles that are detected in the network, rather than the one already set.
          Note that the angle types are determined by 
          :func:`~pylimer_tools_cpp.Universe.hashDihedralAngleType`,
          which serves angle types that should be mapped by you back to smaller numbers, 
          before serving them to :func:`~pylimer_tools_cpp.Universe.addDiheralAngles`.
         )pbdoc")
    .def("has_infinite_strand",
         &Universe::hasInfiniteStrand,
         R"pbdoc(
           Checks whether there is a strand (with cross-linker) in the universe that loops through periodic images without coming back.
           
            CAUTION:
               There are exponentially many paths between two cross-linkers of a network,
               and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
           )pbdoc")
    .def("determine_functionality_per_type",
         &Universe::determineFunctionalityPerType,
         R"pbdoc(
            Find the maximum functionality of each atom type in the network.
            )pbdoc")
    .def("determine_effective_functionality_per_type",
         &Universe::determineEffectiveFunctionalityPerType,
         R"pbdoc(
            Find the average functionality of each atom type in the network.
            )pbdoc")
    .def("compute_mean_strand_length",
         &Universe::getMeanStrandLength,
         R"pbdoc(
              Compute the mean number of beads per strand.
              )pbdoc",
         py::arg("crosslinker_type"))
    .def("compute_total_mass", &Universe::computeTotalMass, R"pbdoc(
          Compute the total mass of this network/universe in whatever mass unit was used when 
          :func:`~pylimer_tools_cpp.Universe.setMasses()` was called.
     )pbdoc")
    .def("compute_number_average_molecular_weight",
         &Universe::computeNumberAverageMolecularWeight,
         R"pbdoc(
              Compute the number average molecular weight.

              NOTE: 
                    Cross-linkers are ignored completely.
              )pbdoc",
         py::arg("crosslinker_type"))
    .def("compute_weight_average_molecular_weight",
         &Universe::computeWeightAverageMolecularWeight,
         R"pbdoc(
              Compute the weight average molecular weight.

              NOTE: 
                    Cross-linkers are ignored completely.
              )pbdoc",
         py::arg("crosslinker_type"))
    .def("compute_polydispersity_index",
         &Universe::computePolydispersityIndex,
         R"pbdoc(
              Compute the polydispersity index: 
              the weight average molecular weight over the number average molecular weight.
              )pbdoc",
         py::arg("crosslinker_type"))
    .def("compute_weight_fractions", &Universe::computeWeightFractions, R"pbdoc(
            Compute the weight fractions of each atom type in the network.

            If no masses are stored, 
            )pbdoc")
    .def("compute_end_to_end_distances",
         &Universe::computeEndToEndDistances,
         R"pbdoc(
          Compute the end-to-end distance of each strand in the network.

          NOTE:
               Internally, this uses the :func:`~pylimer_tools_cpp.Molecule.computeEndToEndDistance`.
               All its cautionary facts apply.
     )pbdoc")
    .def("compute_mean_end_to_end_distance",
         &Universe::computeMeanEndToEndDistance,
         R"pbdoc(
          Computes the mean of the end-to-end distances of each strand in the network.

          NOTE:
               Internally, this uses the :func:`~pylimer_tools_cpp.Molecule.computeEndToEndDistance`.
               All its cautionary facts apply.
               Invalid strands (where said function returns 0.0 or -1.0) are ignored.
     )pbdoc")
    .def("compute_mean_square_end_to_end_distance",
         &Universe::computeMeanSquareEndToEndDistance,
         R"pbdoc(
          Computes the mean square of the end-to-end distances of each strand in the network.

          NOTE:
               Internally, this uses the :func:`~pylimer_tools_cpp.Molecule.computeEndToEndDistance`.
               All its cautionary facts apply.
               Invalid strands (where said function returns 0.0 or -1.0) are ignored.
     )pbdoc",
         py::arg("crosslinker_type"),
         py::arg("only_those_with_two_crosslinkers") = false)
    .def("compute_dxs",
         &Universe::computeDxs,
         "Compute the dx distance for certain bonds (length in x direction).",
         py::arg("atomIdsTo"),
         py::arg("atomIdsFrom"))
    .def("compute_dys",
         &Universe::computeDys,
         "Compute the dy distance for certain bonds (length in y direction).",
         py::arg("atomIdsTo"),
         py::arg("atomIdsFrom"))
    .def("compute_dzs",
         &Universe::computeDzs,
         "Compute the dz distance for certain bonds (length in z direction).",
         py::arg("atomIdsTo"),
         py::arg("atomIdsFrom"))
    .def(
      "compute_temperature",
      &Universe::computeTemperature,
      R"pbdoc(Use the velocities per atom to compute the temperature from the kinetic energy of the system.)pbdoc",
      py::arg("dimensions") = 3,
      py::arg("k_b") = 1.)
    .def("simplify",
         &Universe::simplify,
         "Remove self links and double bonds. This function is called "
         "automatically after adding bonds.")
    .def(py::pickle(
      [](const Universe& u) {
        return py::make_tuple(pylimer_tools::utils::serializeToString(u));
      },
      [](py::tuple t) {
        std::string in = t[0].cast<std::string>();
        Universe u;
        pylimer_tools::utils::deserializeFromString(u, in);
        return u;
      }))
    .def("__copy__",
         [](const Universe& universe) { return Universe(universe); });

  struct LazyUniverseSequenceIterator
  {
    LazyUniverseSequenceIterator(UniverseSequence& us, py::object ref)
      : us(us)
      , ref(ref)
    {
    }

    Universe next()
    {
      if (index == us.getLength()) {
        throw py::stop_iteration();
      }
      Universe toReturn = us.atIndex(index);
      us.forgetAtIndex(index);
      index += 1;
      return toReturn;
    }

    UniverseSequence& us;
    py::object ref;   // keep a reference
    size_t index = 0; // the index to access
  };

  py::class_<LazyUniverseSequenceIterator>(
    m, "LazyUniverseSequenceIterator", R"pbdoc(
       An iterator to iterate throught the universes in :obj:`~pylimer_tools_cpp.UniverseSequence`.
  )pbdoc")
    .def("__iter__",
         [](const LazyUniverseSequenceIterator& it)
           -> const LazyUniverseSequenceIterator& { return it; })
    .def("__next__", &LazyUniverseSequenceIterator::next);

  py::class_<UniverseSequence>(m, "UniverseSequence", R"pbdoc(
     This class represents a sequence of Universes, with the Universe's data
     only being read on request. Dump files are read at once in order
     to know how many timesteps/universes are available in total 
     (but the universes' data is not read on first look through the file).
     This, while it can lead to two (or more) reads of the whole file, 
     is a measure in order to enable low memory useage if needed (i.e. for large dump files).
     Use Python's iterator to have this UniverseSequence only ever retain one universe in memory.
     Alternatively, use :func:`~pylimer_tools_cpp.UniverseSequence.forgetAtIndex`
     to have the UniverseSequence forget about already read universes.
     )pbdoc")
    .def(py::init<>())
    .def("initialize_from_dump_file",
         &UniverseSequence::initializeFromDumpFile,
         R"pbdoc(
          Reset and initialize the Universes from a Lammps :code:`dump` output. 
        
          NOTE:
               If you have not output the id of the atoms in the dump file, they will be assigned sequentially. 
               If you have not output the type of the atoms in the dump file, they will be set to -1 if they cannot be infered from the data file.
        )pbdoc",
         py::arg("initial_data_file"),
         py::arg("dump_file"))
    .def("initialize_from_data_sequence",
         &UniverseSequence::initializeFromDataSequence,
         "Reset and initialize the Universes from an ordered list of Lammps "
         "data (:code:`write_data`) files.",
         py::arg("data_files"))
    .def("set_data_file_atom_style",
         &UniverseSequence::setDataFileAtomStyle,
         R"pbdoc(
          Set the format of the data files to be read. See :obj:`~pylimer_tools_cpp.AtomStyle`.
     )pbdoc",
         py::arg("atom_styles"))
    .def("next",
         &UniverseSequence::next,
         R"pbdoc(Get the Universe that's next in the sequence.)pbdoc")
    .def("at_index",
         &UniverseSequence::atIndex,
         "Get the Universe at the given index (as of in the sequence given "
         "by the dump file).",
         py::arg("index"))
    .def(
      "forget_at_index",
      &UniverseSequence::forgetAtIndex,
      R"pbdoc(Clear the memory of the Universe at the given index (as of in the 
           sequence given by the dump file).)pbdoc",
      py::arg("index"))
    .def("reset_iterator",
         &UniverseSequence::resetIterator,
         R"pbdoc(
          Reset the internal iterator, such that a subsequent call to 
          :func:`~pylimer_tools_cpp.UniverseSequence.next` returns the first one again.
          )pbdoc")
    .def("get_length", &UniverseSequence::getLength, R"pbdoc(
            Get the number of universes in this sequence.
            )pbdoc")
    .def("get_all", &UniverseSequence::getAll, R"pbdoc(
            Get all universes initialized back in a list.
            For big dump files or lots of data files, this might lead to memory issues.
            Use :func:`~pylimer_tools_cpp.UniverseSequence.__iter__`
            to have
            or :func:`~pylimer_tools_cpp.UniverseSequence.atIndex`
            and :func:`~pylimer_tools_cpp.UniverseSequence.forgetAtIndex`
            to craft a more memory-efficient retrieval mechanism.
            )pbdoc")
    // computations
    .def("compute_msd_for_atoms",
         &UniverseSequence::computeMsdForAtoms,
         R"pbdoc(
          Compute the mean square displacement for atoms with the specified ids
     )pbdoc",
         py::arg("atom_ids"),
         py::arg("nr_of_origins") = 25,
         py::arg("reduce_memory") = false)
    .def("compute_msd_for_atom_properties",
         &UniverseSequence::computeMsdForAtomProperties,
         R"pbdoc(
          Compute the mean square displacement for atoms with the specified ids
     )pbdoc",
         py::arg("atom_ids"),
         py::arg("x_property"),
         py::arg("y_property"),
         py::arg("z_property"),
         py::arg("nr_of_origins") = 25,
         py::arg("reduce_memory") = false)
    .def("compute_distance_autocorrelation_from_to",
         &UniverseSequence::computeDistanceAutocorrelationFromToAtoms,
         R"pbdoc(
          Compute the autocorrelation of the dot product of the distance vector from certain to other atoms.

          For example, this can be used to compute Eq. 4.51 from Masao Doi, Introduction to Polymer Physics, p. 74.
         )pbdoc",
         py::arg("atom_ids_from"),
         py::arg("atom_ids_to"),
         py::arg("nr_of_origins") = 25,
         py::arg("reduce_memory") = false)
    .def("compute_distance_from_to_atoms",
         &UniverseSequence::computeDistanceFromToAtoms,
         R"pbdoc(
          Compute the root square norm of all the (unwrapped!) distances for the given pair of atoms.
          
          Can be used to somewhat faster compute e.g. all the end-to-end or bond distances.
          Pay attention that the image flags are correct, otherwise, this data may not be useable.
         )pbdoc",
         py::arg("atom_ids_from"),
         py::arg("atom_ids_to"),
         py::arg("reduce_memory") = false)
    .def("compute_vector_from_to_atoms",
         &UniverseSequence::computeVectorFromToAtoms,
         R"pbdoc(
          Compute the (unwrapped!) distances for the given pair of atoms.
          
          Can be used to somewhat faster compute e.g. all the end-to-end or bond vectors.
          Pay attention that the image flags are correct, otherwise, this data may not be useable.
         )pbdoc",
         py::arg("atom_ids_from"),
         py::arg("atom_ids_to"),
         py::arg("reduce_memory") = false)
    // operators
    .def(
      "__getitem__",
      [](UniverseSequence& us, size_t index) {
        if (index > us.getLength()) {
          throw py::index_error();
        }
        return us.atIndex(index);
      },
      "Get a universe by its index.")
    .def("__len__", &UniverseSequence::getLength, "Get the number of universes")
    .def(
      "__iter__",
      [](py::object us) {
        return LazyUniverseSequenceIterator(us.cast<UniverseSequence&>(), us);
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
