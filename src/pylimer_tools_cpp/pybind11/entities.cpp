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
  py::class_<Box>(m,
                  "Box",
                  R"pbdoc(
        The box that the simulation is run in.

        .. note:: 
          Currently, only rectangular boxes are supported.
        )pbdoc",
                  py::module_local())
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

          .. warning::
            Currently, this is not supported for all operations.

          For shear magnitude, you specify the angle :math:`\gamma`.

          :param shear_magnitude: The shear magnitude (angle :math:`\gamma`)
          :param shear_direction: Direction of shear: 0 for x, 1 for y, 2 for z. 
                                 Use any other integer to disable shear.
         )pbdoc",
         py::arg("shear_magnitude"),
         py::arg("shear_direction") = 0)
    .def("get_volume", &Box::getVolume, R"pbdoc(
            Compute the volume of the box.

            :math:`V = L_x \cdot L_y \cdot L_z`

            :return: The volume of the box
            )pbdoc")
    .def(
      "get_l", [](const Box& box) { return box.getL(); }, R"pbdoc(
          Get the three lengths of the box in an array/list.

          :return: Array containing [Lx, Ly, Lz] box dimensions
      )pbdoc")
    .def("get_lx", &Box::getLx, R"pbdoc(
            Get the length of the box in x direction.

            :return: The x-dimension length of the box
            )pbdoc")
    .def("get_low_x", &Box::getLowX, R"pbdoc(
            Get the lower bound of the box in x direction.

            :return: The lower x-coordinate boundary
            )pbdoc")
    .def("get_high_x", &Box::getHighX, R"pbdoc(
            Get the upper bound of the box in x direction.

            :return: The upper x-coordinate boundary
            )pbdoc")
    .def("get_ly", &Box::getLy, R"pbdoc(
            Get the length of the box in y direction.

            :return: The y-dimension length of the box
            )pbdoc")
    .def("get_low_y", &Box::getLowY, R"pbdoc(
            Get the lower bound of the box in y direction.

            :return: The lower y-coordinate boundary
            )pbdoc")
    .def("get_high_y", &Box::getHighY, R"pbdoc(
            Get the upper bound of the box in y direction.

            :return: The upper y-coordinate boundary
            )pbdoc")
    .def("get_lz", &Box::getLz, R"pbdoc(
            Get the length of the box in z direction.

            :return: The z-dimension length of the box
            )pbdoc")
    .def("get_low_z", &Box::getLowZ, R"pbdoc(
            Get the lower bound of the box in z direction.

            :return: The lower z-coordinate boundary
            )pbdoc")
    .def("get_high_z", &Box::getHighZ, R"pbdoc(
            Get the upper bound of the box in z direction.

            :return: The upper z-coordinate boundary
            )pbdoc")
    .def("get_bounding_box", &Box::getBoundingBox, R"pbdoc(
     Get an orthogonal box that encloses this box.
     
     For non-sheared boxes, the resulting box is identical to the current box.

     :return: A new Box object representing the bounding box
    )pbdoc")
    .def(
      "apply_pbc",
      [](const Box& box, const Eigen::VectorXd& distances) {
        Eigen::VectorXd dist = distances;
        box.handlePBC(dist);
        return dist;
      },
      R"pbdoc(
      Apply periodic boundary conditions (PBC): adjust the specified distances to fit into this box.
      
      :param distances: The distances to adjust
      :return: The adjusted distances
      )pbdoc",
      py::arg("distances"))
    .def("get_offset",
         &Box::getOffset,
         R"pbdoc(
     Compute the offset required to compensate for periodic boundary conditions.

     Useful e.g. if you are using absolute coordinates for distances, but 
     still need an infinite network, 
     e.g., if the bonds need to be able to get longer than half the box.
     
     :param distances: The distances to compute offset for
     :return: The computed offset
    )pbdoc",
         py::arg("distances"))
    .def("is_valid_offset",
         &Box::isValidOffset,
         R"pbdoc(
          Check whether the passed offset is a valid one in this box.
          
          :param potential_offset: The offset to validate
          :param abs_precision: Absolute precision for the validation
          :return: True if the offset is valid, False otherwise
     )pbdoc",
         py::arg("potential_offset"),
         py::arg("abs_precision") = 1e-5)
#ifdef CEREALIZABLE
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
#else
    ;
#endif

  py::class_<Atom>(m,
                   "Atom",
                   R"pbdoc(
       A single bead or atom.
  )pbdoc",
                   py::module_local())
    .def(py::init<const long int,
                  const int,
                  const double,
                  const double,
                  const double,
                  const int,
                  const int,
                  const int>(),
         R"pbdoc(
         Construct this atom.

         :param id: Unique identifier for the atom
         :param type: Type classification of the atom
         :param x: X coordinate position
         :param y: Y coordinate position  
         :param z: Z coordinate position
         :param nx: Periodic image flag in x direction
         :param ny: Periodic image flag in y direction
         :param nz: Periodic image flag in z direction
         )pbdoc",
         py::arg("id"),
         py::arg("type"),
         py::arg("x"),
         py::arg("y"),
         py::arg("z"),
         py::arg("nx"),
         py::arg("ny"),
         py::arg("nz"))
    .def(py::init<std::unordered_map<std::string, double>&>(),
         R"pbdoc(
         Construct this atom from a properties dictionary.

         The dictionary should contain at least the following keys:
         - "id": Unique identifier for the atom
         - "type": Type classification of the atom  
         - "x", "y", "z": Coordinate positions
         - "nx", "ny", "nz": Periodic image flags
         
         Any additional properties will be stored as extra data.

         :param properties: Dictionary containing atom properties
         )pbdoc",
         py::arg("properties"))
    .def("compute_vector_to",
         &Atom::vectorTo,
         R"pbdoc(
            Compute the vector to another atom.
            
            :param to_atom: The target atom
            :param pbc_box: The periodic boundary conditions box
            :return: Vector pointing from this atom to the target atom
            )pbdoc",
         py::arg("to_atom"),
         py::arg("pbc_box"))
    .def("distance_to",
         &Atom::distanceTo,
         R"pbdoc(
            Compute the distance to another atom.
            
            :param to_atom: The target atom
            :param pbc_box: The periodic boundary conditions box
            :return: Euclidean distance between the atoms
         )pbdoc",
         py::arg("to_atom"),
         py::arg("pbc_box"))
    .def("vector_to_unwrapped",
         &Atom::vectorToUnwrapped,
         R"pbdoc(
         Compute the vector to another atom respecting the periodic image flags.
         
         :param to_atom: The target atom
         :return: Unwrapped vector to the target atom
         )pbdoc")
    .def("distance_to_unwrapped",
         &Atom::distanceToUnwrapped,
         R"pbdoc(
         Compute the distance to another atom respecting the periodic image flags.
         
         :param to_atom: The target atom
         :return: Unwrapped distance to the target atom
         )pbdoc")
    .def("mean_position_with",
         &Atom::meanPositionWith,
         R"pbdoc(
         Compute the mean position between this atom and another atom, considering periodic boundaries.
         
         :param other_atom: The other atom
         :param pbc_box: The periodic boundary conditions box
         :return: Vector representing the mean position
         )pbdoc",
         py::arg("other_atom"),
         py::arg("pbc_box"))
    .def("mean_position_with_unwrapped",
         &Atom::meanPositionWithUnwrapped,
         R"pbdoc(
         Compute the mean position between this atom and another atom using unwrapped coordinates.
         
         :param other_atom: The other atom
         :param pbc_box: The periodic boundary conditions box  
         :return: Vector representing the mean position (unwrapped)
         )pbdoc",
         py::arg("other_atom"),
         py::arg("pbc_box"))
    .def("get_id", &Atom::getId, R"pbdoc(
            Get the ID of the atom.
            
            :return: The atom's unique identifier
            )pbdoc")
    .def("get_type", &Atom::getType, R"pbdoc(
            Get the type of the atom.
            
            :return: The atom's type classification
            )pbdoc")
    .def("get_x", &Atom::getX, R"pbdoc(
            Get the x coordinate of the atom.
            
            :return: The x coordinate
            )pbdoc")
    .def("get_y", &Atom::getY, R"pbdoc(
            Get the y coordinate of the atom.
            
            :return: The y coordinate
            )pbdoc")
    .def("get_z", &Atom::getZ, R"pbdoc(
            Get the z coordinate of the atom.
            
            :return: The z coordinate
            )pbdoc")
    .def("get_nx",
         &Atom::getNX,
         R"pbdoc(
         Get the box image that the atom is in in x direction (also known as `ix` or `nx`).
         
         :return: The periodic image flag in x direction
         )pbdoc")
    .def("get_ny",
         &Atom::getNY,
         R"pbdoc(
         Get the box image that the atom is in in y direction (also known as `iy` or `ny`).
         
         :return: The periodic image flag in y direction
         )pbdoc")
    .def("get_nz",
         &Atom::getNZ,
         R"pbdoc(
         Get the box image that the atom is in in z direction (also known as `iz` or `nz`).
         
         :return: The periodic image flag in z direction
         )pbdoc")
    .def("get_unwrapped_x",
         &Atom::getUnwrappedX,
         R"pbdoc(
         Get the unwrapped x coordinate of the atom.
         
         :param box: The simulation box to use for unwrapping
         :return: The unwrapped x coordinate
         )pbdoc",
         py::arg("box"))
    .def("get_unwrapped_y",
         &Atom::getUnwrappedY,
         R"pbdoc(
         Get the unwrapped y coordinate of the atom.
         
         :param box: The simulation box to use for unwrapping
         :return: The unwrapped y coordinate
         )pbdoc",
         py::arg("box"))
    .def("get_unwrapped_z",
         &Atom::getUnwrappedZ,
         R"pbdoc(
         Get the unwrapped z coordinate of the atom.
         
         :param box: The simulation box to use for unwrapping
         :return: The unwrapped z coordinate
         )pbdoc",
         py::arg("box"))
    .def(
      "get_coordinates",
      [](const Atom& a) { return a.getCoordinates(); },
      R"pbdoc(
         Get the coordinates of this atom as a vector.
         
         :return: A vector containing the x, y, z coordinates
         )pbdoc")
    .def(
      "get_unwrapped_coordinates",
      [](const Atom& a, const Box& box) {
        return a.getUnwrappedCoordinates(box);
      },
      R"pbdoc(
         Get the unwrapped coordinates of this atom.
         
         :param box: The simulation box to use for unwrapping
         :return: A vector containing the unwrapped x, y, z coordinates
         )pbdoc")
    .def("get_extra_data",
         &Atom::getExtraData,
         R"pbdoc(
         Get all extra data properties stored with this atom (e.g., charge, dipole, etc.).
         
         :return: Dictionary containing all extra properties
         )pbdoc")
    .def("get_property",
         &Atom::getProperty,
         R"pbdoc(
         Get a specific property value from the extra data.
         
         :param property: The name of the property to retrieve
         :return: The value of the specified property
         :raises: std::out_of_range if the property doesn't exist
         )pbdoc",
         py::arg("property"))
    .def(py::pickle(
           [](const Atom& b) { // __getstate__
             /* Return a tuple that fully encodes the state of the object */
             return py::make_tuple(b.getId(),
                                   b.getType(),
                                   b.getX(),
                                   b.getY(),
                                   b.getZ(),
                                   b.getNX(),
                                   b.getNY(),
                                   b.getNZ(),
                                   b.getExtraData());
           },
           [](py::tuple t) { // __setstate__
             if (t.size() != 9 && t.size() != 8) {
               throw std::runtime_error("Invalid state! Expected 9 elements.");
             }

             /* If there's extra data, create atom with properties instead */
             if (t.size() == 9) {
               auto extra_data =
                 t[8].cast<std::unordered_map<std::string, double>>();
               if (!extra_data.empty()) {
                 // Add standard properties to extra data
                 extra_data["id"] = static_cast<double>(t[0].cast<long int>());
                 extra_data["type"] = static_cast<double>(t[1].cast<int>());
                 extra_data["x"] = t[2].cast<double>();
                 extra_data["y"] = t[3].cast<double>();
                 extra_data["z"] = t[4].cast<double>();
                 extra_data["nx"] = static_cast<double>(t[5].cast<int>());
                 extra_data["ny"] = static_cast<double>(t[6].cast<int>());
                 extra_data["nz"] = static_cast<double>(t[7].cast<int>());

                 return Atom(extra_data);
               }
             }

             return Atom(t[0].cast<long int>(),
                         t[1].cast<int>(),
                         t[2].cast<double>(),
                         t[3].cast<double>(),
                         t[4].cast<double>(),
                         t[5].cast<int>(),
                         t[6].cast<int>(),
                         t[7].cast<int>());
           }),
         "Provides support for pickling")
    .def(pybind11::self == pybind11::self);

  py::class_<MoleculeIterator>(m,
                               "MoleculeIterator",
                               R"pbdoc(
       An iterator to iterate through the atoms in :class:`~pylimer_tools_cpp.Molecule`.
  )pbdoc",
                               py::module_local())
    .def("__iter__",
         [](MoleculeIterator& it) -> MoleculeIterator& { return it; })
    .def("__next__", &MoleculeIterator::next);

  py::enum_<MoleculeType>(
    m,
    "MoleculeType",
    "An enum representing the type of molecule/chain/strand.")
    .value("UNDEFINED",
           MoleculeType::UNDEFINED,
           "This value indicates that either the property was not set or not "
           "discovered.")
    .value("NETWORK_STRAND", MoleculeType::NETWORK_STRAND, R"pbdoc(
           A network strand is a strand in a network.
      )pbdoc")
    .value("PRIMARY_LOOP", MoleculeType::PRIMARY_LOOP, R"pbdoc(
           A primary loop is a network strand looping from and to the same crosslinker.
      )pbdoc")
    .value("DANGLING_CHAIN", MoleculeType::DANGLING_CHAIN, R"pbdoc(
           A dangling chain is a network strand where only one end is attached to a crosslinker.
      )pbdoc")
    .value("FREE_CHAIN", MoleculeType::FREE_CHAIN, R"pbdoc(
           A free chain is a strand not connected to any crosslinker.
      )pbdoc");

  py::class_<Molecule>(m,
                       "Molecule",
                       R"pbdoc(
       An (ideally) connected series of atoms/beads.
  )pbdoc",
                       py::module_local())
    .def(py::init<Box&, igraph_t*, MoleculeType, std::map<int, double>>(),
         R"pbdoc(
         Construct a molecule from a graph structure.
         
         :param box: The simulation box containing this molecule
         :param graph: Pointer to the igraph structure representing connectivity
         :param type: The type of molecule (see MoleculeType enum)
         :param mass_map: Map of atom types to their masses
         )pbdoc")
    // getters
    .def("get_nr_of_bonds",
         &Molecule::getNrOfBonds,
         R"pbdoc(
         Count and return the number of bonds associated with this molecule.
         
         :return: The number of bonds in this molecule
         )pbdoc")
    .def("get_nr_of_atoms",
         &Molecule::getNrOfAtoms,
         R"pbdoc(
         Count and return the number of atoms associated with this molecule.
         
         :return: The number of atoms in this molecule
         )pbdoc")
    .def("get_strand_type", &Molecule::getType, R"pbdoc(
           Get the type of this molecule (see :class:`~pylimer_tools_cpp.MoleculeType` enum).

           .. note:: 
              This type might be unset; currently, only 
              :meth:`~pylimer_tools_cpp.Universe.get_chains_with_crosslinker` assigns them automatically.
      )pbdoc")
    .def("get_strand_ends",
         &Molecule::getChainEnds,
         R"pbdoc(
          Get the ends of the given strand (= molecule).
          In case of a primary loop, the crosslink is returned, if there is one.
          Use the argument `close_loop` to decide, whether this should be returned once or twice.

          .. note:: 
               Currently only works for linear strands.
               
          :param crosslinker_type: The type of crosslinker atoms
          :param close_loop: Whether to return the crosslinker twice for loops
          :return: List of end atoms
     )pbdoc",
         py::arg("crosslinker_type") = 2,
         py::arg("close_loop") = false)
    .def("get_atoms", &Molecule::getAtoms, R"pbdoc(
            Return all atom objects enclosed in this molecule, ordered by vertex id.
            
            :return: List of atoms in vertex order
            )pbdoc")
    .def("get_atoms_lined_up",
         &Molecule::getAtomsLinedUp,
         R"pbdoc(
            Return all atom objects enclosed in this molecule based on the connectivity.

            This method works only for lone chains, atoms and loops, 
            as it throws an error if the molecule does not allow such a "line-up", 
            for example because of crosslinkers.

            Use the `crosslinker_type` parameter to force the atoms in a primary loop 
            to start with the crosslink.
            
            :param crosslinker_type: The type of crosslinker atoms
            :param assumed_coordinates: Whether to assume coordinates are valid
            :param close_loop: Whether to close loops
            :return: List of atoms in connected order
            )pbdoc",
         py::arg("crosslinker_type") = 2,
         py::arg("assumed_coordinates") = false,
         py::arg("close_loop") = false)
    .def("get_atoms_by_type",
         &Molecule::getAtomsOfType,
         R"pbdoc(
            Get the atoms with the specified type.
            
            :param type: The atom type to search for
            :return: List of atoms with the specified type
            )pbdoc",
         py::arg("type"))
    .def("get_atoms_by_degree",
         &Molecule::getAtomsOfDegree,
         R"pbdoc(
            Get the atoms that have the specified number of bonds.
            
            :param degree: The number of bonds (degree/functionality)
            :return: List of atoms with the specified degree
            )pbdoc",
         py::arg("degree"))
    .def("get_atoms_connected_to_vertex",
         &Molecule::getAtomsConnectedTo,
         R"pbdoc(
            Get the atoms connected to a specified vertex id.
            
            :param vertex_idx: The vertex index to query
            :return: List of connected atoms
            )pbdoc",
         py::arg("vertex_idx"))
    .def("get_atoms_connected_to",
         &Molecule::getConnectedAtoms,
         R"pbdoc(
            Get the atoms connected to a specified atom.

            Internally uses :meth:`~pylimer_tools_cpp.Molecule.get_atoms_connected_to`.
            
            :param atom: The atom to query connections for
            :return: List of connected atoms
            )pbdoc",
         py::arg("atom"))
    .def("get_edges", &Molecule::getEdges, R"pbdoc(
            Get all bonds. Returns a dict with three properties: 'edge_from', 'edge_to' and 'edge_type'.
            The order is not necessarily related to any structural property.
            
            .. note::
               The integer values returned refer to the vertex ids, not the atom ids.
               Use :meth:`~pylimer_tools_cpp.Molecule.get_atom_id_by_idx` to translate them to atom ids, or 
               :meth:`~pylimer_tools_cpp.Molecule.get_bonds` to have that done for you.
               
            :return: Dictionary with edge information
            )pbdoc")
    .def("get_bonds", &Molecule::getBonds, R"pbdoc(
            Get all bonds. Returns a dict with three properties: 'bond_from', 'bond_to' and 'bond_type'.
            
            :return: Dictionary with bond information
            )pbdoc")
    .def("get_atom_types",
         &Molecule::getAtomTypes,
         R"pbdoc(
         Query all types (each one for each atom) ordered by atom vertex id.
         
         :return: A vector of atom types in vertex order
         )pbdoc")
    .def("get_atom_by_vertex_idx",
         &Molecule::getAtomByVertexIdx,
         R"pbdoc(
            Get an atom for a specific vertex.
            
            :param vertex_idx: The vertex index to query
            :return: The atom at the specified vertex
            )pbdoc",
         py::arg("vertex_idx"))
    .def(
      "get_atom_by_id",
      [](const Molecule& molecule, const size_t id) {
        return molecule.getAtomByVertexIdx(molecule.getIdxByAtomId(id));
      },
      R"pbdoc(
            Get an atom by its id.
            
            :param atom_id: The atom ID to search for
            :return: The atom with the specified ID
            )pbdoc",
      py::arg("atom_id"))
    .def("get_atom_id_by_vertex_idx",
         &Molecule::getAtomIdByIdx,
         R"pbdoc(
         Get the ID of the atom by the vertex ID of the underlying graph.
         
         :param vertex_id: The vertex index in the underlying graph
         :return: The atom ID corresponding to the vertex
         )pbdoc",
         py::arg("vertex_id"))
    .def("get_vertex_idx_by_atom_id",
         &Molecule::getIdxByAtomId,
         R"pbdoc(
         Get the vertex ID of the underlying graph for an atom with a specified ID.
         
         :param atom_id: The atom ID to look up
         :return: The vertex index corresponding to the atom
         )pbdoc",
         py::arg("atom_id"))
    .def("get_key", &Molecule::getKey, R"pbdoc(
            Get a unique identifier for this molecule.
            
            :return: A unique string identifier for this molecule
            )pbdoc")
    // computations
    .def("compute_total_mass", &Molecule::computeTotalMass, R"pbdoc(
            Compute the total mass of this molecule.
            
            :return: The total mass of all atoms in this molecule
            )pbdoc")
    .def("compute_bond_lengths",
         &Molecule::computeBondLengths,
         R"pbdoc(
         Compute the length of each bond in the molecule, respecting periodic boundaries.
         
         :return: A vector of bond lengths
         )pbdoc")
    .def("compute_radius_of_gyration",
         &Molecule::computeRadiusOfGyration,
         R"pbdoc(
            Compute the radius of gyration, :math:`R_g^2` of this molecule.
            
            :math:`{R_g}^2 = \\frac{1}{M} \sum_i m_i (r_i - r_{cm})^2`,
            where :math:`M` is the total mass of the molecule, :math:`r_{cm}`
            are the coordinates of the center of mass of the molecule and the
            sum is over all contained atoms.
            
            :return: The radius of gyration squared
            )pbdoc")
    .def("compute_radius_of_gyration_with_derived_image_flags",
         &Molecule::computeRadiusOfGyrationWithDerivedImageFlags,
         R"pbdoc(
            Compute the radius of gyration, :math:`R_g^2` of this molecule,
            but ignoring the image flags attached to the atoms.
            This only works for Molecules that can be lined up with 
            :meth:`~pylimer_tools_cpp.Molecule.get_atoms_lined_up`,
            as it needs the atoms sorted such that the periodic box can still be respected somewhat.
            In other words, this function computes the radius of gyration 
            assuming the distance between two lined-up beads 
            is smaller than half the periodic box in each direction.
            
            See also: :meth:`~pylimer_tools_cpp.Molecule.compute_radius_of_gyration`.
            
            :return: The radius of gyration squared with derived image flags
            )pbdoc")
    .def("compute_end_to_end_vector",
         &Molecule::computeEndToEndVector,
         R"pbdoc(
            Compute the end-to-end vector (:math:`\overrightarrow{R}_{ee}`) of this molecule. 

            .. warning::
               Returns 0.0 if the molecule does not have two or more atoms.
               
            :return: The end-to-end vector
            )pbdoc")
    .def("compute_end_to_end_distance",
         &Molecule::computeEndToEndDistance,
         R"pbdoc(
            Compute the end-to-end distance (:math:`R_{ee}`) of this molecule. 

            .. warning::
               Returns 0.0 if the molecule does not have two or more atoms.
               
            :return: The end-to-end distance
            )pbdoc")
    .def("compute_end_to_end_vector_with_derived_image_flags",
         &Molecule::computeEndToEndVectorWithDerivedImageFlags,
         R"pbdoc(
            Compute the end-to-end vector (:math:`\overrightarrow{R}_{ee}`) of this molecule,
            but ignoring the image flags attached to the atoms. 
            This only works for Molecules that can be lined up with 
            :meth:`~pylimer_tools_cpp.Molecule.get_atoms_lined_up`,
            as it needs the atoms sorted such that the periodic box can still be respected somewhat.

            .. warning::
               Returns 0.0 if the molecule does not have two or more atoms.
               Requires bonds to be shorter than half the box length.
               
            :return: The end-to-end vector with derived image flags
            )pbdoc")
    .def("compute_end_to_end_distance_with_derived_image_flags",
         &Molecule::computeEndToEndDistanceWithDerivedImageFlags,
         R"pbdoc(
            Compute the end-to-end distance (:math:`R_{ee}`) of this molecule,
            but ignoring the image flags attached to the atoms. 
            This only works for Molecules that can be lined up with 
            :meth:`~pylimer_tools_cpp.Molecule.get_atoms_lined_up`,
            as it needs the atoms sorted such that the periodic box can still be respected somewhat.

            .. warning::
               Returns 0.0 if the molecule does not have two or more atoms.
               Requires bonds to be shorter than half the box length.
               
            :return: The end-to-end distance with derived image flags
            )pbdoc")
    .def("compute_total_vector",
         &Molecule::getOverallBondSum,
         R"pbdoc(
               Compute the sum of all bond vectors.
               
               :param crosslinker_type: The type of crosslinker atoms
               :param close_loop: Whether to close the loop for calculations
               :return: The vector sum of all bonds
            )pbdoc",
         py::arg("crosslinker_type") = 2,
         py::arg("close_loop") = true)
    .def("compute_vector_from_to",
         &Molecule::getOverallBondSumFromTo,
         R"pbdoc(
               Compute the sum of all bond vectors between two specified atoms.
               
               :param atom_id_from: ID of the starting atom
               :param atom_id_to: ID of the ending atom
               :param crosslinker_type: The type of crosslinker atoms
               :param require_order: Whether to require specific ordering
               :return: The vector sum from start to end atom
            )pbdoc",
         py::arg("atom_id_from"),
         py::arg("atom_id_to"),
         py::arg("crosslinker_type") = 2,
         py::arg("require_order") = true)
    .def("compute_total_length", &Molecule::computeTotalLength, R"pbdoc(
         Compute the sum of the lengths of all bonds.
         In most cases, this is equal to the contour length.
         
         :return: The total contour length of the molecule
         )pbdoc")
    .def("get_edge_ids_from",
         &Molecule::getIncidentEdgeIds,
         R"pbdoc(
         Get the edge IDs incident to a specific vertex.
         
         :param vertex_id: The vertex to query
         :return: A vector of edge IDs connected to the vertex
         )pbdoc",
         py::arg("vertex_id"))
    .def("get_edge_ids_from_to",
         &Molecule::getEdgeIdsFromTo,
         R"pbdoc(
         Get the edge IDs of the edges between two specific vertices.
         
         :param vertex_id_from: Starting vertex ID
         :param vertex_id_to: Ending vertex ID
         :return: Vector of edge IDs between the specified vertices
         )pbdoc",
         py::arg("vertex_id_from"),
         py::arg("vertex_id_to"))
    .def("get_nr_of_edges_from_to",
         &Molecule::getPathLength,
         R"pbdoc(
         Get the number of edges in the shortest path between two specific vertices.

         If `max_length` is provided and positive, it will only consider paths up to that length.
         
         :param vertex_id_from: Starting vertex ID
         :param vertex_id_to: Ending vertex ID  
         :param max_length: Maximum path length to consider (-1 for no limit)
         :return: Number of edges in the shortest path, or -1 if no path exists
         )pbdoc",
         py::arg("vertex_id_from"),
         py::arg("vertex_id_to"),
         py::arg("max_length") = -1)
    // operators
    .def(pybind11::self == pybind11::self)
    .def(
      "__getitem__",
      [](const Molecule& molecule, const size_t index) {
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
    .def(
      "__copy__",
      [](const Molecule& molecule) { return Molecule(molecule); },
      R"pbdoc(
         Create a copy of this molecule.
         
         :return: A new Molecule instance that is a copy of this one
         )pbdoc");

  py::class_<NeighbourList>(m,
                            "NeighbourList",
                            R"pbdoc(
    Gives access to somewhat fast queries on the neighbourhood of atoms.
    
    This class provides efficient spatial queries for finding atoms within
    a specified distance of each other.
    )pbdoc",
                            py::module_local())
    .def(py::init<const std::vector<pylimer_tools::entities::Atom>&,
                  const pylimer_tools::entities::Box&,
                  double>(),
         R"pbdoc(
         Instantiate a new neighbour list.
         
         :param atoms: Vector of atoms to include in the neighbour list
         :param box: The simulation box
         :param cutoff: Maximum distance for neighbour searches
         )pbdoc",
         py::arg("atoms"),
         py::arg("box"),
         py::arg("cutoff"))
    .def("get_atoms_close_to",
         py::overload_cast<const Atom&, double, double, bool, bool>(
           &NeighbourList::getAtomsCloseTo),
         R"pbdoc(
          List all atoms that are close to a given one. 

          It is possible to request it within a new cutoff, 
          though the underlying neighbour list will not be regenerated.
          For performance reasons, it is recommended to initialize a 
          new NeighbourList if you require a different cutoff, depending on your use case.

          You can use a negative value for the upper_cutoff to use the cutoff used for 
          filling the neighbour list buckets.
          
          :param atom: The reference atom
          :param upper_cutoff: Maximum distance for neighbours
          :param lower_cutoff: Minimum distance for neighbours
          :param unwrapped: Whether to use unwrapped coordinates
          :param expect_self: Whether to expect the atom itself in results
          :return: List of neighbouring atoms
         )pbdoc",
         py::arg("atom"),
         py::arg("upper_cutoff") = 1.0,
         py::arg("lower_cutoff") = 0.0,
         py::arg("unwrapped") = false,
         py::arg("expect_self") = false)
    .def("remove_atom",
         &NeighbourList::removeAtom,
         R"pbdoc(
         Remove an atom from this neighbour list.
         It will not show up when querying for neighbours, 
         but its neighbours cannot be queried either.
         
         :param atom: The atom to remove
         :param debug_hint: Optional debug information
         )pbdoc",
         py::arg("atom"),
         py::arg("debug_hint") = "");

  py::class_<Universe>(m,
                       "Universe",
                       R"pbdoc(
    Represents a full Polymer Network structure, a collection of molecules.
    
    This is the main class for representing molecular systems, containing
    atoms, bonds, angles, and the simulation box.
    )pbdoc",
                       py::module_local())
    .def(py::init<const double, const double, const double>(),
         R"pbdoc(
         Instantiate this Universe (Collection of Molecules) providing the box lengths.
         
         :param Lx: Box length in x direction
         :param Ly: Box length in y direction
         :param Lz: Box length in z direction
         )pbdoc",
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
         R"pbdoc(
         Add atoms to the Universe, vertices to the underlying graph.
         
         :param ids: Vector of atom IDs
         :param types: Vector of atom types
         :param x: Vector of x coordinates
         :param y: Vector of y coordinates
         :param z: Vector of z coordinates
         :param nx: Vector of periodic image flags in x direction
         :param ny: Vector of periodic image flags in y direction
         :param nz: Vector of periodic image flags in z direction
         )pbdoc",
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
          Remove atoms and all associated bonds by their atom IDs. 
          
          :param atom_ids: Vector of atom IDs to remove
          )pbdoc",
         py::arg("atom_ids"))
    .def("replace_atom",
         &Universe::replaceAtom,
         R"pbdoc(
          Replace the properties of an atom with the properties of another given atom.
          
          :param atom_id: ID of the atom to replace
          :param replacement_atom: The atom with new properties
          )pbdoc",
         py::arg("atom_id"),
         py::arg("replacement_atom"))
    .def("replace_atom_type",
         &Universe::replaceAtomType,
         R"pbdoc(
          Replace the type of an atom with another type.
          
          :param atom_id: ID of the atom to modify
          :param new_type: The new atom type
          )pbdoc",
         py::arg("atom_id"),
         py::arg("new_type"))
    .def("resample_velocities",
         &Universe::resampleVelocities,
         R"pbdoc(
         Resample velocities for atoms in the universe.
         
         :param mean: Mean velocity
         :param variance: Velocity variance
         :param seed: Random seed for velocity generation
         :param is_2d: Whether to omit sampling the z direction
         )pbdoc",
         py::arg("mean"),
         py::arg("variance"),
         py::arg("seed") = "",
         py::arg("is_2d") = false)
    .def("add_bonds",
         py::overload_cast<const std::vector<long int>&,
                           const std::vector<long int>&>(&Universe::addBonds),
         R"pbdoc(
         Add bonds to the underlying atoms, edges to the underlying graph. 
         If the connected atoms are not found, the bonds are silently skipped.
         
         :param bonds_from: Vector of atom IDs for bond start points
         :param bonds_to: Vector of atom IDs for bond end points
         )pbdoc",
         py::arg("bonds_from"),
         py::arg("bonds_to"))
    .def("add_bonds",
         py::overload_cast<const size_t,
                           const std::vector<long int>&,
                           const std::vector<long int>&,
                           const std::vector<int>&,
                           const bool,
                           const bool>(&Universe::addBonds),
         R"pbdoc(
         Add bonds to the underlying atoms, edges to the underlying graph.
         
         :param nr_of_bonds: Number of bonds to add
         :param bonds_from: Vector of atom IDs for bond start points
         :param bonds_to: Vector of atom IDs for bond end points
         :param bond_types: Vector of bond types
         :param ignore_non_existent_atoms: Whether to skip bonds to non-existent atoms
         :param simplify_universe: Whether to simplify the universe after adding bonds
         )pbdoc",
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
         R"pbdoc(
         Add bonds to the underlying atoms, edges to the underlying graph. 
         If the connected atoms are not found, the bonds are silently skipped.
         
         :param bonds_from: Vector of atom IDs for bond start points
         :param bonds_to: Vector of atom IDs for bond end points
         :param bond_types: Vector of bond types
         )pbdoc",
         py::arg("bonds_from"),
         py::arg("bonds_to"),
         py::arg("bond_types"))
    .def("remove_bonds",
         &Universe::removeBonds,
         R"pbdoc(
          Remove bonds by their connected atom IDs. 
          
          :param bonds_from: Vector of starting atom IDs
          :param bonds_to: Vector of ending atom IDs
          )pbdoc",
         py::arg("bonds_from"),
         py::arg("bonds_to"))
    .def("remove_bonds_by_type",
         &Universe::removeBondsOfType,
         R"pbdoc(
          Remove bonds with a specific type. 
          
          :param bond_type: The bond type to remove
          )pbdoc",
         py::arg("bond_type"))
    .def("add_angles",
         &Universe::addAngles,
         R"pbdoc(
         Add angles to the Universe. No relation to the underlying graph, 
         just a method to preserve read & write capabilities.
         
         :param angles_from: Vector of atom IDs for angle start points
         :param angles_via: Vector of atom IDs for angle middle points
         :param angles_to: Vector of atom IDs for angle end points
         :param angle_types: Vector of angle types
         )pbdoc",
         py::arg("angles_from"),
         py::arg("angles_via"),
         py::arg("angles_to"),
         py::arg("angle_types"))
    .def("add_dihedral_angles",
         &Universe::addDihedralAngles,
         R"pbdoc(
         Add dihedral angles to the Universe. No relation to the underlying graph, 
         just a method to preserve read & write capabilities.
         
         :param angles_from: Vector of atom IDs for dihedral start points
         :param angles_via1: Vector of atom IDs for first middle points
         :param angles_via2: Vector of atom IDs for second middle points
         :param angles_to: Vector of atom IDs for dihedral end points
         :param angle_types: Vector of dihedral angle types
         )pbdoc",
         py::arg("angles_from"),
         py::arg("angles_via1"),
         py::arg("angles_via2"),
         py::arg("angles_to"),
         py::arg("angle_types"))
    .def("remove_all_angles", &Universe::removeAllAngles, R"pbdoc(
          Remove all angles from the Universe. 
          This will not remove the atoms or bonds, just the angles.
          
          .. warning::
            This will not remove dihedral angles, use :meth:`~pylimer_tools_cpp.Universe.remove_all_dihedral_angles` for that.
     )pbdoc")
    .def(
      "remove_all_dihedral_angles", &Universe::removeAllDihedralAngles, R"pbdoc(
          Remove all dihedral angles from the Universe. 
          This will not remove the atoms or bonds, just the stored dihedral angles.
     )pbdoc")
    .def("hash_angle_type",
         &Universe::hashAngleType,
         R"pbdoc(
          Convert the three integers to one long number/hash.
          Used internally for duplicate detection.
          
          :param angle_from: First atom ID in the angle
          :param angle_via: Middle atom ID in the angle
          :param angle_to: Last atom ID in the angle
          :return: Hash value for the angle type
     )pbdoc",
         py::arg("angle_from"),
         py::arg("angle_via"),
         py::arg("angle_to"))
    .def("hash_dihedral_angle_type",
         &Universe::hashDihedralAngleType,
         R"pbdoc(
          Convert the four integers to one long number/hash.
          Used internally for duplicate detection.
          
          :param angle_from: First atom ID in the dihedral
          :param angle_via1: Second atom ID in the dihedral
          :param angle_via2: Third atom ID in the dihedral
          :param angle_to: Fourth atom ID in the dihedral
          :return: Hash value for the dihedral angle type
     )pbdoc",
         py::arg("angle_from"),
         py::arg("angle_via1"),
         py::arg("angle_via2"),
         py::arg("angle_to"))
    .def("set_masses",
         &Universe::setMasses,
         R"pbdoc(
         Set the mass per type of atom.
         
         :param mass_per_type: Map of atom types to their masses
         )pbdoc",
         py::arg("mass_per_type"))
    .def("set_mass",
         &Universe::setMassForType,
         R"pbdoc(
         Set the mass for a specific atom type.
         
         :param atom_type: The atom type to set mass for
         :param mass: The mass value to assign
         )pbdoc",
         py::arg("atom_type"),
         py::arg("mass"))
    .def("set_timestep",
         &Universe::setTimestep,
         R"pbdoc(
         Set the timestep when this Universe was captured.
         
         :param timestep: The timestep value
         )pbdoc",
         py::arg("timestep"))
    .def("set_box_lengths",
         &Universe::setBoxLengths,
         R"pbdoc(
          Override the currently assigned box with one with the side lengths specified.
          
          :param lx: Length in x direction
          :param ly: Length in y direction
          :param lz: Length in z direction
          :param rescale_atoms: Whether to rescale atom positions
          )pbdoc",
         py::arg("lx"),
         py::arg("ly"),
         py::arg("lz"),
         py::arg("rescale_atoms") = false)
    .def("set_box",
         &Universe::setBox,
         R"pbdoc(
          Override the currently assigned box with the one specified.
          
          :param box: The new box to assign
          :param rescale_atoms: Whether to rescale atom positions
          )pbdoc",
         py::arg("box"),
         py::arg("rescale_atoms") = false)
    .def("set_vertex_property",
         &Universe::setPropertyValue<double>,
         R"pbdoc(
          Set a specific property for a specific vertex.
          
          :param vertex_id: The vertex ID to modify
          :param property_name: Name of the property to set
          :param value: Value to assign to the property
          )pbdoc",
         py::arg("vertex_id"),
         py::arg("property_name"),
         py::arg("value"))
    // getters
    .def("get_clusters", &Universe::getClusters, R"pbdoc(
          Get the components of the universe that are not connected to each other.
          Returns a list of :class:`~pylimer_tools_cpp.Universe` objects.
          Unconnected, free atoms/beads become their own :class:`~pylimer_tools_cpp.Universe`.
          
          :return: List of disconnected Universe components
          )pbdoc")
    .def("get_molecules",
         &Universe::getMolecules,
         R"pbdoc(
          Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms.
          
          Reduces the Universe to a list of molecules. 
          Specify the crosslinker_type to an existing type ID, 
          then those atoms will be omitted, and this function returns chains instead.

          :param atom_type_to_omit: The type of atom to omit from the universe to end up with the desired molecules (e.g., the type of the crosslinkers).
          :return: List of molecules
          )pbdoc",
         py::arg("atom_type_to_omit"))
    .def("get_atoms_connected_to_vertex",
         &Universe::getAtomsConnectedTo,
         R"pbdoc(
            Get the atoms connected to a specified vertex ID.
            
            :param vertex_idx: The vertex index to query
            :return: List of connected atoms
            )pbdoc",
         py::arg("vertex_idx"))
    .def("get_atoms_connected_to",
         &Universe::getConnectedAtoms,
         R"pbdoc(
            Get the atoms connected to a specified atom.

            Internally uses :meth:`~pylimer_tools_cpp.Universe.get_atoms_connected_to`.
            
            :param atom: The atom to query connections for
            :return: List of connected atoms
            )pbdoc",
         py::arg("atom"))
    .def("get_atoms_by_degree",
         &Universe::getAtomsOfDegree,
         R"pbdoc(
            Get the atoms that have the specified number of bonds.
            )pbdoc",
         py::arg("functionality"))
    .def("get_vertex_degrees", &Universe::getVertexDegrees, R"pbdoc(
          Get the degree (functionality) of each vertex.
     )pbdoc")
    .def("count_loop_lengths",
         &Universe::countLoopLengths,
         R"pbdoc(
          Find all loops (below a specific length) and count the number of atoms involved in them.
          Returns the count, how many loops per length are found.
         )pbdoc",
         py::arg("max_length") = -1)
    .def("find_loops",
         &Universe::findLoopsOfAtoms,
         R"pbdoc(
            Decompose the Universe into loops.
            The primary index specifies the degree of the loop.

            CAUTION:
               There are exponentially many paths between two crosslinkers of a network,
               and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
               You can use the `max_length` parameter to restrict the algorithm to only search for loops up to a certain length.
               Use a negative value to find all loops and paths.
            )pbdoc",
         py::arg("crosslinker_type"),
         py::arg("max_length") = -1,
         py::arg("skip_self_loops") = false)
    .def("find_minimal_order_loop_from",
         &Universe::findMinimalOrderLoopFrom,
         R"pbdoc(
            Find the loops in the network starting with one connection.

            .. warning::
               There are exponentially many paths between two crosslinkers of a network,
               and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
               You can use the `max_length` parameter to restrict the algorithm to only search for loops up to a certain length.
               Use a negative value to find all loops and paths.

            :param loop_start: The atom ID to start the search from
            :param loop_step1: The first step to take, i.e., the first bond to follow
            :param max_length: The maximum length of the loop to find, or -1 for no limit
            :param skip_self_loops: Whether to skip self-loops (i.e., loops that start and end at the same atom; only relevant if `loop_start` is equal to `loop_step1`)
            :return: List of loops found
            )pbdoc",
         py::arg("loop_start"),
         py::arg("loop_step1"),
         py::arg("max_length") = -1,
         py::arg("skip_self_loops") = false)
    .def("count_atom_types",
         &Universe::countAtomTypes,
         R"pbdoc(
          Count how often each atom type is present.
          
          :return: Dictionary mapping atom types to their counts
     )pbdoc")
    .def("count_atoms_in_skin_distance",
         &Universe::countAtomsInSkinDistance,
         R"pbdoc(
          This is a function that may help you to compute the radial distribution function.
          It loops through all atoms and counts neighbors within specified distance ranges.

          :param distances: The edges of the bins for distance counting
          :param unwrapped: Whether to measure the distance in unwrapped coordinates or as PBC-corrected distance
          :return: Array of counts for each distance bin
     )pbdoc",
         py::arg("distances"),
         py::arg("unwrapped") = false)
    .def("get_chains_with_crosslinker",
         &Universe::getChainsWithCrosslinker,
         R"pbdoc(
            Decompose the Universe into strands (molecules, which could be either chains, or even lonely atoms), without omitting the crosslinkers
            (as in :meth:`~pylimer_tools_cpp.Universe.get_molecules`).
            In turn, e.g. for a tetrafunctional crosslinker, it will be 4 times in the resulting molecules.
            
            .. note::
               Crosslinkers without bonds to non-crosslinkers are not returned 
               (i.e., single crosslinkers, are not counted as strands).
               
            :param crosslinker_type: The type of crosslinker atoms
            :return: List of chains including crosslinkers
           )pbdoc",
         py::arg("crosslinker_type"))
    .def("get_network_of_crosslinker",
         &Universe::getNetworkOfCrosslinker,
         R"pbdoc(
            Reduce the network to contain only crosslinkers, replacing all the strands with a single bond.
            Useful e.g. to reduce the memory usage and runtime of 
            :meth:`~pylimer_tools_cpp.Universe.find_loops` or 
            :meth:`~pylimer_tools_cpp.Universe.has_infinite_strand`.
            
            Further use :meth:`~pylimer_tools_cpp.Universe.simplify` to remove primary loops.
            
            :param crosslinker_type: The type of crosslinker atoms
            :return: Reduced network containing only crosslinkers
          )pbdoc",
         py::arg("crosslinker_type"))
    .def("contract_vertices_along_bond_type",
         &Universe::contractVerticesAlongBondType,
         R"pbdoc(
          Merge vertices along a specific bond type.

          May result in new self-loops; use :meth:`~pylimer_tools_cpp.Universe.simplify` to remove them.
          
          :param bond_type: The bond type to contract along
         )pbdoc",
         py::arg("bond_type"))
    .def("get_atom_types",
         &Universe::getAtomTypes,
         R"pbdoc(
          Get all types (each one for each atom) ordered by atom vertex ID.
          
          :return: Vector of atom types in vertex order
          )pbdoc")
    .def("get_atom",
         &Universe::getAtom,
         R"pbdoc(
         Find an atom by its ID.
         
         :param atom_id: The ID of the atom to find
         :return: The atom with the specified ID
         )pbdoc",
         py::arg("atom_id"))
    .def("get_atom_by_vertex_id",
         &Molecule::getAtomByVertexIdx,
         R"pbdoc(
      Find an atom by the ID of the vertex of the underlying graph.
      
      :param vertex_id: The vertex ID to query
      :return: The atom at the specified vertex
      )pbdoc",
         py::arg("vertex_id"))
    .def("get_atoms", &Universe::getAtoms, R"pbdoc(
            Get all atoms in the universe.
            
            :return: List of all atoms
            )pbdoc")
    .def("get_atoms_by_type",
         &Universe::getAtomsOfType,
         R"pbdoc(
            Query all atoms by their type.
            )pbdoc",
         py::arg("atom_type"))
    .def("get_atom_id_by_vertex_idx",
         &Universe::getAtomIdByIdx,
         R"pbdoc(
         Get the ID of the atom by the vertex ID of the underlying graph.
         
         :param vertex_id: The vertex index in the underlying graph
         :return: The atom ID corresponding to the vertex
         )pbdoc",
         py::arg("vertex_id"))
    .def("get_vertex_idx_by_atom_id",
         &Universe::getIdxByAtomId,
         R"pbdoc(
         Get the vertex ID of the underlying graph for an atom with a specified ID.
         
         :param atom_id: The atom ID to look up
         :return: The vertex index corresponding to the atom
         )pbdoc",
         py::arg("atom_id"))
    .def("get_edges", &Universe::getEdges, R"pbdoc(
            Get all edges. Returns a dict with three properties: 'edge_from', 'edge_to' and 'edge_type'.
            The order is not necessarily related to any structural characteristic.
            
            .. note::
               The integer values returned refer to the vertex IDs, not the atom IDs.
               Use :meth:`~pylimer_tools_cpp.Universe.get_atom_id_by_idx` to translate them to atom IDs, or
               :meth:`~pylimer_tools_cpp.Universe.get_bonds` to have that done for you.
               
            :return: Dictionary with edge information
            )pbdoc")
    .def("interpolate_edges",
         &Universe::interpolateEdges,
         R"pbdoc(
          Get more or less edges than currently present, interpolating between junctions.
          
          :param crosslinker_type: The type of crosslinker atoms
          :param interpolation_factor: Factor for edge interpolation
          :return: Interpolated edge structure
         )pbdoc",
         py::arg("crosslinker_type"),
         py::arg("interpolation_factor"))
    .def("get_bonds", &Universe::getBonds, R"pbdoc(
            Get all bonds. Returns a dict with three properties: 'bond_from', 'bond_to' and 'bond_type'.
            The order is not necessarily related to any structural characteristic.
            
            :return: Dictionary with bond information
            )pbdoc")
    .def("get_angles", &Universe::getAngles, R"pbdoc(
           Get all angles added to this network.

           Returns a dict with three properties: 'angle_from', 'angle_via' and 'angle_to'.

           .. note::
               The integer values returned refer to the atom IDs, not the vertex IDs.
               Use :meth:`~pylimer_tools_cpp.Universe.get_idx_by_atom_id` to translate them to vertex IDs.
               
           :return: Dictionary with angle information
           )pbdoc")
    .def("get_box", &Universe::getBox, R"pbdoc(
            Get the underlying bounding box object.
            
            :return: The simulation box
            )pbdoc")
    .def("get_masses", &Universe::getMasses, R"pbdoc(
            Get the mass of one atom per type.
            
            :return: Dictionary mapping atom types to masses
            )pbdoc")
    .def("get_volume", &Universe::getVolume, R"pbdoc(
            Query the volume of the underlying bounding box.
            
            :return: The volume of the simulation box
            )pbdoc")
    .def("get_nr_of_atoms", &Universe::getNrOfAtoms, R"pbdoc(
            Query the number of atoms in this universe.
            
            :return: Number of atoms
            )pbdoc")
    .def("get_nr_of_bonds", &Universe::getNrOfBonds, R"pbdoc(
            Query the number of bonds associated with this universe.
            )pbdoc")
    .def("get_nr_of_angles", &Universe::getNrOfAngles, R"pbdoc(
            Query the number of angles that have been added to this universe.
            )pbdoc")
    .def("get_nr_of_dihedral_angles", &Universe::getNrOfDihedralAngles, R"pbdoc(
            Query the number of dihedral angles that have been added to this universe.
            )pbdoc")
    .def("get_timestep", &Universe::getTimestep, R"pbdoc(
            Query the timestep when this universe was captured.
            )pbdoc")
    .def("get_nr_of_bonds_of_atom",
         &Universe::computeFunctionalityForAtom,
         R"pbdoc(
       Count the number of immediate neighbors of an atom, specified by its ID.

       :param atom_id: The ID of the atom to query
       :return: Number of bonds connected to the atom
      )pbdoc")
    .def("get_nr_of_bonds_of_vertex",
         &Universe::computeFunctionalityForVertex,
         R"pbdoc(
       Count the number of immediate neighbors of an atom, specified by its vertex ID.

       :param vertex_id: The vertex ID of the atom to query
       :return: Number of bonds connected to the vertex
      )pbdoc")
    .def("get_edge_ids_from",
         &Universe::getIncidentEdgeIds,
         R"pbdoc(
          Get the edge IDs incident to a specific vertex.

          :param vertex_id: The ID of the vertex to query
          :return: List of edge IDs connected to the vertex
         )pbdoc",
         py::arg("vertex_id"))
    .def("get_edge_ids_from_to",
         &Universe::getEdgeIdsFromTo,
         R"pbdoc(
          Get the edge IDs of all edges between two specific vertices.

          :param vertex_id_from: The starting vertex ID
          :param vertex_id_to: The ending vertex ID
          :return: List of edge IDs connecting the two vertices
         )pbdoc",
         py::arg("vertex_id_from"),
         py::arg("vertex_id_to"))
    .def("get_nr_of_edges_from_to",
         &Universe::getPathLength,
         R"pbdoc(
          Get the number of edges in the shortest path between two specific vertices.

          If `max_length` is provided and positive, it will only consider paths up to that length.

          :param vertex_id_from: The starting vertex ID
          :param vertex_id_to: The ending vertex ID  
          :param max_length: Maximum path length to consider (default: -1 for no limit)
          :return: Number of edges in shortest path, or -1 if no path exists
         )pbdoc",
         py::arg("vertex_id_from"),
         py::arg("vertex_id_to"),
         py::arg("max_length") = -1)
    // computations
    .def("compute_bond_vectors",
         &Universe::computeBondVectors,
         R"pbdoc(
         Compute the vectors of each bond in the molecule, respecting periodic boundaries.
         
         :return: A list of bond vectors
         )pbdoc")
    .def("compute_bond_lengths",
         &Universe::computeBondLengths,
         R"pbdoc(
         Compute the length of each bond in the molecule, respecting periodic boundaries.
         
         :return: A list of bond lengths
         )pbdoc")
    .def("compute_angles",
         &Universe::computeAngles,
         R"pbdoc(
         Compute the angle of each angle in the molecule, respecting periodic boundaries.
         
         :return: A list of angle values in radians
         )pbdoc")
    .def("detect_angles",
         &Universe::detectAngles,
         R"pbdoc(
          Detect angles in the network based on the current bonds.
          Return the result in the same format as :meth:`~pylimer_tools_cpp.Universe.get_angles`, 
          but all angles that are detected in the network, rather than the ones already set.
          Note that the angle types are determined by 
          :meth:`~pylimer_tools_cpp.Universe.hash_angle_type`,
          which serves angle types that should be mapped by you back to smaller numbers, 
          before serving them again to :meth:`~pylimer_tools_cpp.Universe.add_angles`,
          if you want to have them written e.g. for LAMMPS.
          
          :return: Dictionary with detected angle information
         )pbdoc")
    .def("detect_dihedral_angles",
         &Universe::detectDihedralAngles,
         R"pbdoc(
          Detect dihedral angles in the network based on the current bonds.
          Return the result in the same format as :meth:`~pylimer_tools_cpp.Universe.get_dihedral_angles`, 
          but all dihedral angles that are detected in the network, rather than the ones already set.
          Note that the angle types are determined by 
          :meth:`~pylimer_tools_cpp.Universe.hash_dihedral_angle_type`,
          which serves angle types that should be mapped by you back to smaller numbers, 
          before serving them to :meth:`~pylimer_tools_cpp.Universe.add_dihedral_angles`,
          if you want to have them written e.g. for LAMMPS.

          :return: Dictionary with detected dihedral angle information
         )pbdoc")
    .def("has_infinite_strand",
         &Universe::hasInfiniteStrand,
         R"pbdoc(
           Check whether there is a strand (with crosslinker) in the universe that loops through periodic images without coming back.
           
            .. warning::
               There are exponentially many paths between two crosslinkers of a network,
               and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
               
           :return: True if infinite strands are detected, False otherwise
           )pbdoc")
    .def("determine_functionality_per_type",
         &Universe::determineFunctionalityPerType,
         R"pbdoc(
            Find the maximum functionality of each atom type in the network.
            
            :return: Dictionary mapping atom types to maximum functionality
            )pbdoc")
    .def("determine_effective_functionality_per_type",
         &Universe::determineEffectiveFunctionalityPerType,
         R"pbdoc(
            Find the average functionality of each atom type in the network.
            
            :return: Dictionary mapping atom types to average functionality
            )pbdoc")
    .def("compute_mean_strand_length",
         &Universe::getMeanStrandLength,
         R"pbdoc(
              Compute the mean number of beads per strand.
              
              :param crosslinker_type: The type of crosslinker atoms
              :return: Mean strand length
              )pbdoc",
         py::arg("crosslinker_type"))
    .def("compute_total_mass", &Universe::computeTotalMass, R"pbdoc(
          Compute the total mass of this network/universe in whatever mass unit was used when 
          :meth:`~pylimer_tools_cpp.Universe.set_masses` was called.
          
          :return: Total mass of the universe
     )pbdoc")
    .def("compute_number_average_molecular_weight",
         &Universe::computeNumberAverageMolecularWeight,
         R"pbdoc(
              Compute the number average molecular weight.

              .. note:: 
                    Crosslinkers are ignored completely.
                    
              :param crosslinker_type: The type of crosslinker atoms
              :return: Number average molecular weight
              )pbdoc",
         py::arg("crosslinker_type"))
    .def("compute_weight_average_molecular_weight",
         &Universe::computeWeightAverageMolecularWeight,
         R"pbdoc(
              Compute the weight average molecular weight.

              .. note:: 
                    Crosslinkers are ignored completely.
                    
              :param crosslinker_type: The type of crosslinker atoms
              :return: Weight average molecular weight
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

            If no masses are stored, assumes a mass of 1 for each atom.

            If the total mass is 0., returns the total mass per atom type.
            )pbdoc")
    .def("compute_end_to_end_distances",
         &Universe::computeEndToEndDistances,
         R"pbdoc(
          Compute the end-to-end distance of each strand in the network.

          .. note::
               Internally, this uses either :meth:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance` 
               or :meth:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags`, 
               depending on `derive_image_flags`.
               Invalid strands (where said function returns 0.0 or -1.0) are ignored.
               
          :param crosslinker_type: The type of crosslinker atoms
          :param derive_image_flags: Whether to derive image flags from connectivity
          :return: List of end-to-end distances
     )pbdoc",
         py::arg("crosslinker_type"),
         py::arg("derive_image_flags") = false)
    .def("compute_mean_end_to_end_distance",
         &Universe::computeMeanEndToEndDistance,
         R"pbdoc(
          Compute the mean of the end-to-end distances of each strand in the network.

          .. note::
               Internally, this uses either :meth:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance` 
               or :meth:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags`, 
               depending on `derive_image_flags`.
               Invalid strands (where said function returns 0.0 or -1.0) are ignored.
               
          :param crosslinker_type: The type of crosslinker atoms
          :param derive_image_flags: Whether to derive image flags from connectivity
          :return: Mean end-to-end distance
     )pbdoc",
         py::arg("crosslinker_type"),
         py::arg("derive_image_flags") = false)
    .def("compute_mean_squared_end_to_end_distance",
         &Universe::computeMeanSquareEndToEndDistance,
         R"pbdoc(
          Compute the mean square of the end-to-end distances of each strand (incl. crosslinks) in the network.

          .. note::
               Internally, this uses either :meth:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance` 
               or :meth:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags`, 
               depending on `derive_image_flags`.
               Invalid strands (where said function returns 0.0 or -1.0) are ignored.
               
          :param crosslinker_type: The type of crosslinker atoms
          :param only_those_with_two_crosslinkers: Whether to only consider strands with two crosslinkers
          :param derive_image_flags: Whether to derive image flags from connectivity
          :return: Mean squared end-to-end distance
     )pbdoc",
         py::arg("crosslinker_type"),
         py::arg("only_those_with_two_crosslinkers") = false,
         py::arg("derive_image_flags") = false)
    .def("compute_dxs",
         &Universe::computeDxs,
         R"pbdoc(
         Compute the dx distance for certain bonds (length in x direction).
         
         :param atom_ids_to: Vector of destination atom IDs
         :param atom_ids_from: Vector of source atom IDs
         :return: Vector of x-direction distances
         )pbdoc",
         py::arg("atom_ids_to"),
         py::arg("atom_ids_from"))
    .def("compute_dys",
         &Universe::computeDys,
         R"pbdoc(
         Compute the dy distance for certain bonds (length in y direction).
         
         :param atom_ids_to: Vector of destination atom IDs
         :param atom_ids_from: Vector of source atom IDs
         :return: Vector of y-direction distances
         )pbdoc",
         py::arg("atom_ids_to"),
         py::arg("atom_ids_from"))
    .def("compute_dzs",
         &Universe::computeDzs,
         R"pbdoc(
         Compute the dz distance for certain bonds (length in z direction).
         
         :param atom_ids_to: Vector of destination atom IDs
         :param atom_ids_from: Vector of source atom IDs
         :return: Vector of z-direction distances
         )pbdoc",
         py::arg("atom_ids_to"),
         py::arg("atom_ids_from"))
    .def("compute_temperature",
         &Universe::computeTemperature,
         R"pbdoc(
      Use the velocities per atom to compute the temperature from the kinetic energy of the system.
      
      :param dimensions: Number of dimensions (typically 3)
      :param k_b: Boltzmann constant value
      :return: Computed temperature
      )pbdoc",
         py::arg("dimensions") = 3,
         py::arg("k_b") = 1.)
    .def("simplify",
         &Universe::simplify,
         R"pbdoc(
         Remove self links and double bonds. This function is called 
         automatically after adding bonds.
         
         This operation cleans up the graph structure by removing
         redundant connections.
         )pbdoc")
    //
    .def("has_atom_with_id",
         &Universe::containsAtomWithId,
         R"pbdoc(
         Check whether this universe contains an atom with the specified ID.

         :param atom_id: The atom ID to check
         :return: True if the atom exists in this universe, False otherwise
         )pbdoc",
         py::arg("atom_id"))
    // operators
    //     .def(pybind11::self == pybind11::self)
    .def(
      "__getitem__",
      [](const Universe& u, const size_t index) {
        if (index > u.getNrOfAtoms()) {
          throw py::index_error();
        }
        return u.getAtom(u.getAtomIdByIdx(index));
      },
      R"pbdoc(
       Access an atom by its vertex index.
  )pbdoc")
    .def("__contains__", &Universe::containsAtom, R"pbdoc(
          Check whether a particular atom is contained in this universe.
     )pbdoc")
    .def("__len__", &Universe::getNrOfAtoms, R"pbdoc(
       Get the number of atoms in this universe.
  )pbdoc")
#ifdef CEREALIZABLE
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
#endif
    .def(
      "__copy__",
      [](const Universe& universe) { return Universe(universe); },
      R"pbdoc(
         Create a copy of this universe.
         
         :return: A new Universe instance that is a copy of this one
         )pbdoc");

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

  py::class_<LazyUniverseSequenceIterator>(m,
                                           "LazyUniverseSequenceIterator",
                                           R"pbdoc(
       An iterator to iterate throught the universes in :obj:`~pylimer_tools_cpp.UniverseSequence`.
  )pbdoc",
                                           py::module_local())
    .def("__iter__",
         [](const LazyUniverseSequenceIterator& it)
           -> const LazyUniverseSequenceIterator& { return it; })
    .def("__next__", &LazyUniverseSequenceIterator::next);

  py::class_<UniverseSequence>(m,
                               "UniverseSequence",
                               R"pbdoc(
     This class represents a sequence of Universes, with the Universe's data
     only being read on request. Dump files are read at once in order
     to know how many timesteps/universes are available in total 
     (but the universes' data is not read on first look through the file).
     This, while it can lead to two (or more) reads of the whole file, 
     is a measure in order to enable low memory useage if needed (i.e. for large dump files).
     Use Python's iterator to have this UniverseSequence only ever retain one universe in memory.
     Alternatively, use :meth:`~pylimer_tools_cpp.UniverseSequence.forget_at_index`
     to have the UniverseSequence forget about already read universes.
     )pbdoc",
                               py::module_local())
    .def(py::init<>(),
         R"pbdoc(
         Construct an empty UniverseSequence.
         
         Use initialization methods to populate it with data.
         )pbdoc")
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
         R"pbdoc(
         Get the Universe that's next in the sequence.
         
         :return: The next Universe in the sequence
         )pbdoc")
    .def("at_index",
         &UniverseSequence::atIndex,
         R"pbdoc(
         Get the Universe at the given index (as of in the sequence given 
         by the dump file).
         
         :param index: The index of the universe to retrieve
         :return: The Universe at the specified index
         )pbdoc",
         py::arg("index"))
    .def("forget_at_index",
         &UniverseSequence::forgetAtIndex,
         R"pbdoc(
      Clear the memory of the Universe at the given index (as of in the
      sequence given by the dump file).
     )pbdoc",
         py::arg("index"))
    .def("reset_iterator",
         &UniverseSequence::resetIterator,
         R"pbdoc(
          Reset the internal iterator, such that a subsequent call to 
          :meth:`~pylimer_tools_cpp.UniverseSequence.next` returns the first one again.
          )pbdoc")
    .def("get_length", &UniverseSequence::getLength, R"pbdoc(
            Get the number of universes in this sequence.
            )pbdoc")
    .def("get_all", &UniverseSequence::getAll, R"pbdoc(
            Get all universes initialized back in a list.
            For big dump files or lots of data files, this might lead to memory issues.
            Use :meth:`~pylimer_tools_cpp.UniverseSequence.__iter__`
            or :meth:`~pylimer_tools_cpp.UniverseSequence.at_index`
            and :meth:`~pylimer_tools_cpp.UniverseSequence.forget_at_index`
            to craft a more memory-efficient retrieval mechanism.
            
            Returns:
                A list of all Universe objects in the sequence
            )pbdoc")
    // computations
    .def("compute_msd_for_atoms",
         &UniverseSequence::computeMsdForAtoms,
         R"pbdoc(
          Compute the mean square displacement for atoms with the specified IDs.
          
          :param atom_ids: List of atom IDs for which to compute the MSD
          :param nr_of_origins: Number of time origins to use for averaging. Higher values provide better statistics but increase computation time (default: 25)
          :param reduce_memory: If True, reduces memory usage by forgetting universes after processing them (default: False)
          :param max_tau: Maximum time lag (tau) to compute. If -1, computes for all possible tau values. For better statistics, consider setting this to approximately half the sequence length (default: -1)
          :return: Dictionary mapping time lag (tau) to mean square displacement values
     )pbdoc",
         py::arg("atom_ids"),
         py::arg("nr_of_origins") = 25,
         py::arg("reduce_memory") = false,
         py::arg("max_tau") = -1)
    .def("compute_msd_for_atom_properties",
         &UniverseSequence::computeMsdForAtomProperties,
         R"pbdoc(
          Compute the mean square displacement for atoms using specified property names.
          
          :param atom_ids: List of atom IDs for which to compute the MSD
          :param x_property: Name of the x-coordinate property in the dump file (e.g., "x", "xu", "xs")
          :param y_property: Name of the y-coordinate property in the dump file (e.g., "y", "yu", "ys")
          :param z_property: Name of the z-coordinate property in the dump file (e.g., "z", "zu", "zs")
          :param nr_of_origins: Number of time origins to use for averaging. Higher values provide better statistics but increase computation time (default: 25)
          :param reduce_memory: If True, reduces memory usage by forgetting universes after processing them (default: False)
          :param max_tau: Maximum time lag (tau) to compute. If -1, computes for all possible tau values. For better statistics, consider setting this to approximately half the sequence length (default: -1)
          :return: Dictionary mapping time lag (tau) to mean square displacement values
     )pbdoc",
         py::arg("atom_ids"),
         py::arg("x_property"),
         py::arg("y_property"),
         py::arg("z_property"),
         py::arg("nr_of_origins") = 25,
         py::arg("reduce_memory") = false,
         py::arg("max_tau") = -1)
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
      [](UniverseSequence& us, const size_t index) {
        if (index > us.getLength()) {
          throw py::index_error();
        }
        return us.atIndex(index);
      },
      R"pbdoc(
      Get a universe by its index.
      
      :param index: The index of the universe to retrieve
      :return: The Universe at the specified index
      )pbdoc")
    .def("__len__", &UniverseSequence::getLength, R"pbdoc(
         Get the number of universes in this sequence.
         
         :return: The total number of universes available
         )pbdoc")
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
