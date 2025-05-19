Nomenclature
============

Here a few words used in this documentation (and/or the library, respectively) with their context:

- :obj:`~pylimer_tools_cpp.Atom`: The atom is used in the same sense as LAMMPS uses it: interchangeably with bead.
- Edge: just like the vertex, the edge refers to the graph's internal connectivity. An edge in this case is used as the bond between the beads/atoms. 
    The distinction between edge and bond is insofar relevant, as that the graph uses a different numbering scheme of the vertices than what might be expected.
    Methods referring to edges require the graph's numbering scheme, whereas methods referring to bonds use the atom's ids as numbering scheme. 
- Junctions: is used identically as cross-linker and corresponds to the beads connecting different strands/molecules in the network/universe.
- `Mass vs. Weight`: in general, we try to follow the official nomenclature. 
    Unfortunately, there are a few inconsistencies where we followed the tradition and kept the word "weight", 
    such as in :func:`~pylimer_tools_cpp.Universe.compute_number_average_molecular_weight()`, 
    :func:`~pylimer_tools_cpp.Universe.compute_weight_average_molecular_weight()`, 
    and  :func:`~pylimer_tools_cpp.Universe.compute_weight_fractions()`.
- :obj:`~pylimer_tools_cpp.Molecule`: a collection of atoms that are necessarily connected. Could be a chain, a strand or a single atom
- :obj:`~pylimer_tools_cpp.Universe`: a concise collection of atoms, not necessarily all connected
- Vertex: this package uses a graph, internally, to monitor the connectivity. A vertex in this case is used as the bead/atom. 

.. _Mass vs. Weight: https://en.wikipedia.org/wiki/Mass_versus_weight

Conventions
===========

In addition to the nomenclature listed above,
we try to follow the following naming conventions:

- all Python expose functionality tries to follow PEP 8
- in C++ codes, in turn, camelCase is used for variable and function names
- whenever there is a configurable parameter on a simulation, the corresponding `config_...()` function can be called with an empty argument to reset the parameter to its default value.
- whenever a subset can be queried by certain properties, the corresponding function uses the naming convention `get_..._by_...` 
