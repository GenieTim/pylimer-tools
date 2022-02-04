Nomenclature
============

Here a few words used in this documentation with their context:

- :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.Atom`: the atom is used in the same sense as LAMMPS uses it: interchangably with bead.
- :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe`: a concise collection of atoms, not necessarily all connected
- :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule`: a collection of atoms that are necessarily connected. Could be a chain, a strand or a single atom
- Vertex: this package uses a graph, internally, to monitor the connectivity. A vertex in this case is used as the bead/atom. 
