Nomenclature
============

Here a few words used in this documentation (and/or the library, respectively) with their context:

- :obj:`~pylimer_tools_cpp.Atom`: the atom is used in the same sense as LAMMPS uses it: interchangably with bead.
- Edge: just like the vertex, the edge refers to the graph's internal connectivity. An edge in this case is used as the bond between the beads/atoms. 
    The distinction between edge and bond is insofar relevant, as that the graph uses a different numbering scheme of the vertices than what might be expected.
    Methods referring to edges require the graph's numbering scheme, whereas methods referring to bonds use the atom's ids as numbering scheme. 
- Junctions: is used identically as cross-linker and corresponds to the beads connecting different strands/molecules in the network/universe.
- `Mass vs. Weight`: in general, we try to follow the official nomenclature. 
    Unfortunately, there are a few inconsistencies where we followed the tradition and kept the word "weight", 
    such as in :func:`~pylimer_tools_cpp.Universe.computeNumberAverageMolecularWeight()`, 
    :func:`~pylimer_tools_cpp.Universe.computeWeightAverageMolecularWeight()`, 
    and  :func:`~pylimer_tools_cpp.Universe.computeWeightFractions()`.
- :obj:`~pylimer_tools_cpp.Molecule`: a collection of atoms that are necessarily connected. Could be a chain, a strand or a single atom
- :obj:`~pylimer_tools_cpp.Universe`: a concise collection of atoms, not necessarily all connected
- Vertex: this package uses a graph, internally, to monitor the connectivity. A vertex in this case is used as the bead/atom. 

.. _Mass vs. Weight: https://en.wikipedia.org/wiki/Mass_versus_weight
