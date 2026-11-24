Structure Analysis
------------------

pylimer-tools offers various tools to manipulate and analyze the topology of polymer networks.

Here you can find examples of some of the functionality provided by the package.

Decomposing Crosslinked Networks to Chains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are a number of ways how a given :class:`~pylimer_tools_cpp.Universe` can be decomposed.

The most common use case is to decompose a crosslinked network into its constituent chains, which can then be used for further analysis or processing.

To find the chains, there are two main approaches:
- Removing the crosslinkers and analysing the remaining polymer strands using :func:`~pylimer_tools_cpp.Universe.get_molecules()`
- Decomposing the network to chains, but keeping the crosslinkers in place, possibly being repeated in all the chains it is associated with. This decomposition can be done using the :func:`~pylimer_tools_cpp.Universe.get_chains_with_crosslinker()` method.

pylimer-tools also provides a more advanced method for decomposing networks:
- If you removed the junctions yourself, you can use :func:`~pylimer_tools_cpp.Universe.get_clusters()`

If you want even more control, you could assemble your own decomposition algorithm using the :class:`~pylimer_tools_cpp.Universe` class and its methods, such as :func:`~pylimer_tools_cpp.Universe.get_atoms_connected_to()`.

Finding Loops
~~~~~~~~~~~~~

You can also find loops in the network using the :func:`~pylimer_tools_cpp.Universe.find_loops()` method.

If you want to count the number of loops present in the network, as a function of their length, you can use the :func:`~pylimer_tools_cpp.Universe.count_loop_lengths()` method.

Finally, to find the shortest loop a specific atom is involved in, you can use the :func:`~pylimer_tools_cpp.Universe.find_minimal_order_loop_from()` method.

.. caution::

   There are exponentially many paths between two crosslinkers of a network,
   and you may run out of memory when using these functions, if your network is lattice-like.

Topological Leaf Pruning & Active Core
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To identify the topologically active 2-core of a network (pruning dangling branches,
1-core trees, and unconnected/free chains), you can use
:func:`~pylimer_tools_cpp.Universe.peel_topological_leaves()` or
:func:`~pylimer_tools_cpp.Universe.identify_active_atoms()`.

.. note::

   Topological leaf pruning operates purely on crosslinker connectivity
   (``strand_from`` / ``strand_to``) and ignores atom coordinates and periodic
   boundary conditions (PBC). As a result, it cannot distinguish a local
   (non-PBC, zero-tension) loop from a PBC-spanning (elastically active) loop.
   (Note that in an entangled network, even a non-PBC primary loop can become
   elastically active if an entanglement or slip-link is associated with it.)
   See :doc:`/auto_examples/structure_analysis/plot_topological_leaf_pruning_vs_pbc`
   for a side-by-side comparison with PBC- and/or entanglement-aware force-balance relaxation
   (:class:`~pylimer_tools_cpp.MEHPSimulator`).

Topological Resistance Distance & Commute Times
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For cyclic and cross-linked networks where many parallel pathways connect pairs of nodes,
:class:`~pylimer_tools_cpp.ResistanceDistanceAnalysis` computes effective resistances
:math:`R_{ij}`, the Kirchhoff index :math:`K(G)`, the Moore-Penrose pseudo-inverse
:math:`\mathbf{L}^+`, random-walk commute times :math:`C_{ij}`, and node resistance
centralities directly from a :class:`~pylimer_tools_cpp.Universe`, sparse graph
Laplacian matrix, or weighted edge list.

See :doc:`/auto_examples/structure_analysis/plot_resistance_distance_analysis` for
analytical verifications, matrix diagnostics, and performance benchmarks on crosslinked networks.

