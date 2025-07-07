Structure Analysis
==================

pylimer-tools offers various tools to manipulate and analyze the topology of polymer networks.

Here you can find examples of some of the functionality provided by the package.

Decomposing Crosslinked Networks to Chains
------------------------------------------

There are a number of ways how a given :class:`~pylimer_tools_cpp.Universe` can be decomposed.

The most common use case is to decompose a crosslinked network into its constituent chains, which can then be used for further analysis or processing.

To find the chains, there are two main approaches:
- Removing the crosslinkers and analysing the remaining polymer strands using :func:`~pylimer_tools_cpp.Universe.get_molecules()`
- Decomposing the network to chains, but keeping the crosslinkers in place, possibly being repeated in all the chains it is associated with.
   This decomposition can be done using the :func:`~pylimer_tools_cpp.Universe.get_chains_with_crosslinker()` method.

pylimer-tools provides two more advanced methods for decomposing networks:
- If you removed the junctions yourself, you can use :func:`~pylimer_tools_cpp.Universe.get_clusters()`

If you want even more control, you could assemble your own decomposition algorithm using the :class:`~pylimer_tools_cpp.Universe` class and its methods, such as :func:`~pylimer_tools_cpp.Universe.get_atoms_connected_to()`.

Finding Loops
-------------

You can also find loops in the network using the :func:`~pylimer_tools_cpp.Universe.find_loops()` method.

If you want to count the number of loops present in the network, as a function of their length, you can use the :func:`~pylimer_tools_cpp.Universe.count_loop_lengths()` method.

Finally, to find the shortest loop a specific atom is involved in, you can use the :func:`~pylimer_tools_cpp.Universe.find_minimal_order_loop_from()` method.

.. caution::

   There are exponentially many paths between two crosslinkers of a network,
   and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
   You can use the `max_length` parameter to restrict the algorithm to only search for loops up to a certain length.
   Use a negative value to find all loops and paths.

