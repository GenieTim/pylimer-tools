# AtomPairEntanglements

### *class* pylimer_tools_cpp.AtomPairEntanglements(self: [pylimer_tools_cpp.AtomPairEntanglements](#pylimer_tools_cpp.AtomPairEntanglements))

Bases: `pybind11_object`

A struct to store pairs of atoms that are close together and could be entanglements.

Get an instance of this struct

### Attributes Summary

| [`pair_of_atom`](#pylimer_tools_cpp.AtomPairEntanglements.pair_of_atom)     | An index in the pairs_of_atoms if the atom is part of a pair, -1 else.         |
|-----------------------------------------------------------------------------|--------------------------------------------------------------------------------|
| [`pairs_of_atoms`](#pylimer_tools_cpp.AtomPairEntanglements.pairs_of_atoms) | A list of pairs of atom ids that are close together and could be entanglements |

### Attributes Documentation

#### pair_of_atom

An index in the pairs_of_atoms if the atom is part of a pair, -1 else.

#### pairs_of_atoms

A list of pairs of atom ids that are close together and could be entanglements
