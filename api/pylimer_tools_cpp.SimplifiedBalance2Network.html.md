# SimplifiedBalance2Network

### *class* pylimer_tools_cpp.SimplifiedBalance2Network

Bases: `pybind11_object`

A more efficient structure of the network for use in
[`MEHPForceBalance2`](pylimer_tools_cpp.MEHPForceBalance2.md#pylimer_tools_cpp.MEHPForceBalance2).
Consists usually only of the cross- and slip-links (and their connectivity),
i.e., no “normal strand beads” in between, in order to reduce the degrees of freedom
and therewith improve performance of the solver.

A note on the terminology: a spring is the connection between two links (crosslink, entanglement-link/slip-link).
A strand is a chain of connected links between two crosslinks.

### Attributes Summary

| [`box_lengths`](#pylimer_tools_cpp.SimplifiedBalance2Network.box_lengths)                             |    |
|-------------------------------------------------------------------------------------------------------|----|
| [`coordinates`](#pylimer_tools_cpp.SimplifiedBalance2Network.coordinates)                             |    |
| [`link_is_entanglement`](#pylimer_tools_cpp.SimplifiedBalance2Network.link_is_entanglement)           |    |
| [`links_of_strand`](#pylimer_tools_cpp.SimplifiedBalance2Network.links_of_strand)                     |    |
| [`nr_of_crosslinks`](#pylimer_tools_cpp.SimplifiedBalance2Network.nr_of_crosslinks)                   |    |
| [`nr_of_links`](#pylimer_tools_cpp.SimplifiedBalance2Network.nr_of_links)                             |    |
| [`nr_of_springs`](#pylimer_tools_cpp.SimplifiedBalance2Network.nr_of_springs)                         |    |
| [`nr_of_strands`](#pylimer_tools_cpp.SimplifiedBalance2Network.nr_of_strands)                         |    |
| [`old_atom_ids`](#pylimer_tools_cpp.SimplifiedBalance2Network.old_atom_ids)                           |    |
| [`old_atom_types`](#pylimer_tools_cpp.SimplifiedBalance2Network.old_atom_types)                       |    |
| [`spring_box_offset`](#pylimer_tools_cpp.SimplifiedBalance2Network.spring_box_offset)                 |    |
| [`spring_contour_length`](#pylimer_tools_cpp.SimplifiedBalance2Network.spring_contour_length)         |    |
| [`spring_coordinate_index_a`](#pylimer_tools_cpp.SimplifiedBalance2Network.spring_coordinate_index_a) |    |
| [`spring_coordinate_index_b`](#pylimer_tools_cpp.SimplifiedBalance2Network.spring_coordinate_index_b) |    |
| [`spring_index_a`](#pylimer_tools_cpp.SimplifiedBalance2Network.spring_index_a)                       |    |
| [`spring_index_b`](#pylimer_tools_cpp.SimplifiedBalance2Network.spring_index_b)                       |    |
| [`spring_is_entanglement`](#pylimer_tools_cpp.SimplifiedBalance2Network.spring_is_entanglement)       |    |
| [`springs_of_strand`](#pylimer_tools_cpp.SimplifiedBalance2Network.springs_of_strand)                 |    |
| [`strand_of_spring`](#pylimer_tools_cpp.SimplifiedBalance2Network.strand_of_spring)                   |    |
| [`strands_of_link`](#pylimer_tools_cpp.SimplifiedBalance2Network.strands_of_link)                     |    |

### Attributes Documentation

#### box_lengths

#### coordinates

#### link_is_entanglement

#### links_of_strand

#### nr_of_crosslinks

#### nr_of_links

#### nr_of_springs

#### nr_of_strands

#### old_atom_ids

#### old_atom_types

#### spring_box_offset

#### spring_contour_length

#### spring_coordinate_index_a

#### spring_coordinate_index_b

#### spring_index_a

#### spring_index_b

#### spring_is_entanglement

#### springs_of_strand

#### strand_of_spring

#### strands_of_link
