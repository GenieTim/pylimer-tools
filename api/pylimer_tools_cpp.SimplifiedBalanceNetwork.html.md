# SimplifiedBalanceNetwork

### *class* pylimer_tools_cpp.SimplifiedBalanceNetwork

Bases: `pybind11_object`

A more efficient structure of the network for use in MEHP force balance,
namely [`MEHPForceBalance`](pylimer_tools_cpp.MEHPForceBalance.md#pylimer_tools_cpp.MEHPForceBalance), though also passable to
namely [`MEHPForceBalance2`](pylimer_tools_cpp.MEHPForceBalance2.md#pylimer_tools_cpp.MEHPForceBalance2).
Consists usually only of the cross- and slip-links (and their connectivity),
i.e., no “normal strand beads” in between, in order to reduce the degrees of freedom
and therewith improve performance of the solver.

A note on the terminology: a spring is the connection between two links (crosslink, entanglement-link/slip-link).
A strand is a chain of connected links between two crosslinks.

### Attributes Summary

| [`box_lengths`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.box_lengths)                                     |    |
|--------------------------------------------------------------------------------------------------------------|----|
| [`coordinates`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.coordinates)                                     |    |
| [`link_is_sliplink`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.link_is_sliplink)                           |    |
| [`links_of_strand`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.links_of_strand)                             |    |
| [`nr_of_crosslink_swaps_endured`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.nr_of_crosslink_swaps_endured) |    |
| [`nr_of_crosslinks`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.nr_of_crosslinks)                           |    |
| [`nr_of_links`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.nr_of_links)                                     |    |
| [`nr_of_springs`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.nr_of_springs)                                 |    |
| [`nr_of_strands`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.nr_of_strands)                                 |    |
| [`old_atom_ids`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.old_atom_ids)                                   |    |
| [`spring_box_offset`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.spring_box_offset)                         |    |
| [`spring_coordinate_index_a`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.spring_coordinate_index_a)         |    |
| [`spring_coordinate_index_b`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.spring_coordinate_index_b)         |    |
| [`spring_index_a`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.spring_index_a)                               |    |
| [`spring_index_b`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.spring_index_b)                               |    |
| [`springs_of_strand`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.springs_of_strand)                         |    |
| [`strand_contour_length`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.strand_contour_length)                 |    |
| [`strand_coordinate_index_a`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.strand_coordinate_index_a)         |    |
| [`strand_coordinate_index_b`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.strand_coordinate_index_b)         |    |
| [`strand_index_a`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.strand_index_a)                               |    |
| [`strand_index_b`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.strand_index_b)                               |    |
| [`strand_of_spring`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.strand_of_spring)                           |    |
| [`strands_of_link`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.strands_of_link)                             |    |
| [`volume`](#pylimer_tools_cpp.SimplifiedBalanceNetwork.volume)                                               |    |

### Attributes Documentation

#### box_lengths

#### coordinates

#### link_is_sliplink

#### links_of_strand

#### nr_of_crosslink_swaps_endured

#### nr_of_crosslinks

#### nr_of_links

#### nr_of_springs

#### nr_of_strands

#### old_atom_ids

#### spring_box_offset

#### spring_coordinate_index_a

#### spring_coordinate_index_b

#### spring_index_a

#### spring_index_b

#### springs_of_strand

#### strand_contour_length

#### strand_coordinate_index_a

#### strand_coordinate_index_b

#### strand_index_a

#### strand_index_b

#### strand_of_spring

#### strands_of_link

#### volume
