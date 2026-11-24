# MEHPSlipLinkSolverConfig

### *class* pylimer_tools_cpp.MEHPSlipLinkSolverConfig(self: [pylimer_tools_cpp.MEHPSlipLinkSolverConfig](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig), max_nr_of_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 50000, xtol: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1e-09, initial_residual_to_use: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1.0, simplification_mode: [pylimer_tools_cpp.StructureSimplificationMode](pylimer_tools_cpp.StructureSimplificationMode.html.md#pylimer_tools_cpp.StructureSimplificationMode) = StructureSimplificationMode.NO_SIMPLIFICATION, inactive_removal_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.001, do_inner_iterations: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, link_swapping_mode: [pylimer_tools_cpp.LinkSwappingMode](pylimer_tools_cpp.LinkSwappingMode.html.md#pylimer_tools_cpp.LinkSwappingMode) = LinkSwappingMode.NO_SWAPPING, swapping_frequency: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 10, one_over_spring_partition_upper_limit: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, nr_of_crosslink_swaps_allowed_per_sliplink: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1, slip_mode: [pylimer_tools_cpp.SlipMode](pylimer_tools_cpp.SlipMode.html.md#pylimer_tools_cpp.SlipMode) = SlipMode.SLIDING, spring_breaking_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1.0)

Bases: `pybind11_object`

Configuration for alternating coordinate/partition slip-link relaxation.

### Attributes Summary

| [`do_inner_iterations`](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig.do_inner_iterations)                        |    |
|---------------------------------------------------------------------------------------------|----|
| [`inactive_removal_cutoff`](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig.inactive_removal_cutoff)                    |    |
| [`initial_residual_to_use`](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig.initial_residual_to_use)                    |    |
| [`link_swapping_mode`](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig.link_swapping_mode)                         |    |
| [`max_nr_of_steps`](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig.max_nr_of_steps)                            |    |
| [`nr_of_crosslink_swaps_allowed_per_sliplink`](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig.nr_of_crosslink_swaps_allowed_per_sliplink) |    |
| [`one_over_spring_partition_upper_limit`](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig.one_over_spring_partition_upper_limit)      |    |
| [`simplification_mode`](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig.simplification_mode)                        |    |
| [`slip_mode`](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig.slip_mode)                                  |    |
| [`spring_breaking_length`](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig.spring_breaking_length)                     |    |
| [`swapping_frequency`](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig.swapping_frequency)                         |    |
| [`xtol`](#pylimer_tools_cpp.MEHPSlipLinkSolverConfig.xtol)                                       |    |

### Attributes Documentation

#### do_inner_iterations

#### inactive_removal_cutoff

#### initial_residual_to_use

#### link_swapping_mode

#### max_nr_of_steps

#### nr_of_crosslink_swaps_allowed_per_sliplink

#### one_over_spring_partition_upper_limit

#### simplification_mode

#### slip_mode

#### spring_breaking_length

#### swapping_frequency

#### xtol
