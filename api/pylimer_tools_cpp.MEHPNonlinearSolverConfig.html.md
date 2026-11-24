# MEHPNonlinearSolverConfig

### *class* pylimer_tools_cpp.MEHPNonlinearSolverConfig(self: [pylimer_tools_cpp.MEHPNonlinearSolverConfig](#pylimer_tools_cpp.MEHPNonlinearSolverConfig), algorithm: [pylimer_tools_cpp.NLoptAlgorithm](pylimer_tools_cpp.NLoptAlgorithm.html.md#pylimer_tools_cpp.NLoptAlgorithm) = NLoptAlgorithm.LD_MMA, max_nr_of_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 250000, xtol: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1e-12, ftol: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1e-09, simplification_mode: [pylimer_tools_cpp.StructureSimplificationMode](pylimer_tools_cpp.StructureSimplificationMode.html.md#pylimer_tools_cpp.StructureSimplificationMode) = StructureSimplificationMode.NO_SIMPLIFICATION, inactive_removal_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1e-06, spring_breaking_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1.0, slip_mode: [pylimer_tools_cpp.SlipMode](pylimer_tools_cpp.SlipMode.html.md#pylimer_tools_cpp.SlipMode) = SlipMode.FIXED)

Bases: `pybind11_object`

Configuration for general nonlinear optimization MEHP equilibrium solving via NLopt.

### Attributes Summary

| [`algorithm`](#pylimer_tools_cpp.MEHPNonlinearSolverConfig.algorithm)               |    |
|--------------------------------------------------------------------------|----|
| [`ftol`](#pylimer_tools_cpp.MEHPNonlinearSolverConfig.ftol)                    |    |
| [`inactive_removal_cutoff`](#pylimer_tools_cpp.MEHPNonlinearSolverConfig.inactive_removal_cutoff) |    |
| [`max_nr_of_steps`](#pylimer_tools_cpp.MEHPNonlinearSolverConfig.max_nr_of_steps)         |    |
| [`simplification_mode`](#pylimer_tools_cpp.MEHPNonlinearSolverConfig.simplification_mode)     |    |
| [`slip_mode`](#pylimer_tools_cpp.MEHPNonlinearSolverConfig.slip_mode)               |    |
| [`spring_breaking_length`](#pylimer_tools_cpp.MEHPNonlinearSolverConfig.spring_breaking_length)  |    |
| [`xtol`](#pylimer_tools_cpp.MEHPNonlinearSolverConfig.xtol)                    |    |

### Attributes Documentation

#### algorithm

#### ftol

#### inactive_removal_cutoff

#### max_nr_of_steps

#### simplification_mode

#### slip_mode

#### spring_breaking_length

#### xtol
