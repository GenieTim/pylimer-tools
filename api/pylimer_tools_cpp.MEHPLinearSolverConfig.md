# MEHPLinearSolverConfig

### *class* pylimer_tools_cpp.MEHPLinearSolverConfig(self: [pylimer_tools_cpp.MEHPLinearSolverConfig](#pylimer_tools_cpp.MEHPLinearSolverConfig), solver: [pylimer_tools_cpp.SLESolver](pylimer_tools_cpp.SLESolver.md#pylimer_tools_cpp.SLESolver) = SLESolver.DEFAULT, simplification_mode: [pylimer_tools_cpp.StructureSimplificationMode](pylimer_tools_cpp.StructureSimplificationMode.md#pylimer_tools_cpp.StructureSimplificationMode) = StructureSimplificationMode.NO_SIMPLIFICATION, inactive_removal_cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1e-06, residual_reduction: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1e-12, max_iterations: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 50000, spring_breaking_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = -1.0, in_place_rebuild_threshold: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.15, cg_fallback_max_iterations: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 20000, gradient_descent_learning_rate: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.01, slip_mode: [pylimer_tools_cpp.SlipMode](pylimer_tools_cpp.SlipMode.md#pylimer_tools_cpp.SlipMode) = SlipMode.FIXED)

Bases: `pybind11_object`

Configuration for the sparse linear-system MEHP equilibrium solver.

### Attributes Summary

| [`cg_fallback_max_iterations`](#pylimer_tools_cpp.MEHPLinearSolverConfig.cg_fallback_max_iterations)     |    |
|---------------------------------------------------------------------------------|----|
| [`gradient_descent_learning_rate`](#pylimer_tools_cpp.MEHPLinearSolverConfig.gradient_descent_learning_rate) |    |
| [`in_place_rebuild_threshold`](#pylimer_tools_cpp.MEHPLinearSolverConfig.in_place_rebuild_threshold)     |    |
| [`inactive_removal_cutoff`](#pylimer_tools_cpp.MEHPLinearSolverConfig.inactive_removal_cutoff)        |    |
| [`max_iterations`](#pylimer_tools_cpp.MEHPLinearSolverConfig.max_iterations)                 |    |
| [`residual_reduction`](#pylimer_tools_cpp.MEHPLinearSolverConfig.residual_reduction)             |    |
| [`simplification_mode`](#pylimer_tools_cpp.MEHPLinearSolverConfig.simplification_mode)            |    |
| [`slip_mode`](#pylimer_tools_cpp.MEHPLinearSolverConfig.slip_mode)                      |    |
| [`solver`](#pylimer_tools_cpp.MEHPLinearSolverConfig.solver)                         |    |
| [`spring_breaking_length`](#pylimer_tools_cpp.MEHPLinearSolverConfig.spring_breaking_length)         |    |

### Attributes Documentation

#### cg_fallback_max_iterations

#### gradient_descent_learning_rate

#### in_place_rebuild_threshold

#### inactive_removal_cutoff

#### max_iterations

#### residual_reduction

#### simplification_mode

#### slip_mode

#### solver

#### spring_breaking_length
