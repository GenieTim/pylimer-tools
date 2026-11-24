# LinearSystemMEHPSolver

### *class* pylimer_tools_cpp.LinearSystemMEHPSolver(self: [pylimer_tools_cpp.LinearSystemMEHPSolver](#pylimer_tools_cpp.LinearSystemMEHPSolver), config: [pylimer_tools_cpp.MEHPLinearSolverConfig](pylimer_tools_cpp.MEHPLinearSolverConfig.md#pylimer_tools_cpp.MEHPLinearSolverConfig) = MEHPLinearSolverConfig())

Bases: [`MEHPSolver`](pylimer_tools_cpp.MEHPSolver.md#pylimer_tools_cpp.MEHPSolver)

Equilibrium solver for linear (Hookean) MEHP networks via sparse linear system assembly.

### Attributes Summary

| [`config`](#pylimer_tools_cpp.LinearSystemMEHPSolver.config)   |    |
|-----------------------------------------------------------|----|

### Methods Summary

| [`get_config`](#pylimer_tools_cpp.LinearSystemMEHPSolver.get_config)(self)                              |                                                                  |
|------------------------------------------------------------------------------------------------|------------------------------------------------------------------|
| [`set_config`](#pylimer_tools_cpp.LinearSystemMEHPSolver.set_config)(self, config)                      |                                                                  |
| [`solve`](#pylimer_tools_cpp.LinearSystemMEHPSolver.solve)(self, net, u, force_evaluator[, is_2d]) | Solve for equilibrium link displacements on net starting from u. |

### Attributes Documentation

#### config

### Methods Documentation

#### get_config(self: [pylimer_tools_cpp.LinearSystemMEHPSolver](#pylimer_tools_cpp.LinearSystemMEHPSolver)) → [pylimer_tools_cpp.MEHPLinearSolverConfig](pylimer_tools_cpp.MEHPLinearSolverConfig.md#pylimer_tools_cpp.MEHPLinearSolverConfig)

#### set_config(self: [pylimer_tools_cpp.LinearSystemMEHPSolver](#pylimer_tools_cpp.LinearSystemMEHPSolver), config: [pylimer_tools_cpp.MEHPLinearSolverConfig](pylimer_tools_cpp.MEHPLinearSolverConfig.md#pylimer_tools_cpp.MEHPLinearSolverConfig)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### solve(self: [pylimer_tools_cpp.MEHPSolver](pylimer_tools_cpp.MEHPSolver.md#pylimer_tools_cpp.MEHPSolver), net: [pylimer_tools_cpp.MEHPNetwork](pylimer_tools_cpp.MEHPNetwork.md#pylimer_tools_cpp.MEHPNetwork), u: [object](https://docs.python.org/3/builtins/functions.html#object), force_evaluator: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.md#pylimer_tools_cpp.MEHPForceEvaluator), is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [pylimer_tools_cpp.MEHPSolveResult](pylimer_tools_cpp.MEHPSolveResult.md#pylimer_tools_cpp.MEHPSolveResult)

Solve for equilibrium link displacements on net starting from u.
