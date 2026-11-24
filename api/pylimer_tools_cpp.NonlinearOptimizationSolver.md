# NonlinearOptimizationSolver

### *class* pylimer_tools_cpp.NonlinearOptimizationSolver(self: [pylimer_tools_cpp.NonlinearOptimizationSolver](#pylimer_tools_cpp.NonlinearOptimizationSolver), config: [pylimer_tools_cpp.MEHPNonlinearSolverConfig](pylimer_tools_cpp.MEHPNonlinearSolverConfig.md#pylimer_tools_cpp.MEHPNonlinearSolverConfig) = MEHPNonlinearSolverConfig())

Bases: [`MEHPSolver`](pylimer_tools_cpp.MEHPSolver.md#pylimer_tools_cpp.MEHPSolver)

General nonlinear energy minimization solver for MEHP networks using NLopt.

### Attributes Summary

| [`config`](#pylimer_tools_cpp.NonlinearOptimizationSolver.config)   |    |
|-----------------------------------------------------------|----|

### Methods Summary

| [`get_config`](#pylimer_tools_cpp.NonlinearOptimizationSolver.get_config)(self)                              |                                                                  |
|------------------------------------------------------------------------------------------------|------------------------------------------------------------------|
| [`set_config`](#pylimer_tools_cpp.NonlinearOptimizationSolver.set_config)(self, config)                      |                                                                  |
| [`solve`](#pylimer_tools_cpp.NonlinearOptimizationSolver.solve)(self, net, u, force_evaluator[, is_2d]) | Solve for equilibrium link displacements on net starting from u. |

### Attributes Documentation

#### config

### Methods Documentation

#### get_config(self: [pylimer_tools_cpp.NonlinearOptimizationSolver](#pylimer_tools_cpp.NonlinearOptimizationSolver)) → [pylimer_tools_cpp.MEHPNonlinearSolverConfig](pylimer_tools_cpp.MEHPNonlinearSolverConfig.md#pylimer_tools_cpp.MEHPNonlinearSolverConfig)

#### set_config(self: [pylimer_tools_cpp.NonlinearOptimizationSolver](#pylimer_tools_cpp.NonlinearOptimizationSolver), config: [pylimer_tools_cpp.MEHPNonlinearSolverConfig](pylimer_tools_cpp.MEHPNonlinearSolverConfig.md#pylimer_tools_cpp.MEHPNonlinearSolverConfig)) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### solve(self: [pylimer_tools_cpp.MEHPSolver](pylimer_tools_cpp.MEHPSolver.md#pylimer_tools_cpp.MEHPSolver), net: [pylimer_tools_cpp.MEHPNetwork](pylimer_tools_cpp.MEHPNetwork.md#pylimer_tools_cpp.MEHPNetwork), u: [object](https://docs.python.org/3/builtins/functions.html#object), force_evaluator: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.md#pylimer_tools_cpp.MEHPForceEvaluator), is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [pylimer_tools_cpp.MEHPSolveResult](pylimer_tools_cpp.MEHPSolveResult.md#pylimer_tools_cpp.MEHPSolveResult)

Solve for equilibrium link displacements on net starting from u.
