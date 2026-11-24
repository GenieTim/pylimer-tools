# MEHPSolver

### *class* pylimer_tools_cpp.MEHPSolver

Bases: `pybind11_object`

Abstract strategy interface for solving mechanical equilibrium on an MEHPNetwork.

### Methods Summary

| [`solve`](#pylimer_tools_cpp.MEHPSolver.solve)(self, net, u, force_evaluator[, is_2d])   | Solve for equilibrium link displacements on net starting from u.   |
|--------------------------------------------------------------------------------------------------|--------------------------------------------------------------------|

### Methods Documentation

#### solve(self: [pylimer_tools_cpp.MEHPSolver](#pylimer_tools_cpp.MEHPSolver), net: [pylimer_tools_cpp.MEHPNetwork](pylimer_tools_cpp.MEHPNetwork.md#pylimer_tools_cpp.MEHPNetwork), u: [object](https://docs.python.org/3/builtins/functions.html#object), force_evaluator: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.md#pylimer_tools_cpp.MEHPForceEvaluator), is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [pylimer_tools_cpp.MEHPSolveResult](pylimer_tools_cpp.MEHPSolveResult.md#pylimer_tools_cpp.MEHPSolveResult)

Solve for equilibrium link displacements on net starting from u.
