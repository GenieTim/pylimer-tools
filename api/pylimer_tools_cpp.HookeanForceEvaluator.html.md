# HookeanForceEvaluator

### *class* pylimer_tools_cpp.HookeanForceEvaluator(self: [pylimer_tools_cpp.HookeanForceEvaluator](#pylimer_tools_cpp.HookeanForceEvaluator), kappa: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0)

Bases: [`MEHPForceEvaluator`](pylimer_tools_cpp.MEHPForceEvaluator.html.md#pylimer_tools_cpp.MEHPForceEvaluator)

Linear Gaussian (Hookean) spring constitutive evaluator:
$E_s = \frac{\kappa}{2 N_s} \|\Delta\mathbf{r}_s\|^2$.

* **Parameters:**
  **kappa** – Spring stiffness prefactor (>= 0).

### Methods Summary

| [`evaluate_network_energy_and_gradient`](#pylimer_tools_cpp.HookeanForceEvaluator.evaluate_network_energy_and_gradient)(self, ...)   | Evaluate total network elastic energy and optional link gradient.           |
|----------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------|
| [`evaluate_spring_energy`](#pylimer_tools_cpp.HookeanForceEvaluator.evaluate_spring_energy)(self, spring_vec, ...)     | Evaluate strain energy of a single spring.                                  |
| [`evaluate_spring_force`](#pylimer_tools_cpp.HookeanForceEvaluator.evaluate_spring_force)(self, spring_vec, ...)      | Evaluate spring force vector dE/d(spring_vec) for a single spring.          |
| [`evaluate_spring_stiffness`](#pylimer_tools_cpp.HookeanForceEvaluator.evaluate_spring_stiffness)(self, spring_vec, ...)  | Evaluate 3x3 tangent stiffness matrix dF/d(spring_vec) for a single spring. |
| [`evaluate_stress_contribution`](#pylimer_tools_cpp.HookeanForceEvaluator.evaluate_stress_contribution)(self, ...)           | Evaluate the (i, j) virial stress contribution of a single spring.          |
| [`evaluate_stress_tensor`](#pylimer_tools_cpp.HookeanForceEvaluator.evaluate_stress_tensor)(self, net, ...[, is_2d])   | Evaluate the 3x3 virial stress tensor for the network.                      |
| [`get_kappa`](#pylimer_tools_cpp.HookeanForceEvaluator.get_kappa)(self)                                   | Return the spring stiffness parameter kappa.                                |
| [`is_linear`](#pylimer_tools_cpp.HookeanForceEvaluator.is_linear)(self)                                   | Indicate whether the constitutive force law is linear in displacement.      |

### Methods Documentation

#### evaluate_network_energy_and_gradient(self: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.html.md#pylimer_tools_cpp.MEHPForceEvaluator), net: [pylimer_tools_cpp.MEHPNetwork](pylimer_tools_cpp.MEHPNetwork.html.md#pylimer_tools_cpp.MEHPNetwork), spring_vectors: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, compute_gradient: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)

Evaluate total network elastic energy and optional link gradient.

* **Parameters:**
  * **net** – Network topology and contour lengths.
  * **spring_vectors** – Flattened spring vectors (size 3 \* nr_of_springs).
  * **is_2d** – If True, treat the system as planar (ignore z components).
  * **compute_gradient** – If True, compute and return link gradient ndarray.
* **Returns:**
  Tuple (energy, gradient_or_none).

#### evaluate_spring_energy(self: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.html.md#pylimer_tools_cpp.MEHPForceEvaluator), spring_vec: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[3, 1]'], contour_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Evaluate strain energy of a single spring.

#### evaluate_spring_force(self: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.html.md#pylimer_tools_cpp.MEHPForceEvaluator), spring_vec: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[3, 1]'], contour_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 1]']

Evaluate spring force vector dE/d(spring_vec) for a single spring.

#### evaluate_spring_stiffness(self: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.html.md#pylimer_tools_cpp.MEHPForceEvaluator), spring_vec: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[3, 1]'], contour_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 3]']

Evaluate 3x3 tangent stiffness matrix dF/d(spring_vec) for a single spring.

#### evaluate_stress_contribution(self: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.html.md#pylimer_tools_cpp.MEHPForceEvaluator), spring_vec: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[3, 1]'], i: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), j: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), contour_length: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Evaluate the (i, j) virial stress contribution of a single spring.

#### evaluate_stress_tensor(self: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.html.md#pylimer_tools_cpp.MEHPForceEvaluator), net: [pylimer_tools_cpp.MEHPNetwork](pylimer_tools_cpp.MEHPNetwork.html.md#pylimer_tools_cpp.MEHPNetwork), spring_vectors: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], volume: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 3]']

Evaluate the 3x3 virial stress tensor for the network.

#### get_kappa(self: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.html.md#pylimer_tools_cpp.MEHPForceEvaluator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Return the spring stiffness parameter kappa.

#### is_linear(self: [pylimer_tools_cpp.MEHPForceEvaluator](pylimer_tools_cpp.MEHPForceEvaluator.html.md#pylimer_tools_cpp.MEHPForceEvaluator)) → [bool](https://docs.python.org/3/builtins/functions.html#bool)

Indicate whether the constitutive force law is linear in displacement.
