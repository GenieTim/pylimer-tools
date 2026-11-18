# NormalModeSpectrum

### *class* pylimer_tools_cpp.NormalModeSpectrum(self: [pylimer_tools_cpp.NormalModeSpectrum](#pylimer_tools_cpp.NormalModeSpectrum), eigenvalues: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], weights: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, eigenvectors: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, n]'] | [None](https://docs.python.org/3/builtins/constants.html#None) = None, is_complete: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True)

Bases: `pybind11_object`

Immutable normal mode spectrum and viscoelastic observable evaluator.

Represents any discrete or quadrature-weighted spectrum $\{(\lambda_i, w_i)\}$.

Construct a NormalModeSpectrum from eigenvalues and optional weights/eigenvectors.

* **Parameters:**
  * **eigenvalues** – 1D array of eigenvalues
  * **weights** – Optional 1D array of quadrature weights (defaults to all 1.0)
  * **eigenvectors** – Optional 2D array of eigenvectors (columns correspond to modes)
  * **is_complete** – Whether this spectrum covers all modes of the system

### Attributes Summary

| [`eigenvalues`](#pylimer_tools_cpp.NormalModeSpectrum.eigenvalues)    | Eigenvalues of the spectrum.                                      |
|-----------------------------------------------------------------|-------------------------------------------------------------------|
| [`eigenvectors`](#pylimer_tools_cpp.NormalModeSpectrum.eigenvectors)   | Eigenvectors matrix (columns are modes), or None if not computed. |
| [`is_complete`](#pylimer_tools_cpp.NormalModeSpectrum.is_complete)    | Whether this spectrum is complete.                                |
| [`num_modes`](#pylimer_tools_cpp.NormalModeSpectrum.num_modes)      | Number of modes/quadrature nodes in the spectrum.                 |
| [`num_zero_modes`](#pylimer_tools_cpp.NormalModeSpectrum.num_zero_modes) | Number of zero modes in the spectrum.                             |
| [`weights`](#pylimer_tools_cpp.NormalModeSpectrum.weights)        | Quadrature weights of the spectrum.                               |

### Methods Summary

| [`loss_modulus`](#pylimer_tools_cpp.NormalModeSpectrum.loss_modulus)(self, omega)       | Evaluate the loss modulus $G''(\omega)$ (without prefactor).                        |
|----------------------------------------------------------------------------------|-------------------------------------------------------------------------------------|
| [`moduli`](#pylimer_tools_cpp.NormalModeSpectrum.moduli)(self, omega)             | Evaluate both storage and loss moduli $(G'(\omega), G''(\omega))$ in a single pass. |
| [`storage_modulus`](#pylimer_tools_cpp.NormalModeSpectrum.storage_modulus)(self, omega)    | Evaluate the storage modulus $G'(\omega)$ (without prefactor).                      |
| [`stress_autocorrelation`](#pylimer_tools_cpp.NormalModeSpectrum.stress_autocorrelation)(self, t) | Evaluate stress autocorrelation $G(t)$.                                             |

### Attributes Documentation

#### eigenvalues

Eigenvalues of the spectrum.

#### eigenvectors

Eigenvectors matrix (columns are modes), or None if not computed.

#### is_complete

Whether this spectrum is complete.

#### num_modes

Number of modes/quadrature nodes in the spectrum.

#### num_zero_modes

Number of zero modes in the spectrum.

#### weights

Quadrature weights of the spectrum.

### Methods Documentation

#### loss_modulus(self: [pylimer_tools_cpp.NormalModeSpectrum](#pylimer_tools_cpp.NormalModeSpectrum), omega: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate the loss modulus $G''(\omega)$ (without prefactor).

* **Parameters:**
  **omega** – Angular frequencies
* **Returns:**
  Loss modulus values

#### moduli(self: [pylimer_tools_cpp.NormalModeSpectrum](#pylimer_tools_cpp.NormalModeSpectrum), omega: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [tuple](https://docs.python.org/3/builtins/stdtypes.html#tuple)[[Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]'], [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']]

Evaluate both storage and loss moduli $(G'(\omega), G''(\omega))$
in a single pass.

* **Parameters:**
  **omega** – Angular frequencies
* **Returns:**
  Tuple of (storage_modulus, loss_modulus)

#### storage_modulus(self: [pylimer_tools_cpp.NormalModeSpectrum](#pylimer_tools_cpp.NormalModeSpectrum), omega: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate the storage modulus $G'(\omega)$ (without prefactor).

* **Parameters:**
  **omega** – Angular frequencies
* **Returns:**
  Storage modulus values

#### stress_autocorrelation(self: [pylimer_tools_cpp.NormalModeSpectrum](#pylimer_tools_cpp.NormalModeSpectrum), t: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]']) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Evaluate stress autocorrelation $G(t)$.

* **Parameters:**
  **t** – Time points at which to evaluate
* **Returns:**
  Stress autocorrelation values
