# Correlator

### *class* pylimer_tools_cpp.Correlator(self: [pylimer_tools_cpp.Correlator](#pylimer_tools_cpp.Correlator), num_correlators: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 32, p: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 16, m: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2)

Bases: `pybind11_object`

Multi-tau scalar autocorrelator computing $f(\tau) = \langle A(t) A(t + \tau) \rangle$.

Uses a hierarchy of coarse-grained shift registers (Likhtman/Ramirez multi-tau algorithm)
to evaluate time correlations across exponentially spanning lag times in O(1) amortized time per sample.

Construct a multi-tau scalar Correlator.

* **Parameters:**
  * **num_correlators** – Number of hierarchical correlator levels (default: 32)
  * **p** – Number of lag points per correlator level (default: 16)
  * **m** – Averaging factor between consecutive levels; recommended `p % m == 0` (default: 2)

### Attributes Summary

| [`accval`](#pylimer_tools_cpp.Correlator.accval)          | Accumulated sum of level-0 samples.                      |
|------------------------------------------------------------------|----------------------------------------------------------|
| [`correlation`](#pylimer_tools_cpp.Correlator.correlation)     | Evaluated correlation values (call evaluate() first).    |
| [`m`](#pylimer_tools_cpp.Correlator.m)               | Averaging factor between consecutive levels.             |
| [`mean`](#pylimer_tools_cpp.Correlator.mean)            | Running sample mean at level 0.                          |
| [`npcorr`](#pylimer_tools_cpp.Correlator.npcorr)          | Number of populated correlation points after evaluate(). |
| [`num_correlators`](#pylimer_tools_cpp.Correlator.num_correlators) | Number of correlator levels configured.                  |
| [`num_samples`](#pylimer_tools_cpp.Correlator.num_samples)     | Total number of samples added at level 0.                |
| [`p`](#pylimer_tools_cpp.Correlator.p)               | Number of lag points per correlator level.               |
| [`raw_f`](#pylimer_tools_cpp.Correlator.raw_f)           | Full pre-allocated correlation buffer.                   |
| [`raw_t`](#pylimer_tools_cpp.Correlator.raw_t)           | Full pre-allocated lag time buffer.                      |
| [`times`](#pylimer_tools_cpp.Correlator.times)           | Evaluated lag times (call evaluate() first).             |

### Methods Summary

| [`add`](#pylimer_tools_cpp.Correlator.add)(self, value[, level])               | Add a scalar observation to correlator level `level`.                    |
|------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|
| [`add_batch`](#pylimer_tools_cpp.Correlator.add_batch)(self, values)                 | Add a 1D array of scalar observations at level 0.                        |
| [`evaluate`](#pylimer_tools_cpp.Correlator.evaluate)(self[, norm])                  | Evaluate the autocorrelation function at all populated lag times.        |
| [`reset`](#pylimer_tools_cpp.Correlator.reset)(self)                             | Reset accumulated state while keeping the current correlator dimensions. |
| [`set_size`](#pylimer_tools_cpp.Correlator.set_size)(self[, num_correlators, p, m]) | Reconfigure correlator dimensions and reset accumulated state.           |

### Attributes Documentation

#### accval

Accumulated sum of level-0 samples.

#### correlation

Evaluated correlation values (call evaluate() first).

#### m

Averaging factor between consecutive levels.

#### mean

Running sample mean at level 0.

#### npcorr

Number of populated correlation points after evaluate().

#### num_correlators

Number of correlator levels configured.

#### num_samples

Total number of samples added at level 0.

#### p

Number of lag points per correlator level.

#### raw_f

Full pre-allocated correlation buffer.

#### raw_t

Full pre-allocated lag time buffer.

#### times

Evaluated lag times (call evaluate() first).

### Methods Documentation

#### add(self: [pylimer_tools_cpp.Correlator](#pylimer_tools_cpp.Correlator), value: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), level: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0) → [None](https://docs.python.org/3/builtins/constants.html#None)

Add a scalar observation to correlator level `level`.

* **Parameters:**
  * **value** – Finite scalar observation
  * **level** – Correlator level index (default: 0)

#### add_batch(self: [pylimer_tools_cpp.Correlator](#pylimer_tools_cpp.Correlator), values: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']) → [None](https://docs.python.org/3/builtins/constants.html#None)

Add a 1D array of scalar observations at level 0.

* **Parameters:**
  **values** – 1D array of finite scalar observations

#### evaluate(self: [pylimer_tools_cpp.Correlator](#pylimer_tools_cpp.Correlator), norm: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [None](https://docs.python.org/3/builtins/constants.html#None)

Evaluate the autocorrelation function at all populated lag times.

* **Parameters:**
  **norm** – If True, subtract the squared sample mean $\langle A \rangle^2$ (default: False)

#### reset(self: [pylimer_tools_cpp.Correlator](#pylimer_tools_cpp.Correlator)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Reset accumulated state while keeping the current correlator dimensions.

#### set_size(self: [pylimer_tools_cpp.Correlator](#pylimer_tools_cpp.Correlator), num_correlators: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 32, p: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 16, m: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2) → [None](https://docs.python.org/3/builtins/constants.html#None)

Reconfigure correlator dimensions and reset accumulated state.

* **Parameters:**
  * **num_correlators** – Number of hierarchical correlator levels (default: 32)
  * **p** – Number of lag points per correlator level (default: 16)
  * **m** – Averaging factor between consecutive levels (default: 2)
