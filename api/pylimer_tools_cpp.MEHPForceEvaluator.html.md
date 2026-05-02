# MEHPForceEvaluator

### *class* pylimer_tools_cpp.MEHPForceEvaluator(self: [pylimer_tools_cpp.MEHPForceEvaluator](#pylimer_tools_cpp.MEHPForceEvaluator))

Bases: `pybind11_object`

The base interface to change the way the force is evaluated during a MEHP run.

To create a custom force evaluator in Python, subclass this class and implement:

- [`evaluate_force_set_gradient()`](#pylimer_tools_cpp.MEHPForceEvaluator.evaluate_force_set_gradient): Compute the force and optionally its gradient
- [`evaluate_stress_contribution()`](#pylimer_tools_cpp.MEHPForceEvaluator.evaluate_stress_contribution): Compute stress tensor contributions
- [`prepare_for_evaluations()`](#pylimer_tools_cpp.MEHPForceEvaluator.prepare_for_evaluations): Prepare for a batch of evaluations

### Attributes Summary

| [`is_2d`](#pylimer_tools_cpp.MEHPForceEvaluator.is_2d)     |    |
|------------------------------------------------------------|----|
| [`network`](#pylimer_tools_cpp.MEHPForceEvaluator.network) |    |

### Methods Summary

| [`evaluate_force_set_gradient`](#pylimer_tools_cpp.MEHPForceEvaluator.evaluate_force_set_gradient)(self, n, ...[, ...])   | Evaluate the force and optionally compute its gradient.   |
|---------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------|
| [`evaluate_stress_contribution`](#pylimer_tools_cpp.MEHPForceEvaluator.evaluate_stress_contribution)(self, ...)           | Evaluate the stress-contribution for a single spring.     |
| [`prepare_for_evaluations`](#pylimer_tools_cpp.MEHPForceEvaluator.prepare_for_evaluations)(self)                          | Prepare the evaluator for a series of evaluations.        |

### Attributes Documentation

#### is_2d

#### network

### Methods Documentation

#### evaluate_force_set_gradient(self: [pylimer_tools_cpp.MEHPForceEvaluator](#pylimer_tools_cpp.MEHPForceEvaluator), n: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), spring_distances: [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.ArrayLike, numpy.float64, '[m, 1]'], compute_gradient: [bool](https://docs.python.org/3/library/functions.html#bool) = False) → [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)

Evaluate the force and optionally compute its gradient.

This is one of the primary methods to override when creating a
custom force evaluator. Override this in Python by defining a method
that returns a tuple of (force, gradient).

* **Parameters:**
  * **n** – The dimensionality of the problem (the number of coordinates)
  * **spring_distances** – The sequential (x, y, z) spring distances as a vector
  * **compute_gradient** – If True, the gradient should be computed and returned
* **Returns:**
  A tuple (force, gradient) where:
  - force: The computed force value (float)
  - gradient: The gradient as a list (or None if compute_gradient was False)

Example implementation in Python:

```default
def evaluate_force_set_gradient(self, n, spring_distances, compute_gradient):
    force = 0.0
    # ... compute force ...

    gradient = None
    if compute_gradient:
        gradient = [0.0] * n
        # ... compute gradient ...

    return (force, gradient)
```

#### evaluate_stress_contribution(self: [pylimer_tools_cpp.MEHPForceEvaluator](#pylimer_tools_cpp.MEHPForceEvaluator), spring_distances: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat), i: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), j: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt), spring_index: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt)) → [float](https://docs.python.org/3/library/functions.html#float)

Evaluate the stress-contribution for a single spring.

* **Parameters:**
  * **spring_distances** – The three coordinate differences for one spring (list of 3 floats)
  * **i** – The row index of the stress tensor
  * **j** – The column index of the stress tensor
  * **spring_index** – The index of the spring being evaluated
* **Returns:**
  The stress contribution value

#### prepare_for_evaluations(self: [pylimer_tools_cpp.MEHPForceEvaluator](#pylimer_tools_cpp.MEHPForceEvaluator)) → [None](https://docs.python.org/3/library/constants.html#None)

Prepare the evaluator for a series of evaluations.
This method is called before a batch of force evaluations,
allowing for any necessary setup or precomputation.
