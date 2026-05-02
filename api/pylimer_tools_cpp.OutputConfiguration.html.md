# OutputConfiguration

### *class* pylimer_tools_cpp.OutputConfiguration(self: [pylimer_tools_cpp.OutputConfiguration](#pylimer_tools_cpp.OutputConfiguration))

Bases: `pybind11_object`

A configuration object to configure the output values and frequency
for simulation classes in this package.

This class specifies which quantities to output and how often to write them
during simulations.

Create a new OutputConfiguration instance.

* **Returns:**
  A new OutputConfiguration object with default settings

### Attributes Summary

| [`append`](#pylimer_tools_cpp.OutputConfiguration.append)               | Whether to append to the file or truncate it                    |
|-------------------------------------------------------------------------|-----------------------------------------------------------------|
| [`double_values`](#pylimer_tools_cpp.OutputConfiguration.double_values) | List of double-valued quantities to output.                     |
| [`filename`](#pylimer_tools_cpp.OutputConfiguration.filename)           | The path and name of the file to write to.                      |
| [`int_values`](#pylimer_tools_cpp.OutputConfiguration.int_values)       | List of integer-valued quantities to output.                    |
| [`output_every`](#pylimer_tools_cpp.OutputConfiguration.output_every)   | How often to write the values to the output.                    |
| [`use_every`](#pylimer_tools_cpp.OutputConfiguration.use_every)         | For autocorrelation and averaging, how often to include values. |

### Attributes Documentation

#### append

Whether to append to the file or truncate it

#### double_values

List of double-valued quantities to output.

Use [`ComputedDoubleValues`](pylimer_tools_cpp.ComputedDoubleValues.md#pylimer_tools_cpp.ComputedDoubleValues) enum to specify which floating-point quantities
should be computed and written to output.

#### filename

The path and name of the file to write to.
An empty string (“”) means standard output (console).

#### int_values

List of integer-valued quantities to output.

Use [`ComputedIntValues`](pylimer_tools_cpp.ComputedIntValues.md#pylimer_tools_cpp.ComputedIntValues) enum to specify which integer quantities
should be computed and written to output.

#### output_every

How often to write the values to the output.
For averages, this value also says how many values will be averaged.

#### use_every

For autocorrelation and averaging, how often to include values.

Use a value of 1 to take average of or autocorrelate, respectively,
all values encountered during the simulation or optimization procedure.
