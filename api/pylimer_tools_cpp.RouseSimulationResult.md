# RouseSimulationResult

### *class* pylimer_tools_cpp.RouseSimulationResult

Bases: `pybind11_object`

Output container for Rouse dynamics simulation runs.

### Attributes Summary

| [`end_to_end_correlation`](#pylimer_tools_cpp.RouseSimulationResult.end_to_end_correlation)   | Normalized end-to-end vector autocorrelation at sampled times.   |
|---------------------------------------------------------------------------|------------------------------------------------------------------|
| [`msd`](#pylimer_tools_cpp.RouseSimulationResult.msd)                      | Bead mean-squared displacement at sampled times.                 |
| [`msd_center_of_mass`](#pylimer_tools_cpp.RouseSimulationResult.msd_center_of_mass)       | Center-of-mass mean-squared displacement at sampled times.       |
| [`times`](#pylimer_tools_cpp.RouseSimulationResult.times)                    | Sampled simulation times.                                        |
| [`trajectory`](#pylimer_tools_cpp.RouseSimulationResult.trajectory)               | List of (N, 3) position matrices if save_trajectory was enabled. |

### Attributes Documentation

#### end_to_end_correlation

Normalized end-to-end vector autocorrelation at sampled times.

#### msd

Bead mean-squared displacement at sampled times.

#### msd_center_of_mass

Center-of-mass mean-squared displacement at sampled times.

#### times

Sampled simulation times.

#### trajectory

List of (N, 3) position matrices if save_trajectory was enabled.
