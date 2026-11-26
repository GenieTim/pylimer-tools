# DPDSimulator

### *class* pylimer_tools_cpp.DPDSimulator(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), universe: [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe), crosslinker_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2, slipspring_bond_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 9, is_2d: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False, seed: [str](https://docs.python.org/3/builtins/stdtypes.html#str) = '', assume_box_large_enough: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False)

Bases: `pybind11_object`

A quick-and-dirty implementation of the dissipative particle dynamics (DPD) simulation
with slip-springs as presented by Langeloth *et al.* [[LMBohmMullerPlathe13](../acknowledgements.html.md#id11)]
and Schneider *et al.* [[SFKarimiVarzanehMullerPlathe21](../acknowledgements.html.md#id10)].

Get an instance of this class

### Methods Summary

| [`assume_box_large_enough`](#pylimer_tools_cpp.DPDSimulator.assume_box_large_enough)(self[, assume])           | Configure whether to assume minimum-image PBC is sufficient for all bonds.                                                |
|----------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------|
| [`assumes_box_large_enough`](#pylimer_tools_cpp.DPDSimulator.assumes_box_large_enough)(self)                    | Return `True` if minimum-image PBC is assumed sufficient for bonds, `False` if per-bond periodic box offsets are tracked. |
| [`config_a`](#pylimer_tools_cpp.DPDSimulator.config_a)(self[, A])                               | Configure the force-field (pair-style) parameter A.                                                                       |
| [`config_allow_relocation_in_network`](#pylimer_tools_cpp.DPDSimulator.config_allow_relocation_in_network)(self[, ...])   | Configure whether a relocation step may happen when a slip-spring has ended at a crosslink.                               |
| [`config_assume_box_large_enough`](#pylimer_tools_cpp.DPDSimulator.config_assume_box_large_enough)(self[, assume])    | Configure whether to assume minimum-image PBC is sufficient for all bonds.                                                |
| [`config_auto_correlator_output`](#pylimer_tools_cpp.DPDSimulator.config_auto_correlator_output)(self, values)       | Set which values to compute multiple-tau autocorrelation for.                                                             |
| [`config_average_output`](#pylimer_tools_cpp.DPDSimulator.config_average_output)(self, values)               | Set which values to compute averages for.                                                                                 |
| [`config_bond_formation`](#pylimer_tools_cpp.DPDSimulator.config_bond_formation)(self, ...[, ...])           | Configure how to do bond formation during the run.                                                                        |
| [`config_box_deformation`](#pylimer_tools_cpp.DPDSimulator.config_box_deformation)(self, target_box)          | Configure where to (incrementally) deform the box to during the next simulation run.                                      |
| [`config_gamma`](#pylimer_tools_cpp.DPDSimulator.config_gamma)(self[, gamma])                       | Configure the DPD dissipative friction parameter gamma independently of sigma.                                            |
| [`config_lambda`](#pylimer_tools_cpp.DPDSimulator.config_lambda)(self[, l])                          | Configure the modified velocity verlet integration parameter lambda.                                                      |
| [`config_num_steps_dpd`](#pylimer_tools_cpp.DPDSimulator.config_num_steps_dpd)(self[, num_steps])           | Configure the number of steps to do in one DPD sequence.                                                                  |
| [`config_num_steps_mc`](#pylimer_tools_cpp.DPDSimulator.config_num_steps_mc)(self[, num_steps])            | Configure the number of steps to do in one MC sequence.                                                                   |
| [`config_restart_output`](#pylimer_tools_cpp.DPDSimulator.config_restart_output)(self, file[, output_every]) | Set when to output a restart file.                                                                                        |
| [`config_shift_one_at_a_time`](#pylimer_tools_cpp.DPDSimulator.config_shift_one_at_a_time)(self[, ...])           | Configure whether to shift atoms one at a time.                                                                           |
| [`config_shift_possibility_empty`](#pylimer_tools_cpp.DPDSimulator.config_shift_possibility_empty)(self[, ...])       | Configure the possibility of shifting to empty positions.                                                                 |
| [`config_sigma`](#pylimer_tools_cpp.DPDSimulator.config_sigma)(self[, sigma])                       | Configure the force-field (pair-style) parameter sigma (and set gamma = 0.5 sigma^2).                                     |
| [`config_slipspring_high_cutoff`](#pylimer_tools_cpp.DPDSimulator.config_slipspring_high_cutoff)(self[, cutoff])     | Configure the higher cut-off of how far a pair may be distanced for a slip-spring to be created.                          |
| [`config_slipspring_low_cutoff`](#pylimer_tools_cpp.DPDSimulator.config_slipspring_low_cutoff)(self[, cutoff])      | Configure the lower cut-off of how far a pair may be distanced for a slip-spring to be created.                           |
| [`config_spring_constant`](#pylimer_tools_cpp.DPDSimulator.config_spring_constant)(self[, k])                 | Configure the force-field (bond-style) parameter k, the spring constant.                                                  |
| [`config_step_output`](#pylimer_tools_cpp.DPDSimulator.config_step_output)(self, values)                  | Set which values to log.                                                                                                  |
| [`create_slip_springs`](#pylimer_tools_cpp.DPDSimulator.create_slip_springs)(self, num[, bond_type])       | Randomly add up to the specified number of slip-springs to eligible neighbours within the configured cut-offs.            |
| [`get_bond_lengths`](#pylimer_tools_cpp.DPDSimulator.get_bond_lengths)(self)                            | Get the lengths of all bonds in the system.                                                                               |
| [`get_coordinates`](#pylimer_tools_cpp.DPDSimulator.get_coordinates)(self)                             | Get the current particle coordinates.                                                                                     |
| [`get_current_timestep`](#pylimer_tools_cpp.DPDSimulator.get_current_timestep)(self)                        | Get the current timestep number.                                                                                          |
| [`get_dpd_gamma`](#pylimer_tools_cpp.DPDSimulator.get_dpd_gamma)(self)                               | Get the configured DPD dissipative friction coefficient gamma.                                                            |
| [`get_nr_of_atoms`](#pylimer_tools_cpp.DPDSimulator.get_nr_of_atoms)(self)                             | Get the total number of atoms in the system.                                                                              |
| [`get_nr_of_bonds`](#pylimer_tools_cpp.DPDSimulator.get_nr_of_bonds)(self)                             | Get the number of regular bonds (excluding slip-springs).                                                                 |
| [`get_nr_of_bonds_to_form`](#pylimer_tools_cpp.DPDSimulator.get_nr_of_bonds_to_form)(self)                     | Get the number of bonds that are configured to have to be formed.                                                         |
| [`get_nr_of_extra_atoms`](#pylimer_tools_cpp.DPDSimulator.get_nr_of_extra_atoms)(self)                       | Get the number of extra atoms (always 0 for DPD simulations).                                                             |
| [`get_nr_of_extra_bonds`](#pylimer_tools_cpp.DPDSimulator.get_nr_of_extra_bonds)(self)                       | Get the number of extra bonds (slip-springs).                                                                             |
| [`get_nr_of_slip_springs`](#pylimer_tools_cpp.DPDSimulator.get_nr_of_slip_springs)(self)                      | Get the current number of slip-springs in the system.                                                                     |
| [`get_nr_of_steps_dpd`](#pylimer_tools_cpp.DPDSimulator.get_nr_of_steps_dpd)(self)                         | Get the configured number of DPD steps per sequence.                                                                      |
| [`get_nr_of_steps_mc`](#pylimer_tools_cpp.DPDSimulator.get_nr_of_steps_mc)(self)                          | Get the configured number of Monte Carlo steps per sequence.                                                              |
| [`get_num_relocations`](#pylimer_tools_cpp.DPDSimulator.get_num_relocations)(self)                         | Get the number of accepted slip-spring relocation moves in the most recent MC sequence.                                   |
| [`get_num_shifts`](#pylimer_tools_cpp.DPDSimulator.get_num_shifts)(self)                              | Get the number of accepted slip-spring shift moves in the most recent MC sequence.                                        |
| [`get_shift_one_at_a_time`](#pylimer_tools_cpp.DPDSimulator.get_shift_one_at_a_time)(self)                     | Get whether slip-springs are shifted one at a time.                                                                       |
| [`get_shift_possibility_empty`](#pylimer_tools_cpp.DPDSimulator.get_shift_possibility_empty)(self)                 | Get whether shifting to empty positions is allowed.                                                                       |
| [`get_sigma`](#pylimer_tools_cpp.DPDSimulator.get_sigma)(self)                                   | Get the configured DPD random force noise amplitude sigma.                                                                |
| [`get_slip_spring_bond_type`](#pylimer_tools_cpp.DPDSimulator.get_slip_spring_bond_type)(self)                   | Get the bond type identifier used for slip-springs.                                                                       |
| [`get_spring_constant`](#pylimer_tools_cpp.DPDSimulator.get_spring_constant)(self)                         | Get the current spring constant value.                                                                                    |
| [`get_stress_tensor`](#pylimer_tools_cpp.DPDSimulator.get_stress_tensor)(self)                           | Get the current stress tensor.                                                                                            |
| [`get_target_kbt`](#pylimer_tools_cpp.DPDSimulator.get_target_kbt)(self)                              | Get the target thermal energy $k_B T = \sigma^2 / (2 \gamma)$ configured on the DPD thermostat.                           |
| [`get_temperature`](#pylimer_tools_cpp.DPDSimulator.get_temperature)(self)                             | Get the current system temperature.                                                                                       |
| [`get_timestep`](#pylimer_tools_cpp.DPDSimulator.get_timestep)(self)                                | Get the timestep used in the simulation.                                                                                  |
| [`get_universe`](#pylimer_tools_cpp.DPDSimulator.get_universe)(self[, with_slipsprings])            | Get a universe instance from the current coordinates (and connectivity).                                                  |
| [`get_velocities`](#pylimer_tools_cpp.DPDSimulator.get_velocities)(self)                              | Get the current particle velocities.                                                                                      |
| [`get_volume`](#pylimer_tools_cpp.DPDSimulator.get_volume)(self)                                  | Get the current system volume.                                                                                            |
| [`read_restart_file`](#pylimer_tools_cpp.DPDSimulator.read_restart_file)(file)                           | Read a restart file in order to continue a simulation.                                                                    |
| [`refresh_current_state`](#pylimer_tools_cpp.DPDSimulator.refresh_current_state)(self)                       | After re-configuring the force-field parameters, this method should be called to update the current stress tensor etc.    |
| [`relocate_slip_springs`](#pylimer_tools_cpp.DPDSimulator.relocate_slip_springs)(self[, kbt])                | Perform one sweep of Monte Carlo slip-spring chain-end relocation moves.                                                  |
| [`run_simulation`](#pylimer_tools_cpp.DPDSimulator.run_simulation)(self, n_steps[, dt, with_MC])      |                                                                                                                           |
| [`shift_slip_springs`](#pylimer_tools_cpp.DPDSimulator.shift_slip_springs)(self[, kbt])                   | Perform one sweep of Monte Carlo slip-spring shift moves along polymer chains.                                            |
| [`start_measuring_msd_for_atoms`](#pylimer_tools_cpp.DPDSimulator.start_measuring_msd_for_atoms)(self, atom_ids)     | Set a new origin for measuing the mean square displacement for a specified set of atoms                                   |
| [`validate_neighbour_list`](#pylimer_tools_cpp.DPDSimulator.validate_neighbour_list)(self, cutoff)             | Validate the neighbor list consistency for debugging purposes.                                                            |
| [`validate_state`](#pylimer_tools_cpp.DPDSimulator.validate_state)(self)                              | Validate the current simulation state for debugging purposes.                                                             |
| [`write_restart_file`](#pylimer_tools_cpp.DPDSimulator.write_restart_file)(self, file)                    | Explicitly force the writing of a restart file, now!                                                                      |

### Methods Documentation

#### assume_box_large_enough(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), assume: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure whether to assume minimum-image PBC is sufficient for all bonds.

If your bonds could exceed half the box length, keep this `False` (default) so per-bond periodic box offsets are tracked.
Set to `True` to apply minimum-image PBC on every bond evaluation.

#### assumes_box_large_enough(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [bool](https://docs.python.org/3/builtins/functions.html#bool)

Return `True` if minimum-image PBC is assumed sufficient for bonds, `False` if per-bond periodic box offsets are tracked.

#### config_a(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), A: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 25.0) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure the force-field (pair-style) parameter A.

#### config_allow_relocation_in_network(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), allow_relocation_in_network: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure whether a relocation step may happen when a slip-spring has ended at a crosslink.

Side-effect: if true, the relocations may also happen *to* a slip-spring next to a crosslink.

* **Parameters:**
  **config_allow_relocation_in_network.****(****bool****)** (*allow_relocation_in_network*) – Whether to allow relocation in the network or not.

#### config_assume_box_large_enough(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), assume: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure whether to assume minimum-image PBC is sufficient for all bonds.

If your bonds could exceed half the box length, keep this `False` (default) so per-bond periodic box offsets are tracked.
Set to `True` to apply minimum-image PBC on every bond evaluation.

#### config_auto_correlator_output(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), values: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[pylimer_tools_cpp.OutputConfiguration](pylimer_tools_cpp.OutputConfiguration.html.md#pylimer_tools_cpp.OutputConfiguration)], num_corr_in: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 32, p: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 16, m: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set which values to compute multiple-tau autocorrelation for.
If you use this, you should cite Ramírez *et al.* [[RamirezSVL10](../acknowledgements.html.md#id4)].

* **Parameters:**
  * **values** – a list of OutputConfiguration structs
  * **num_corr_in** – Number of correlations in
  * **p** – Parameter p for the autocorrelator
  * **m** – Parameter m for the autocorrelator

#### config_average_output(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), values: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[pylimer_tools_cpp.OutputConfiguration](pylimer_tools_cpp.OutputConfiguration.html.md#pylimer_tools_cpp.OutputConfiguration)]) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set which values to compute averages for.

* **Parameters:**
  **values** – A list of OutputConfiguration structs specifying what to average

#### config_bond_formation(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), num_bonds_to_form: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), max_bonds_per_atom_type: [collections.abc.Mapping](https://docs.python.org/3/library/collections.abc.html#collections.abc.Mapping)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)], bond_formation_dist: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0, attempt_bond_formation_every: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 50, atom_type_form_from: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2, atom_type_form_to: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure how to do bond formation during the run.

* **Parameters:**
  * **config_bond_formation.****(****int****)** (*atom_type_form_to*) – The nr of bonds to form in total. Use 0 to stop bond formation.
  * **config_bond_formation.****(****dict****)** (*num_bonds_per_atom_type*) – The nr of bonds each atom type may have at most (e.g., 2 for strand atoms, 4 for a tertiary crosslinkers)
  * **config_bond_formation.****(****float****)** (*bond_formation_dist*) – The maximum distance allowed to form bonds
  * **config_bond_formation.****(****int****)** – attempt to form bonds every this many steps during the simulation run
  * **config_bond_formation.****(****int****)** – The atom type to start forming bonds from.
  * **config_bond_formation.****(****int****)** – The atom type to start forming bonds to.

#### config_box_deformation(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), target_box: [pylimer_tools_cpp.Box](pylimer_tools_cpp.Box.html.md#pylimer_tools_cpp.Box)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure where to (incrementally) deform the box to during the next simulation run.

#### config_gamma(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), gamma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 4.5) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure the DPD dissipative friction parameter gamma independently of sigma.

#### config_lambda(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), l: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.65) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure the modified velocity verlet integration parameter lambda.

#### config_num_steps_dpd(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), num_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 500) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure the number of steps to do in one DPD sequence.

* **Parameters:**
  **num_steps** – Number of DPD steps per sequence

#### config_num_steps_mc(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), num_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 500) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure the number of steps to do in one MC sequence.

* **Parameters:**
  **num_steps** – Number of Monte Carlo steps per sequence

#### config_restart_output(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), file: [str](https://docs.python.org/3/builtins/stdtypes.html#str), output_every: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 50000) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set when to output a restart file.

#### NOTE
The filename determines the type of serialization:
.json, .xml are supported; other file endings will lead to binary serialization (fastest!).

#### CAUTION
This method may not be backwards- nor forward-compatible.
Use the same version of pylimer-tools if you want to be sure that things work.

* **Parameters:**
  * **file** – The file path to the restart file to write
  * **output_every** – How often to write the restart file (default: 50000)

#### config_shift_one_at_a_time(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), shift_one_at_a_time: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure whether to shift atoms one at a time.

This setting affects Monte Carlo move behavior in the simulation.

#### config_shift_possibility_empty(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), shift_possibility_empty: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure the possibility of shifting to empty positions.

This setting affects Monte Carlo moves in the simulation.

#### config_sigma(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), sigma: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 3.0) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure the force-field (pair-style) parameter sigma (and set gamma = 0.5 sigma^2).

#### config_slipspring_high_cutoff(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2.0) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure the higher cut-off of how far a pair may be distanced for a slip-spring to be created.

#### config_slipspring_low_cutoff(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.5) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure the lower cut-off of how far a pair may be distanced for a slip-spring to be created.

#### config_spring_constant(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), k: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 2.0) → [None](https://docs.python.org/3/builtins/constants.html#None)

Configure the force-field (bond-style) parameter k, the spring constant.

#### config_step_output(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), values: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[pylimer_tools_cpp.OutputConfiguration](pylimer_tools_cpp.OutputConfiguration.html.md#pylimer_tools_cpp.OutputConfiguration)]) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set which values to log.

* **Parameters:**
  **values** – a list of OutputConfiguration structs

#### create_slip_springs(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), num: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), bond_type: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0) → [int](https://docs.python.org/3/builtins/functions.html#int)

Randomly add up to the specified number of slip-springs to eligible neighbours within the configured cut-offs.

* **Parameters:**
  * **num** – Maximum number of slip-springs to create
  * **bond_type** – Bond type for created slip-springs (`0` inherits `slipspring_bond_type`)
* **Returns:**
  Number of slip-springs actually created

#### get_bond_lengths(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the lengths of all bonds in the system.

* **Returns:**
  Vector containing the length of each bond

#### get_coordinates(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the current particle coordinates.

* **Returns:**
  Vector of particle coordinates (x1,y1,z1,x2,y2,z2,…)

#### get_current_timestep(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the current timestep number.

* **Returns:**
  The current timestep index

#### get_dpd_gamma(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the configured DPD dissipative friction coefficient gamma.

#### get_nr_of_atoms(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the total number of atoms in the system.

* **Returns:**
  Total number of atoms

#### get_nr_of_bonds(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the number of regular bonds (excluding slip-springs).

* **Returns:**
  Number of regular bonds

#### get_nr_of_bonds_to_form(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the number of bonds that are configured to have to be formed.

#### get_nr_of_extra_atoms(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the number of extra atoms (always 0 for DPD simulations).

* **Returns:**
  Number of extra atoms

#### get_nr_of_extra_bonds(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the number of extra bonds (slip-springs).

* **Returns:**
  Number of slip-springs

#### get_nr_of_slip_springs(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the current number of slip-springs in the system.

* **Returns:**
  Number of slip-springs

#### get_nr_of_steps_dpd(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the configured number of DPD steps per sequence.

* **Returns:**
  Number of DPD steps

#### get_nr_of_steps_mc(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the configured number of Monte Carlo steps per sequence.

* **Returns:**
  Number of MC steps

#### get_num_relocations(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the number of accepted slip-spring relocation moves in the most recent MC sequence.

#### get_num_shifts(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the number of accepted slip-spring shift moves in the most recent MC sequence.

#### get_shift_one_at_a_time(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [bool](https://docs.python.org/3/builtins/functions.html#bool)

Get whether slip-springs are shifted one at a time.

* **Returns:**
  True if shifting one at a time, False otherwise

#### get_shift_possibility_empty(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [bool](https://docs.python.org/3/builtins/functions.html#bool)

Get whether shifting to empty positions is allowed.

* **Returns:**
  True if shifting to empty positions is allowed

#### get_sigma(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the configured DPD random force noise amplitude sigma.

#### get_slip_spring_bond_type(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [int](https://docs.python.org/3/builtins/functions.html#int)

Get the bond type identifier used for slip-springs.

* **Returns:**
  Bond type for slip-springs

#### get_spring_constant(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the current spring constant value.

* **Returns:**
  The current spring constant for bond interactions

#### get_stress_tensor(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[3, 3]']

Get the current stress tensor.

* **Returns:**
  3x3 stress tensor matrix

#### get_target_kbt(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the target thermal energy $k_B T = \sigma^2 / (2 \gamma)$ configured on the DPD thermostat.

* **Returns:**
  Target thermal energy $k_B T$

#### get_temperature(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the current system temperature.

* **Returns:**
  The current temperature of the system

#### get_timestep(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the timestep used in the simulation.

* **Returns:**
  The simulation timestep value

#### get_universe(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), with_slipsprings: [bool](https://docs.python.org/3/builtins/functions.html#bool) = True) → [pylimer_tools_cpp.Universe](pylimer_tools_cpp.Universe.html.md#pylimer_tools_cpp.Universe)

Get a universe instance from the current coordinates (and connectivity).

* **Parameters:**
  **with_slipsprings** – Whether to include slip-springs in the returned universe (default: True)

#### get_velocities(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [Annotated](https://docs.python.org/3/library/typing.html#typing.Annotated)[numpy.typing.NDArray[numpy.float64], '[m, 1]']

Get the current particle velocities.

* **Returns:**
  Vector of particle velocities (vx1,vy1,vz1,vx2,vy2,vz2,…)

#### get_volume(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [float](https://docs.python.org/3/builtins/functions.html#float)

Get the current system volume.

* **Returns:**
  Current simulation box volume

#### *static* read_restart_file(file: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)

Read a restart file in order to continue a simulation.

* **Parameters:**
  **file** – The file path to the restart file to read

#### refresh_current_state(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [None](https://docs.python.org/3/builtins/constants.html#None)

After re-configuring the force-field parameters,
this method should be called to update the current stress tensor etc.

#### relocate_slip_springs(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), kbt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0) → [int](https://docs.python.org/3/builtins/functions.html#int)

Perform one sweep of Monte Carlo slip-spring chain-end relocation moves.

* **Parameters:**
  **kbt** – Thermal energy scale $k_B T$ for Metropolis acceptance (must be `>= 0`)
* **Returns:**
  Number of accepted relocation moves

#### run_simulation(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), n_steps: [SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex), dt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 0.06, with_MC: [bool](https://docs.python.org/3/builtins/functions.html#bool) = False) → [None](https://docs.python.org/3/builtins/constants.html#None)

#### shift_slip_springs(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), kbt: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex) = 1.0) → [int](https://docs.python.org/3/builtins/functions.html#int)

Perform one sweep of Monte Carlo slip-spring shift moves along polymer chains.

* **Parameters:**
  **kbt** – Thermal energy scale $k_B T$ for Metropolis acceptance (must be `>= 0`)
* **Returns:**
  Number of accepted shift moves

#### start_measuring_msd_for_atoms(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), atom_ids: [collections.abc.Sequence](https://docs.python.org/3/library/collections.abc.html#collections.abc.Sequence)[[SupportsInt](https://docs.python.org/3/library/typing.html#typing.SupportsInt) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)]) → [None](https://docs.python.org/3/builtins/constants.html#None)

Set a new origin for measuing the mean square displacement for a specified set of atoms

#### validate_neighbour_list(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), cutoff: [SupportsFloat](https://docs.python.org/3/library/typing.html#typing.SupportsFloat) | [SupportsIndex](https://docs.python.org/3/library/typing.html#typing.SupportsIndex)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Validate the neighbor list consistency for debugging purposes.

* **Parameters:**
  **cutoff** – Cutoff distance for validation

#### validate_state(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Validate the current simulation state for debugging purposes.

Checks internal data structure consistency and throws exceptions if issues are found.

#### write_restart_file(self: [pylimer_tools_cpp.DPDSimulator](#pylimer_tools_cpp.DPDSimulator), file: [str](https://docs.python.org/3/builtins/stdtypes.html#str)) → [None](https://docs.python.org/3/builtins/constants.html#None)

Explicitly force the writing of a restart file, now!

* **Parameters:**
  **file** – The file path and name of the restart file to be written.
  Can end in .xml, .json or anything else (-> binary)
