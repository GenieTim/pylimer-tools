Dissipative Particle Dynamics (DPD) Simulations
===============================================

pylimer-tools provides a powerful framework for simulating coarse-grained polymer systems using Dissipative Particle Dynamics (DPD) with slip-springs. 
This allows for realistic modeling of polymer behavior at a mesoscopic scale.
The implementation is based on `doi:10.1021/acs.macromol.1c00567 <https://doi.org/10.1021/acs.macromol.1c00567>`_.

Refer to the full documentation of the DPD simulator in the :class:`~pylimer_tools_cpp.DPDSimulator` class for more details.

Example Usage
-------------

Here's a minimal example of setting up and running a DPD simulation:

.. code-block:: python

    from pylimer_tools_cpp import DPDSimulator, AtomStyle, OutputConfiguration, ComputedIntValues, ComputedDoubleValues
    from pylimer_tools.io.read_lammps_output_file import read_data_file

    # Load a LAMMPS data file
    filePath = "your_lammps_data_file.structure.out"
    universe = read_data_file(
      filePath,
      atom_style=AtomStyle.FULL
    )

    # Create the DPD simulator
    dpd_simulator = DPDSimulator(
        universe=universe,
        seed=12345,  # Set a random seed for reproducibility
    )

    # Randomly sample slip-springs
    dpd_simulator.create_slip_springs(
      num=100,  # Number of slip-springs to create
    )

    step_output_config = OutputConfiguration()
    step_output_config.int_values = [
      ComputedIntValues.STEP,  # Output the simulation step
    ]
    step_output_config.double_values = [
      ComputedDoubleValues.TEMPERATURE,  # Output the temperature
      ComputedDoubleValues.PRESSURE,  # Output the pressure
      ComputedDoubleValues.STRESS_XX,  # Output the xx component of the stress tensor
      ComputedDoubleValues.STRESS_YY,  # Output the yy component of the stress tensor
      ComputedDoubleValues.STRESS_ZZ,  # Output the zz component of the stress tensor
    ]
    step_output_config.output_every = 10  # Output every 10 steps
    dpd_simulator.config_step_output(
      step_output_config
    )

    dpd_simulator.run_simulation(
      n_steps=1000,  # Number of simulation steps
      dt=0.01,  # Time step size
      with_MC=True,  # Enable Monte Carlo moves of slip-springs
    )

    print("DPD simulation completed.")
