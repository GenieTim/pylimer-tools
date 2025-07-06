Trajectory Analysis
===================

pylimer-tools provides advanced tools for analyzing molecular dynamics trajectories, particularly for polymer systems. This includes reading trajectory data and computing various properties.


Example Usage
-------------

Here's a minimal example of analyzing a polymer trajectory:

.. code-block:: python
  
    import numpy as np
    from pylimer_tools_cpp import UniverseSequence
    from pylimer_tools.io.read_lammps_output_file import read_data_file

    # Load a LAMMPS dump trjectory file
    # This file should contain multiple frames of the polymer system
    filePath = "your_lammps_dump_file.lammpstrj"
    sequence = UniverseSequence()
    sequence.initialize_from_dump_file(
        initial_data_file="initial_complete_frame.data",
        dump_file=filePath
    )

    # Compute the mean square displacement for all atoms
    msd = sequence.compute_msd_for_atoms(
      atom_ids=[
        a.get_id() for a in sequence[0].get_atoms()
      ],
      nr_of_origins = 10,  # Number of origins to use for MSD calculation
    )

    # plot the results
    pass

    # You can also compute other properties like radius of gyration, end-to-end distance, etc.
    # by using the "lazy" iterator over the sequence
    # This allows you to process large trajectories without loading everything into memory at once.
    rgs = []
    timesteps = []
    for universe in sequence:
        rgs.append(np.mean([m.compute_radius_of_gyration() for m in universe.get_molecules()]))
        timesteps.append(universe.get_timestep())

    # plot the radius of gyration over time
    pass

For more information, refer to the documentation of the :class:`~pylimer_tools_cpp.UniverseSequence` class and the :class:`~pylimer_tools_cpp.Universe` class for a list of methods available for property computations.
