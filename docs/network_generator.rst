Network Generation
==================

pylimer-tools includes a powerful Monte Carlo-based network generator for creating realistic polymer systems. 
This tools is particularly useful for generating quick initial configurations for molecular dynamics simulations,
or for doing fast predictions of the viscoelasticity using the :class:`~pylimer_tools_cpp.NormalModeAnalyzer` or the 
:class:`pylimer_tools_cpp.MEHPForceBalance` and its variants.

Overview
--------

The network generation system allows you to:

- Create end-linked and vulcanized polymer networks with controllable topology
- Generate solvent chains and background molecules
- Control crosslinking degree and functionality
- Set realistic bead distances and box sizes

The :class:`~pylimer_tools_cpp.MCUniverseGenerator` class is the main interface for network generation.

Basic Network Generation
-----------------------

Simple Crosslinked Network
~~~~~~~~~~~~~~~~~~~~~~~~~~

Here's a minimal example of creating a crosslinked polymer network:

.. code-block:: python

   from pylimer_tools_cpp import MCUniverseGenerator, DataFileWriter
   
   # Create generator for a 50x50x50 simulation box
   generator = MCUniverseGenerator(50.0, 50.0, 50.0)
   
   # Set random seed for reproducibility
   generator.set_seed(12345)
   
   # Set mean bead-to-bead distance
   generator.set_bead_distance(1.0)
   
   # Add crosslinking sites
   generator.add_crosslinkers(100, 4)  # 100 crosslinkers of functionality 4
   
   # Add polymer chains
   generator.add_strands(
       nr_of_strands=200,                      # 200 polymer chains
       strand_lengths=[20 for _ in range(200)],# 20 beads per chain
   )

   # Crosslink polymer strands
   generator.link_strands_to_conversion(
      crosslinker_conversion=0.925,  # 92.5% of crosslinker sites used
   )
   
   # Get the generated network
   universe = generator.get_universe()
   
   # Set masses for LAMMPS
   universe.set_masses({1: 1.0, 2: 2.0})
   
   print(f"Generated network with {universe.get_nr_of_atoms()} atoms")

Advanced Configuration
----------------------

Multi-Component Systems
~~~~~~~~~~~~~~~~~~~~~~~

Create complex systems with multiple components:

.. code-block:: python

   import datetime
   from pylimer_tools_cpp import MCUniverseGenerator, DataFileWriter
   
   # System parameters
   box_size = 60.0
   seed = 42
   bead_distance = 1.0
   
   # Network parameters
   n_crosslinkers = 150

   n_polymer_chains = 80
   beads_per_chain = 25

   n_solvent_chains = 200
   beads_per_solvent = 5

   n_monofunctional_chains = 50
   beads_per_monofunctional_chain = 10

   crosslinker_conversion = 0.85
   crosslinker_functionality = 4
   
   print(f"Starting network generation at {datetime.datetime.now()}")
   
   # Initialize generator
   generator = MCUniverseGenerator(box_size, box_size, box_size)
   generator.set_seed(seed)
   generator.set_bead_distance(bead_distance)
   
   # Add components in order
   print("Adding crosslinkers...")
   generator.add_crosslinkers(
    n_crosslinkers, 
    crosslinker_functionality=crosslinker_functionality,
    crosslinker_type=2)  # Type 2
   
   print("Adding solvent chains...")
   generator.add_solvent_chains(n_solvent_chains, beads_per_solvent, 3)  # Type 3
   
   print("Adding and linking polymer strands...")
   # Add bifunctional polymer chains
   generator.add_strands(
       nr_of_strands=n_polymer_chains,                   
       strand_lengths=[beads_per_chain for _ in range(n_polymer_chains)],
       strand_atom_type=1  # Type 1 for "normal" beads
   )

   generator.add_monofunctional_strands(
       nr_of_monofunctional_strands=n_monofunctional_chains,
       monofunctional_strand_length=[beads_per_monofunctional_chain for _ in range(n_monofunctional_chains)],
       strand_atom_type=4  # Type 4 to recognize these as "monofunctional" strands later
   )

   # Crosslink polymer strands
   generator.link_strands_to_conversion(
      crosslinker_conversion=crosslinker_conversion,  # 92.5% of crosslinker sites used
   )
   
   # Get final system
   universe = generator.get_universe()
   universe.set_masses({1: 1.0, 2: 2.0, 3: 0.5, 4: 1.})  # Different masses
   
   print(f"Network generation completed at {datetime.datetime.now()}")

Monte Carlo Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~

Fine-tune the Monte Carlo generation process:

.. code-block:: python

   # Configure the main options of the generator
   generator.set_seed(12345)  # Set random seed for reproducibility
   generator.set_bead_distance(1.0)  # Set mean bead-to-bead distance
   
   # Configure whether to generate loops or not
   generator.config_primary_loop_probability(0.0)  # Primary loops prohibited
   generator.config_secondary_loop_probability(0.0)  # Secondary loops prohibited

   # Improve performance of sampling pairs of cross-links for large systems
   # Does not sample all pairs of cross-links anymore, 
   # only 99.9994% of those that are selected, but is ca. 60 times faster
   generator.use_zscore_max_distance(3.29) 

   # Configure MC parameters
   # If you do not trust the Brownian Bridge process we use,
   # you can add additional Monte Carlo steps
   generator.config_nr_of_mc_steps(5000)  # More MC steps for better equilibration


   
Adding Topology Information
--------------------------

The generated networks can include topological information:

.. code-block:: python

   # Generate basic network
   generator = MCUniverseGenerator(40, 40, 40)
   generator.set_seed(54321)
   
   # Add crosslinking sites
   generator.add_crosslinkers(100, 4)  # 100 crosslinkers of functionality 4
   
   # Add polymer chains
   generator.add_strands(
       nr_of_strands=200,                      # 200 polymer chains
       strand_lengths=[20 for _ in range(200)],# 20 beads per chain
   )

   # Crosslink polymer strands
   generator.link_strands_to_conversion(
      crosslinker_conversion=0.925,  # 92.5% of crosslinker sites used
   )
   
   # Get the generated universe
   universe = generator.get_universe()
   
   # Detect and add angles
   angles = universe.detect_angles()
   if len(angles["angle_from"]) > 0:
       universe.add_angles(
           angles["angle_from"],
           angles["angle_via"],
           angles["angle_to"],
           angle_types=[1 for _ in range(len(angles["angle_from"]))]
       )
       print(f"Added {len(angles['angle_from'])} angles to the system")
   
   # Detect dihedral angles (if needed)
   dihedrals = universe.detect_dihedral_angles()
   if len(dihedrals["dihedral_from"]) > 0:
       # Add dihedral information if your force field requires it
       universe.add_dihedral_angles(
           dihedrals["dihedral_from"],
           dihedrals["dihedral_via1"],
           dihedrals["dihedral_via2"],
           dihedrals["dihedral_to"],
           angle_types=[1 for _ in range(len(dihedrals["dihedral_from"]))]
       )

Output and Visualization
------------------------

Save Generated Networks
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from pathlib import Path
   
   # Generate network
   generator = MCUniverseGenerator(45, 45, 45)
   # ... configure and generate ...
   universe = generator.get_universe()
   
   # Create output directory
   output_dir = Path("generated_networks")
   output_dir.mkdir(exist_ok=True)
   
   # Write LAMMPS data file
   writer = DataFileWriter(universe)
   writer.config_include_angles(True)
   writer.config_reindex_atoms(True)
   writer.config_atom_style(pylimer_tools_cpp.AtomStyle.ANGLE)
   
   output_file = output_dir / "polymer_network.data"
   writer.write_to_file(str(output_file))
   
   print(f"Network saved to: {output_file}")


For more options, examples and advanced usage, see the `test suite <https://github.com/GenieTim/pylimer-tools/tree/main/tests>`_ 
and the :doc:`modules` API reference, particularly of the :class:`~pylimer_tools_cpp.MCUniverseGenerator` class.
