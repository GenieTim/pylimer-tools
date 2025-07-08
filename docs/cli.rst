Command Line Interface
======================

Apart from the Python API, `pylimer-tools` provides a command line interface (CLI) for various functionalities.
This allows you to perform tasks directly from the terminal without writing Python scripts.

Network Generation
------------------

Generate crosslinked polymer networks using Monte Carlo procedures.

.. click:: pylimer_tools.generate_network:cli
   :prog: pylimer-generate-network
   :nested: full

Network Analysis
----------------

Analyze existing polymer networks to compute structural properties.

.. click:: pylimer_tools.analyse_networks:cli
   :prog: pylimer-analyse-networks
   :nested: full

LAMMPS Data Statistics
----------------------

Compute basic statistics from LAMMPS data files.

.. click:: pylimer_tools.basic_lammps_structure_stats:cli
   :prog: pylimer-basic-lammps-stats
   :nested: full

Random Displacement
-------------------

Randomly displace atoms in a structure file.

.. click:: pylimer_tools.displace_randomly:cli
   :prog: pylimer-displace-randomly
   :nested: full
