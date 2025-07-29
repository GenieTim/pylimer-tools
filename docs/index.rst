pylimer-tools Documentation
===========================

.. image:: https://codecov.io/gh/GenieTim/pylimer-tools/branch/main/graph/badge.svg?token=5ZE1VSDXJQ 
   :target: https://codecov.io/gh/GenieTim/pylimer-tools
   :alt: Test Coverage

.. image:: https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml/badge.svg
   :target: https://github.com/GenieTim/pylimer-tools/actions/workflows/run-tests.yml
   :alt: Tests

.. image:: https://badge.fury.io/py/pylimer-tools.svg
   :target: https://badge.fury.io/py/pylimer-tools
   :alt: PyPI version

.. image:: https://img.shields.io/pypi/dm/pylimer-tools.svg
   :target: https://pypi.python.org/pypi/pylimer-tools/
   :alt: Downloads

Welcome to pylimer-tools, a comprehensive Python library for working with bead-spring polymer systems and for analyzing LAMMPS molecular dynamics simulation output. 
This library provides powerful tools for reading, processing, and analyzing polymer networks with a focus on performance and ease of use.

🚀 **Quick Start**
------------------

Install pylimer-tools with pip:

.. code-block:: bash

   pip install pylimer-tools

Basic usage example:

.. code-block:: python

   import numpy as np
   from pylimer_tools_cpp import UniverseSequence

   # Load a LAMMPS data file
   filePath = "your_lammps_data_file.structure.out"
   universeSequence = UniverseSequence()
   universeSequence.initialize_from_data_sequence([filePath])
   universe = universeSequence.at_index(0)
   
   # Analyze the system
   print(f"System size: {universe.get_size()}")
   print(f"Volume: {universe.get_volume()} u³")
   print(f"Mean bond length: {np.mean(universe.compute_bond_lengths())} u")

📚 **Documentation Sections**
-----------------------------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started
   :hidden:

   installation
   usage
   
.. toctree::
   :maxdepth: 2
   :caption: User Guide
   :hidden:

   nomenclature
   assumptions
   cli
   auto_examples/index
   acknowledgements

.. toctree::
   :maxdepth: 3
   :caption: API Reference
   :hidden:

   modules

🔧 **Key Features**
-------------------

- **Polymer Analysis**: Calculate radius of gyration, end-to-end distances, molecular weights, bond and loop statistics
- **Network Generation**: Monte Carlo network generators for polymer systems with highly flexible parameters
- **Force Balance**: Reduce polymer networks to their minimum energy, maximum entropy homogenized state, predict the equilibrium shear modulus
- **Normal Mode Analysis**: Predict the loss and storage modulus of polymer networks
- **Dissipative Particle Dynamics**: Simulate coarse-grained polymer systems with slip-springs and DPD
- **High Performance**: C++ backend with Python bindings for optimal performance on large polymer structures
- **LAMMPS Integration**: Seamlessly read and write LAMMPS data and dump files with memory-efficient streaming
- **Comprehensive**: Support for angles, dihedrals, crosslinkers, velocities, charges, and custom polymer properties
- **Trajectory Analysis**: Advanced tools for analyzing molecular dynamics trajectories
- **Batch Processing**: Efficient tools for processing multiple files and time series data

💡 **Quick Navigation**
-----------------------

- **New to pylimer-tools?** Start with :doc:`installation` and :doc:`usage`
- **Working with files?** Check out :doc:`auto_examples/readers_writers/index`
- **Generating networks?** See :doc:`auto_examples/network_generator/index`
- **Need examples?** Browse :doc:`auto_examples/index` for practical code samples
- **Understanding limitations?** Review :doc:`assumptions` and :doc:`nomenclature`

🔗 **Useful Links**
-------------------

- `GitHub Repository <https://github.com/GenieTim/pylimer-tools>`_
- `PyPI Package <https://pypi.org/project/pylimer-tools/>`_
- `Issue Tracker <https://github.com/GenieTim/pylimer-tools/issues>`_
- `Tests <https://github.com/GenieTim/pylimer-tools/tree/main/tests>`_
- `Examples <https://github.com/GenieTim/pylimer-tools/tree/main/examples>`_

📖 **Indices and Tables**
-------------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
* :ref:`references`
