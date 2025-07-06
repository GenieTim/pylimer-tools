Normal Mode Analysis
====================

pylimer-tools provides a normal mode analysis (NMA) implementation to predict the loss and storage modulus of polymer networks as published in `doi:10.1021/acs.macromol.4c01429 <https://doi.org/10.1021/acs.macromol.4c01429>`_. 
This is useful for understanding the viscoelastic properties of the network.
The corresponding class is :class:`~pylimer_tools_cpp.NormalModeAnalyzer`.

Example Usage
-------------

Here's a minimal example of using the normal mode analysis:

.. code-block:: python

    from pylimer_tools_cpp import NormalModeAnalyzer
    from pylimer_tools.io.read_lammps_output_file import read_data_file

    # Load a LAMMPS data file
    filePath = "your_lammps_data_file.structure.out"
    universe = read_data_file(filePath)
    edges = universe.get_edges()

    # Create the normal mode analyzer
    nma = NormalModeAnalyzer(
      edges["edge_from"],
      edges["edge_to"],
    )

    # Perform the analysis
    nma.find_all_eigenvalues(compute_eigenvectors=False)

    # Get the storage and loss modulus at a given frequency
    frequency = 1.0  # Example frequency in appropriate units
    storage_modulus = nma.evaluate_storage_modulus(
      omega=[frequency],
    )
    loss_modulus = nma.evaluate_loss_modulus(
      omega=[frequency],
    )

    print(f"Storage Modulus: {storage_modulus}")
    print(f"Loss Modulus: {loss_modulus}")
