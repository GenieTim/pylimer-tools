Network Generators
==================

This package also provides a Monte-Carlo type network generator.

It could be used e.g. like this:

.. code:: python

  import os
  import random
  import sys
  from datetime import datetime

  import numpy as np
  import pandas as pd
  from numpy.random import default_rng
  from pylimer_tools.calc import doMMTAnalysis
  from pylimer_tools_cpp import (DataFileWriter, MCUniverseGenerator,)


  # BEGIN TODO: define the parameters you want to generate the network with
  side_len = 0
  seed = 0
  bead_distance = 0
  nr_of_crosslinkers = 0
  number_of_chains = 0
  number_of_solvent_chains = 0
  beads_per_solvent_chain = 0
  crosslinker_conversion = 0
  crosslinker_functionality = 0
  # END TODO

  # Start generating the network
  generator = MCUniverseGenerator(side_len, side_len, side_len)
  generator.set_seed(seed)
  generator.set_bead_distance(bead_distance)
  print("Starting MC Generator {}".format(datetime.now().strftime("%H:%M:%S")))
  generator.add_crosslinkers(nr_of_crosslinkers, 2)
  print("Added cross-linkers {}".format(datetime.now().strftime("%H:%M:%S")))
  generator.add_solvent_chains(numberOfSolventChains, beadsPerSolventChain, 3)
  print("Added solvent chains {}".format(datetime.now().strftime("%H:%M:%S")))
  generator.add_and_link_strands(number_of_chains, nr_of_beads_per_chain,
                              crosslinker_conversion, crosslinker_functionality, 1)
  print("Added and linked strands {}".format(
      datetime.now().strftime("%H:%M:%S")))
  universe = generator.get_universe()
  universe.set_masses({1: 1, 2: 1, 3: 1})

  angles = universe.detect_angles()
  universe.add_angles(angles["angle_from"],
                    angles["angle_via"], angles["angle_to"])
  print("Added angles {}".format(datetime.now().strftime("%H:%M:%S")))

  # output new system
  dataWriter = DataFileWriter(universe)
  dataWriter.config_include_angles(True)
  dataWriter.config_molecule_idx_for_swap(mode_swappable)
  dataWriter.write_to_file(fileToWrite)
