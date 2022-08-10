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
  sideLen = 0
  seed = 0
  beadDistance = 0
  nrOfCrosslinkers = 0
  numberOfChains = 0
  numberOfSolventChains = 0
  beadsPerSolventChain = 0
  crosslinkerConversion = 0
  crosslinkerFunctionality = 0
  # END TODO

  # Start generating the network
  generator = MCUniverseGenerator(sideLen, sideLen, sideLen)
  generator.setSeed(seed)
  generator.setBeadDistance(beadDistance)
  print("Starting MC Generator {}".format(datetime.now().strftime("%H:%M:%S")))
  generator.addCrosslinkers(nrOfCrosslinkers, 2)
  print("Added cross-linkers {}".format(datetime.now().strftime("%H:%M:%S")))
  generator.addSolventChains(numberOfSolventChains, beadsPerSolventChain, 3)
  print("Added solvent chains {}".format(datetime.now().strftime("%H:%M:%S")))
  generator.addAndLinkStrands(numberOfChains, nrOfBeadsPerChain,
                              crosslinkerConversion, crosslinkerFunctionality, 1)
  print("Added and linked strands {}".format(
      datetime.now().strftime("%H:%M:%S")))
  universe = generator.getUniverse()
  universe.setMasses({1: 1, 2: 1, 3: 1})

  angles = universe.detectAngles()
  universe.addAngles(angles["angle_from"],
                    angles["angle_via"], angles["angle_to"])
  print("Added angles {}".format(datetime.now().strftime("%H:%M:%S")))

  # output new system
  dataWriter = DataFileWriter(universe)
  dataWriter.configIncludeAngles(True)
  dataWriter.configMoleculeIdxForSwap(mode_swappable)
  dataWriter.writeToFile(fileToWrite)
