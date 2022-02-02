Readers & Writers 
=================

Readers
-------

pylimer_tools provides an interface to write LAMMPS data and dump files.


Writers
-------

pylimer_tools provides an interface to write LAMMPS data files.


Examples
--------

Converting a system of polymer chains to have the "molecule index" in the LAMMPS file 
set in such a way that a `fix bond/swap` would not change the distribution of chain lengths:

.. code:: python

  import os

  from pylimer_tools_cpp import Universe, UniverseSequence, DataFileWriter

  fileToRead = "yourLammpsDataFileToConvert.structure.out"
  filename, file_extension = os.path.splitext(fileToRead)
  fileToWrite = filename + "_converted" + file_extension

  universeSeq = UniverseSequence()
  universeSeq.initializeFromDataSequence([fileToRead])

  universe = universeSeq.atIndex(0)

  writer = DataFileWriter(universe)
  writer.configIncludeAngles(True)
  writer.configReindexAtoms(True)
  writer.configCrosslinkerType(2)
  writer.configMoleculeIdxForSwap(False)

  writer.writeToFile(fileToWrite)
  print("Written file: {}".format(fileToWrite))
