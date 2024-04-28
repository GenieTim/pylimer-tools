Readers & Writers 
=================

Readers
-------

pylimer_tools provides an interface to read LAMMPS data 
(:obj:`~pylimer_tools_cpp.DataFileReader`) 
and dump 
(:obj:`~pylimer_tools_cpp.DumpFileReader`) 
files.

Additionally, in :obj:`~pylimer_tools_cpp.UniverseSequence`,
an additional level of abstraction is provided: the UniverseSequence automatically 
translates the data from an (internal) :obj:`~pylimer_tools_cpp.DataFileReader`
or :obj:`~pylimer_tools_cpp.DumpFileReader` to a 
:obj:`~pylimer_tools_cpp.Universe`, while also handling this whole operation 
of reading the file into :obj:`~pylimer_tools_cpp.Universe`s memory-efficient
(if used correctly).

Writers
-------

In :obj:`~pylimer_tools_cpp.DataFileWriter`, 
pylimer_tools provides an interface to write LAMMPS data files.

Some numbering schemes in the written file can be configured as seen in the examples below,
allowing for conversions from one numbering scheme to another.

Examples
--------

Converting a system of polymer chains to have the "molecule index" in the LAMMPS file 
set in such a way that a `fix bond/swap` would not change the distribution of chain lengths:

.. code:: python

  import os

  from pylimer_tools_cpp import Universe, UniverseSequence, DataFileWriter

  # TODO: change the following parameters to your liking
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
