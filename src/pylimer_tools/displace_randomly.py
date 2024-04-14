#!/usr/bin/env python
# cli.py
import os
import click
import numpy as np
import random

from pylimer_tools.io.readLammpsOutputFile import readDataFile
from pylimer_tools_cpp.pylimer_tools_cpp import UniverseSequence, DataFileWriter, Atom


@click.command()
@click.argument('file', type=click.Path(exists=True), required=True)
@click.argument('max_displacement', default=0.5, type=click.FLOAT)
def cli(file, max_displacement):
  """
  Basic CLI application iterating all atoms in a file, displacing them by a bit.
  
  Arguments:
    - file: the file to read (and write, with prefix "random-displaced-")
    - max_displacement: the maximum displacement
  """
  universe = readDataFile(file)

  atoms = universe.getAtoms()
  for atom in atoms:
    new_atom = Atom(
      atom.getId(),
      atom.getType(),
      atom.getX() + (random.random() - 0.5)*max_displacement,
      atom.getY() + (random.random() - 0.5)*max_displacement,
      atom.getZ() + (random.random() - 0.5)*max_displacement,
      atom.getNX(),
      atom.getNY(),
      atom.getNZ()
    )
    universe.replaceAtom(atom.getId(), new_atom)
  
  writer = DataFileWriter(universe)

  target_file = os.path.join(
    os.path.dirname(file),
    "random-displaced-" + os.path.basename(file)
  )
  writer.writeToFile(
    target_file
  )
  click.echo("Written file '{}'".format(target_file))


if __name__ == "__main__":
    cli()

