#!/usr/bin/env python
# cli.py
import os
import random

import click

from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import Atom, DataFileWriter


@click.command()
@click.argument("file", type=click.Path(exists=True), required=True)
@click.argument("max_displacement", default=0.5, type=click.FLOAT)
def cli(file, max_displacement):
    """
    Basic CLI application iterating all atoms in a file, displacing them by a bit.

    Arguments:
      - file: The file to read (and write, with prefix "random-displaced-")
      - max_displacement: The maximum displacement
    """
    universe = read_data_file(file)

    atoms = universe.get_atoms()
    for atom in atoms:
        new_atom = Atom(
            atom.get_id(),
            atom.get_type(),
            atom.get_x() + (random.random() - 0.5) * max_displacement,
            atom.get_y() + (random.random() - 0.5) * max_displacement,
            atom.get_z() + (random.random() - 0.5) * max_displacement,
            atom.get_nx(),
            atom.get_ny(),
            atom.get_nz(),
        )
        universe.replace_atom(atom.get_id(), new_atom)

    writer = DataFileWriter(universe)

    target_file = os.path.join(
        os.path.dirname(file), "random-displaced-" + os.path.basename(file)
    )
    writer.write_to_file(target_file)
    click.echo("Written file '{}'".format(target_file))


if __name__ == "__main__":
    cli()
