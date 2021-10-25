#!/usr/bin/env python
# cli.py
import click
import numpy as np

from pylimer_tools.entities.universum import Universum
from pylimer_tools.io.readLammpData import readLammpData


@click.command()
@click.argument('files', nargs=-1, type=click.Path(exists=True))
def cli(files):
    """
    Basic CLI application reading all passed files, outputting some stats on the structures therein

    Arguments:
      - files: list of files to read
    """
    click.echo("Processing {} files".format(len(files)))
    for filePath in files:
        click.echo("Analysing File " + filePath)
        allData = readLammpData(filePath)
        universe = Universum(boxSizes=[
                             allData["Lx"], allData["Ly"], allData["Lz"]])
        universe.addAtomBondData(allData['atom_data'], allData['bond_data'])
        click.echo("Size: {}. Volume: {} u^3".format(
            universe.getSize(), universe.getVolume()))
        bondLengths = [m.computeBondLengths().mean() for m in universe]
        click.echo("Mean bond length: {} u".format(
            np.mean([l for l in bondLengths if l is not None])))
        endToEndDistances = [m.computeEndToEndDistance() for m in universe]
        click.echo("Mean end to end distance: {} u".format(
            np.mean([e for e in endToEndDistances if e is not None])))
    click.echo("Arbitrary units used. E.g.: Length: u")


if __name__ == "__main__":
    cli()
