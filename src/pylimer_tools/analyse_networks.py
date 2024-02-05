#!/usr/bin/env python
# cli.py
import click
import numpy as np

import pylimer_tools.calc.doMEHPAnalysis as mehp
import pylimer_tools.calc.doMMTAnalysis as mmt
from pylimer_tools_cpp import UniverseSequence


@click.command()
@click.argument('files', nargs=-1, type=click.Path(exists=True))
@click.option('--crosslinkertype', type=int, default=2)
def cli(files, crosslinkertype):
    """
    Basic CLI application reading all passed files, outputting some stats on the structures therein

    Arguments:
      - files: list of files to read
    """
    click.echo("Processing {} files".format(len(files)))
    crosslinkerType = crosslinkertype
    for filePath in files:
        click.echo("\nAnalysing File " + filePath)

        universeSequence = UniverseSequence()
        universeSequence.initializeFromDataSequence([filePath])
        universe = universeSequence.atIndex(0)
        click.echo("Size: {}. Volume: {} u^3".format(
            universe.getNrOfAtoms(), universe.getVolume()))
        click.echo("{} atoms and {} bonds, {} angles, {} dihedrals".format(universe.getNrOfAtoms(
        ), universe.getNrOfBonds(), universe.getNrOfAngles(), universe.getNrOfDihedralAngles()))
        molecules = universe.getMolecules(crosslinkerType)
        bondLengths = [np.mean(m.computeBondLengths()) for m in molecules]
        nonNoneBondLengths = [
            l for l in bondLengths if l is not None and l > 0]
        click.echo("Mean bond length: {} u, (min: {}, max: {}, median: {}) u".format(
            np.mean(nonNoneBondLengths), np.min(nonNoneBondLengths), np.max(nonNoneBondLengths), np.median(nonNoneBondLengths)))
        endToEndDistances = [m.computeEndToEndDistance() for m in molecules]
        click.echo("Mean end to end distance: {} u".format(
            np.mean([e for e in endToEndDistances if e is not None and e > 0])))
        click.echo("For {} molecules of mean length of {} atoms".format(
            len(molecules), np.mean([m.getNrOfAtoms() for m in molecules])))
        click.echo("r = {}, p = {} ({}), D = {}".format(
            mmt.computeStoichiometricInbalance(
                universe, crosslinkerType),
            mmt.computeExtentOfReaction(universe, crosslinkerType),
            # mehp.calculateEffectiveCrosslinkerFunctionality(
            #     universe, crosslinkerType),
            mehp.computeCrosslinkerConversion(universe, crosslinkerType),
            universe.computePolydispersityIndex(crosslinkerType)
        ))
        click.echo("")
    click.echo("Arbitrary units used. E.g.: Length: u")


if __name__ == "__main__":
    cli()
