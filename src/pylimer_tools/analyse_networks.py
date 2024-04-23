#!/usr/bin/env python
# cli.py
from pylimer_tools_cpp.pylimer_tools_cpp import UniverseSequence
import click
import numpy as np

from pylimer_tools.calc.structure_analysis import (compute_stoichiometric_imbalance,
                                                   compute_crosslinker_conversion, compute_extent_of_reaction
                                                   )


@click.command()
@click.argument('files', nargs=-1, type=click.Path(exists=True))
@click.option('--crosslinker_type', type=int, default=2)
def cli(files, crosslinker_type):
    """
    Basic CLI application reading all passed files, outputting some stats on the structures therein

    Arguments:
      - files: list of files to read
    """
    click.echo("Processing {} files".format(len(files)))
    crosslinker_type = crosslinker_type
    for file_path in files:
        click.echo("\nAnalysing File " + file_path)

        universe_sequence = UniverseSequence()
        universe_sequence.initializeFromDataSequence([file_path])
        universe = universe_sequence.atIndex(0)
        click.echo("Size: {}. Volume: {} u^3 (ρ = {})".format(
            universe.getNrOfAtoms(), universe.getVolume(), universe.getNrOfAtoms() / universe.getVolume()))
        click.echo("{} atoms and {} bonds, {} angles, {} dihedrals".format(universe.getNrOfAtoms(
        ), universe.getNrOfBonds(), universe.getNrOfAngles(), universe.getNrOfDihedralAngles()))
        molecules = universe.getMolecules(crosslinker_type)
        bond_lengths = [np.mean(m.computeBondLengths()) for m in molecules]
        non_none_bond_lengths = [
            bl for bl in bond_lengths if bl is not None and bl > 0]
        click.echo("Mean bond length: {} u, (min: {}, max: {}, median: {}) u".format(
            np.mean(non_none_bond_lengths), np.min(non_none_bond_lengths),
            np.max(non_none_bond_lengths), np.median(non_none_bond_lengths)))
        end_to_end_distances = [m.computeEndToEndDistance() for m in molecules]
        click.echo("Mean end to end distance: {} u".format(
            np.mean([e for e in end_to_end_distances if e is not None and e > 0])))
        click.echo("For {} molecules of mean length of {} atoms".format(
            len(molecules), np.mean([m.getNrOfAtoms() for m in molecules])))
        click.echo("r = {}, p = {} ({}), D = {}".format(
            compute_stoichiometric_imbalance(
                universe, crosslinker_type),
            compute_extent_of_reaction(universe, crosslinker_type),
            # mehp.calculateEffectiveCrosslinkerFunctionality(
            #     universe, crosslinker_type),
            compute_crosslinker_conversion(universe, crosslinker_type),
            universe.computePolydispersityIndex(crosslinker_type)
        ))
        click.echo("")
    click.echo("Arbitrary units used. E.g.: Length: u")


if __name__ == "__main__":
    cli()
