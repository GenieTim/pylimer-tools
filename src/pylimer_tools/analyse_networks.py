#!/usr/bin/env python
# cli.py
import click
import numpy as np

from pylimer_tools.calc.structure_analysis import (
    compute_crosslinker_conversion,
    compute_extent_of_reaction,
    compute_stoichiometric_imbalance,
)
from pylimer_tools_cpp import UniverseSequence


@click.command()
@click.argument("files", nargs=-1, type=click.Path(exists=True))
@click.option("--crosslinker_type", type=int, default=2)
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
        universe_sequence.initialize_from_data_sequence([file_path])
        universe = universe_sequence.at_index(0)
        click.echo(
            "Size: {}. Volume: {} u^3 (ρ = {})".format(
                universe.get_nr_of_atoms(),
                universe.get_volume(),
                universe.get_nr_of_atoms() / universe.get_volume(),
            )
        )
        click.echo(
            "{} atoms and {} bonds, {} angles, {} dihedrals".format(
                universe.get_nr_of_atoms(),
                universe.get_nr_of_bonds(),
                universe.get_nr_of_angles(),
                universe.get_nr_of_dihedral_angles(),
            )
        )
        molecules = universe.get_molecules(crosslinker_type)
        bond_lengths = universe.compute_bond_lengths()
        non_none_bond_lengths = [
            bl for bl in bond_lengths if bl is not None and bl > 0]
        click.echo(
            "Bond length b: <b> = {} u, (min: {}, max: {}, median: {}) u, <b^2> = {} u^2".format(
                np.mean(non_none_bond_lengths),
                np.min(non_none_bond_lengths),
                np.max(non_none_bond_lengths),
                np.median(non_none_bond_lengths),
                np.mean(np.square(bond_lengths)),
            )
        )
        end_to_end_distances = universe.compute_end_to_end_distances(
            crosslinker_type=crosslinker_type, derive_image_flags=True
        )
        non_none_end_to_end_distances = [
            e for e in end_to_end_distances if e is not None and e > 0
        ]
        click.echo(
            "End to end distance R_ee: <R_ee> = {} u, <R_ee^2> = {}".format(
                np.mean(non_none_end_to_end_distances),
                np.mean(np.square(non_none_end_to_end_distances)),
            )
        )
        click.echo(
            "For {} molecules of mean length of {} atoms".format(
                len(molecules), np.mean(
                    [m.get_nr_of_atoms() for m in molecules])
            )
        )
        click.echo(
            "f = {}, r = {}, p = {} ({}), D = {}".format(
                universe.determine_functionality_per_type()[crosslinker_type],
                compute_stoichiometric_imbalance(universe, crosslinker_type),
                compute_extent_of_reaction(universe, crosslinker_type),
                # mehp.calculateEffectiveCrosslinkerFunctionality(
                #     universe, crosslinker_type),
                compute_crosslinker_conversion(universe, crosslinker_type),
                universe.compute_polydispersity_index(crosslinker_type),
            )
        )
        click.echo("")
    click.echo("Arbitrary units used. E.g.: Length: u")


if __name__ == "__main__":
    cli()
