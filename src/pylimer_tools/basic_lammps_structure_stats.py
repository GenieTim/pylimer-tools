#!/usr/bin/env python
# cli.py
import click
import numpy as np

from pylimer_tools.io.read_lammps_output_file import read_data_file


@click.command()
@click.argument("files", nargs=-1, type=click.Path(exists=True))
def cli(files):
    """
    Basic CLI application reading all passed files, outputting some stats on the structures therein

    Arguments:
      - files: list of files to read
    """
    click.echo("Processing {} files".format(len(files)))
    for file_path in files:
        click.echo("\nAnalysing File " + file_path)

        universe = read_data_file(file_path)
        click.echo(
            "Size: {}. Volume: {} u^3".format(
                universe.get_nr_of_atoms(), universe.get_volume()
            )
        )
        molecules = universe.get_molecules(2)
        bond_lengths = [np.mean(m.compute_bond_lengths()) for m in molecules]
        non_none_bond_lengths = [
            bl for bl in bond_lengths if bl is not None and bl > 0]
        click.echo(
            "Mean bond length: {} u, (min: {}, max: {}, median: {}) u".format(
                np.mean(non_none_bond_lengths),
                np.min(non_none_bond_lengths),
                np.max(non_none_bond_lengths),
                np.median(non_none_bond_lengths),
            )
        )
        end_to_end_distances = [m.compute_end_to_end_distance()
                                for m in molecules]
        click.echo(
            "Mean end to end distance: {} u".format(
                np.mean(
                    [e for e in end_to_end_distances if e is not None and e > 0])
            )
        )
        click.echo(
            "For {} molecules of mean length of {} atoms".format(
                len(molecules), np.mean(
                    [m.get_nr_of_atoms() for m in molecules])
            )
        )
    click.echo("Arbitrary units used. E.g.: Length: u")


if __name__ == "__main__":
    cli()
