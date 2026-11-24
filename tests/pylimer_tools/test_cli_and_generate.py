#!/usr/bin/env python
import os
import statistics
import tempfile
import unittest

from click.testing import CliRunner

from pylimer_tools.analyse_networks import cli as analyse_networks_cli
from pylimer_tools.basic_lammps_structure_stats import cli as basic_stats_cli
from pylimer_tools.displace_randomly import cli as displace_randomly_cli
from pylimer_tools.generate_network import cli as generate_network_cli
from pylimer_tools.generate_network import generate_structure
from pylimer_tools.io.bead_spring_parameter_provider import (
    ParameterType,
    get_parameters_for_polymer,
)
from pylimer_tools.io.read_lammps_output_file import read_data_file


class TestCliAndGenerate(unittest.TestCase):
    def test_generate_structure_and_bond_length_stats_equivalence(self) -> None:
        params = get_parameters_for_polymer("PDMS", ParameterType.GAUSSIAN)
        universe = generate_structure(
            params=params,
            target_p=0.8,
            n_chains_crosslinkers=20,
            target_f=4,
            n_chains_1=40,
            n_beads_per_chain_1=10,
        )
        self.assertEqual(universe.get_nr_of_atoms(), 40 * 10 + 20)

        bond_stats = universe.compute_bond_length_stats()
        bond_lengths = universe.compute_bond_lengths()
        self.assertEqual(bond_stats.nr_of_bonds, len(bond_lengths))
        self.assertAlmostEqual(
            bond_stats.mean_squared_bond_length,
            statistics.mean(bl**2 for bl in bond_lengths),
            places=6,
        )
        self.assertAlmostEqual(
            bond_stats.mean_bond_length,
            statistics.mean(bond_lengths),
            places=6,
        )
        self.assertAlmostEqual(bond_stats.min_bond_length, min(bond_lengths), places=6)
        self.assertAlmostEqual(bond_stats.max_bond_length, max(bond_lengths), places=6)

    def test_generate_network_cli(self) -> None:
        runner = CliRunner()
        help_res = runner.invoke(generate_network_cli, ["--help"])
        self.assertEqual(help_res.exit_code, 0)
        self.assertIn("--n-chains-crosslinkers", help_res.output)

        with tempfile.TemporaryDirectory() as tmpdir:
            out_file = os.path.join(tmpdir, "network.dat")

            # Missing target-p / target-wsol should fail cleanly
            err_res = runner.invoke(generate_network_cli, ["--target-file", out_file])
            self.assertNotEqual(err_res.exit_code, 0)

            res = runner.invoke(
                generate_network_cli,
                [
                    "--polymer-name",
                    "PDMS",
                    "--parameter-type",
                    "GAUSSIAN",
                    "--target-p",
                    "0.8",
                    "--n-chains-crosslinkers",
                    "20",
                    "--n-chains-1",
                    "40",
                    "--n-beads-per-chain-1",
                    "10",
                    "--target-file",
                    out_file,
                ],
            )
            self.assertEqual(res.exit_code, 0, msg=res.output)
            self.assertTrue(os.path.exists(out_file))

            universe = read_data_file(out_file)
            self.assertEqual(universe.get_nr_of_atoms(), 40 * 10 + 20)

    def test_analyse_networks_cli(self) -> None:
        runner = CliRunner()
        help_res = runner.invoke(analyse_networks_cli, ["--help"])
        self.assertEqual(help_res.exit_code, 0)
        self.assertIn("--parameter-type", help_res.output)

        with tempfile.TemporaryDirectory() as tmpdir:
            net_file = os.path.join(tmpdir, "network.dat")
            gen_res = runner.invoke(
                generate_network_cli,
                [
                    "--polymer-name",
                    "PDMS",
                    "--parameter-type",
                    "GAUSSIAN",
                    "--target-p",
                    "0.8",
                    "--n-chains-crosslinkers",
                    "20",
                    "--n-chains-1",
                    "40",
                    "--n-beads-per-chain-1",
                    "10",
                    "--target-file",
                    net_file,
                ],
            )
            self.assertEqual(gen_res.exit_code, 0, msg=gen_res.output)

            res = runner.invoke(
                analyse_networks_cli,
                [
                    "--polymer-name",
                    "PDMS",
                    "--parameter-type",
                    "GAUSSIAN",
                    net_file,
                ],
            )
            self.assertEqual(res.exit_code, 0, msg=res.output)
            self.assertIn("Processing 1 files", res.output)
            self.assertIn("Bond length b:", res.output)
            self.assertIn("Shear Modulus [MPa]:", res.output)
            self.assertIn("w_sol =", res.output)

    def test_basic_lammps_structure_stats_cli(self) -> None:
        runner = CliRunner()
        help_res = runner.invoke(basic_stats_cli, ["--help"])
        self.assertEqual(help_res.exit_code, 0)

        with tempfile.TemporaryDirectory() as tmpdir:
            net_file = os.path.join(tmpdir, "network.dat")
            gen_res = runner.invoke(
                generate_network_cli,
                [
                    "--polymer-name",
                    "PDMS",
                    "--parameter-type",
                    "GAUSSIAN",
                    "--target-p",
                    "0.8",
                    "--n-chains-crosslinkers",
                    "20",
                    "--n-chains-1",
                    "40",
                    "--n-beads-per-chain-1",
                    "10",
                    "--target-file",
                    net_file,
                ],
            )
            self.assertEqual(gen_res.exit_code, 0, msg=gen_res.output)

            res = runner.invoke(basic_stats_cli, [net_file])
            self.assertEqual(res.exit_code, 0, msg=res.output)
            self.assertIn("Processing 1 files", res.output)
            self.assertIn("Mean bond length:", res.output)
            self.assertIn("Mean end to end distance:", res.output)

    def test_displace_randomly_cli(self) -> None:
        runner = CliRunner()
        help_res = runner.invoke(displace_randomly_cli, ["--help"])
        self.assertEqual(help_res.exit_code, 0)

        with tempfile.TemporaryDirectory() as tmpdir:
            net_file = os.path.join(tmpdir, "network.dat")
            gen_res = runner.invoke(
                generate_network_cli,
                [
                    "--polymer-name",
                    "PDMS",
                    "--parameter-type",
                    "GAUSSIAN",
                    "--target-p",
                    "0.8",
                    "--n-chains-crosslinkers",
                    "10",
                    "--n-chains-1",
                    "20",
                    "--n-beads-per-chain-1",
                    "5",
                    "--target-file",
                    net_file,
                ],
            )
            self.assertEqual(gen_res.exit_code, 0, msg=gen_res.output)

            res = runner.invoke(displace_randomly_cli, [net_file, "0.25"])
            self.assertEqual(res.exit_code, 0, msg=res.output)

            displaced_file = os.path.join(tmpdir, "random-displaced-network.dat")
            self.assertTrue(os.path.exists(displaced_file))
            displaced_universe = read_data_file(displaced_file)
            self.assertEqual(displaced_universe.get_nr_of_atoms(), 20 * 5 + 10)


if __name__ == "__main__":
    unittest.main()
