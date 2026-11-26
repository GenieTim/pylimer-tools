"""
Tests for structure_analysis module.

This module contains tests for functions in pylimer_tools.calc.structure_analysis
that analyze polymer network structures, including:
- Strand number density calculations
- Crosslinker functionality and conversion
- Weight fractions of network components (backbone, dangling chains, soluble material)
- End-to-end distance and vector computations
"""

import os
import sys
import unittest

import numpy as np

from pylimer_tools.calc.structure_analysis import (
    compute_crosslinker_conversion,
    compute_effective_crosslinker_functionalities,
    compute_effective_crosslinker_functionality,
    compute_extent_of_reaction,
    compute_fraction_of_bifunctional_reactive_sites,
    compute_mean_end_to_end_distances,
    compute_mean_end_to_end_vectors,
    compute_stoichiometric_imbalance,
    determine_chain_functionality_per_type,
    measure_weight_fraction_of_backbone,
    measure_weight_fraction_of_dangling_chains,
    measure_weight_fraction_of_soluble_material,
)
from pylimer_tools_cpp import MoleculeType, Universe

if __name__ == "__main__":
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../"))

from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestStructureAnalysisWithUniverse(UniverseUsingTestCase):
    """Tests for structure_analysis functions that require Universe fixtures."""

    def test_weight_fraction_calculations(self):
        """Test weight fraction calculations for dangling chains, backbone, and soluble material."""
        self.assertEqual(
            (0.0, 0.0),
            measure_weight_fraction_of_dangling_chains(self.emptyUniverse, 2),
        )
        # empty weight -> empty weight fraction
        self.testUniverse.set_masses({1: 0, 2: 0})
        self.assertEqual(
            0.0, measure_weight_fraction_of_soluble_material(self.testUniverse, 2)
        )
        self.assertEqual(1.0, measure_weight_fraction_of_backbone(self.testUniverse, 2))
        # non-empty weights
        self.testUniverse.set_masses({1: 1, 2: 0})
        self.assertTrue(self.testUniverse.get_nr_of_atoms() > 0)
        self.assertEqual(self.testUniverse.get_masses(), {1: 1, 2: 0})
        all_chains = self.testUniverse.get_chains_with_crosslinker(2)
        self.assertEqual(all_chains[2].get_strand_type(), MoleculeType.DANGLING_CHAIN)
        self.assertEqual(
            (0.2, 0.25),
            measure_weight_fraction_of_dangling_chains(
                self.testUniverse, crosslinker_type=2
            ),
        )

    def test_crosslinker_functionality_calculation(self):
        """Test crosslinker functionality and conversion calculations."""
        self.assertCountEqual(
            [], compute_effective_crosslinker_functionalities(self.emptyUniverse, 2)
        )
        self.assertSequenceEqual(
            [0, 2, 3],
            compute_effective_crosslinker_functionalities(self.testUniverse, 2),
        )
        self.assertEqual(
            5.0 / 3.0, compute_effective_crosslinker_functionality(self.testUniverse, 2)
        )
        self.assertEqual(
            5.0 / 3.0 / 3.0, compute_crosslinker_conversion(self.testUniverse, 2, f=3)
        )
        self.assertRaises(
            ValueError,
            lambda: compute_crosslinker_conversion(self.testUniverse, 2, np.inf),  # type: ignore
        )

    def test_mean_end_to_end_computation(self):
        """Test end-to-end distance and vector calculations."""
        self.assertCountEqual([], compute_mean_end_to_end_vectors([], 2))
        self.assertDictEqual(
            {"1-2-3-6-7": 1.0, "5-6-7": 1.0},
            compute_mean_end_to_end_distances([self.testUniverse], 2),
        )

    def test_compute_fraction_of_bifunctional_reactive_sites(self):
        """Test b2 calculation with monofunctional chains, empty/crosslinker-only, and solvent."""
        # Empty and all-crosslinker networks should not raise ZeroDivisionError
        self.assertEqual(
            1.0, compute_fraction_of_bifunctional_reactive_sites(self.emptyUniverse, 2)
        )
        xlink_only = Universe(10, 10, 10)
        xlink_only.add_atoms(
            [1, 2], [2, 2], [0.0, 1.0], [0.0, 0.0], [0.0, 0.0], [0, 0], [0, 0], [0, 0]
        )
        self.assertEqual(
            1.0, compute_fraction_of_bifunctional_reactive_sites(xlink_only, 2)
        )

        # Low-conversion network with solvent (type 4) and bifunctional chain (type 1)
        low_conv = Universe(10, 10, 10)
        low_conv.add_atoms(
            [1, 2, 3, 4],
            [2, 1, 1, 4],
            [0.0, 1.0, 2.0, 5.0],
            [0.0, 0.0, 0.0, 5.0],
            [0.0, 0.0, 0.0, 5.0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
        )
        low_conv.add_bonds([1, 2], [2, 3])
        # With explicit functionality_per_type={1: 2, 2: 4, 4: 0}:
        # solvent (4) and crosslinker (2) must not be marked monofunctional
        self.assertEqual(
            1.0,
            compute_fraction_of_bifunctional_reactive_sites(
                low_conv, 2, functionality_per_type={1: 2, 2: 4, 4: 0}
            ),
        )

        # Network with 1 bifunctional chain (type 1) and 1 monofunctional chain (type 3)
        # b2 = (2 * 1) / (1 + 2 * 1) = 2/3
        mixed = Universe(10, 10, 10)
        mixed.add_atoms(
            [1, 2, 3, 4, 5, 6],
            [2, 1, 1, 2, 3, 3],
            [0.0, 1.0, 2.0, 3.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 2.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
        )
        mixed.add_bonds([1, 2, 3, 1, 5], [2, 3, 4, 5, 6])
        self.assertAlmostEqual(
            2.0 / 3.0,
            compute_fraction_of_bifunctional_reactive_sites(
                mixed, 2, functionality_per_type={1: 2, 2: 4, 3: 1}
            ),
        )

    def test_stoichiometric_imbalance_does_not_mutate_ignore_types(self):
        """Test that compute_stoichiometric_imbalance does not mutate ignore_types."""
        ignore = [3]
        compute_stoichiometric_imbalance(
            self.testUniverse, crosslinker_type=2, ignore_types=ignore
        )
        self.assertEqual([3], ignore)

    def test_measure_weight_fraction_of_dangling_chains_unset_masses(self):
        """Test measure_weight_fraction_of_dangling_chains with unset masses."""
        universe = Universe(10, 10, 10)
        universe.add_atoms(
            [1, 2, 3, 4, 5],
            [2, 1, 2, 1, 1],
            [0.0, 1.0, 2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 2.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0],
        )
        universe.add_bonds([1, 2, 1, 4], [2, 3, 4, 5])
        # Masses are unset; default mass 1.0 per atom => dangling chain [1, 4, 5] has 3 atoms out of 5
        w_dang, n_dang = measure_weight_fraction_of_dangling_chains(universe, 2)
        self.assertAlmostEqual(3.0 / 5.0, w_dang)
        self.assertAlmostEqual(3.0 / 5.0, n_dang)

    def test_compute_mean_end_to_end_vectors_missing_keys(self):
        """Test that missing molecules in a subset of networks raises ValueError."""
        with self.assertRaises(ValueError):
            compute_mean_end_to_end_vectors(
                [self.testUniverse, self.testUniverse, self.saturatedTestUniverse],
                crosslinker_type=2,
            )

    def test_determine_chain_functionality_per_type(self):
        """Test chain functionality inference for B2, B1, B0, primary loops, and unreacted networks."""
        self.assertEqual(
            {}, determine_chain_functionality_per_type(self.emptyUniverse, 2)
        )

        # Network with A4 crosslinkers (2), multi-bead B2 (1), multi-bead B1 (3), multi-bead B0 (4)
        net = Universe(10, 10, 10)
        # xlinks: 1, 2; B2: 10-11-12 (primary loop on 1), 13-14-15 (between 1 & 2)
        # B1: 20-21-22 (attached to 2); B0: 30-31-32 (free solvent)
        ids = [1, 2, 10, 11, 12, 13, 14, 15, 20, 21, 22, 30, 31, 32]
        types = [2, 2, 1, 1, 1, 1, 1, 1, 3, 3, 3, 4, 4, 4]
        n = len(ids)
        net.add_atoms(
            ids,
            types,
            [float(i) for i in range(n)],
            [0.0] * n,
            [0.0] * n,
            [0] * n,
            [0] * n,
            [0] * n,
        )
        net.add_bonds(
            [10, 11, 13, 14, 20, 21, 30, 31, 1, 1, 1, 2, 2],
            [11, 12, 14, 15, 21, 22, 31, 32, 10, 12, 13, 15, 20],
        )
        inferred = determine_chain_functionality_per_type(net, crosslinker_type=2)
        self.assertEqual(inferred, {1: 2, 2: 3, 3: 1, 4: 0})

        # Primary-loop-only B2 network: single xlink (1) with a 3-bead primary loop (10-11-12)
        loop_net = Universe(10, 10, 10)
        loop_net.add_atoms(
            [1, 10, 11, 12],
            [2, 1, 1, 1],
            [0.0, 1.0, 2.0, 3.0],
            [0.0] * 4,
            [0.0] * 4,
            [0] * 4,
            [0] * 4,
            [0] * 4,
        )
        loop_net.add_bonds([10, 11, 1, 1], [11, 12, 10, 12])
        self.assertEqual(
            determine_chain_functionality_per_type(loop_net, crosslinker_type=2),
            {1: 2, 2: 2},
        )

        # Unreacted network (0 crosslinker bonds): non-crosslinker types default to 2,
        # and conversion / extent of reaction must be 0.0 (not 1.0 or ValueError)
        unreacted = Universe(10, 10, 10)
        unreacted.add_atoms(
            [1, 10, 11, 12],
            [2, 1, 1, 1],
            [0.0, 1.0, 2.0, 3.0],
            [0.0] * 4,
            [0.0] * 4,
            [0] * 4,
            [0] * 4,
            [0] * 4,
        )
        unreacted.add_bonds([10, 11], [11, 12])
        self.assertEqual(
            determine_chain_functionality_per_type(unreacted, crosslinker_type=2),
            {1: 2, 2: 0},
        )
        self.assertAlmostEqual(
            0.0, compute_extent_of_reaction(unreacted, crosslinker_type=2)
        )
        self.assertAlmostEqual(
            0.0, compute_crosslinker_conversion(unreacted, crosslinker_type=2)
        )

        # Direct crosslinker-crosslinker bond (1-2) with unreacted B2 chain (10-11-12):
        # unreacted type-1 chain must still default to 2 (not 0)
        direct_xlink = Universe(10, 10, 10)
        direct_xlink.add_atoms(
            [1, 2, 10, 11, 12],
            [2, 2, 1, 1, 1],
            [0.0, 0.5, 1.0, 2.0, 3.0],
            [0.0] * 5,
            [0.0] * 5,
            [0] * 5,
            [0] * 5,
            [0] * 5,
        )
        direct_xlink.add_bonds([1, 10, 11], [2, 11, 12])
        self.assertEqual(
            determine_chain_functionality_per_type(direct_xlink, crosslinker_type=2),
            {1: 2, 2: 1},
        )

    def test_stoichiometric_imbalance_and_extent_of_reaction_with_b1_and_b0(self):
        """Test r, p, and b2 when monofunctional (B1) and solvent (B0) chains are present."""
        # 2 crosslinkers (type 2, target f=4 => 8 A sites)
        # 1 B2 chain (type 1, 3 beads => 2 B sites)
        # 2 B1 chains (type 3, 3 beads each => 2 B sites)
        # 2 B0 solvent chains (type 4, 3 beads each => 0 B sites)
        # Total B sites = 2 + 2 + 0 = 4
        # With f=4: r = 8 / 4 = 2.0; max_formable_bonds = min(4, 8) = 4
        # Formed bonds: B2 bonded to xlink 1 & xlink 2 (2 bonds),
        #               both B1 bonded to xlink 1 (2 bonds) => 4 bonds total
        # Crosslinker max degree observed: xlink 1 has 3 bonds, xlink 2 has 1 bond => inferred f = 3
        net = Universe(10, 10, 10)
        ids = (
            [1, 2] + [10, 11, 12] + [20, 21, 22, 23, 24, 25] + [30, 31, 32, 33, 34, 35]
        )
        types = [2, 2] + [1] * 3 + [3] * 6 + [4] * 6
        n = len(ids)
        net.add_atoms(
            ids,
            types,
            [float(i) for i in range(n)],
            [0.0] * n,
            [0.0] * n,
            [0] * n,
            [0] * n,
            [0] * n,
        )
        bonds_from = [10, 11, 20, 21, 23, 24, 30, 31, 33, 34, 1, 2, 1, 1]
        bonds_to = [11, 12, 21, 22, 24, 25, 31, 32, 34, 35, 10, 12, 20, 23]
        net.add_bonds(bonds_from, bonds_to)

        explicit_func = {1: 2, 2: 4, 3: 1, 4: 0}

        # With explicit functionality_per_type
        self.assertAlmostEqual(
            2.0,
            compute_stoichiometric_imbalance(
                net, crosslinker_type=2, functionality_per_type=explicit_func
            ),
        )
        self.assertAlmostEqual(
            1.0,
            compute_extent_of_reaction(
                net, crosslinker_type=2, functionality_per_type=explicit_func
            ),
        )
        self.assertAlmostEqual(
            0.5,
            compute_fraction_of_bifunctional_reactive_sites(
                net, crosslinker_type=2, functionality_per_type=explicit_func
            ),
        )

        # Without explicit functionality_per_type (inferred f=3 => 2*3=6 A sites, 4 B sites)
        self.assertAlmostEqual(
            6.0 / 4.0,
            compute_stoichiometric_imbalance(net, crosslinker_type=2),
        )
        self.assertAlmostEqual(
            1.0,
            compute_extent_of_reaction(net, crosslinker_type=2),
        )
        self.assertAlmostEqual(
            0.5,
            compute_fraction_of_bifunctional_reactive_sites(net, crosslinker_type=2),
        )

        # Partial functionality_per_type={2: 4} still infers B1=1 and B0=0
        self.assertAlmostEqual(
            2.0,
            compute_stoichiometric_imbalance(
                net, crosslinker_type=2, functionality_per_type={2: 4}
            ),
        )
        self.assertAlmostEqual(
            1.0,
            compute_extent_of_reaction(
                net, crosslinker_type=2, functionality_per_type={2: 4}
            ),
        )

        # ignore_types on compute_extent_of_reaction (including ignoring bonded B1 type 3)
        self.assertAlmostEqual(
            1.0,
            compute_extent_of_reaction(
                net, crosslinker_type=2, functionality_per_type={2: 4}, ignore_types=[4]
            ),
        )
        self.assertAlmostEqual(
            1.0,
            compute_extent_of_reaction(
                net, crosslinker_type=2, functionality_per_type={2: 4}, ignore_types=[3]
            ),
        )

        # Omitting crosslinker_type from functionality_per_type warns but preserves chain entries
        with self.assertWarns(UserWarning):
            r_warn = compute_stoichiometric_imbalance(
                net,
                crosslinker_type=2,
                functionality_per_type={1: 2, 3: 1, 4: 0},
            )
        self.assertAlmostEqual(6.0 / 4.0, r_warn)


if __name__ == "__main__":
    unittest.main()
