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
    compute_fraction_of_bifunctional_reactive_sites,
    compute_mean_end_to_end_distances,
    compute_mean_end_to_end_vectors,
    compute_stoichiometric_imbalance,
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


if __name__ == "__main__":
    unittest.main()
