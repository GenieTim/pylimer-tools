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
    compute_mean_end_to_end_distances, compute_mean_end_to_end_vectors,
    measure_weight_fraction_of_backbone,
    measure_weight_fraction_of_dangling_chains,
    measure_weight_fraction_of_soluble_material)
from pylimer_tools_cpp import MoleculeType

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
            0.0, measure_weight_fraction_of_soluble_material(
                self.testUniverse, 2)
        )
        self.assertEqual(
            1.0, measure_weight_fraction_of_backbone(
                self.testUniverse, 2))
        # non-empty weights
        self.testUniverse.set_masses({1: 1, 2: 0})
        self.assertTrue(self.testUniverse.get_nr_of_atoms() > 0)
        self.assertEqual(self.testUniverse.get_masses(), {1: 1, 2: 0})
        all_chains = self.testUniverse.get_chains_with_crosslinker(2)
        self.assertEqual(
            all_chains[2].get_strand_type(),
            MoleculeType.DANGLING_CHAIN)
        self.assertEqual(
            (0.2, 0.25),
            measure_weight_fraction_of_dangling_chains(
                self.testUniverse, crosslinker_type=2
            ),
        )

    def test_crosslinker_functionality_calculation(self):
        """Test crosslinker functionality and conversion calculations."""
        self.assertCountEqual(
            [], compute_effective_crosslinker_functionalities(
                self.emptyUniverse, 2)
        )
        self.assertSequenceEqual(
            [0, 2, 3],
            compute_effective_crosslinker_functionalities(
                self.testUniverse, 2),
        )
        self.assertEqual(
            5.0 /
            3.0, compute_effective_crosslinker_functionality(
                self.testUniverse, 2)
        )
        self.assertEqual(
            5.0 / 3.0 /
            3.0, compute_crosslinker_conversion(self.testUniverse, 2, f=3)
        )
        self.assertRaises(
            ValueError,
            lambda: compute_crosslinker_conversion(
                self.testUniverse, 2, np.inf),  # type: ignore
        )

    def test_mean_end_to_end_computation(self):
        """Test end-to-end distance and vector calculations."""
        self.assertCountEqual([], compute_mean_end_to_end_vectors([], 2))
        self.assertDictEqual(
            {"1-2-3-6-7": 1.0, "5-6-7": 1.0},
            compute_mean_end_to_end_distances([self.testUniverse], 2),
        )


if __name__ == "__main__":
    unittest.main()
