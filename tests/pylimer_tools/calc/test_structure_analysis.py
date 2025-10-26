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
from pint import UnitRegistry

from pylimer_tools.calc.structure_analysis import (
    compute_crosslinker_conversion,
    compute_effective_crosslinker_functionalities,
    compute_effective_crosslinker_functionality,
    compute_mean_end_to_end_distances, compute_mean_end_to_end_vectors,
    compute_strand_number_density, measure_weight_fraction_of_backbone,
    measure_weight_fraction_of_dangling_chains,
    measure_weight_fraction_of_soluble_material)
from pylimer_tools_cpp import MoleculeType

if __name__ == "__main__":
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../"))

from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestStrandNumberDensity(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures."""
        self.ureg = UnitRegistry()

    def test_compute_strand_number_density_bifunctional_only(self):
        """Test strand number density with only bifunctional strands (b2=1.0)."""
        mw_bifunctional = 1000 * self.ureg("g/mol")
        density = 1.2 * self.ureg("g/cm^3")
        
        result = compute_strand_number_density(
            mw_bifunctional=mw_bifunctional,
            density=density,
            mw_monofunctional=None,
            b2=1.0
        )
        
        # Expected: density / mw_bifunctional = 1.2 g/cm^3 / 1000 g/mol
        expected = density / mw_bifunctional
        self.assertAlmostEqual(
            result.to("mol/cm^3").magnitude,
            expected.to("mol/cm^3").magnitude,
            places=10
        )

    def test_compute_strand_number_density_bifunctional_only_default_b2(self):
        """Test strand number density with only bifunctional strands (default b2)."""
        mw_bifunctional = 2000 * self.ureg("g/mol")
        density = 1.0 * self.ureg("g/cm^3")
        
        result = compute_strand_number_density(
            mw_bifunctional=mw_bifunctional,
            density=density
        )
        
        expected = density / mw_bifunctional
        self.assertAlmostEqual(
            result.to("mol/cm^3").magnitude,
            expected.to("mol/cm^3").magnitude,
            places=10
        )

    def test_compute_strand_number_density_mixed_strands_b2_equals_1(self):
        """Test that b2=1.0 with monofunctional strands gives same result as bifunctional only."""
        mw_bifunctional = 1500 * self.ureg("g/mol")
        mw_monofunctional = 800 * self.ureg("g/mol")
        density = 1.1 * self.ureg("g/cm^3")
        
        result_bifunctional_only = compute_strand_number_density(
            mw_bifunctional=mw_bifunctional,
            density=density,
            mw_monofunctional=None,
            b2=1.0
        )
        
        result_with_mono = compute_strand_number_density(
            mw_bifunctional=mw_bifunctional,
            density=density,
            mw_monofunctional=mw_monofunctional,
            b2=1.0
        )
        
        # When b2=1.0, all reactive sites are in B2 (bifunctional)
        # So monofunctional fraction should be 0, making results equal
        self.assertAlmostEqual(
            result_bifunctional_only.to("mol/cm^3").magnitude,
            result_with_mono.to("mol/cm^3").magnitude,
            places=10
        )

    def test_compute_strand_number_density_mixed_strands_b2_equals_0(self):
        """Test strand number density with b2=0 (all monofunctional)."""
        mw_bifunctional = 1500 * self.ureg("g/mol")
        mw_monofunctional = 800 * self.ureg("g/mol")
        density = 1.1 * self.ureg("g/cm^3")
        
        result = compute_strand_number_density(
            mw_bifunctional=mw_bifunctional,
            density=density,
            mw_monofunctional=mw_monofunctional,
            b2=0.0
        )
        
        # When b2=0, all reactive sites are in B1 (monofunctional)
        # So bifunctional fraction should be 0
        expected = density / mw_monofunctional
        self.assertAlmostEqual(
            result.to("mol/cm^3").magnitude,
            expected.to("mol/cm^3").magnitude,
            places=10
        )

    def test_compute_strand_number_density_mixed_strands_b2_half(self):
        """Test strand number density with b2=0.5 (equal reactive sites)."""
        mw_bifunctional = 2000 * self.ureg("g/mol")
        mw_monofunctional = 1000 * self.ureg("g/mol")
        density = 1.0 * self.ureg("g/cm^3")
        
        result = compute_strand_number_density(
            mw_bifunctional=mw_bifunctional,
            density=density,
            mw_monofunctional=mw_monofunctional,
            b2=0.5
        )
        
        # When b2=0.5: (2 * n_bifunctional) / (n_monofunctional + 2 * n_bifunctional) = 0.5
        # This means: 2 * n_bifunctional = 0.5 * (n_monofunctional + 2 * n_bifunctional)
        # Solving: n_bifunctional = 0.5 * n_monofunctional
        # If we have 1000 total strands, n_bifunctional = 333.33..., n_monofunctional = 666.66...
        # bifunctional_fraction = 1/3, monofunctional_fraction = 2/3
        expected_num_density = (
            (1/3) * density / mw_bifunctional +
            (2/3) * density / mw_monofunctional
        )
        
        self.assertAlmostEqual(
            result.to("mol/cm^3").magnitude,
            expected_num_density.to("mol/cm^3").magnitude,
            places=10
        )

    def test_compute_strand_number_density_different_b2_values(self):
        """Test strand number density with various b2 values."""
        mw_bifunctional = 1800 * self.ureg("g/mol")
        mw_monofunctional = 900 * self.ureg("g/mol")
        density = 1.15 * self.ureg("g/cm^3")
        
        b2_values = [0.1, 0.25, 0.5, 0.75, 0.9]
        
        for b2 in b2_values:
            with self.subTest(b2=b2):
                result = compute_strand_number_density(
                    mw_bifunctional=mw_bifunctional,
                    density=density,
                    mw_monofunctional=mw_monofunctional,
                    b2=b2
                )
                
                # Result should be a valid pint.Quantity with proper units
                self.assertTrue(hasattr(result, 'magnitude'))
                self.assertTrue(hasattr(result, 'units'))
                
                # Result should be positive
                self.assertGreater(result.magnitude, 0)

    def test_compute_strand_number_density_assertion_no_mono_with_b2_not_1(self):
        """Test that assertion fails when no monofunctional strands but b2 != 1.0."""
        mw_bifunctional = 1000 * self.ureg("g/mol")
        density = 1.2 * self.ureg("g/cm^3")
        
        with self.assertRaises(AssertionError) as context:
            compute_strand_number_density(
                mw_bifunctional=mw_bifunctional,
                density=density,
                mw_monofunctional=None,
                b2=0.8
            )
        
        self.assertIn("b2 must be 1.0", str(context.exception))

    def test_compute_strand_number_density_assertion_b2_out_of_range_low(self):
        """Test that assertion fails when b2 < 0."""
        mw_bifunctional = 1000 * self.ureg("g/mol")
        mw_monofunctional = 500 * self.ureg("g/mol")
        density = 1.2 * self.ureg("g/cm^3")
        
        with self.assertRaises(AssertionError) as context:
            compute_strand_number_density(
                mw_bifunctional=mw_bifunctional,
                density=density,
                mw_monofunctional=mw_monofunctional,
                b2=-0.1
            )
        
        self.assertIn("b2 must be between 0 and 1", str(context.exception))

    def test_compute_strand_number_density_assertion_b2_out_of_range_high(self):
        """Test that assertion fails when b2 > 1."""
        mw_bifunctional = 1000 * self.ureg("g/mol")
        mw_monofunctional = 500 * self.ureg("g/mol")
        density = 1.2 * self.ureg("g/cm^3")
        
        with self.assertRaises(AssertionError) as context:
            compute_strand_number_density(
                mw_bifunctional=mw_bifunctional,
                density=density,
                mw_monofunctional=mw_monofunctional,
                b2=1.1
            )
        
        self.assertIn("b2 must be between 0 and 1", str(context.exception))

    def test_compute_strand_number_density_unit_consistency(self):
        """Test that different but compatible units work correctly."""
        # Test with different density units
        mw_bifunctional = 1000 * self.ureg("g/mol")
        density_g_cm3 = 1.2 * self.ureg("g/cm^3")
        density_kg_m3 = 1200 * self.ureg("kg/m^3")
        
        result1 = compute_strand_number_density(
            mw_bifunctional=mw_bifunctional,
            density=density_g_cm3
        )
        
        result2 = compute_strand_number_density(
            mw_bifunctional=mw_bifunctional,
            density=density_kg_m3
        )
        
        # Both should give the same result when converted to the same units
        self.assertAlmostEqual(
            result1.to("mol/cm^3").magnitude,
            result2.to("mol/cm^3").magnitude,
            places=10
        )

    def test_compute_strand_number_density_realistic_polymer_values(self):
        """Test with realistic polymer network values."""
        # Typical PDMS network parameters
        mw_bifunctional = 5000 * self.ureg("g/mol")  # 5 kg/mol
        mw_monofunctional = 2500 * self.ureg("g/mol")  # 2.5 kg/mol
        density = 0.97 * self.ureg("g/cm^3")  # Typical PDMS density
        b2 = 0.8  # 80% of reactive sites are bifunctional
        
        result = compute_strand_number_density(
            mw_bifunctional=mw_bifunctional,
            density=density,
            mw_monofunctional=mw_monofunctional,
            b2=b2
        )
        
        # Check that result is in a reasonable range for polymer networks
        # Typically between 10^-5 and 10^-3 mol/cm^3
        result_mol_cm3 = result.to("mol/cm^3").magnitude
        self.assertGreater(result_mol_cm3, 1e-6)
        self.assertLess(result_mol_cm3, 1e-2)

    def test_compute_strand_number_density_edge_case_b2_boundary(self):
        """Test boundary conditions for b2 (exactly 0 and exactly 1)."""
        mw_bifunctional = 1000 * self.ureg("g/mol")
        mw_monofunctional = 500 * self.ureg("g/mol")
        density = 1.0 * self.ureg("g/cm^3")
        
        # Test b2 = 0.0 (boundary)
        result_b2_0 = compute_strand_number_density(
            mw_bifunctional=mw_bifunctional,
            density=density,
            mw_monofunctional=mw_monofunctional,
            b2=0.0
        )
        
        # Test b2 = 1.0 (boundary)
        result_b2_1 = compute_strand_number_density(
            mw_bifunctional=mw_bifunctional,
            density=density,
            mw_monofunctional=mw_monofunctional,
            b2=1.0
        )
        
        # Both should produce valid results
        self.assertIsNotNone(result_b2_0)
        self.assertIsNotNone(result_b2_1)
        self.assertGreater(result_b2_0.magnitude, 0)
        self.assertGreater(result_b2_1.magnitude, 0)


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
