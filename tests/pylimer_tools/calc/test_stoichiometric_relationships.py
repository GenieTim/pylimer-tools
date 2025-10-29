"""
Tests for stoichiometric_relationships module.

This module contains tests for functions in pylimer_tools.calc.stoichiometric_relationships
that calculate stoichiometric relationships in polymer networks, including:
- Number fractions of bifunctional and monofunctional strands
- Weight fractions of bifunctional and monofunctional strands
- Strand number density calculations
"""
import os
import sys
import unittest

from pint import UnitRegistry

if __name__ == "__main__":
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../"))

from pylimer_tools.calc.stoichiometric_relationships import (
    compute_number_fractions,
    compute_strand_number_density,
    compute_weight_fractions,
)


class TestComputeNumberFractions(unittest.TestCase):
    """Test the compute_number_fractions function."""

    def test_compute_number_fractions_b2_equals_1(self):
        """Test number fractions when b2=1.0 (all bifunctional)."""
        bifunctional_fraction, monofunctional_fraction = compute_number_fractions(b2=1.0)
        
        self.assertAlmostEqual(bifunctional_fraction, 1.0, places=10)
        self.assertAlmostEqual(monofunctional_fraction, 0.0, places=10)

    def test_compute_number_fractions_b2_equals_0(self):
        """Test number fractions when b2=0.0 (all monofunctional)."""
        bifunctional_fraction, monofunctional_fraction = compute_number_fractions(b2=0.0)
        
        self.assertAlmostEqual(bifunctional_fraction, 0.0, places=10)
        self.assertAlmostEqual(monofunctional_fraction, 1.0, places=10)

    def test_compute_number_fractions_b2_equals_half(self):
        """Test number fractions when b2=0.5."""
        # When b2=0.5: (2 * n_bifunctional) / (n_monofunctional + 2 * n_bifunctional) = 0.5
        # This means: 2 * n_bifunctional = 0.5 * (n_monofunctional + 2 * n_bifunctional)
        # Solving: n_bifunctional = 0.5 * n_monofunctional
        # bifunctional_fraction = 1/3, monofunctional_fraction = 2/3
        bifunctional_fraction, monofunctional_fraction = compute_number_fractions(b2=0.5)
        
        self.assertAlmostEqual(bifunctional_fraction, 1/3, places=10)
        self.assertAlmostEqual(monofunctional_fraction, 2/3, places=10)

    def test_compute_number_fractions_various_b2_values(self):
        """Test number fractions with various b2 values."""
        b2_values = [0.1, 0.25, 0.5, 0.75, 0.9]
        
        for b2 in b2_values:
            with self.subTest(b2=b2):
                bifunctional_fraction, monofunctional_fraction = compute_number_fractions(b2=b2)
                
                # Check that fractions sum to 1
                self.assertAlmostEqual(
                    bifunctional_fraction + monofunctional_fraction,
                    1.0,
                    places=10,
                    msg=f"Fractions should sum to 1 for b2={b2}"
                )
                
                # Check that fractions are non-negative
                self.assertGreaterEqual(bifunctional_fraction, 0.0)
                self.assertGreaterEqual(monofunctional_fraction, 0.0)
                
                # Check that fractions are at most 1
                self.assertLessEqual(bifunctional_fraction, 1.0)
                self.assertLessEqual(monofunctional_fraction, 1.0)

    def test_compute_number_fractions_b2_relationship(self):
        """Test that the b2 relationship holds: b2 = (2*n_bi)/(n_mono + 2*n_bi)."""
        b2_values = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        
        for b2 in b2_values:
            with self.subTest(b2=b2):
                bifunctional_fraction, monofunctional_fraction = compute_number_fractions(b2=b2)
                
                # Verify the relationship: b2 = (2 * x_bi) / (x_mono + 2 * x_bi)
                if monofunctional_fraction + 2 * bifunctional_fraction > 0:
                    computed_b2 = (2 * bifunctional_fraction) / (
                        monofunctional_fraction + 2 * bifunctional_fraction
                    )
                    self.assertAlmostEqual(
                        computed_b2,
                        b2,
                        places=10,
                        msg=f"b2 relationship should hold for b2={b2}"
                    )

    def test_compute_number_fractions_assertion_b2_too_low(self):
        """Test that assertion fails when b2 < 0."""
        with self.assertRaises(AssertionError) as context:
            compute_number_fractions(b2=-0.1)
        
        self.assertIn("b2 must be between 0 and 1", str(context.exception))

    def test_compute_number_fractions_assertion_b2_too_high(self):
        """Test that assertion fails when b2 > 1."""
        with self.assertRaises(AssertionError) as context:
            compute_number_fractions(b2=1.1)
        
        self.assertIn("b2 must be between 0 and 1", str(context.exception))


class TestComputeWeightFractions(unittest.TestCase):
    """Test the compute_weight_fractions function."""

    def setUp(self):
        """Set up test fixtures."""
        self.ureg = UnitRegistry()

    def test_compute_weight_fractions_no_monofunctional(self):
        """Test weight fractions when no monofunctional strands are present."""
        mw_bifunctional = 1000 * self.ureg("g/mol")
        
        w_bifunctional, w_monofunctional = compute_weight_fractions(
            mw_bifunctional=mw_bifunctional,
            mw_monofunctional=None,
            b2=1.0
        )
        
        self.assertAlmostEqual(w_bifunctional, 1.0, places=10)
        self.assertAlmostEqual(w_monofunctional, 0.0, places=10)

    def test_compute_weight_fractions_no_monofunctional_infinite_mw(self):
        """Test weight fractions when monofunctional MW is infinite."""
        mw_bifunctional = 1000 * self.ureg("g/mol")
        mw_monofunctional = float('inf') * self.ureg("g/mol")
        
        w_bifunctional, w_monofunctional = compute_weight_fractions(
            mw_bifunctional=mw_bifunctional,
            mw_monofunctional=mw_monofunctional,
            b2=1.0
        )
        
        self.assertAlmostEqual(w_bifunctional, 1.0, places=10)
        self.assertAlmostEqual(w_monofunctional, 0.0, places=10)

    def test_compute_weight_fractions_b2_equals_1(self):
        """Test weight fractions when b2=1.0 (all bifunctional)."""
        mw_bifunctional = 2000 * self.ureg("g/mol")
        mw_monofunctional = 1000 * self.ureg("g/mol")
        
        w_bifunctional, w_monofunctional = compute_weight_fractions(
            mw_bifunctional=mw_bifunctional,
            mw_monofunctional=mw_monofunctional,
            b2=1.0
        )
        
        # When b2=1.0, all strands are bifunctional (number-wise)
        self.assertAlmostEqual(w_bifunctional, 1.0, places=10)
        self.assertAlmostEqual(w_monofunctional, 0.0, places=10)

    def test_compute_weight_fractions_b2_equals_0(self):
        """Test weight fractions when b2=0.0 (all monofunctional)."""
        mw_bifunctional = 2000 * self.ureg("g/mol")
        mw_monofunctional = 1000 * self.ureg("g/mol")
        
        w_bifunctional, w_monofunctional = compute_weight_fractions(
            mw_bifunctional=mw_bifunctional,
            mw_monofunctional=mw_monofunctional,
            b2=0.0
        )
        
        # When b2=0.0, all strands are monofunctional (number-wise)
        self.assertAlmostEqual(w_bifunctional, 0.0, places=10)
        self.assertAlmostEqual(w_monofunctional, 1.0, places=10)

    def test_compute_weight_fractions_equal_molecular_weights(self):
        """Test weight fractions when molecular weights are equal."""
        mw = 1500 * self.ureg("g/mol")
        
        for b2 in [0.0, 0.5, 1.0]:
            with self.subTest(b2=b2):
                w_bifunctional, w_monofunctional = compute_weight_fractions(
                    mw_bifunctional=mw,
                    mw_monofunctional=mw,
                    b2=b2
                )
                
                # When molecular weights are equal, weight fractions should equal number fractions
                bifunctional_fraction, monofunctional_fraction = compute_number_fractions(b2=b2)
                
                self.assertAlmostEqual(w_bifunctional, bifunctional_fraction, places=10)
                self.assertAlmostEqual(w_monofunctional, monofunctional_fraction, places=10)

    def test_compute_weight_fractions_sum_to_one(self):
        """Test that weight fractions always sum to 1."""
        mw_bifunctional = 2000 * self.ureg("g/mol")
        mw_monofunctional = 1000 * self.ureg("g/mol")
        b2_values = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        
        for b2 in b2_values:
            with self.subTest(b2=b2):
                w_bifunctional, w_monofunctional = compute_weight_fractions(
                    mw_bifunctional=mw_bifunctional,
                    mw_monofunctional=mw_monofunctional,
                    b2=b2
                )
                
                self.assertAlmostEqual(
                    w_bifunctional + w_monofunctional,
                    1.0,
                    places=10,
                    msg=f"Weight fractions should sum to 1 for b2={b2}"
                )

    def test_compute_weight_fractions_heavier_bifunctional(self):
        """Test weight fractions when bifunctional chains are heavier."""
        mw_bifunctional = 3000 * self.ureg("g/mol")
        mw_monofunctional = 1000 * self.ureg("g/mol")
        
        w_bifunctional, w_monofunctional = compute_weight_fractions(
            mw_bifunctional=mw_bifunctional,
            mw_monofunctional=mw_monofunctional,
            b2=0.5
        )
        
        # For b2=0.5, number fractions are 1/3 bifunctional, 2/3 monofunctional
        # Weight: (1/3 * 3000) + (2/3 * 1000) = 1000 + 666.67 = 1666.67
        # w_bifunctional = 1000 / 1666.67 = 0.6
        # w_monofunctional = 666.67 / 1666.67 = 0.4
        expected_w_bifunctional = 1000 / 1666.666667
        expected_w_monofunctional = 666.666667 / 1666.666667
        
        self.assertAlmostEqual(w_bifunctional, expected_w_bifunctional, places=5)
        self.assertAlmostEqual(w_monofunctional, expected_w_monofunctional, places=5)

    def test_compute_weight_fractions_realistic_values(self):
        """Test with realistic polymer values."""
        # PDMS-like values
        mw_bifunctional = 5000 * self.ureg("g/mol")
        mw_monofunctional = 2500 * self.ureg("g/mol")
        
        w_bifunctional, w_monofunctional = compute_weight_fractions(
            mw_bifunctional=mw_bifunctional,
            mw_monofunctional=mw_monofunctional,
            b2=0.8
        )
        
        # Check that values are reasonable
        self.assertGreaterEqual(w_bifunctional, 0.0)
        self.assertLessEqual(w_bifunctional, 1.0)
        self.assertGreaterEqual(w_monofunctional, 0.0)
        self.assertLessEqual(w_monofunctional, 1.0)
        
        # Since b2 is high (0.8), we expect more bifunctional strands
        self.assertGreater(w_bifunctional, w_monofunctional)

    def test_compute_weight_fractions_assertion_no_mono_with_b2_not_1(self):
        """Test that assertion fails when no monofunctional strands but b2 != 1.0."""
        mw_bifunctional = 1000 * self.ureg("g/mol")
        
        with self.assertRaises(AssertionError) as context:
            compute_weight_fractions(
                mw_bifunctional=mw_bifunctional,
                mw_monofunctional=None,
                b2=0.8
            )
        
        self.assertIn("b2 must be 1.0", str(context.exception))

    def test_compute_weight_fractions_assertion_b2_out_of_range_low(self):
        """Test that assertion fails when b2 < 0."""
        mw_bifunctional = 1000 * self.ureg("g/mol")
        mw_monofunctional = 500 * self.ureg("g/mol")
        
        with self.assertRaises(AssertionError) as context:
            compute_weight_fractions(
                mw_bifunctional=mw_bifunctional,
                mw_monofunctional=mw_monofunctional,
                b2=-0.1
            )
        
        self.assertIn("b2 must be between 0 and 1", str(context.exception))

    def test_compute_weight_fractions_assertion_b2_out_of_range_high(self):
        """Test that assertion fails when b2 > 1."""
        mw_bifunctional = 1000 * self.ureg("g/mol")
        mw_monofunctional = 500 * self.ureg("g/mol")
        
        with self.assertRaises(AssertionError) as context:
            compute_weight_fractions(
                mw_bifunctional=mw_bifunctional,
                mw_monofunctional=mw_monofunctional,
                b2=1.1
            )
        
        self.assertIn("b2 must be between 0 and 1", str(context.exception))


class TestStrandNumberDensity(unittest.TestCase):
    """Test the compute_strand_number_density function."""

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


if __name__ == "__main__":
    unittest.main()
