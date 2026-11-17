#!/usr/bin/env python

import unittest

from pylimer_tools.io.bead_spring_parameter_provider import (
    ParameterType,
    Parameters,
    get_parameters_for_polymer,
    get_supported_polymer_names,
)


class TestBeadSpringParameterProvider(unittest.TestCase):
    def test_different_parameter_types_returned_without_issue(self):
        """Test that different parameter types are returned without errors."""
        polymer_names = get_supported_polymer_names()
        self.assertGreater(len(polymer_names), 0)

        # Use the first polymer for testing
        polymer_name = polymer_names[0]

        # Test all parameter types
        gaussian_params = get_parameters_for_polymer(
            polymer_name, ParameterType.GAUSSIAN
        )
        kg_lj_params = get_parameters_for_polymer(polymer_name, ParameterType.KG_LJ)
        kuhn_params = get_parameters_for_polymer(polymer_name, ParameterType.KUHN)

        # Assert they are all Parameters objects
        self.assertIsNotNone(gaussian_params)
        self.assertIsInstance(gaussian_params, Parameters)
        self.assertIsNotNone(kg_lj_params)
        self.assertIsInstance(kg_lj_params, Parameters)
        self.assertIsNotNone(kuhn_params)
        self.assertIsInstance(kuhn_params, Parameters)

    def test_parameter_types_are_different(self):
        """Test that parameters from different types are different."""
        polymer_names = get_supported_polymer_names()
        polymer_name = polymer_names[0]

        gaussian_params = get_parameters_for_polymer(
            polymer_name, ParameterType.GAUSSIAN
        )
        kg_lj_params = get_parameters_for_polymer(polymer_name, ParameterType.KG_LJ)
        kuhn_params = get_parameters_for_polymer(polymer_name, ParameterType.KUHN)

        # Compare some key values that should differ
        # For example, distance_units should be different
        gaussian_distance = gaussian_params.get("distance_units")
        kg_lj_distance = kg_lj_params.get("distance_units")
        kuhn_distance = kuhn_params.get("distance_units")

        # At least one pair should be different
        distances = [gaussian_distance, kg_lj_distance, kuhn_distance]
        self.assertTrue(
            any(d1 != d2 for d1 in distances for d2 in distances if d1 is not d2)
        )

        # Check names are different
        names = [
            gaussian_params.get_name(),
            kg_lj_params.get_name(),
            kuhn_params.get_name(),
        ]
        self.assertEqual(len(set(names)), 3)  # All unique

    def test_dynamic_methods_work(self):
        """Test that dynamic methods like get_characteristic_ratio() work."""
        polymer_names = get_supported_polymer_names()
        polymer_name = polymer_names[0]

        params = get_parameters_for_polymer(polymer_name, ParameterType.GAUSSIAN)

        # Test dynamic method for characteristic ratio
        c_inf_direct = params.get("characteristic_ratio")
        c_inf_dynamic = params.get_characteristic_ratio()
        c_inf_alias = params.get("C_inf")
        self.assertEqual(c_inf_direct, c_inf_dynamic)
        self.assertEqual(c_inf_direct, c_inf_alias)

        # Test another dynamic method
        mw_direct = params.get("Mw")
        mw_dynamic = params.get_Mw()
        self.assertEqual(mw_direct, mw_dynamic)

        # Test temperature
        t_direct = params.get("T")
        t_dynamic = params.get_T()
        self.assertEqual(t_direct, t_dynamic)


if __name__ == "__main__":
    unittest.main()
