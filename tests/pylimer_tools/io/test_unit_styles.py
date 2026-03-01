import unittest
import warnings

from pint import Quantity, UnitRegistry

from pylimer_tools.io.bead_spring_parameter_provider import (
    Parameters,
    ParameterType,
    get_parameters_for_polymer,
    get_supported_polymer_names,
)
from pylimer_tools.io.unit_styles import UnitStyle, UnitStyleFactory


class UnitStyleTest(unittest.TestCase):
    def compare_unit_only(
        self, unit_style1: UnitStyle, unit_style2: UnitStyle, unit: str
    ):
        base_unit = 1 * unit_style1.get_base_unit_of(unit)
        units_to_compare = 1 * unit_style2.get_base_unit_of(unit)
        self.assertEqual(
            base_unit.to_root_units().units, units_to_compare.to_root_units().units
        )

    def test_all_styles_are_sensible(self):
        unit_style_factory = UnitStyleFactory()
        self.assertIsInstance(unit_style_factory.get_unit_registry(), UnitRegistry)
        base_style = unit_style_factory.get_unit_style("si")
        other_styles = [
            unit_style_factory.get_unit_style("nano"),
            unit_style_factory.get_unit_style("real"),
        ]
        for style in other_styles:
            self.assertIsInstance(style, UnitStyle)
            units_to_check = [
                "mass",
                "distance",
                "time",
                "energy",
                "velocity",
                "force",
                "torque",
                "temperature",
                "pressure",
                "viscosity",
                "charge",
                "dipole",
                "electric field",
                "density",
                "dt",
                "skin",
            ]
            for unit in units_to_check:
                self.compare_unit_only(base_style, style, unit)

        # LJ is a special kind due to its own type of time
        lj_style = unit_style_factory.get_unit_style(
            "lj", polymer="pdms", warning=False
        )
        less_units_to_check = [
            "mass",
            "distance",
            "energy",
            "force",
            "torque",
            "temperature",
            "pressure",
            "density",
            "charge",
            "dipole",
            "electric field",
        ]
        for unit in less_units_to_check:
            self.compare_unit_only(base_style, lj_style, unit)

    def test_errors_are_thrown(self):
        unit_style_factory = UnitStyleFactory()
        self.assertRaises(ValueError, lambda: unit_style_factory.get_unit_style("lj"))

    def test_get_attr_equivalence(self):
        unit_style_factory = UnitStyleFactory()
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            unit_style = unit_style_factory.get_unit_style(
                "lj", polymer="pdms", warning=True
            )
            # Assert that the warning has been triggered
            self.assertTrue(len(w) == 1)
        self.assertEqual(1 * unit_style.mass, 1 * unit_style.get_base_unit_of("mass"))

    def test_get_parameters_for_polymer(self):
        self.assertTrue(len(get_supported_polymer_names()) > 10)

        def check_parameter_types(params: Parameters):
            self.assertIsInstance(params.get_unit_registry(), UnitRegistry)
            self.assertIsInstance(params.get("distance_units"), Quantity)
            self.assertIsInstance(params.get_bead_density(), float)
            self.assertIsInstance(params.get_entanglement_density(), float)
            self.assertIsInstance(params.get_sampling_cutoff(), float)
            self.assertIsInstance(params.get_fb_stress_conversion(), float)
            self.assertIsInstance(params.get_kappa(), Quantity)
            self.assertIsInstance(params.get("T"), Quantity)
            self.assertIsInstance(params.get_gamma_conversion_factor(), Quantity)

        # test the Kremer-Grest/Lennard-Jones parameters
        for polymer_name in get_supported_polymer_names():
            params = get_parameters_for_polymer(
                polymer_name=polymer_name, parameter_type=ParameterType.KG_LJ
            )
            self.assertIsInstance(params, Parameters)
            self.assertEqual(params.get_name(), "kg-lj-" + polymer_name)
            check_parameter_types(params)

        # test the kuhn parameters
        for polymer_name in get_supported_polymer_names():
            params = get_parameters_for_polymer(
                polymer_name=polymer_name, parameter_type=ParameterType.KUHN
            )
            self.assertIsInstance(params, Parameters)
            self.assertEqual(params.get_name(), "kuhn-si-" + polymer_name)
            check_parameter_types(params)

        # test the gaussian parameters
        for polymer_name in get_supported_polymer_names():
            params = get_parameters_for_polymer(
                polymer_name, parameter_type=ParameterType.GAUSSIAN
            )
            self.assertIsInstance(params, Parameters)
            self.assertEqual(params.get_name(), "si-" + polymer_name)
            check_parameter_types(params)

    def test_get_polymer_names(self):
        polymer_names = get_supported_polymer_names()
        self.assertGreater(len(polymer_names), 0)
        for name in polymer_names:
            self.assertIsInstance(name, str)
            self.assertGreater(len(name), 0)
            params = get_parameters_for_polymer(
                name, parameter_type=ParameterType.GAUSSIAN
            )
            self.assertIsInstance(params, Parameters)
            self.assertEqual(params.get_name(), "si-" + name)


if __name__ == "__main__":
    unittest.main()
