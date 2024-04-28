
import unittest
import warnings

from pint import UnitRegistry

from pylimer_tools.io.unit_styles import UnitStyle, UnitStyleFactory


class UnitStyleTest(unittest.TestCase):
    def compare_unit_only(self, unit_style1: UnitStyle, unit_style2: UnitStyle, unit: str):
        base_unit = 1 * unit_style1.get_base_unit_of(unit)
        units_to_compare = 1 * unit_style2.get_base_unit_of(unit)
        self.assertEqual(base_unit.to_root_units().units,
                         units_to_compare.to_root_units().units)

    def test_all_styles_are_sensible(self):
        unit_style_factory = UnitStyleFactory()
        self.assertIsInstance(
            unit_style_factory.get_unit_registry(), UnitRegistry)
        base_style = unit_style_factory.get_unit_style("si")
        other_styles = [unit_style_factory.get_unit_style(
            "nano"), unit_style_factory.get_unit_style("real")]
        for style in other_styles:
            self.assertIsInstance(style, UnitStyle)
            units_to_check = [
                'mass', 'distance', 'time', 'energy', 'velocity', 'force', 'torque', 'temperature',
                'pressure', 'viscosity', 'charge', 'dipole', 'electric field', 'density', 'dt', 'skin'
            ]
            for unit in units_to_check:
                self.compare_unit_only(base_style, style, unit)

        # LJ is a special kind due to its own type of time
        lj_style = unit_style_factory.get_unit_style(
            "lj", polymer="pdms", warning=False)
        less_units_to_check = [
            'mass', 'distance', 'energy', 'force', 'torque', 'temperature', 'pressure', 'density',
            'charge', 'dipole', 'electric field'
        ]
        for unit in less_units_to_check:
            self.compare_unit_only(base_style, lj_style, unit)

    def test_errors_are_thrown(self):
        unit_style_factory = UnitStyleFactory()
        self.assertRaises(
            ValueError, lambda: unit_style_factory.get_unit_style("lj"))

    def test_get_attr_equivalence(self):
        unit_style_factory = UnitStyleFactory()
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            unit_style = unit_style_factory.get_unit_style(
                "lj", polymer="pdms", warning=True)
            # Assert that the warning has been triggered
            self.assertTrue(len(w) == 1)
        self.assertEqual(1 * unit_style.mass, 1 *
                         unit_style.get_base_unit_of("mass"))
