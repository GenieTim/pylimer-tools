
import unittest
import warnings

from pylimer_tools.io.unitStyles import UnitStyle, UnitStyleFactory


class UnitStyleTest(unittest.TestCase):
    def compareUnitOnly(self, unitStyle1: UnitStyle, unitStyle2: UnitStyle, unit: str):
        baseUnit = 1*unitStyle1.getBaseUnitOf(unit)
        unitsToCompare = 1*unitStyle2.getBaseUnitOf(unit)
        self.assertEqual(baseUnit.to_root_units().units,
                         unitsToCompare.to_root_units().units)

    def test_allStylesAreSensible(self):
        unitStyleFactory = UnitStyleFactory()
        baseStyle = unitStyleFactory.getUnitStyle("si")
        otherStyles = [unitStyleFactory.getUnitStyle(
            "nano"), unitStyleFactory.getUnitStyle("real")]
        for style in otherStyles:
            self.assertIsInstance(style, UnitStyle)
            unitsToCheck = [
                'mass', 'distance', 'time', 'energy', 'velocity', 'force', 'torque', 'temperature', 'pressure', 'viscosity', 'charge', 'dipole', 'electric field', 'density', 'dt', 'skin'
            ]
            for unit in unitsToCheck:
                self.compareUnitOnly(baseStyle, style, unit)

        # LJ is a special kind due to its own type of time
        ljStyle = unitStyleFactory.getUnitStyle("lj", polymer="pdms", warning=False)
        lessUnitsToCheck = [
            'mass', 'distance', 'energy', 'force', 'torque', 'temperature', 'pressure', 'density', 'charge', 'dipole', 'electric field'
        ]
        for unit in lessUnitsToCheck:
            self.compareUnitOnly(baseStyle, ljStyle, unit)
    
    def test_errorsAreThrown(self):
        unitStyleFactory = UnitStyleFactory()
        self.assertRaises(ValueError, lambda: unitStyleFactory.getUnitStyle("lj"))

    def test_getAttrEquivalence(self):
        unitStyleFactory = UnitStyleFactory()
        with warnings.catch_warnings(record=True) as w:
          # Cause all warnings to always be triggered.
          warnings.simplefilter("always")
          # Trigger a warning.
          unitStyle = unitStyleFactory.getUnitStyle("lj", polymer="pdms", warning=True)
          # Assert that the warning has been triggered
          self.assertTrue(len(w) == 1)
        self.assertEqual(1*unitStyle.mass, 1*unitStyle.getBaseUnitOf("mass"))



