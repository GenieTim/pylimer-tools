
import unittest

import pandas as pd
from pylimer_tools.calc.calculateBondLen import (calculateBondLen,
                                                 calculateMeanBondLen)
from pylimer_tools.entities.universum import Universum


class TestEntityCalculations(unittest.TestCase):
    def test_universe(self):
        universe = Universum(boxSizes=[10, 10, 10])
        self.assertIsInstance(universe, Universum)
        self.assertEqual(universe.getVolume(), 10*10*10)
        universe.setBoxSizes([1, 1, 1])
        self.assertEqual(universe.getVolume(), 1)

    def test_calculateMeanBondLen(self):
        baseAtom = {
            "id": 1,
            "x": 0,
            "y": 0,
            "z": 0,
            "nx": 0,
            "ny": 0,
            "nz": 0,
            "type": 1
        }

        bondsDf = pd.DataFrame([{"to": 2, "from": 1}, {"to": 3, "from": 2}])

        for dir in ["x", "y", "z"]:
            secondAtom = baseAtom.copy()
            secondAtom[dir] = 1
            secondAtom["id"] = 2

            thirdAtom = baseAtom.copy()
            thirdAtom[dir] = 2
            thirdAtom["id"] = 3

            coordsDf = pd.DataFrame([thirdAtom, secondAtom, baseAtom])
            universe = Universum(boxSizes=[1, 1, 1])
            universe.addAtomBondData(coordsDf, bondsDf)
            self.assertEqual(len(universe.getMolecules()), 1)
            molecule = universe.getMolecules()[0]
            self.assertEqual(molecule.computeBondLengths().mean(), 1)
