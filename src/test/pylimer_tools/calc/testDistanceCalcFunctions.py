
import unittest

import pandas as pd
from pylimer_tools.calc.calculateBondLen import (calculateBondLen,
                                                 calculateMeanBondLen)
from pylimer_tools.calc.calculateDistanceBetweenAtoms import \
    calculateDistanceBetweenAtoms


class TestDistanceCalcFunctions(unittest.TestCase):

    def test_calculateMeanBondLen(self):
        baseAtom = {
            "id": 1,
            "xsu": 0,
            "ysu": 0,
            "zsu": 0
        }

        for dir in ["xsu", "ysu", "zsu"]:
            secondAtom = baseAtom.copy()
            secondAtom[dir] = 1
            secondAtom["id"] = 2

            thirdAtom = baseAtom.copy()
            thirdAtom[dir] = 2
            thirdAtom["id"] = 3

            coordsDf = pd.DataFrame([thirdAtom, secondAtom, baseAtom])
            self.assertEqual(1, calculateMeanBondLen(coordsDf, [1, 1, 1]))

    def test_calculateBondLen(self):
        baseAtom = {
            "id": 0,
            "x": 0,
            "y": 0,
            "z": 0,
            "type": 0
        }

        for dir in ["x", "y", "z"]:
            secondAtom = baseAtom.copy()
            secondAtom[dir] = 1
            secondAtom["id"] = 1

            thirdAtom = baseAtom.copy()
            thirdAtom[dir] = 2
            thirdAtom["id"] = 2

            coordsDf = pd.DataFrame([baseAtom, secondAtom, thirdAtom])
            bondsDf = pd.DataFrame(
                [{"to": 1, "bondFrom": 0}, {"to": 2, "bondFrom": 1}])
            self.assertEqual(1, calculateBondLen(
                coordsDf, bondsDf, [10, 10, 10]).mean())

    def test_calculateDistanceBetweenAtoms(self):
        self.assertEqual(1, calculateDistanceBetweenAtoms({
            "x": 0,
            "y": 0,
            "z": 0
        }, {
            "x": 0,
            "y": 0,
            "z": 1
        }))


if __name__ == '__main__':
    unittest.main()
