
import unittest

import pandas as pd
from pylimer_tools_cpp import Atom, Box


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

            atom1 = Atom(
                0, 0, secondAtom["xsu"], secondAtom["ysu"], secondAtom["zsu"], 0, 0, 0)
            atom2 = Atom(
                0, 0, thirdAtom["xsu"], thirdAtom["ysu"], thirdAtom["zsu"], 0, 0, 0)
            self.assertEqual(1, atom1.distanceTo(
                atom2, Box(10, 10, 10)))


if __name__ == '__main__':
    unittest.main()
