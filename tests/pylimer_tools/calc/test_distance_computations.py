#!/usr/bin/env python
import unittest

from pylimer_tools_cpp import Atom, Box


class TestDistanceCalcFunctions(unittest.TestCase):

    def test_compute_mean_bond_len(self):
        base_atom = {
            "id": 1,
            "xsu": 0,
            "ysu": 0,
            "zsu": 0
        }

        for dir in ["xsu", "ysu", "zsu"]:
            second_atom = base_atom.copy()
            second_atom[dir] = 1
            second_atom["id"] = 2

            third_atom = base_atom.copy()
            third_atom[dir] = 2
            third_atom["id"] = 3

            atom1 = Atom(
                0, 0, second_atom["xsu"], second_atom["ysu"], second_atom["zsu"], 0, 0, 0)
            atom2 = Atom(
                0, 0, third_atom["xsu"], third_atom["ysu"], third_atom["zsu"], 0, 0, 0)
            self.assertEqual(1, atom1.distanceTo(
                atom2, Box(10, 10, 10)))


if __name__ == '__main__':
    unittest.main()
