#!/usr/bin/env python

import unittest

import numpy as np
import pandas as pd

from pylimer_tools_cpp import Box, Universe


class TestEntityCalculations(unittest.TestCase):
    def test_universe(self):
        universe = Universe(10, 10, 10)
        self.assertIsInstance(universe, Universe)
        self.assertEqual(universe.get_volume(), 10*10*10)
        universe.setBox(Box(1, 1, 1))
        self.assertEqual(universe.get_volume(), 1)
        universe.setBoxLengths(100, 1, 1)
        self.assertEqual(universe.get_volume(), 100*1*1)

    def test_compute_mean_bond_len(self):
        base_atom = {
            "id": 1,
            "x": 0,
            "y": 0,
            "z": 0,
            "nx": 0,
            "ny": 0,
            "nz": 0,
            "type": 1
        }

        bonds_df = pd.DataFrame([{"to": 2, "from": 1}, {"to": 3, "from": 2}])

        for dir in ["x", "y", "z"]:
            second_atom = base_atom.copy()
            second_atom[dir] = 1
            second_atom["id"] = 2

            third_atom = base_atom.copy()
            third_atom[dir] = 2
            third_atom["id"] = 3

            coords_df = pd.DataFrame([third_atom, second_atom, base_atom])
            universe = Universe(2, 2, 2)
            universe.addAtoms(coords_df["id"].tolist(), coords_df["type"].tolist(),
                              coords_df["x"].tolist(),
                              coords_df["y"].tolist(),
                              coords_df["z"].tolist(),
                              coords_df["nx"].tolist(), coords_df["ny"].tolist(), coords_df["nz"].tolist())
            universe.addBonds(
                bonds_df["from"].tolist(), bonds_df["to"].tolist())
            self.assertEqual(len(universe.get_molecules(-1)), 1)
            molecule = universe.get_molecules(-1)[0]
            self.assertEqual(molecule.getLength(), 3)
            self.assertEqual(molecule.get_nr_of_bonds(), 2)
            self.assertEqual(np.mean(molecule.compute_bond_lengths()), 1)

    def test_compute_end_to_end_distance(self):
        base_atom = {
            "id": 1,
            "x": 0,
            "y": 0,
            "z": 0,
            "nx": 0,
            "ny": 0,
            "nz": 0,
            "type": 1
        }
        universe = Universe(1, 1, 1)
        coords_df = pd.DataFrame([base_atom])
        bonds_df = pd.DataFrame([], columns=["to", "bondFrom"])
        universe.addAtoms(coords_df["id"].tolist(), coords_df["type"].tolist(),
                          coords_df["x"].tolist(), coords_df["y"].tolist(),
                          coords_df["z"].tolist(),
                          coords_df["nx"].tolist(), coords_df["ny"].tolist(), coords_df["nz"].tolist())
        universe.addBonds(
            bonds_df["bondFrom"].tolist(), bonds_df["to"].tolist())
        molecules = universe.get_molecules(-1)
        self.assertEqual(len(molecules), 1)
        self.assertEqual(molecules[0].computeEndToEndDistance(), 0)
        atoms = []
        for i in range(3):
            new_atom = base_atom.copy()
            new_atom["id"] = i
            atoms.append(new_atom)
        universe = Universe(1, 1, 1)
        coords_df = pd.DataFrame(atoms)
        bonds_df = pd.DataFrame(
            [{"to": 1, "bondFrom": 0}, {"to": 2, "bondFrom": 1}, {"to": 0, "bondFrom": 2}])
        universe.addAtoms(coords_df["id"].tolist(), coords_df["type"].tolist(),
                          coords_df["x"].tolist(), coords_df["y"].tolist(),
                          coords_df["z"].tolist(),
                          coords_df["nx"].tolist(), coords_df["ny"].tolist(), coords_df["nz"].tolist())
        universe.addBonds(
            bonds_df["bondFrom"].tolist(), bonds_df["to"].tolist())
        molecules = universe.get_molecules(-1)
        self.assertEqual(len(molecules), 1)
        self.assertEqual(-1, molecules[0].computeEndToEndDistance())

    def test_compute_distance_through_periodic_image(self):
        base_atom = {
            "id": 1,
            "x": 0,
            "y": 0,
            "z": 0,
            "nx": 0,
            "ny": 0,
            "nz": 0,
            "type": 1
        }
        for dir in ["x", "y", "z"]:
            second_atom = base_atom.copy()
            second_atom[dir] = 1
            second_atom["id"] = 2
            second_atom["n" + dir] = 1
            universe = Universe(1, 1, 1)
            coords_df = pd.DataFrame([base_atom, second_atom])
            bonds_df = pd.DataFrame([{
                "to": 1, "bondFrom": 2
            }])
            universe.addAtoms(coords_df["id"].tolist(), coords_df["type"].tolist(),
                              coords_df["x"].tolist(
            ), coords_df["y"].tolist(),
                coords_df["z"].tolist(),
                coords_df["nx"].tolist(), coords_df["ny"].tolist(), coords_df["nz"].tolist())
            universe.addBonds(
                bonds_df["bondFrom"].tolist(), bonds_df["to"].tolist())
            self.assertEqual(len(universe.get_molecules(-1)), 1)
            self.assertEqual(0.0, np.mean(universe.get_molecules(-1)[
                0].compute_bond_lengths()))


if __name__ == '__main__':
    unittest.main()
