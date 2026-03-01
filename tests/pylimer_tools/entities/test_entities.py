#!/usr/bin/env python

import copy
import os
import sys
import unittest

import numpy as np
import pandas as pd
import pandas.testing as pd_testing

from pylimer_tools_cpp import Atom, Molecule, MoleculeType, Universe

if __name__ == "__main__":
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../.."))
from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestEntities(UniverseUsingTestCase):
    def assert_series_equal(self, a, b, msg):
        try:
            pd_testing.assert_series_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.Series, self.assert_series_equal)
        super().setUp()

    def test_universe(self):
        self.assertIsInstance(self.emptyUniverse, Universe)
        # self.assertIsInstance(universe.getUnderlyingGraph(), igraph.Graph)
        # check that the except paths work too: non-existant atom ids & type
        self.assertEqual([], self.emptyUniverse.get_atoms_by_type(1))
        self.assertRaises(IndexError, lambda: self.emptyUniverse.get_atom(1))

        self.assertCountEqual([], self.emptyUniverse.get_molecules(0))
        self.assertCountEqual([], self.emptyUniverse.get_chains_with_crosslinker(0))
        self.assertEqual(0, self.emptyUniverse.get_nr_of_atoms())

        atom = self.testUniverseSmall.get_atom(1)
        self.assertEqual(
            atom,
            Atom(
                self.testAtomsSmall.iloc[0]["id"],
                self.testAtomsSmall.iloc[0]["type"],
                self.testAtomsSmall.iloc[0]["x"],
                self.testAtomsSmall.iloc[0]["y"],
                self.testAtomsSmall.iloc[0]["z"],
                self.testAtomsSmall.iloc[0]["nx"],
                self.testAtomsSmall.iloc[0]["ny"],
                self.testAtomsSmall.iloc[0]["nz"],
            ),
        )
        self.assertIsInstance(atom, Atom)

    def test_molecule_entity(self):
        universe = self.testUniverseSmall
        self.assertEqual(4, len(universe.get_atoms_by_type(1)))
        self.assertEqual(2, len(universe.get_atoms_by_type(2)))
        molecules = universe.get_molecules(0)
        self.assertEqual(len(molecules), 2)
        self.assertEqual(molecules[0].get_nr_of_atoms(), 3)
        self.assertEqual(
            np.sum([m.get_nr_of_atoms() for m in molecules]), len(self.testAtomsSmall)
        )
        molecules = universe.get_molecules(2)
        self.assertEqual(len(molecules), 2)
        self.assertEqual(len(universe.get_chains_with_crosslinker(0)), 2)
        for molecule in molecules:
            self.assertIsInstance(molecule, Molecule)
        for molecule in universe.get_molecules(0):
            self.assertIsInstance(molecule, Molecule)
        for molecule in universe.get_chains_with_crosslinker(0):
            self.assertIsInstance(molecule, Molecule)
            self.assertEqual(molecule.get_strand_type(), MoleculeType.FREE_CHAIN)

        chains_with_crosslinker = universe.get_chains_with_crosslinker(2)
        self.assertEqual(
            chains_with_crosslinker[0].get_strand_type(), MoleculeType.FREE_CHAIN
        )
        self.assertEqual(
            chains_with_crosslinker[1].get_strand_type(), MoleculeType.DANGLING_CHAIN
        )
        universe_clone = copy.copy(universe)
        self.assertEqual(universe.get_nr_of_atoms(), universe_clone.get_nr_of_atoms())

    def test_molecule_entity_iterations(self):
        molecules = self.testUniverseSmall.get_molecules(0)
        # test iteration & return type
        for molecule in molecules:
            self.assertIsInstance(molecule, Molecule)
        # test calculations
        self.assertEqual(molecules[0].compute_end_to_end_distance(), 2)
        self.assertEqual(np.mean(molecules[0].compute_bond_lengths()), 1.0)

    def test_atom_entity(self):
        atom1 = Atom(1, 1, 0.0, 0.0, 0.0, 0, 0, 0)
        self.assertIsInstance(atom1, Atom)
        self.assertEqual(atom1.get_type(), 1)
        self.assertEqual(atom1.get_id(), 1)
        self.assertEqual(atom1.get_x(), 0)
        self.assertEqual(atom1.get_y(), 0)
        self.assertEqual(atom1.get_z(), 0)
        self.assertEqual(atom1.get_nx(), 0)
        self.assertEqual(atom1.get_ny(), 0)
        self.assertEqual(atom1.get_nz(), 0)


if __name__ == "__main__":
    unittest.main()
