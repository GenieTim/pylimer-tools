
import os
import sys
import unittest

import numpy as np
import pandas as pd
import pandas.testing as pd_testing
from pylimer_tools_cpp import Atom, Molecule, Universe, MoleculeType

if __name__ == '__main__':
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../.."))
from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestEntities(UniverseUsingTestCase):

    def assertSeriesEqual(self, a, b, msg):
        try:
            pd_testing.assert_series_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.Series, self.assertSeriesEqual)
        super().setUp()

    def test_universe(self):
        self.assertIsInstance(self.emptyUniverse, Universe)
        # self.assertIsInstance(universe.getUnderlyingGraph(), igraph.Graph)
        # check that the except paths work too: non-existant atom ids & type
        self.assertEqual([], self.emptyUniverse.getAtomsWithType(1))
        self.assertRaises(ValueError, lambda: self.emptyUniverse.getAtom(1))

        self.assertCountEqual([], self.emptyUniverse.getMolecules(0))
        self.assertCountEqual(
            [], self.emptyUniverse.getChainsWithCrosslinker(0))
        self.assertEqual(0, self.emptyUniverse.getNrOfAtoms())

        atom = self.testUniverseSmall.getAtom(1)
        # self.assertEqual(atom.getUnderlyingData(), self.testAtomsSmall.iloc[0])
        self.assertIsInstance(atom, Atom)

    def test_moleculeEntity(self):
        universe = self.testUniverseSmall
        self.assertEqual(4, len(universe.getAtomsWithType(1)))
        self.assertEqual(2, len(universe.getAtomsWithType(2)))
        molecules = universe.getMolecules(0)
        self.assertEqual(len(molecules), 2)
        self.assertEqual(molecules[0].getLength(), 3)
        self.assertEqual(np.sum([m.getLength()
                                 for m in molecules]), len(self.testAtomsSmall))
        molecules = universe.getMolecules(2)
        self.assertEqual(len(molecules), 2)
        self.assertEqual(len(universe.getChainsWithCrosslinker(0)), 2)
        for molecule in molecules:
            self.assertIsInstance(molecule, Molecule)
        for molecule in universe.getMolecules(0):
            self.assertIsInstance(molecule, Molecule)
        for molecule in universe.getChainsWithCrosslinker(0):
            self.assertIsInstance(molecule, Molecule)
            self.assertEqual(molecule.getType(),
                             MoleculeType.FREE_CHAIN)

        chainsWithCrosslinker = universe.getChainsWithCrosslinker(2)
        self.assertEqual(chainsWithCrosslinker[0].getType(
        ), MoleculeType.FREE_CHAIN)
        self.assertEqual(
            chainsWithCrosslinker[1].getType(), MoleculeType.DANGLING_CHAIN)

    def test_moleculeEntityIterations(self):
        molecules = self.testUniverseSmall.getMolecules(0)
        # test iteration & return type
        for molecule in molecules:
            self.assertIsInstance(molecule, Molecule)
        # test calculations
        self.assertEqual(molecules[0].computeEndToEndDistance(), 2)
        self.assertEqual(np.mean(molecules[0].computeBondLengths()), 1.0)

    def test_atomEntity(self):
        atom1 = Atom(1, 1, 0.0, 0.0, 0.0, 0, 0, 0)
        self.assertIsInstance(atom1, Atom)
        self.assertEqual(atom1.getType(), 1)
        self.assertEqual(atom1.getId(), 1)
        self.assertEqual(atom1.getX(), 0)
        self.assertEqual(atom1.getY(), 0)
        self.assertEqual(atom1.getZ(), 0)
        self.assertEqual(atom1.getNX(), 0)
        self.assertEqual(atom1.getNY(), 0)
        self.assertEqual(atom1.getNZ(), 0)


if __name__ == '__main__':
    unittest.main()
