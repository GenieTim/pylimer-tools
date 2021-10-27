
import unittest
from test.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase

import igraph
import numpy as np
import pandas as pd
import pandas.testing as pd_testing
from pylimer_tools.calc.calculateBondLen import (calculateBondLen,
                                                 calculateMeanBondLen)
from pylimer_tools.entities.atom import Atom
from pylimer_tools.entities.molecule import Molecule
from pylimer_tools.entities.universum import Universum


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
        universe = Universum(boxSizes=[10, 10, 10])
        self.assertIsInstance(universe, Universum)
        universe.addAtomBondData(self.testAtomsSmall, self.testBondsSmall)
        atom = universe.getAtom(1)
        self.assertEqual(atom.getUnderlyingData(), self.testAtomsSmall.iloc[0])
        universe.reset()
        self.assertCountEqual([], universe.getMolecules())
        self.assertCountEqual([], universe.getChainsWithCrosslinker(0))
        self.assertEqual(0, universe.getSize())
        self.assertIsInstance(universe.getUnderlyingGraph(), igraph.Graph)
        # check that the except paths work too: non-existant atom ids & type
        self.assertEqual(None, universe.getAtomsWithType(1))
        self.assertEqual(None, universe.getAtom(1))

    def test_moleculeEntity(self):
        universe = self.testUniverseSmall
        self.assertEqual(4, len(universe.getAtomsWithType(1)))
        self.assertEqual(2, len(universe.getAtomsWithType(2)))
        molecules = universe.getMolecules()
        self.assertEqual(len(molecules), 2)
        self.assertEqual(molecules[0].getLength(), 3)
        self.assertEqual(np.sum([m.getLength()
                                 for m in molecules]), len(self.testAtomsSmall))
        molecules = universe.getMolecules(ignoreAtomType=2)
        self.assertEqual(len(molecules), 2)
        self.assertEqual(len(universe.getChainsWithCrosslinker(0)), 2)
        for molecule in molecules:
            self.assertIsInstance(molecule, Molecule)
        for molecule in universe:
            self.assertIsInstance(molecule, Molecule)
        for molecule in universe.getChainsWithCrosslinker(0):
            self.assertIsInstance(molecule, Molecule)
            self.assertEqual(molecule.getType(),
                             Molecule.MoleculeType.FREE_CHAIN)

        chainsWithCrosslinker = universe.getChainsWithCrosslinker(
            crosslinkerType=2)
        self.assertEqual(chainsWithCrosslinker[0].getType(
        ), Molecule.MoleculeType.FREE_CHAIN)
        self.assertEqual(
            chainsWithCrosslinker[1].getType(), Molecule.MoleculeType.DANGLING_CHAIN)

    def test_moleculeEntityIterations(self):
        molecules = self.testUniverseSmall.getMolecules()
        # test iteration & return type
        for molecule in molecules:
            self.assertIsInstance(molecule, Molecule)
        # test calculations
        self.assertEqual(molecules[0].computeEndToEndDistance(), 2)
        self.assertEqual(molecules[0].computeBondLengths().mean(), 1.0)

        # test iteration for empty and 1 atom molecules
        universe = Universum(boxSizes=[10, 10, 10])
        universe.addAtomBondData(pd.DataFrame([{
            "id": 1,
            "x": 0,
            "y": 0,
            "z": 0,
            "nx": 0,
            "ny": 0,
            "nz": 0,
            "type": 1
        }]), pd.DataFrame([], columns=["to", "bondFrom"]))
        molecules = universe.getMolecules()
        self.assertEqual(1, len(molecules))
        self.assertEqual(molecules[0].getLength(), 1)
        iterations = 0
        for atom in molecules[0]:
            self.assertIsInstance(atom, Atom)
            iterations += 1
        self.assertEqual(1, iterations)

    def test_atomEntity(self):
        atom1 = Atom(data=pd.Series({
            "id": 1,
            "x": 0,
            "y": 0,
            "z": 0,
            "nx": 0,
            "ny": 0,
            "nz": 0,
            "type": 1
        }), boxSizes=[1, 1, 1])
        self.assertIsInstance(atom1, Atom)
        self.assertEqual(atom1.name, "atom1")
