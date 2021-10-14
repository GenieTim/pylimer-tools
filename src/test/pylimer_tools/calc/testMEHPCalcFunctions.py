import unittest

import pandas as pd
from pylimer_tools.calc.doMEHPAnalysis import *
from pylimer_tools.entities.universum import Universum


class TestDistanceCalcFunctions(unittest.TestCase):

    # The system looks like this (in terms of bonds, not 3D placement):
    # 1-2-3-*6
    # |      |
    # *7-5---|
    # 8
    #
    # *4
    testAtoms = pd.DataFrame([
        {"id": 1, "nx": 1, "ny": 1, "nz": 1,
         "type": 1, "x": 1, "y": 1, "z": 1},
        {"id": 2, "nx": 1, "ny": 1, "nz": 1,
         "type": 1, "x": 2, "y": 1, "z": 1},
        {"id": 3, "nx": 1, "ny": 1, "nz": 1,
         "type": 1, "x": 3, "y": 1, "z": 1},
        {"id": 4, "nx": 1, "ny": 1, "nz": 1,
         "type": 2, "x": 2, "y": 2, "z": 1},
        {"id": 5, "nx": 1, "ny": 1, "nz": 1,
         "type": 1, "x": 1, "y": 3, "z": 1},
        {"id": 6, "nx": 1, "ny": 1, "nz": 1,
         "type": 2, "x": 1, "y": 1, "z": 2},
        {"id": 7, "nx": 1, "ny": 1, "nz": 1,
         "type": 2, "x": 1, "y": 1, "z": 3},
        {"id": 8, "nx": 1, "ny": 1, "nz": 1,
         "type": 1, "x": 2, "y": 2, "z": 2},
    ])
    testBonds = pd.DataFrame([
        {"to": 1, "bondFrom": 2},
        {"to": 3, "bondFrom": 2},
        {"to": 5, "bondFrom": 6},
        {"to": 1, "bondFrom": 7},
        {"to": 5, "bondFrom": 7},
        {"to": 3, "bondFrom": 6},
        {"to": 7, "bondFrom": 8}
    ])

    def testCycleRankCalculation(self):
        self.assertEqual(1, calculateCycleRank(None, 1, 0))
        self.assertEqual(0, calculateCycleRank(None, 1, 1))
        self.assertEqual(-1, calculateCycleRank(None, 0, 1))
        universe = Universum([10, 10, 10])
        universe.addAtomBondData(self.testAtoms, self.testBonds)
        # test basic exception thrown when specifiying the wrong arguments
        self.assertRaises(ValueError, lambda: calculateCycleRank(universe))
        self.assertRaises(
            ValueError, lambda: calculateCycleRank(universe, nu=1))
        self.assertEqual(
            3.0/8.0, calculateCycleRank(universe, None, None, 1, 1, 2))

    def testEffectiveNrDensityOfJunctionCalculation(self):
        universe = Universum([10, 10, 10])
        self.assertIsNone(calculateEffectiveNrDensityOfJunctions(universe))
        universe.addAtomBondData(self.testAtoms, self.testBonds)
        # Border cases
        self.assertRaises(
            ValueError, lambda: calculateEffectiveNrDensityOfJunctions(universe, None, 2))
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfJunctions(universe, 0, 0))
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfJunctions(universe, 1000, 2))
        # Other border
        self.assertEqual(
            1.0, calculateEffectiveNrDensityOfJunctions(universe, 0, 2, 0))
        # actual calc: 6 & 7 are active, 4 not
        self.assertEqual(
            2.0/3.0, calculateEffectiveNrDensityOfJunctions(universe, 0, 2, 2))

    def testEffectiveNrDensityOfNetworkCalculation(self):
        universe = Universum([10, 10, 10])
        self.assertIsNone(calculateEffectiveNrDensityOfNetwork(universe))
        universe.addAtomBondData(self.testAtoms, self.testBonds)
        self.assertEqual(3, len(universe.getMolecules(2)))
        # Border cases
        self.assertEqual(0.0, calculateEffectiveNrDensityOfNetwork(
            universe, None, 10, junctionType=2))
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfNetwork(universe, 100, 100, junctionType=2))
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfNetwork(universe, 1000, 1, junctionType=2))
        # actual calc
        self.assertEqual(
            3.0/8.0, calculateEffectiveNrDensityOfNetwork(universe, 0, 2, junctionType=2))

    def testWeightFractionCalculations(self):
        universe = Universum([10, 10, 10])
        self.assertEqual(
            (0.0, 0.0), calculateWeightFractionOfDanglingChains(universe, 2, 1))
        universe.addAtomBondData(self.testAtoms, self.testBonds)
        # empty weight -> empty weight fraction
        self.assertEqual(
            (0.0, 0.25), calculateWeightFractionOfDanglingChains(universe, 2, 0))
        self.assertEqual(
            1.0, calculateWeightFractionOfBackbone(universe, 2, 0))
        # non-empty weights
        self.assertEqual(
            (0.2, 0.25), calculateWeightFractionOfDanglingChains(universe, crosslinkerType=2, weights={1: 1, 2: 0}))

    def testGelationPointPrediction(self):
        universe = Universum([10, 10, 10])
        universe.addAtomBondData(self.testAtoms, self.testBonds)
        self.assertEqual(1, predictGelationPoint(1, 2))
        self.assertEqual(1, predictGelationPoint(1, 2, 2))

    def testCrosslinkerFunctionalityCalculation(self):
        universe = Universum([10, 10, 10])
        self.assertCountEqual(
            [], calculateEffectiveCrosslinkerFunctionalities(universe, 2))
        universe.addAtomBondData(self.testAtoms, self.testBonds)
        self.assertSequenceEqual(
            [0, 2, 3], calculateEffectiveCrosslinkerFunctionalities(universe, 2))
        self.assertEqual(
            5.0/3.0, calculateEffectiveCrosslinkerFunctionality(universe, 2))
        self.assertEqual(
            5.0/3.0/3.0, computeCrosslinkerConversion(universe, 2, 3))

    def testTopologicalFactorComputation(self):
        universe = Universum([10, 10, 10])
        universe.addAtomBondData(self.testAtoms, self.testBonds)
        self.assertEqual(1 + 1.0/3.0, calculateTopologicalFactor(universe, 2))
