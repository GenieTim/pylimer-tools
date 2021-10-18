import unittest

import pandas as pd
from pylimer_tools.calc.doMEHPAnalysis import *
from pylimer_tools.entities.universum import Universum
from test.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestMEHPAnalysisFunctions(UniverseUsingTestCase):

    def testWeightFractionCalculations(self):
        self.assertEqual(
            (0.0, 0.0), calculateWeightFractionOfDanglingChains(self.emptyUniverse, 2, 1))
        # empty weight -> empty weight fraction
        self.assertEqual(
            (0.0, 0.25), calculateWeightFractionOfDanglingChains(self.testUniverse, 2, 0))
        self.assertEqual(
            1.0, calculateWeightFractionOfBackbone(self.testUniverse, 2, 0))
        # non-empty weights
        self.assertEqual(
            (0.2, 0.25), calculateWeightFractionOfDanglingChains(self.testUniverse, crosslinkerType=2, weights={1: 1, 2: 0}))

    def testCrosslinkerFunctionalityCalculation(self):
        self.assertCountEqual(
            [], calculateEffectiveCrosslinkerFunctionalities(self.emptyUniverse, 2))
        self.assertSequenceEqual(
            [0, 2, 3], calculateEffectiveCrosslinkerFunctionalities(self.testUniverse, 2))
        self.assertEqual(
            5.0/3.0, calculateEffectiveCrosslinkerFunctionality(self.testUniverse, 2))
        self.assertEqual(
            5.0/3.0/3.0, computeCrosslinkerConversion(self.testUniverse, 2, 3))

    def testMeanEndToEndComputation(self):
        self.assertCountEqual([], computeMeanEndToEndVectors([], 2))
        self.assertDictEqual({
            'atom6+atom7+atom1atom2atom3atom6atom7': 1.0,
            'atom6+atom7+atom5atom6atom7': 1.0
        }, computeMeanEndToEndDistances([self.testUniverse], 2))

    def testEffectiveNrDensityOfJunctionCalculation(self):
        universe = Universum([1, 1, 1])
        self.assertIsNone(calculateEffectiveNrDensityOfJunctions([]))
        universe.addAtomBondData(self.testAtoms, self.testBonds)
        # Border cases
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfJunctions([universe], 0, 0))
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfJunctions([universe], 1000, junctionType=2))
        # Other border
        # 3 junctions, volume of 1
        self.assertEqual(
            3.0/universe.getVolume(), calculateEffectiveNrDensityOfJunctions([universe], 0, junctionType=2, minNumEffectiveStrands=0))
        # actual calc: 6 & 7 are active, 4 not
        self.assertEqual(
            2.0/universe.getVolume(), calculateEffectiveNrDensityOfJunctions([universe], 0, junctionType=2, minNumEffectiveStrands=2))

    def testEffectiveNrDensityOfNetworkCalculation(self):
        universe = Universum([1, 1, 1])
        self.assertIsNone(calculateEffectiveNrDensityOfNetwork([]))
        universe.addAtomBondData(self.testAtoms, self.testBonds)
        self.assertEqual(3, len(universe.getMolecules(2)))
        # Border cases
        self.assertEqual(0.0, calculateEffectiveNrDensityOfNetwork(
            [universe], None, 10, junctionType=2))
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfNetwork([universe], 100, 100, junctionType=2))
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfNetwork([universe], 1000, 1, junctionType=2))
        # actual calc: we got 2 active strands in a Volume of 1
        self.assertEqual(
            2.0, calculateEffectiveNrDensityOfNetwork([universe], 0, 2, junctionType=2))

    def testCycleRankCalculation(self):
        self.assertEqual(1, calculateCycleRank(None, 1, 0))
        self.assertEqual(0, calculateCycleRank(None, 1, 1))
        self.assertEqual(-1, calculateCycleRank(None, 0, 1))
        universe = Universum([10, 10, 10])
        universe.addAtomBondData(self.testAtoms, self.testBonds)
        # test basic exception thrown when specifiying the wrong arguments
        self.assertRaises(ValueError, lambda: calculateCycleRank([universe]))
        self.assertRaises(
            ValueError, lambda: calculateCycleRank([universe], nu=1))
        # same nr of active strands as junctions
        self.assertEqual(
            0.0, calculateCycleRank([universe], None, None, 1, 1, 2))

    def testTopologicalFactorComputation(self):
        universe = Universum([10, 10, 10])
        universe.addAtomBondData(self.testAtoms, self.testBonds)
        self.assertEqual(
            1 + 1.0/3.0, calculateTopologicalFactor([universe], 2, b=1))
        self.assertEqual(
            0.5485762961986437, calculateTopologicalFactor([universe], 2))
