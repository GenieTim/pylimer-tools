import unittest
import sys
import os

import pandas as pd
from pylimer_tools.calc.doMEHPAnalysis import *
from pylimer_tools_cpp import Universe

if __name__ == '__main__':
    sys.path.append(os.path.join(os.path.dirname(__file__), "../../../"))

from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestMEHPAnalysisFunctions(UniverseUsingTestCase):

    def test_weightFractionCalculations(self):
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

    def test_crosslinkerFunctionalityCalculation(self):
        self.assertCountEqual(
            [], calculateEffectiveCrosslinkerFunctionalities(self.emptyUniverse, 2))
        self.assertSequenceEqual(
            [0, 2, 3], calculateEffectiveCrosslinkerFunctionalities(self.testUniverse, 2))
        self.assertEqual(
            5.0/3.0, calculateEffectiveCrosslinkerFunctionality(self.testUniverse, 2))
        self.assertEqual(
            5.0/3.0/3.0, computeCrosslinkerConversion(self.testUniverse, 2, 3))

    def test_meanEndToEndComputation(self):
        self.assertCountEqual([], computeMeanEndToEndVectors([], 2))
        self.assertDictEqual({
            'atom6+atom7+atom1atom2atom3atom6atom7': 1.0,
            'atom6+atom7+atom5atom6atom7': 1.0
        }, computeMeanEndToEndDistances([self.testUniverse], 2))

    def test_meanUniverseVolume(self):
        self.assertRaises(NotImplementedError, lambda: calculateMeanUniverseVolume(
            [self.testUniverse, self.testUniverseSmall]))

    def test_effectiveNrDensityOfJunctionCalculation(self):
        self.assertIsNone(calculateEffectiveNrDensityOfJunctions([]))
        # Border cases
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfJunctions([self.testUniverse], 0, 0))
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfJunctions([self.testUniverse], 1000, junctionType=2))
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfJunctions([self.emptyUniverse], 1000, junctionType=2))
        # Other border
        # 3 junctions, volume of 1
        self.assertEqual(
            3.0/self.testUniverse.getVolume(), calculateEffectiveNrDensityOfJunctions([self.testUniverse], 0, junctionType=2, minNumEffectiveStrands=0))
        # actual calc: 6 & 7 are active, 4 not
        self.assertEqual(
            2.0/self.testUniverse.getVolume(), calculateEffectiveNrDensityOfJunctions([self.testUniverse], absTol=None, relTol=0, junctionType=2, minNumEffectiveStrands=2))
        self.assertEqual(
            2.0/self.testUniverse.getVolume(), calculateEffectiveNrDensityOfJunctions([self.testUniverse], 0, junctionType=2, minNumEffectiveStrands=2))

    def test_effectiveNrDensityOfNetworkCalculation(self):
        self.assertIsNone(calculateEffectiveNrDensityOfNetwork([]))
        self.assertEqual(3, len(self.testUniverse.getMolecules(2)))
        # Border cases
        self.assertEqual(0.0, calculateEffectiveNrDensityOfNetwork(
            [self.testUniverse], None, 10, junctionType=2))
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfNetwork([self.testUniverse], 100, 100, junctionType=2))
        self.assertEqual(
            0.0, calculateEffectiveNrDensityOfNetwork([self.testUniverse], 1000, 1, junctionType=2))
        # actual calc: we got 2 active strands in a Volume of 1
        self.assertEqual(
            2.0, calculateEffectiveNrDensityOfNetwork([self.testUniverse], 0, 2, junctionType=2))

    def test_cycleRankCalculation(self):
        self.assertEqual(1, calculateCycleRank(None, 1, 0))
        self.assertEqual(0, calculateCycleRank(None, 1, 1))
        self.assertEqual(-1, calculateCycleRank(None, 0, 1))
        universe = Universe(10, 10, 10)
        universe = self.addAtomBondData(universe, self.testAtoms, self.testBonds)
        # test basic exception thrown when specifiying the wrong arguments
        self.assertRaises(ValueError, lambda: calculateCycleRank([universe]))
        self.assertRaises(
            ValueError, lambda: calculateCycleRank([universe], nu=1))
        # same nr of active strands as junctions
        self.assertEqual(
            0.0, calculateCycleRank([universe], None, None, 1, 1, 2))
        # other system
        self.assertEqual(
            1/(10*10*10), calculateCycleRank([self.saturatedTestUniverse], None, None, 1, 1, 2))

    def test_topologicalFactorComputation(self):
        self.assertEqual(
            1 + 1.0/3.0, calculateTopologicalFactor([self.testUniverse], 2, b=1))
        self.assertEqual(
            0.5485762961986437, calculateTopologicalFactor([self.testUniverse], 2))
        # larger system
        # g = self.saturatedTestUniverse.getUnderlyingGraph()
        # igraph.plot(g, vertex_label=g.vs["name"], vertex_color=["green" if n["type"] == 2 else "red" for n in g.vs], target="large_test.png", vertex_label_dist=1)
        self.assertEqual(0.7249043914053506, calculateTopologicalFactor(
            [self.saturatedTestUniverse], 2))

    def test_shearModulusPrediction(self):
        self.assertEqual(0.0, predictShearModulus(
            [self.emptyUniverse], foreignAtomType=2))
        self.assertEqual(0.003624521957026753, predictShearModulus(
            [self.saturatedTestUniverse], foreignAtomType=2))


if __name__ == '__main__':
    unittest.main()
