import os
import sys
import unittest

from pylimer_tools.calc.doMMTAnalysis import *

if __name__ == '__main__':
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../"))
from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestMMTAnalysisFunctions(UniverseUsingTestCase):

    def testStoichiometricInbalance(self):
        self.assertEqual(
            0, computeStoichiometricInbalance(self.emptyUniverse, 2))
        self.assertEqual(
            (3*3)/(5*2), computeStoichiometricInbalance(self.testUniverse, 2, strandLength=1))
        self.assertEqual(
            (3*3)/(2), computeStoichiometricInbalance(self.testUniverse, 2, strandLength=5))

    def testExtentOfReaction(self):
        self.assertEqual(1.0, computeExtentOfReaction(self.emptyUniverse, 2))
        self.assertEqual(5.0/6.0, computeExtentOfReaction(self.testUniverse, 2))

    def testGelationPointPrediction(self):
        self.assertEqual(1, predictGelationPoint(1, 2))
        self.assertEqual(1, predictGelationPoint(1, 2, 2))

    def testShearModulusPrediction(self):
        self.assertIsNone(predictShearModulus(self.emptyUniverse, 2, None))
        self.saturatedTestUniverse.setMasses({1: 1, 2: 1})
        self.assertEqual(1.6138142499798294e-07, predictShearModulus(
            self.saturatedTestUniverse, 2, strandLength=2))

    def testWeightFractionCalculations(self):
        self.assertDictEqual(
            {}, computeWeightFractions(self.emptyUniverse))
        self.testUniverse.setMasses({1: 1, 2: 1})
        weightFractions = computeWeightFractions(
            self.testUniverse)
        self.assertDictEqual(weightFractions, {1: 1-3./8., 2: 3./8.})

    def testSolubleMaterialWeightFractionCalculation(self):
        self.testUniverse.setMasses({1: 1, 2: 1})
        self.assertRaises(NotImplementedError, lambda: computeWeightFractionOfSolubleMaterial(
            self.testUniverse, 2, functionalityPerType={1: 2, 2: 2}))
        self.assertRaises(NotImplementedError, lambda: computeWeightFractionOfSolubleMaterial(
            self.testUniverse, 2, functionalityPerType={1: 1, 2: 3}))
        self.saturatedTestUniverse.setMasses({1: 1, 2: 1})
        self.assertEqual((0.25407891551682393, {1: 0.85, 2: 0.15}, 0.8888888888888888, 0.4183006535947712), computeWeightFractionOfSolubleMaterial(
            self.saturatedTestUniverse, 2, strandLength=2))

    def testProbabilityCalculations(self):
        self.assertRaises(
            ValueError, lambda: computeMMsProbabilities(0.9, 2, 2))
        self.assertRaises(
            ValueError, lambda: computeMMsProbabilities(0.1, 0.9, 3))

    def testBackboneWeightFractionCalculations(self):
        self.assertEqual(0, calculateWeightFractionOfBackbone(
            self.emptyUniverse, 2, {}))
        self.assertEqual(1, calculateWeightFractionOfDanglingChains(
            self.emptyUniverse, 2, {}))
        self.saturatedTestUniverse.setMasses({1: 1, 2: 0})
        bb = calculateWeightFractionOfBackbone(
            self.saturatedTestUniverse, junctionType=2, strandLength=2)
        self.assertEqual(0.41013824884792627, bb)
        self.saturatedTestUniverse.setMasses({1: 1, 2: 0})
        self.assertEqual(1-bb, calculateWeightFractionOfDanglingChains(
            self.saturatedTestUniverse, 2, strandLength=2))

        self.saturatedTestUniverse.setMasses({1: 1, 2: 1})
        # test also as if the functionality was 4
        self.assertRaises(ValueError, lambda: calculateWeightFractionOfBackbone(self.saturatedTestUniverse, junctionType=2, functionalityPerType={
            1: 2, 2: 4
        }))
        # NOTE: requires a short strand length with these systems, as otherwise, r > 1 which is not supported by the formulas implemented
        self.assertEqual(0.4101065117216994, calculateWeightFractionOfBackbone(self.saturatedTestUniverse, junctionType=2, strandLength=2, functionalityPerType={
            1: 2, 2: 4
        }))


if __name__ == '__main__':
    unittest.main()
