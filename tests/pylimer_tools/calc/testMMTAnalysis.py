
from test.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase

from pylimer_tools.calc.doMMTAnalysis import *
from pylimer_tools.entities.universum import Universum


class TestMMTAnalysisFunctions(UniverseUsingTestCase):

    def testStoichiometricInbalance(self):
        self.assertEqual(
            0, computeStoichiometricInbalance(self.emptyUniverse, 2))
        self.assertEqual(
            (3*3)/(5*2), computeStoichiometricInbalance(self.testUniverse, 2, strandLength=1))
        self.assertEqual(
            (3*3)/(2), computeStoichiometricInbalance(self.testUniverse, 2, strandLength=5))

    def testExtentOfReaction(self):
        self.assertEqual(1.0, computeExtentOfReaction(self.emptyUniverse))
        self.assertEqual(14.0/19.0, computeExtentOfReaction(self.testUniverse))

    def testGelationPointPrediction(self):
        self.assertEqual(1, predictGelationPoint(1, 2))
        self.assertEqual(1, predictGelationPoint(1, 2, 2))

    def testShearModulusPrediction(self):
        self.assertIsNone(predictShearModulus(self.emptyUniverse, 2, None))
        self.assertEqual(1.0521245069791487e-09, predictShearModulus(
            self.saturatedTestUniverse, 2, {1: 1, 2: 1}, strandLength=2))

    def testWeightFractionCalculations(self):
        self.assertDictEqual(
            {}, computeWeightFractions(self.emptyUniverse, {}))
        weightFractions = computeWeightFractions(
            self.testUniverse, {1: 1, 2: 1})
        self.assertDictEqual(weightFractions, {1: 1-3./8., 2: 3./8.})

    def testSolubleMaterialWeightFractionCalculation(self):
        self.assertRaises(NotImplementedError, lambda: computeWeightFractionOfSolubleMaterial(
            self.testUniverse, 2, {1: 1, 2: 1}, None, {1: 2, 2: 2}))
        self.assertRaises(NotImplementedError, lambda: computeWeightFractionOfSolubleMaterial(
            self.testUniverse, 2, {1: 1, 2: 1}, None, {1: 1, 2: 3}))
        self.assertEqual((0.3506977562162289, {1: 0.85, 2: 0.15}, 0.9799067775258254, 0.49652823066801094), computeWeightFractionOfSolubleMaterial(
            self.saturatedTestUniverse, 2, {1: 1, 2: 1}, strandLength=2))

    def testProbabilityCalculations(self):
        self.assertRaises(
            ValueError, lambda: computeMMsProbabilities(0.9, 2, 2))
        self.assertRaises(
            ValueError, lambda: computeMMsProbabilities(0.1, 0.9, 2))

    def testBackboneWeightFractionCalculations(self):
        self.assertEqual(0, calculateWeightFractionOfBackbone(
            self.emptyUniverse, 2, {}))
        self.assertEqual(1, calculateWeightFractionOfDanglingChains(
            self.emptyUniverse, 2, {}))
        bb = calculateWeightFractionOfBackbone(
            self.saturatedTestUniverse, 2, {1: 1, 2: 0}, strandLength=2)
        self.assertEqual(0.33642650971392124, bb)
        self.assertEqual(1-bb, calculateWeightFractionOfDanglingChains(
            self.saturatedTestUniverse, 2, {1: 1, 2: 0}, strandLength=2))

        # test also as if the functionality was 4
        self.assertRaises(ValueError, lambda: calculateWeightFractionOfBackbone(self.saturatedTestUniverse, junctionType=2, weightPerType={1: 1, 2: 1}, functionalityPerType={
            1: 2, 2: 4
        }))
        # NOTE: requires a short strand length with these systems, as otherwise, r > 1 which is not supported by the formulas implemented
        self.assertEqual(0.43094694702110886, calculateWeightFractionOfBackbone(self.saturatedTestUniverse, junctionType=2, strandLength=2, weightPerType={1: 1, 2: 1}, functionalityPerType={
            1: 2, 2: 4
        }))
