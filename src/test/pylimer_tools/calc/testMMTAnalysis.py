
import unittest

import pandas as pd
from pylimer_tools.calc.doMMTAnalysis import *
from pylimer_tools.entities.universum import Universum


class TestMMTAnalysisFunctions(unittest.TestCase):

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

    def setUp(self):
        self.testUniverse = Universum([10, 10, 10])
        self.testUniverse.addAtomBondData(self.testAtoms, self.testBonds)
        self.emptyUniverse = Universum([10, 10, 10])
        # an additional larget test universe where the stoichiometric inbalance is < 1
        # even when imposing a crosslinker functionality of 1
        # in essence, it is on loop around 4 plus a connction to 6.
        self.largeTestUniverse = Universum([10, 10, 10])
        self.largeTestUniverse.addAtomBondData(self.testAtoms.append([
            {"id": 9, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 10, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 11, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 12, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 13, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 14, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 15, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 16, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 17, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 18, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 19, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 20, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
        ]), self.testBonds.append([
            {"to": 9, "bondFrom": 4},
            {"to": 10, "bondFrom": 9},
            {"to": 11, "bondFrom": 10},
            {"to": 12, "bondFrom": 11},
            {"to": 13, "bondFrom": 12},
            {"to": 4, "bondFrom": 13},
            {"to": 14, "bondFrom": 4},
            {"to": 15, "bondFrom": 14},
            {"to": 16, "bondFrom": 15},
            {"to": 17, "bondFrom": 16},
            {"to": 18, "bondFrom": 17},
            {"to": 19, "bondFrom": 18},
            {"to": 20, "bondFrom": 19},
            {"to": 6, "bondFrom": 20},
        ]))

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
            self.largeTestUniverse, 2, {1: 1, 2: 1}, strandLength=2))

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
            self.largeTestUniverse, 2, {1: 1, 2: 1}, strandLength=2))

    def testProbabilityCalculations(self):
        self.assertRaises(ValueError, lambda: computeMMsProbabilities(0.9, 2, 2))
        self.assertRaises(ValueError, lambda: computeMMsProbabilities(0.1, 0.9, 2))

    def testBackboneWeightFractionCalculations(self):
        self.assertEqual(0, calculateWeightFractionOfBackbone(
            self.emptyUniverse, 2, {}))
        self.assertEqual(1, calculateWeightFractionOfDanglingChains(
            self.emptyUniverse, 2, {}))
        bb = calculateWeightFractionOfBackbone(
            self.largeTestUniverse, 2, {1: 1, 2: 0}, strandLength=2)
        self.assertEqual(0.33642650971392124, bb)
        self.assertEqual(1-bb, calculateWeightFractionOfDanglingChains(
            self.largeTestUniverse, 2, {1: 1, 2: 0}, strandLength=2))

        # test also as if the functionality was 4
        self.assertRaises(ValueError, lambda: calculateWeightFractionOfBackbone(self.largeTestUniverse, junctionType=2, weightPerType={1: 1, 2: 1}, functionalityPerType={
            1: 2, 2: 4
        }))
        # NOTE: requires a short strand length with these systems, as otherwise, r > 1 which is not supported by the formulas implemented
        self.assertEqual(0.43094694702110886, calculateWeightFractionOfBackbone(self.largeTestUniverse, junctionType=2, strandLength=2, weightPerType={1: 1, 2: 1}, functionalityPerType={
            1: 2, 2: 4
        }))
