import copy
import os
import sys
import unittest

from pylimer_tools.calc.doMMTAnalysis import *
from pylimer_tools.io.unitStyles import UnitStyleFactory

if __name__ == '__main__':
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../"))
from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestMMTAnalysisFunctions(UniverseUsingTestCase):

    def testStoichiometricInbalance(self):
        self.assertAlmostEqual(
            0, computeStoichiometricInbalance(self.emptyUniverse, 2))
        self.assertAlmostEqual(
            (1*2 + 1*3 + 0)/(4*2 + 1*1), computeStoichiometricInbalance(self.testUniverse, 2, strandLength=1, effective=True))
        self.assertAlmostEqual(
            (3*3)/(5*2), computeStoichiometricInbalance(self.testUniverse, 2, strandLength=1, functionalityPerType={
                1: 2, 2: 3
            }))
        self.assertAlmostEqual(
            (3*3)/(5*2), computeStoichiometricInbalance(self.testUniverse, 2, strandLength=1))
        self.assertAlmostEqual(
            (3*3)/(2), computeStoichiometricInbalance(self.testUniverse, 2, strandLength=5, functionalityPerType={
                1: 2, 2: 3
            }))
        self.assertRaises(
            ValueError, computeStoichiometricInbalance, self.testUniverse, 7)
        self.assertAlmostEqual(computeStoichiometricInbalance(
            self.testUniverse, 2, functionalityPerType={
                1: 2, 2: 3
            }), ((3*3)/((5*2)/(5/3))))
        self.testUniverse.addAtoms([100], [3], [0], [0], [0], [0], [0], [0])
        self.assertAlmostEqual(computeStoichiometricInbalance(
            self.testUniverse, 2, ignoreTypes=[3]), ((3*3)/((5*2)/(5/3))))

    def testExtentOfReaction(self):
        self.assertAlmostEqual(
            1.0, computeExtentOfReaction(self.emptyUniverse, 2))
        self.assertAlmostEqual(
            5.0/6.0, computeExtentOfReaction(self.testUniverse, 2))

    def testGelationPointPrediction(self):
        self.assertAlmostEqual(1, predictGelationPoint(1, 2))
        self.assertAlmostEqual(1, predictGelationPoint(1, 2, 2))

    def testShearModulusPrediction(self):
        self.assertRaises(ValueError, lambda: predictShearModulus(
            self.emptyUniverse, 2, None))
        self.assertRaises(ValueError, lambda: predictShearModulus(
            self.emptyUniverse, None))
        self.saturatedTestUniverse.setMasses({1: 1, 2: 1})

        unitStyleFactory = UnitStyleFactory()
        unitStyle = unitStyleFactory.getUnitStyle("si")
        # TODO: find literature motiviation for results fo the functions
        self.assertAlmostEqual(0.13734491693339432, predictShearModulus(
            self.saturatedTestUniverse, unitStyle, crosslinkerType=2, strandLength=2).to('MPa').magnitude)
        self.assertAlmostEqual(0.13734491693339432, predictShearModulus(
            self.saturatedTestUniverse, unitStyle, crosslinkerType=2, strandLength=2, functionalityPerType={2: 4}).to('MPa').magnitude)
        self.saturatedTestUniverse.setMasses({1: 1, 2: 1})
        self.assertAlmostEqual(0.13734491693339432, predictShearModulus(
            self.saturatedTestUniverse, unitStyle, crosslinkerType=2).to('MPa').magnitude)

    def testPredictNumberDensityOfJunctionPoints(self):
        self.testUniverse.setMasses({1: 1, 2: 1})
        self.assertAlmostEqual(
            1.5, computeStoichiometricInbalance(self.testUniverse, 2))
        # TODO: find literature motiviation for results fo the functions
        self.assertAlmostEqual(predictNumberDensityOfJunctionPoints(
            self.testUniverse, 2), 0.0)
        self.assertRaises(NotImplementedError, lambda: predictNumberDensityOfJunctionPoints(
            self.testUniverse, 2, functionalityPerType={1: 2, 2: 5}))
        self.assertAlmostEqual(predictNumberDensityOfJunctionPoints(
            self.testUniverse, 2, functionalityPerType={1: 2, 2: 3}), 0.0)
        self.saturatedTestUniverse.setMasses({1: 1, 2: 1})
        self.assertAlmostEqual(predictNumberDensityOfJunctionPoints(
            self.saturatedTestUniverse, 2, functionalityPerType={1: 2, 2: 4}), 0.08304616298035744)

    def testPredictNumberDensityOfNetworkStrands(self):
        self.testUniverse.setMasses({1: 1, 2: 1})
        # TODO: find literature motiviation for results fo the functions
        self.assertAlmostEqual(predictNumberDensityOfNetworkStrands(
            self.testUniverse, 2), 0.0)
        self.assertRaises(NotImplementedError, lambda: predictNumberDensityOfNetworkStrands(
            self.testUniverse, 2, functionalityPerType={1: 2, 2: 5}))
        self.assertAlmostEqual(predictNumberDensityOfNetworkStrands(
            self.testUniverse, 2, functionalityPerType={1: 2, 2: 3}), 0.0)
        self.saturatedTestUniverse.setMasses({1: 1, 2: 1})
        self.assertAlmostEqual(predictNumberDensityOfNetworkStrands(
            self.saturatedTestUniverse, 2, functionalityPerType={1: 2, 2: 4}), 0.13752853164771697)

    def testWeightFractionCalculations(self):
        self.assertDictEqual(
            {}, computeWeightFractions(self.emptyUniverse))

        self.assertRaises(IndexError, lambda: computeWeightFractions(
            self.testUniverse))

        self.testUniverse.setMasses({1: 1, 2: 1})
        weightFractions = computeWeightFractions(
            self.testUniverse)
        self.assertDictEqual(weightFractions, {1: 1-3./8., 2: 3./8.})
        testUniverseCopy = copy.copy(self.testUniverse)
        testUniverseCopy.removeAtoms([1, 2, 3, 4, 5, 6])
        self.assertTrue(testUniverseCopy.getNrOfAtoms() == 2)
        self.assertDictEqual(computeWeightFractions(
            testUniverseCopy), {1: 1./2., 2: 1./2.})

    def testSolubleWeightFractionMeasurement(self):
        self.testUniverse.setMasses({1: 1, 2: 1})
        self.assertEqual(measureWeightFractioOfSolubleMaterial(
            self.emptyUniverse), None)
        self.assertEqual(measureWeightFractioOfSolubleMaterial(
            self.testUniverse, relTol=0), 0.0)
        self.assertEqual(measureWeightFractioOfSolubleMaterial(
            self.testUniverse, absTol=0), 0.0)
        self.assertEqual(measureWeightFractioOfSolubleMaterial(
            self.testUniverse, absTol=10000), 1.0)
        self.assertEqual(measureWeightFractioOfSolubleMaterial(
            self.testUniverse), 1/8)

    def testSolubleMaterialWeightFractionCalculation(self):
        self.testUniverse.setMasses({1: 1, 2: 1})
        self.assertRaises(NotImplementedError, lambda: computeWeightFractionOfSolubleMaterial(
            self.testUniverse, 2, functionalityPerType={1: 2, 2: 2}))
        self.assertRaises(NotImplementedError, lambda: computeWeightFractionOfSolubleMaterial(
            self.testUniverse, 2, functionalityPerType={1: 1, 2: 3}))
        self.saturatedTestUniverse.setMasses({1: 1, 2: 1})
        resTuple = computeWeightFractionOfSolubleMaterial(
            self.saturatedTestUniverse, 2)
        expectedTuple = (0.010699588477366243, {
                         1: 0.85, 2: 0.15}, 0.111111111111111, 0.111111111111111)
        self.assertEqual(len(resTuple), len(expectedTuple))
        for i in range(len(resTuple)):
            if (isinstance(resTuple[i], int) or isinstance(resTuple[i], float)):
                self.assertAlmostEqual(resTuple[i], expectedTuple[i])
            elif (isinstance(resTuple[i], dict)):
                for key in resTuple[i].keys():
                    self.assertAlmostEqual(
                        resTuple[i][key], expectedTuple[i][key])
            else:
                raise ValueError(
                    "Expected integer, float or dict for comparison")

    # def testProbabilityCalculations(self):
    #     self.assertRaises(
    #         ValueError, lambda: computeMMsProbabilities(0.9, 2, 2))
    #     self.assertRaises(
    #         ValueError, lambda: computeMMsProbabilities(0.1, 0.9, 3))

    def testBackboneWeightFractionCalculations(self):
        self.assertEqual(0, calculateWeightFractionOfBackbone(
            self.emptyUniverse, 2, {}))
        self.assertEqual(1, calculateWeightFractionOfDanglingChains(
            self.emptyUniverse, 2, {}))
        self.saturatedTestUniverse.setMasses({1: 1, 2: 0})
        self.assertAlmostEqual(1.0, computeExtentOfReaction(
            self.saturatedTestUniverse, 2))
        # TODO: get some literature backed values to test for
        bb = calculateWeightFractionOfBackbone(
            self.saturatedTestUniverse, crosslinkerType=2)
        self.assertAlmostEqual(0.8, bb)
        self.saturatedTestUniverse.setMasses({1: 1, 2: 0})
        self.assertEqual(1-bb, calculateWeightFractionOfDanglingChains(
            self.saturatedTestUniverse, 2))

        self.saturatedTestUniverse.setMasses({1: 1, 2: 1})
        # test also as if the functionality was 4
        # self.assertRaises(ValueError, lambda: calculateWeightFractionOfBackbone(self.saturatedTestUniverse, crosslinkerType=2, functionalityPerType={
        #     1: 2, 2: 4
        # }))
        # NOTE: requires a short strand length with these systems, as otherwise, r > 1 which is not supported by the formulas implemented
        self.assertEqual(0.14931407018789813, calculateWeightFractionOfBackbone(self.saturatedTestUniverse, crosslinkerType=2, strandLength=2, functionalityPerType={
            1: 2, 2: 4
        }))

    def testModulusDecompositions(self):
        unitStyleFactory = UnitStyleFactory()
        unitStyle = unitStyleFactory.getUnitStyle("si")
        self.assertAlmostEqual(unitStyle.kb.to('J/K').magnitude, 1.381e-23)

        # these results are pretty certain, align with experimental results, confirmed
        G_MMT_phantom, G_MMT_entanglement, G_ANM, G_PNM = computeModulusDecomposition(
            network=None, unitStyle=unitStyle, crosslinkerType=2, r=1., p=0.95, f=4, nu=4.69218e25*(unitStyle.getUnderlyingUnitRegistry()(
                'meter')**-3), T=298*unitStyle.getUnderlyingUnitRegistry()('kelvin')
        )

        alpha, beta = computeMMsProbabilities(r=1., p=0.95, f=4.)
        self.assertAlmostEqual(alpha, 0.0983588, places=5)

        self.assertAlmostEqual(G_ANM.to('MPa').magnitude, 0.193101, places=5)
        self.assertAlmostEqual(G_PNM.to('MPa').magnitude, 0.0965506, places=5)
        self.assertAlmostEqual(G_MMT_entanglement.to(
            'MPa').magnitude, 0.17851, places=5)
        self.assertAlmostEqual(G_MMT_phantom.to(
            'MPa').magnitude, 0.0777321, places=5)

        # these in turn require further investigation into the involvement of r
        G_MMT_phantom, G_MMT_entanglement, G_ANM, G_PNM = computeModulusDecomposition(
            network=None, unitStyle=unitStyle, crosslinkerType=2, r=1.3, p=0.6465, f=4, nu=1.25981e25*(unitStyle.getUnderlyingUnitRegistry()(
                'meter')**-3), T=298*unitStyle.getUnderlyingUnitRegistry()('kelvin')
        )

        self.assertAlmostEqual(G_ANM.to('MPa').magnitude, 0.051846, places=5)
        self.assertAlmostEqual(G_PNM.to('MPa').magnitude, 0.025923, places=5)
        self.assertAlmostEqual(G_MMT_entanglement.to(
            'MPa').magnitude, 0.05433807, places=5)
        self.assertAlmostEqual(G_MMT_phantom.to(
            'MPa').magnitude, 0.00492674, places=5)

    def testProbabilityCalculations(self):
        # most already in other tests
        alpha, beta = computeMMsProbabilities(r=1.0, p=0.85, f=5)
        self.assertAlmostEqual(0.282074, alpha, places=5)
        alpha, beta = computeMMsProbabilities(r=1.0, p=0.85, f=6)
        self.assertAlmostEqual(0.278715, alpha, places=5)
        self.assertRaises(
            ValueError, lambda: computeMMsProbabilities(r=1.0, p=2, f=7))


if __name__ == '__main__':
    unittest.main()
