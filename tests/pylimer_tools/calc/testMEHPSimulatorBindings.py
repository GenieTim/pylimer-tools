import os
import sys

from pylimer_tools.calc.doMEHPAnalysis import *
from pylimer_tools_cpp.pylimer_tools_cpp import (
    MEHPForceRelaxation, NonGaussianSpringForceEvaluator,
    SimpleSpringMEHPForceEvaluator, Universe)

if __name__ == '__main__':
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../"))

from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestMEHPSimulatorBindings(UniverseUsingTestCase):
    def test_simulatorInstantiations(self):
        # test to make sure that none of the following raise an exception
        forceRelaxer = MEHPForceRelaxation(self.testUniverse)
        # gaussianForceEvaluator = forceRelaxer.getForceEvaluator()
        gaussianForceEvaluator = SimpleSpringMEHPForceEvaluator(1.0)
        self.assertIsInstance(gaussianForceEvaluator,
                              SimpleSpringMEHPForceEvaluator)
        self.assertFalse(gaussianForceEvaluator.is2D)
        forceRelaxer.setForceEvaluator(gaussianForceEvaluator)
        nonGaussianForceEvaluator = NonGaussianSpringForceEvaluator(
            1.0, 20.0, 0.98)
        self.assertFalse(nonGaussianForceEvaluator.is2D)
        forceRelaxer2 = MEHPForceRelaxation(
            self.testUniverse, 2, False, nonGaussianForceEvaluator)
