import os
import sys

from pylimer_tools_cpp import (
    MEHPForceRelaxation,
    NonGaussianSpringForceEvaluator,
    SimpleSpringMEHPForceEvaluator,
)

if __name__ == "__main__":
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../"))

from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestMEHPSimulatorBindings(UniverseUsingTestCase):
    def test_simulator_instantiations(self):
        # test to make sure that none of the following raise an exception
        force_relaxer = MEHPForceRelaxation(self.testUniverse)
        # gaussianForceEvaluator = forceRelaxer.getForceEvaluator()
        gaussian_force_evaluator = SimpleSpringMEHPForceEvaluator(1.0)
        self.assertIsInstance(
            gaussian_force_evaluator,
            SimpleSpringMEHPForceEvaluator)
        self.assertFalse(gaussian_force_evaluator.is_2d)
        force_relaxer.set_force_evaluator(gaussian_force_evaluator)
        non_gaussian_force_evaluator = NonGaussianSpringForceEvaluator(
            1.0, 20.0, 0.98)
        self.assertFalse(non_gaussian_force_evaluator.is_2d)
        force_relaxer_2 = MEHPForceRelaxation(
            self.testUniverse, 2, False, non_gaussian_force_evaluator
        )
        self.assertEqual(force_relaxer_2.get_nr_of_iterations(), 0)
