import gc
import math
import os
import sys

import numpy as np

from pylimer_tools_cpp import (
    Box,
    MEHPForceRelaxation,
    MEHPForceBalance2,
    SLESolver,
    NonGaussianSpringForceEvaluator,
    SimpleSpringMEHPForceEvaluator,
    inverse_langevin,
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

    def test_inverse_langevin_symmetry(self):
        """Test langevin_inv symmetry: langevin_inv(-x) == -langevin_inv(x)."""
        test_values = [0.0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.999, 1.5]
        for val in test_values:
            pos = inverse_langevin(val)
            neg = inverse_langevin(-val)
            self.assertAlmostEqual(neg, -pos, places=7)
            if val > 0:
                self.assertGreater(pos, 0.0)
                self.assertLess(neg, 0.0)

    def test_box_is_valid_offset_negative(self):
        """Test Box::isValidOffset with negative offsets."""
        box = Box(10.0, 10.0, 10.0)
        # Periodic offset of -1 box length in x should be valid
        valid_offset = np.array([-10.0, 0.0, 0.0])
        self.assertTrue(box.is_valid_offset(valid_offset))

        # Offset of -1.5 box lengths in x should be invalid
        invalid_offset = np.array([-15.0, 0.0, 0.0])
        self.assertFalse(box.is_valid_offset(invalid_offset))

        # Mixed signs
        mixed_valid = np.array([-20.0, 10.0, -30.0])
        self.assertTrue(box.is_valid_offset(mixed_valid))

    def test_keep_alive_evaluator_lifetime(self):
        """Test keep_alive<1, 5> ensures evaluator stays alive with MEHPForceRelaxation."""

        def create_relaxer_with_temp_evaluator():
            evaluator = SimpleSpringMEHPForceEvaluator(2.5)
            relaxer = MEHPForceRelaxation(
                self.testUniverse, 2, False, evaluator)
            return relaxer

        relaxer = create_relaxer_with_temp_evaluator()
        # Force garbage collection to collect the temporary evaluator if not
        # kept alive
        gc.collect()
        # Running force relaxation must not segfault or crash
        relaxer.run_force_relaxation()
        self.assertGreater(relaxer.get_nr_of_iterations(), 0)

    def test_non_gaussian_stress_contribution_at_zero(self):
        """Test stress contribution at r=0 returns 0.0 without NaN."""
        evaluator = NonGaussianSpringForceEvaluator(1.0, 10.0, 1.0)
        stress = evaluator.evaluate_stress_contribution(
            [0.0, 0.0, 0.0], 0, 0, 0)
        self.assertEqual(stress, 0.0)
        self.assertFalse(math.isnan(stress))

    def test_metal_cg_force_balance_relaxation(self):
        """Test that MEHPForceBalance2 with METAL_CG completes force relaxation."""
        try:
            fb = MEHPForceBalance2(self.testUniverseSmall, 2, False)
            fb.run_force_relaxation(sle_solver=SLESolver.METAL_CG)
            self.assertGreater(fb.get_nr_of_iterations(), 0)
        except RuntimeError as e:
            if (
                "Metal acceleration is not enabled" in str(e)
                or "Metal GPU backend is not available" in str(e)
            ):
                self.skipTest(f"Metal not available: {e}")
            else:
                raise

    def test_metal_bb_force_balance_relaxation(self):
        """Test that MEHPForceBalance2 with METAL_BARZILAI_BORWEIN completes force relaxation."""
        try:
            fb = MEHPForceBalance2(self.testUniverseSmall, 2, False)
            fb.run_force_relaxation(sle_solver=SLESolver.METAL_BARZILAI_BORWEIN)
            self.assertGreater(fb.get_nr_of_iterations(), 0)
        except RuntimeError as e:
            if (
                "Metal acceleration is not enabled" in str(e)
                or "Metal GPU backend is not available" in str(e)
            ):
                self.skipTest(f"Metal not available: {e}")
            else:
                raise

