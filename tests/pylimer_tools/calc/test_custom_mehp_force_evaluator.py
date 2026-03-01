import os
import sys
import unittest

import numpy as np

from pylimer_tools_cpp import (
    MEHPForceEvaluator,
    MEHPForceRelaxation,
    SimpleSpringMEHPForceEvaluator,
)

if __name__ == "__main__":
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../"))

from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class CustomForceEvaluator(MEHPForceEvaluator):
    """A custom force evaluator that implements a simple harmonic spring."""

    def __init__(self, kappa=1.0):
        super().__init__()
        self.kappa = kappa

    def evaluate_force_set_gradient(
            self, n, spring_distances, compute_gradient):
        """
        Evaluate force as 0.5 * kappa * sum(r_i^2 / N_i)
        where r_i is the spring length and N_i is the contour length.

        Returns a tuple (force, gradient) where gradient is a list or None.
        """
        network = self.network
        nr_springs = network.nr_of_springs

        # Compute force
        force = 0.0
        for i in range(nr_springs):
            spring_vec = spring_distances[3 * i: 3 * i + 3]
            r_squared = np.sum(spring_vec**2)
            contour_length = network.spring_contour_length[i]
            force += r_squared / contour_length
        force *= 0.5 * self.kappa

        # Compute gradient if required
        gradient = None
        if compute_gradient:
            gradient = [0.0] * n
            nrOfDim = 2 if self.is_2d else 3

            for j in range(nr_springs):
                a = network.spring_index_a[j]
                b = network.spring_index_b[j]
                contour_length = network.spring_contour_length[j]

                for dir_idx in range(nrOfDim):
                    spring_dist = spring_distances[3 * j + dir_idx]
                    grad_term = spring_dist * self.kappa / contour_length
                    gradient[3 * a + dir_idx] += grad_term
                    gradient[3 * b + dir_idx] -= grad_term

        return (force, gradient)

    def evaluate_stress_contribution(
            self, spring_distances, i, j, spring_index):
        """
        Evaluate stress contribution for spring at given indices.
        """
        network = self.network
        contour_length = network.spring_contour_length[spring_index]
        return self.kappa * spring_distances[i] * \
            spring_distances[j] / contour_length

    def prepare_for_evaluations(self):
        """Prepare any necessary data before evaluations."""
        pass


class TestCustomMEHPForceEvaluator(UniverseUsingTestCase):
    def test_custom_force_evaluator_instantiation(self):
        """Test that we can create a custom force evaluator."""
        custom_evaluator = CustomForceEvaluator(kappa=1.5)
        self.assertIsInstance(custom_evaluator, MEHPForceEvaluator)
        self.assertFalse(custom_evaluator.is_2d)

    def test_custom_force_evaluator_with_relaxation(self):
        """Test that custom force evaluator works with MEHPForceRelaxation."""
        custom_evaluator = CustomForceEvaluator(kappa=1.0)

        # Create force relaxation with custom evaluator
        force_relaxer = MEHPForceRelaxation(
            self.testUniverse,
            crosslinker_type=2,
            is_2d=False,
            force_evaluator=custom_evaluator,
        )

        # Run the relaxation
        force_relaxer.run_force_relaxation(
            algorithm="LD_LBFGS",
            max_nr_of_steps=1000,
            x_tolerance=1e-12,
            f_tolerance=1e-9,
        )

        # Check that simulation ran
        self.assertGreater(force_relaxer.get_nr_of_iterations(), 0)

        # Check that we can get results
        nr_springs = force_relaxer.get_nr_of_springs()
        self.assertGreater(nr_springs, 0)

        spring_lengths = force_relaxer.get_spring_lengths()
        self.assertEqual(len(spring_lengths), nr_springs)

    def test_custom_vs_builtin_evaluator(self):
        """Test that custom evaluator gives same results as built-in for same force law."""
        kappa = 1.0

        # Run with custom evaluator
        custom_evaluator = CustomForceEvaluator(kappa=kappa)
        force_relaxer_custom = MEHPForceRelaxation(
            self.testUniverse,
            crosslinker_type=2,
            is_2d=False,
            force_evaluator=custom_evaluator,
        )
        force_relaxer_custom.run_force_relaxation(
            algorithm="LD_LBFGS",
            max_nr_of_steps=5000,
            x_tolerance=1e-12,
            f_tolerance=1e-9,
        )

        # Run with built-in evaluator
        builtin_evaluator = SimpleSpringMEHPForceEvaluator(kappa=kappa)
        force_relaxer_builtin = MEHPForceRelaxation(
            self.testUniverse,
            crosslinker_type=2,
            is_2d=False,
            force_evaluator=builtin_evaluator,
        )
        force_relaxer_builtin.run_force_relaxation(
            algorithm="LD_LBFGS",
            max_nr_of_steps=5000,
            x_tolerance=1e-12,
            f_tolerance=1e-9,
        )

        # Compare results - spring lengths should be very close
        custom_lengths = np.array(force_relaxer_custom.get_spring_lengths())
        builtin_lengths = np.array(force_relaxer_builtin.get_spring_lengths())

        # They should be nearly identical
        np.testing.assert_allclose(
            custom_lengths, builtin_lengths, rtol=1e-4, atol=1e-6
        )

        # Compare gamma factors
        custom_gamma = force_relaxer_custom.get_gamma_factor()
        builtin_gamma = force_relaxer_builtin.get_gamma_factor()

        self.assertAlmostEqual(custom_gamma, builtin_gamma, places=4)

    def test_set_force_evaluator(self):
        """Test that we can set a force evaluator after creation."""
        # Create with default evaluator
        force_relaxer = MEHPForceRelaxation(self.testUniverse)

        # Set custom evaluator
        custom_evaluator = CustomForceEvaluator(kappa=2.0)
        force_relaxer.set_force_evaluator(custom_evaluator)

        # Run simulation
        force_relaxer.run_force_relaxation(
            algorithm="LD_LBFGS",
            max_nr_of_steps=1000,
            x_tolerance=1e-12,
            f_tolerance=1e-9,
        )

        # Should complete successfully
        self.assertGreater(force_relaxer.get_nr_of_iterations(), 0)


if __name__ == "__main__":
    unittest.main()
