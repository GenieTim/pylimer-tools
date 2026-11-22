import gc
import os
import sys
import unittest

import numpy as np

from pylimer_tools_cpp import (
    Box,
    HookeanForceEvaluator,
    MEHPForceEvaluator,
    MEHPNetwork,
    MEHPNonlinearSolverConfig,
    MEHPSimulator,
    NLoptAlgorithm,
    NonlinearOptimizationSolver,
)

if __name__ == "__main__":
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../"))

from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class CustomHookeanEvaluator(MEHPForceEvaluator):
    """Stateless custom evaluator implementing only per-spring energy and force."""

    def __init__(self, kappa: float = 1.0) -> None:
        super().__init__()
        self.kappa = float(kappa)

    def get_kappa(self) -> float:
        return self.kappa

    def evaluate_spring_energy(
        self, spring_vec: np.ndarray, contour_length: float
    ) -> float:
        vec = np.asarray(spring_vec, dtype=np.float64)
        return float(0.5 * self.kappa * np.dot(vec, vec) / contour_length)

    def evaluate_spring_force(
        self, spring_vec: np.ndarray, contour_length: float
    ) -> np.ndarray:
        vec = np.asarray(spring_vec, dtype=np.float64)
        return (self.kappa / contour_length) * vec


class VectorizedHookeanEvaluator(MEHPForceEvaluator):
    """Stateless custom evaluator overriding vectorized network energy/gradient."""

    def __init__(self, kappa: float = 1.0) -> None:
        super().__init__()
        self.kappa = float(kappa)
        self.vectorized_call_count = 0

    def get_kappa(self) -> float:
        return self.kappa

    def evaluate_spring_energy(
        self, spring_vec: np.ndarray, contour_length: float
    ) -> float:
        vec = np.asarray(spring_vec, dtype=np.float64)
        return float(0.5 * self.kappa * np.dot(vec, vec) / contour_length)

    def evaluate_spring_force(
        self, spring_vec: np.ndarray, contour_length: float
    ) -> np.ndarray:
        vec = np.asarray(spring_vec, dtype=np.float64)
        return (self.kappa / contour_length) * vec

    def evaluate_network_energy_and_gradient(
        self,
        net: MEHPNetwork,
        spring_vectors: np.ndarray,
        is_2d: bool,
        compute_gradient: bool,
    ) -> tuple[float, np.ndarray | None]:
        self.vectorized_call_count += 1
        nr_springs = net.nr_of_springs
        vecs = (
            np.asarray(spring_vectors, dtype=np.float64).reshape(nr_springs, 3).copy()
        )
        if is_2d:
            vecs[:, 2] = 0.0
        inv_n = 1.0 / np.asarray(net.spring_contour_length, dtype=np.float64)
        energy = float(0.5 * self.kappa * np.sum(np.sum(vecs * vecs, axis=1) * inv_n))
        if not compute_gradient:
            return (energy, None)

        forces = (self.kappa * inv_n)[:, None] * vecs
        grad = np.zeros((net.nr_of_links, 3), dtype=np.float64)
        idx_a = np.asarray(net.spring_index_a, dtype=np.int64)
        idx_b = np.asarray(net.spring_index_b, dtype=np.int64)
        valid = (idx_a >= 0) & (idx_b >= 0)
        np.add.at(grad, idx_a[valid], -forces[valid])
        np.add.at(grad, idx_b[valid], forces[valid])
        if is_2d:
            grad[:, 2] = 0.0
        return (energy, grad.reshape(-1))


class QuarticNonlinearEvaluator(MEHPForceEvaluator):
    """Stateless anharmonic spring evaluator E = 0.5*k/N*r^2 + 0.25*a/N^3*r^4."""

    def __init__(self, kappa: float = 1.0, alpha: float = 0.25) -> None:
        super().__init__()
        self.kappa = float(kappa)
        self.alpha = float(alpha)

    def get_kappa(self) -> float:
        return self.kappa

    def evaluate_spring_energy(
        self, spring_vec: np.ndarray, contour_length: float
    ) -> float:
        vec = np.asarray(spring_vec, dtype=np.float64)
        r2 = float(np.dot(vec, vec))
        return 0.5 * self.kappa * r2 / contour_length + 0.25 * self.alpha * (
            r2 * r2
        ) / (contour_length**3)

    def evaluate_spring_force(
        self, spring_vec: np.ndarray, contour_length: float
    ) -> np.ndarray:
        vec = np.asarray(spring_vec, dtype=np.float64)
        r2 = float(np.dot(vec, vec))
        coeff = (self.kappa / contour_length) + (self.alpha * r2 / (contour_length**3))
        return coeff * vec


class TestCustomMEHPForceEvaluator(UniverseUsingTestCase):
    def test_custom_force_evaluator_instantiation(self) -> None:
        """Test stateless custom force evaluator properties and single-spring methods."""
        custom_evaluator = CustomHookeanEvaluator(kappa=1.5)
        self.assertIsInstance(custom_evaluator, MEHPForceEvaluator)
        self.assertFalse(custom_evaluator.is_linear())
        self.assertAlmostEqual(custom_evaluator.get_kappa(), 1.5)

        spring_vec = np.array([1.2, -0.4, 0.8])
        contour_length = 4.0
        builtin = HookeanForceEvaluator(1.5)

        self.assertAlmostEqual(
            custom_evaluator.evaluate_spring_energy(spring_vec, contour_length),
            builtin.evaluate_spring_energy(spring_vec, contour_length),
        )
        np.testing.assert_allclose(
            custom_evaluator.evaluate_spring_force(spring_vec, contour_length),
            builtin.evaluate_spring_force(spring_vec, contour_length),
        )
        np.testing.assert_allclose(
            custom_evaluator.evaluate_spring_stiffness(spring_vec, contour_length),
            builtin.evaluate_spring_stiffness(spring_vec, contour_length),
            atol=1e-7,
        )
        self.assertAlmostEqual(
            custom_evaluator.evaluate_stress_contribution(
                spring_vec, 0, 2, contour_length
            ),
            builtin.evaluate_stress_contribution(spring_vec, 0, 2, contour_length),
        )

    def test_custom_force_evaluator_with_simulator(self) -> None:
        """Test stateless custom force evaluator works with MEHPSimulator."""
        custom_evaluator = CustomHookeanEvaluator(kappa=1.0)
        net = MEHPNetwork.from_universe(self.testUniverse)
        sim = MEHPSimulator(net, custom_evaluator, False)

        config = MEHPNonlinearSolverConfig()
        config.algorithm = NLoptAlgorithm.LD_LBFGS
        config.max_nr_of_steps = 1000
        config.xtol = 1e-12
        config.ftol = 1e-9

        result = sim.run(config)

        self.assertTrue(sim.get_simulation_has_run())
        self.assertGreater(result.nr_of_iterations, 0)
        self.assertGreater(sim.get_nr_of_iterations(), 0)

        nr_springs = sim.get_nr_of_springs()
        self.assertGreater(nr_springs, 0)
        self.assertEqual(len(sim.get_spring_lengths()), nr_springs)

    def test_custom_vs_builtin_evaluator(self) -> None:
        """Test that stateless custom evaluator matches HookeanForceEvaluator."""
        kappa = 1.0
        config = MEHPNonlinearSolverConfig()
        config.algorithm = NLoptAlgorithm.LD_LBFGS
        config.max_nr_of_steps = 5000
        config.xtol = 1e-12
        config.ftol = 1e-12

        net_custom = MEHPNetwork.from_universe(self.testUniverse)
        sim_custom = MEHPSimulator(
            net_custom, CustomHookeanEvaluator(kappa=kappa), False
        )
        sim_custom.run(config)

        net_builtin = MEHPNetwork.from_universe(self.testUniverse)
        sim_builtin = MEHPSimulator(net_builtin, HookeanForceEvaluator(kappa), False)
        sim_builtin.run(config)

        np.testing.assert_allclose(
            np.array(sim_custom.get_spring_lengths()),
            np.array(sim_builtin.get_spring_lengths()),
            rtol=1e-4,
            atol=1e-6,
        )
        self.assertAlmostEqual(
            sim_custom.get_gamma_factor(), sim_builtin.get_gamma_factor(), places=4
        )
        self.assertAlmostEqual(
            sim_custom.get_energy(), sim_builtin.get_energy(), places=6
        )
        np.testing.assert_allclose(
            sim_custom.get_stress_tensor(),
            sim_builtin.get_stress_tensor(),
            rtol=1e-4,
            atol=1e-6,
        )

    def test_vectorized_network_energy_and_gradient_override(self) -> None:
        """Test overriding evaluate_network_energy_and_gradient in Python."""
        kappa = 1.5
        vec_eval = VectorizedHookeanEvaluator(kappa=kappa)
        builtin_eval = HookeanForceEvaluator(kappa)

        net = MEHPNetwork.from_universe(self.saturatedTestUniverse)
        sim_ref = MEHPSimulator(net, builtin_eval, False)
        sim_ref.deform_to(Box(11.5, 9.5, 10.2))
        spring_vecs = sim_ref.get_spring_vectors()
        deformed_net = sim_ref.get_network()

        for is_2d in (False, True):
            e_vec, g_vec = vec_eval.evaluate_network_energy_and_gradient(
                deformed_net, spring_vecs, is_2d, True
            )
            e_ref, g_ref = builtin_eval.evaluate_network_energy_and_gradient(
                deformed_net, spring_vecs, is_2d, True
            )
            self.assertAlmostEqual(e_vec, e_ref, places=10)
            self.assertIsNotNone(g_vec)
            self.assertIsNotNone(g_ref)
            np.testing.assert_allclose(g_vec, g_ref, rtol=1e-10, atol=1e-10)

            e_vec_ng, g_vec_ng = vec_eval.evaluate_network_energy_and_gradient(
                deformed_net, spring_vecs, is_2d, False
            )
            self.assertAlmostEqual(e_vec_ng, e_ref, places=10)
            self.assertIsNone(g_vec_ng)

        vec_eval.vectorized_call_count = 0
        sim_vec = MEHPSimulator(deformed_net, vec_eval, False)
        config = MEHPNonlinearSolverConfig()
        config.algorithm = NLoptAlgorithm.LD_LBFGS
        config.max_nr_of_steps = 5000
        config.xtol = 1e-12
        config.ftol = 1e-12
        sim_vec.run(NonlinearOptimizationSolver(config))

        self.assertGreater(vec_eval.vectorized_call_count, 0)

        sim_ref.run()
        self.assertAlmostEqual(sim_vec.get_energy(), sim_ref.get_energy(), places=6)
        self.assertAlmostEqual(
            sim_vec.get_gamma_factor(), sim_ref.get_gamma_factor(), places=4
        )
        np.testing.assert_allclose(
            sim_vec.get_stress_tensor(),
            sim_ref.get_stress_tensor(),
            rtol=1e-4,
            atol=1e-6,
        )

    def test_set_force_evaluator_and_stateless_reuse(self) -> None:
        """Test set_force_evaluator and reusing one stateless evaluator across simulators."""
        sim1 = MEHPSimulator(self.testUniverse)
        shared_evaluator = CustomHookeanEvaluator(kappa=2.0)
        sim1.set_force_evaluator(shared_evaluator)
        self.assertAlmostEqual(sim1.get_force_evaluator().get_kappa(), 2.0)

        sim2 = MEHPSimulator(self.saturatedTestUniverse)
        sim2.set_force_evaluator(shared_evaluator)

        sim1.run()
        sim2.run()

        self.assertGreater(sim1.get_nr_of_iterations(), 0)
        self.assertGreater(sim2.get_nr_of_iterations(), 0)

    def test_nonlinear_custom_evaluator_and_lifetime(self) -> None:
        """Test genuine nonlinear evaluator, stiffness finite differences, and lifetime."""
        evaluator = QuarticNonlinearEvaluator(kappa=1.2, alpha=0.5)
        spring_vec = np.array([1.5, -0.5, 0.25])
        contour_length = 2.0
        stiffness = evaluator.evaluate_spring_stiffness(spring_vec, contour_length)
        r2 = float(np.dot(spring_vec, spring_vec))
        expected_stiffness = (
            1.2 / contour_length + 0.5 * r2 / (contour_length**3)
        ) * np.eye(3) + (2.0 * 0.5 / (contour_length**3)) * np.outer(
            spring_vec, spring_vec
        )
        np.testing.assert_allclose(stiffness, expected_stiffness, rtol=1e-5, atol=1e-6)

        def build_sim() -> MEHPSimulator:
            net = MEHPNetwork.from_universe(self.saturatedTestUniverse)
            return MEHPSimulator(net, QuarticNonlinearEvaluator(kappa=1.2, alpha=0.5))

        sim = build_sim()
        gc.collect()
        sim.deform_to(Box(11.0, 10.0, 10.0))
        res = sim.run()
        self.assertTrue(sim.get_simulation_has_run())
        self.assertGreater(res.nr_of_iterations, 0)
        self.assertGreater(sim.get_energy(), 0.0)
        stress = sim.get_stress_tensor()
        np.testing.assert_allclose(stress, stress.T, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
