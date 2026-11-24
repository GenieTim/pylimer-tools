import math
import unittest

import numpy as np

from pylimer_tools_cpp import (
    ResistanceDistanceAnalysis,
    RouseBrownianSimulator,
    RouseSolverType,
)


class TestResistanceDistanceAnalysis(unittest.TestCase):
    def test_linear_chain_and_laplacian_queries(self) -> None:
        n = 6
        rda = ResistanceDistanceAnalysis.create_linear_chain(n)
        self.assertEqual(rda.get_nr_of_nodes(), n)
        self.assertEqual(rda.get_nr_of_edges(), n - 1)
        self.assertAlmostEqual(rda.get_total_edge_weight(), float(n - 1))
        self.assertTrue(rda.is_connected())
        self.assertEqual(rda.get_nr_of_connected_components(), 1)
        self.assertEqual(len(rda.get_connected_components()), 1)

        lap = rda.get_laplacian()
        pinned_lap = rda.get_pinned_laplacian()
        self.assertEqual(lap.shape, (n, n))
        self.assertEqual(pinned_lap.shape, (n, n))

        for i in range(n):
            for j in range(n):
                self.assertAlmostEqual(
                    rda.compute_effective_resistance(i, j), float(abs(i - j))
                )

        pairs = [(0, 1), (0, 5), (2, 4), (3, 3)]
        batched = rda.compute_effective_resistance_batched(pairs)
        np.testing.assert_allclose(batched, [1.0, 5.0, 2.0, 0.0], atol=1e-9)

        expected_k = n * (n * n - 1) / 6.0
        self.assertAlmostEqual(rda.compute_kirchhoff_index(), expected_k)
        self.assertAlmostEqual(rda.compute_kirchhoff_index_spectral(), expected_k)

    def test_factories_and_matrices(self) -> None:
        ring = ResistanceDistanceAnalysis.create_ring(8)
        self.assertAlmostEqual(ring.compute_kirchhoff_index(), (8**3 - 8) / 12.0)

        complete = ResistanceDistanceAnalysis.create_complete(5)
        self.assertAlmostEqual(complete.compute_kirchhoff_index(), 4.0)

        star = ResistanceDistanceAnalysis.create_star(4, 1)
        self.assertAlmostEqual(star.compute_kirchhoff_index(), 16.0)

        all_r = star.compute_all_pairs_effective_resistance()
        self.assertEqual(all_r.shape, (5, 5))
        tot_r = star.compute_total_resistance_per_node()
        np.testing.assert_allclose(tot_r, all_r.sum(axis=1), atol=1e-9)

        centrality = star.compute_resistance_centrality()
        self.assertGreater(centrality[0], centrality[1])

        l_plus = star.compute_moore_penrose_pseudo_inverse()
        self.assertAlmostEqual(
            5.0 * float(np.trace(l_plus)), star.compute_kirchhoff_index()
        )

        commute = star.compute_all_pairs_commute_time()
        self.assertAlmostEqual(commute[0, 1], star.compute_commute_time(0, 1))

    def test_disconnected_and_validation(self) -> None:
        rda = ResistanceDistanceAnalysis([0, 2], [1, 3], [1.0, 2.0], 4)
        self.assertFalse(rda.is_connected())
        self.assertEqual(rda.get_nr_of_connected_components(), 2)
        self.assertTrue(math.isinf(rda.compute_effective_resistance(0, 2)))
        self.assertTrue(math.isinf(rda.compute_kirchhoff_index()))
        tot = rda.compute_total_resistance_per_node()
        self.assertTrue(np.all(np.isinf(tot)))

        with self.assertRaises(ValueError):
            ResistanceDistanceAnalysis([], [], [], 0)


class TestRouseBrownianSimulator(unittest.TestCase):
    def test_parameters_positions_and_solvers(self) -> None:
        sim = RouseBrownianSimulator.create_linear_chain(
            5, gamma=2.0, k=1.5, T=1.0, dt=0.02
        )
        self.assertEqual(sim.get_nr_of_nodes(), 5)
        self.assertAlmostEqual(sim.get_gamma(), 2.0)
        self.assertAlmostEqual(sim.get_k(), 1.5)
        self.assertAlmostEqual(sim.get_temperature(), 1.0)
        self.assertAlmostEqual(sim.get_dt(), 0.02)
        self.assertAlmostEqual(sim.get_kb(), 1.0)

        sim.set_kb(2.0)
        self.assertAlmostEqual(sim.get_kb(), 2.0)
        self.assertAlmostEqual(sim.get_diffusion_coefficient(), 1.0)
        sim.set_kb(1.0)

        sim.set_solver_type(RouseSolverType.CONJUGATE_GRADIENT)
        self.assertEqual(sim.get_solver_type(), RouseSolverType.CONJUGATE_GRADIENT)
        sim.set_solver_type(RouseSolverType.DIRECT_LDLT)

        sim.initialize_linear_positions(1.0, 0)
        self.assertAlmostEqual(sim.get_end_to_end_distance(), 4.0)
        np.testing.assert_allclose(
            sim.get_end_to_end_vector(), [4.0, 0.0, 0.0], atol=1e-12
        )

        flat = sim.get_positions_flattened(interleaved=True)
        self.assertEqual(flat.shape, (15,))
        sim.initialize_random_walk(1.0, seed=7)
        sim.set_positions_flattened(flat, interleaved=True)
        self.assertAlmostEqual(sim.get_end_to_end_distance(), 4.0)

    def test_stepping_and_run(self) -> None:
        sim = RouseBrownianSimulator.create_linear_chain(
            6, gamma=1.0, k=1.0, T=0.0, dt=0.05
        )
        sim.initialize_linear_positions(1.0, 0)
        e0 = sim.get_potential_energy()
        sim.step_deterministic(10)
        self.assertLess(sim.get_potential_energy(), e0)

        forces = np.zeros((6, 3))
        forces[0, 0] = 6.0
        cm0 = sim.get_center_of_mass().copy()
        sim.step_deterministic_with_external_forces(forces, 20)
        cm1 = sim.get_center_of_mass()
        self.assertAlmostEqual(cm1[0] - cm0[0], 20 * 0.05 * 1.0)

        sim.set_temperature(1.0)
        sim.set_seed(123)
        res = sim.run(total_steps=40, sample_interval=10, save_trajectory=True)
        self.assertEqual(len(res.times), 5)
        self.assertEqual(len(res.trajectory), 5)

        msd = RouseBrownianSimulator.compute_msd(res.trajectory, time_averaged=False)
        np.testing.assert_allclose(msd, res.msd, atol=1e-10)

        cm_msd = RouseBrownianSimulator.compute_center_of_mass_msd(
            res.trajectory, time_averaged=False
        )
        np.testing.assert_allclose(cm_msd, res.msd_center_of_mass, atol=1e-10)

        ee_corr = RouseBrownianSimulator.compute_end_to_end_autocorrelation(
            res.trajectory, time_averaged=True
        )
        self.assertAlmostEqual(float(ee_corr[0]), 1.0)


if __name__ == "__main__":
    unittest.main()
