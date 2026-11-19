#!/usr/bin/env python
import pickle
import unittest

import numpy as np
from pylimer_tools_cpp import (
    Correlator,
    ResistanceDistanceAnalysis,
    RouseBrownianSimulator,
    RouseSolverType,
    Universe,
)
from scipy import sparse


def _make_chain_universe(num_atoms: int = 4) -> Universe:
    u = Universe(10.0, 10.0, 10.0)
    ids = list(range(1, num_atoms + 1))
    types = [1] * num_atoms
    xs = [float(i) for i in range(num_atoms)]
    zeros = [0.0] * num_atoms
    zeros_i = [0] * num_atoms
    u.add_atoms(ids, types, xs, zeros, zeros, zeros_i, zeros_i, zeros_i)
    if num_atoms > 1:
        u.add_bonds(ids[:-1], ids[1:])
    return u


class TestResistanceDistanceAnalysis(unittest.TestCase):
    def test_canonical_topologies_and_factories(self) -> None:
        chain = ResistanceDistanceAnalysis.create_linear_chain(5)
        self.assertEqual(chain.get_nr_of_nodes(), 5)
        self.assertEqual(chain.get_nr_of_edges(), 4)
        self.assertTrue(chain.is_connected())
        self.assertAlmostEqual(chain.compute_effective_resistance(0, 4), 4.0)
        self.assertAlmostEqual(
            chain.compute_kirchhoff_index(chunk_size=2),
            5 * (25 - 1) / 6.0,
            places=8,
        )

        # Batched & all-pairs & pseudo-inverse
        pairs = [(0, 4), (1, 3), (2, 2)]
        batched_default = chain.compute_effective_resistance_batched(pairs)
        batched_chunked = chain.compute_effective_resistance_batched(
            pairs, chunk_size=2
        )
        np.testing.assert_allclose(batched_default, [4.0, 2.0, 0.0], atol=1e-9)
        np.testing.assert_allclose(batched_chunked, batched_default, atol=1e-9)

        all_r = chain.compute_all_pairs_effective_resistance(chunk_size=2)
        self.assertEqual(all_r.shape, (5, 5))
        np.testing.assert_allclose(
            chain.compute_total_resistance_per_node(chunk_size=2),
            all_r.sum(axis=1),
            atol=1e-9,
        )
        np.testing.assert_allclose(
            chain.compute_resistance_centrality(chunk_size=2),
            chain.compute_resistance_centrality(),
            atol=1e-9,
        )
        np.testing.assert_allclose(
            chain.compute_all_pairs_commute_time(chunk_size=2),
            chain.compute_all_pairs_commute_time(),
            atol=1e-9,
        )
        lplus = chain.compute_moore_penrose_pseudo_inverse(chunk_size=2)
        np.testing.assert_allclose(lplus @ np.ones(5), np.zeros(5), atol=1e-10)

        with self.assertRaises(ValueError):
            chain.compute_effective_resistance_batched(pairs, chunk_size=0)
        with self.assertRaises(ValueError):
            chain.compute_resistance_centrality(chunk_size=0)
        with self.assertRaises(ValueError):
            chain.compute_all_pairs_commute_time(chunk_size=0)

        # From Universe and from_matrix
        u = _make_chain_universe(4)
        rda_u = ResistanceDistanceAnalysis(u)
        rda_from_u = ResistanceDistanceAnalysis.from_universe(u)
        rda_from_m = ResistanceDistanceAnalysis.from_matrix(rda_u.get_laplacian())
        for obj in (rda_u, rda_from_u, rda_from_m):
            self.assertAlmostEqual(obj.compute_effective_resistance(0, 3), 3.0)
            self.assertAlmostEqual(obj.compute_kirchhoff_index(), 10.0)

    def test_validation_rejects_invalid_laplacians(self) -> None:
        adj = sparse.csc_matrix(np.array([[0.0, 1.0], [1.0, 0.0]]))
        with self.assertRaises(ValueError):
            ResistanceDistanceAnalysis(adj)

        shifted = sparse.csc_matrix(np.array([[2.0, -1.0], [-1.0, 2.0]]))
        with self.assertRaises(ValueError):
            ResistanceDistanceAnalysis(shifted)


class TestRouseBrownianSimulator(unittest.TestCase):
    def test_factories_solvers_and_run(self) -> None:
        u = _make_chain_universe(4)
        sim_u = RouseBrownianSimulator(u, gamma=2.0, k=1.5, T=0.0, dt=0.05)
        sim_from_u = RouseBrownianSimulator.from_universe(
            u, gamma=2.0, k=1.5, T=0.0, dt=0.05
        )
        sim_from_m = RouseBrownianSimulator.from_matrix(
            sim_u.get_laplacian(), gamma=2.0, k=1.5, T=0.0, dt=0.05
        )
        self.assertEqual(sim_u.get_nr_of_nodes(), 4)
        self.assertAlmostEqual(sim_u.get_end_to_end_distance(), 3.0)
        self.assertAlmostEqual(
            sim_u.get_end_to_end_distance(node_a=0, node_b=None), 3.0
        )
        np.testing.assert_allclose(
            sim_u.get_end_to_end_vector(node_a=0, node_b=None),
            [3.0, 0.0, 0.0],
            atol=1e-12,
        )
        np.testing.assert_allclose(sim_u.get_positions(), sim_from_u.get_positions())
        self.assertEqual(sim_from_m.get_nr_of_nodes(), 4)

        sim_complete = RouseBrownianSimulator.create_complete(
            4, gamma=1.0, k=2.0, T=0.0, dt=0.01
        )
        sim_complete.initialize_linear_positions(1.0, 0)
        self.assertAlmostEqual(sim_complete.get_potential_energy(), 20.0)

        # Weighted springs
        sim_weighted = RouseBrownianSimulator(
            [0], [1], gamma=1.0, k=1.0, T=0.0, dt=0.01, num_nodes=2, weights=[3.0]
        )
        sim_weighted.initialize_linear_positions(2.0, 0)
        self.assertAlmostEqual(sim_weighted.get_potential_energy(), 6.0)

        # Stochastic run + trajectory post-processing with end_to_end_node_b=None
        sim_run = RouseBrownianSimulator.create_linear_chain(
            6, gamma=1.0, k=1.0, T=1.0, dt=0.01
        )
        sim_run.set_seed(7)
        sim_run.initialize_random_walk(1.0, seed=11)
        sim_run.set_solver_type(RouseSolverType.DIRECT_LLT)
        res = sim_run.run(
            total_steps=20,
            sample_interval=5,
            save_trajectory=True,
            end_to_end_node_a=0,
            end_to_end_node_b=None,
        )
        self.assertEqual(len(res.times), 5)
        self.assertEqual(len(res.trajectory), 5)
        msd = RouseBrownianSimulator.compute_msd(res.trajectory, time_averaged=False)
        np.testing.assert_allclose(msd, res.msd, atol=1e-10)
        cm_msd = RouseBrownianSimulator.compute_center_of_mass_msd(
            res.trajectory, time_averaged=False
        )
        np.testing.assert_allclose(cm_msd, res.msd_center_of_mass, atol=1e-10)
        ee_ac = RouseBrownianSimulator.compute_end_to_end_autocorrelation(
            res.trajectory, node_a=0, node_b=None, time_averaged=False
        )
        self.assertEqual(len(ee_ac), 5)
        self.assertAlmostEqual(ee_ac[0], 1.0, places=9)

        # DIRECT_LDLT and CONJUGATE_GRADIENT solvers match DIRECT_LLT
        for solver_type in (
            RouseSolverType.DIRECT_LDLT,
            RouseSolverType.CONJUGATE_GRADIENT,
        ):
            sim_alt = RouseBrownianSimulator.create_linear_chain(
                6, gamma=1.0, k=1.0, T=1.0, dt=0.01
            )
            sim_alt.set_seed(7)
            sim_alt.initialize_random_walk(1.0, seed=11)
            sim_alt.set_solver_type(solver_type)
            res_alt = sim_alt.run(
                total_steps=20,
                sample_interval=5,
                save_trajectory=True,
                end_to_end_node_b=None,
            )
            np.testing.assert_allclose(
                res_alt.trajectory[-1], res.trajectory[-1], atol=1e-6
            )

    def test_validation_rejects_invalid_laplacians(self) -> None:
        adj = sparse.csc_matrix(np.array([[0.0, 1.0], [1.0, 0.0]]))
        with self.assertRaises(ValueError):
            RouseBrownianSimulator(adj)

        shifted = sparse.csc_matrix(np.array([[2.0, -1.0], [-1.0, 2.0]]))
        with self.assertRaises(ValueError):
            RouseBrownianSimulator(shifted)


class TestCorrelator(unittest.TestCase):
    def test_correlator_workflow_and_pickle(self) -> None:
        corr = Correlator(num_correlators=3, p=4, m=2)
        self.assertEqual(corr.num_correlators, 3)
        self.assertEqual(corr.p, 4)
        self.assertEqual(corr.m, 2)
        self.assertEqual(corr.num_samples, 0)

        # Empty evaluate
        corr.evaluate(norm=True)
        self.assertEqual(corr.npcorr, 0)
        self.assertEqual(len(corr.times), 0)
        self.assertEqual(len(corr.correlation), 0)

        # Single-sample add(value, level=0) and readonly properties
        corr.add(2.0, level=0)
        corr.add(4.0, level=0)
        corr.add(6.0, level=0)
        self.assertEqual(corr.num_samples, 3)
        self.assertAlmostEqual(corr.accval, 12.0)
        corr.evaluate(norm=False)
        self.assertEqual(corr.npcorr, 3)
        self.assertGreaterEqual(len(corr.raw_t), corr.npcorr)
        self.assertGreaterEqual(len(corr.raw_f), corr.npcorr)
        np.testing.assert_allclose(corr.raw_t[: corr.npcorr], corr.times)
        np.testing.assert_allclose(corr.raw_f[: corr.npcorr], corr.correlation)

        # Constant negative large signal
        corr.reset()
        val = -5e10
        corr.add_batch(np.full(12, val))
        self.assertEqual(corr.num_samples, 12)
        self.assertAlmostEqual(corr.mean, val)
        corr.evaluate(norm=False)
        self.assertGreater(corr.npcorr, 4)
        np.testing.assert_allclose(
            corr.correlation, np.full_like(corr.correlation, val * val), rtol=1e-12
        )

        # Pickle round-trip
        restored = pickle.loads(pickle.dumps(corr))
        self.assertEqual(restored.num_samples, corr.num_samples)
        np.testing.assert_allclose(restored.times, corr.times)
        np.testing.assert_allclose(restored.correlation, corr.correlation)

        # Reset & cosine signal
        corr.set_size(num_correlators=6, p=16, m=2)
        self.assertEqual(corr.num_samples, 0)
        t = np.arange(1024, dtype=float)
        omega = 0.25
        corr.add_batch(np.cos(omega * t))
        corr.evaluate(norm=True)
        np.testing.assert_allclose(
            corr.correlation[:16], 0.5 * np.cos(omega * corr.times[:16]), atol=1.5e-2
        )

        corr.reset()
        self.assertEqual(corr.num_samples, 0)

        with self.assertRaises(ValueError):
            Correlator(0, 16, 2)
        with self.assertRaises(ValueError):
            Correlator(4, 2, 4)
        with self.assertRaises(ValueError):
            corr.add(float("nan"))
        with self.assertRaises(ValueError):
            corr.add_batch(np.array([1.0, float("nan"), 2.0]))
        self.assertEqual(corr.num_samples, 0)


if __name__ == "__main__":
    unittest.main()
