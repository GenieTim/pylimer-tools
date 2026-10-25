#!/usr/bin/env python
import pickle
import unittest

import numpy as np
from pylimer_tools_cpp import NormalModeAnalyzer, NormalModeSolver


class TestNormalModeAnalyzer(unittest.TestCase):
    def setUp(self):
        # 3-bead ring
        self.ring_from = [0, 1, 2]
        self.ring_to = [1, 2, 0]
        self.nma_ring = NormalModeAnalyzer(self.ring_from, self.ring_to)

        # 10-bead network with chain, ring, isolated beads
        self.net_from = [0, 1, 2, 3, 5, 6, 7]
        self.net_to = [1, 2, 3, 4, 6, 7, 5]
        self.nma_net = NormalModeAnalyzer(self.net_from, self.net_to)

    def test_matrix_properties(self):
        self.assertEqual(self.nma_ring.get_matrix_size(), 3)
        mat = self.nma_ring.get_matrix().toarray()
        self.assertEqual(mat.shape, (3, 3))
        np.testing.assert_array_equal(np.diag(mat), [2.0, 2.0, 2.0])

    def test_solvers_and_eigenvalues(self):
        solvers = [
            NormalModeSolver.AUTO,
            NormalModeSolver.DENSE_LAPACK,
            NormalModeSolver.DENSE_EIGEN,
            NormalModeSolver.SPARSE_SLICING,
            "auto",
            "dense_lapack",
            "dense_eigen",
            "sparse_slicing",
        ]
        for solver in solvers:
            with self.subTest(solver=solver):
                nma = NormalModeAnalyzer(self.ring_from, self.ring_to)
                nma.find_all_eigenvalues(
                    solver=solver, decompose_components=False)
                evals = nma.get_eigenvalues()
                self.assertEqual(len(evals), 3)
                self.assertAlmostEqual(evals[0], 0.0, places=5)
                self.assertAlmostEqual(evals[1], 3.0, places=3)
                self.assertAlmostEqual(evals[2], 3.0, places=3)

    def test_legacy_boolean_signature(self):
        # Old signature was find_all_eigenvalues(compute_eigenvectors: bool = False)
        # Passing True as first positional argument should compute eigenvectors
        nma = NormalModeAnalyzer(self.ring_from, self.ring_to)
        nma.find_all_eigenvalues(True)
        evecs = nma.get_eigenvectors()
        self.assertEqual(evecs.shape, (3, 3))
        # Orthonormality check
        np.testing.assert_allclose(evecs.T @ evecs, np.eye(3), atol=1e-8)

    def test_invalid_solver_argument(self):
        nma = NormalModeAnalyzer(self.ring_from, self.ring_to)
        with self.assertRaises(ValueError):
            nma.find_all_eigenvalues(solver="nonexistent_solver")
        with self.assertRaises(TypeError):
            nma.find_all_eigenvalues(solver=12345)
        with self.assertRaises(TypeError):
            nma.find_all_eigenvalues(solver=None)

    def test_component_decomposition(self):
        nma_decomp = NormalModeAnalyzer(self.net_from, self.net_to)
        nma_decomp.find_all_eigenvalues(decompose_components=True)
        evals_decomp = nma_decomp.get_eigenvalues()

        nma_dense = NormalModeAnalyzer(self.net_from, self.net_to)
        nma_dense.find_all_eigenvalues(
            solver=NormalModeSolver.DENSE_LAPACK, decompose_components=False
        )
        evals_dense = nma_dense.get_eigenvalues()

        np.testing.assert_allclose(evals_decomp, evals_dense, atol=1e-5)
        self.assertEqual(nma_decomp.get_nr_of_soluble_clusters(), 2)

    def test_exact_moduli_and_stress_autocorrelation(self):
        self.nma_ring.find_all_eigenvalues()
        omega = np.array([0.1, 1.0, 10.0])
        g_prime = self.nma_ring.evaluate_storage_modulus(omega)
        g_double_prime = self.nma_ring.evaluate_loss_modulus(omega)
        self.assertEqual(len(g_prime), 3)
        self.assertEqual(len(g_double_prime), 3)
        self.assertTrue(g_prime[0] <= g_prime[1] <= g_prime[2])

        t = np.array([0.0, 1.0, 5.0])
        gt = self.nma_ring.evaluate_stress_autocorrelation(t)
        self.assertEqual(len(gt), 3)
        self.assertTrue(gt[0] >= gt[1] >= gt[2])

    def test_kpm_methods(self):
        nma = NormalModeAnalyzer(self.net_from, self.net_to)

        # Chebyshev moments
        moments = nma.compute_chebyshev_moments_kpm(
            num_moments=64, num_random_vectors=10
        )
        self.assertEqual(len(moments), 64)
        self.assertAlmostEqual(moments[0], 1.0, delta=0.2)

        # Density of states
        grid = np.linspace(0.5, 4.0, 10)
        dos = nma.compute_density_of_states_kpm(
            grid, num_moments=64, num_random_vectors=10
        )
        self.assertEqual(len(dos), 10)
        self.assertTrue(np.all(dos >= 0.0))

        # Moduli via KPM
        omega = np.logspace(-1, 1, 5)
        g_prime_kpm = nma.evaluate_storage_modulus_kpm(
            omega, num_moments=128, num_random_vectors=15, num_quadrature_points=256
        )
        g_double_prime_kpm = nma.evaluate_loss_modulus_kpm(
            omega, num_moments=128, num_random_vectors=15, num_quadrature_points=256
        )
        self.assertEqual(len(g_prime_kpm), 5)
        self.assertEqual(len(g_double_prime_kpm), 5)
        self.assertTrue(np.all(np.isfinite(g_prime_kpm)))
        self.assertTrue(np.all(np.isfinite(g_double_prime_kpm)))

        # Stress autocorrelation via KPM
        t = np.array([0.0, 0.5, 2.0])
        gt_kpm = nma.evaluate_stress_autocorrelation_kpm(
            t, num_moments=128, num_random_vectors=15, num_quadrature_points=256
        )
        self.assertEqual(len(gt_kpm), 3)
        self.assertTrue(gt_kpm[0] >= gt_kpm[1] >= gt_kpm[2])

    def test_serialization(self):
        self.nma_ring.find_all_eigenvalues()
        evals_before = self.nma_ring.get_eigenvalues()
        pickled = pickle.dumps(self.nma_ring)
        unpickled = pickle.loads(pickled)
        evals_after = unpickled.get_eigenvalues()
        np.testing.assert_array_equal(evals_before, evals_after)


if __name__ == "__main__":
    unittest.main()
