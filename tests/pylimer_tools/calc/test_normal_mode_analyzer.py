#!/usr/bin/env python
import pickle
import unittest

import numpy as np
from pylimer_tools_cpp import (
    ComputeBackend,
    HybridConfig,
    KPMConfig,
    NormalModeAnalyzer,
    NormalModeSolver,
    NormalModeSpectrum,
    SLQConfig,
    Universe,
)
from scipy import sparse


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

    def test_matrix_properties_and_factories(self):
        self.assertEqual(self.nma_ring.size, 3)
        mat = self.nma_ring.matrix.toarray()
        self.assertEqual(mat.shape, (3, 3))
        np.testing.assert_array_equal(np.diag(mat), [2.0, 2.0, 2.0])
        self.assertEqual(self.nma_ring.num_soluble_clusters, 1)
        self.assertEqual(self.nma_net.num_soluble_clusters, 2)

        # from_matrix factory
        nma_from_mat = NormalModeAnalyzer.from_matrix(self.nma_ring.matrix)
        self.assertEqual(nma_from_mat.size, 3)
        np.testing.assert_allclose(
            nma_from_mat.solve_exact().eigenvalues,
            self.nma_ring.solve_exact().eigenvalues,
        )

        # Weighted springs
        nma_weighted = NormalModeAnalyzer([0], [1], 2, weights=[2.5])
        evals_w = nma_weighted.solve_exact().eigenvalues
        np.testing.assert_allclose(evals_w, [0.0, 5.0], atol=1e-9)

        # from_universe factory
        u = Universe(10.0, 10.0, 10.0)
        u.add_atoms(
            [1, 2, 3],
            [1, 1, 1],
            [0.0, 1.0, 2.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
        )
        u.add_bonds([1, 2], [2, 3])
        nma_u = NormalModeAnalyzer.from_universe(u)
        self.assertEqual(nma_u.size, 3)
        self.assertEqual(nma_u.num_soluble_clusters, 1)
        self.assertEqual(nma_u.solve_exact().num_zero_modes, 1)

        # Non-square matrix raises
        with self.assertRaises(ValueError):
            NormalModeAnalyzer.from_matrix(sparse.csc_matrix((2, 3)))

    def test_solvers_and_eigenvalues(self):
        solvers = [
            NormalModeSolver.AUTO,
            NormalModeSolver.DENSE_LAPACK,
            NormalModeSolver.DENSE_EIGEN,
            NormalModeSolver.SPARSE_SLICING,
        ]
        for solver in solvers:
            with self.subTest(solver=solver):
                spec = self.nma_ring.solve_exact(
                    solver=solver, decompose_components=False
                )
                self.assertTrue(spec.is_complete)
                self.assertIsNone(spec.eigenvectors)
                evals = spec.eigenvalues
                self.assertEqual(len(evals), 3)
                self.assertEqual(len(spec), 3)
                self.assertAlmostEqual(evals[0], 0.0, places=5)
                self.assertAlmostEqual(evals[1], 3.0, places=3)
                self.assertAlmostEqual(evals[2], 3.0, places=3)

    def test_auto_solver_eigenvectors_and_sparse_modes(self):
        spec = self.nma_ring.solve_exact(
            solver=NormalModeSolver.AUTO, compute_eigenvectors=True
        )
        evecs = spec.eigenvectors
        self.assertIsNotNone(evecs)
        self.assertEqual(evecs.shape, (3, 3))
        np.testing.assert_allclose(evecs.T @ evecs, np.eye(3), atol=1e-8)

        # Partial spectrum via num_modes
        sparse_spec = self.nma_net.solve_exact(num_modes=4)
        self.assertFalse(sparse_spec.is_complete)
        self.assertEqual(sparse_spec.num_modes, 4)

    def test_invalid_solver_argument(self):
        with self.assertRaises(TypeError):
            self.nma_ring.solve_exact(solver="nonexistent_solver")
        with self.assertRaises(TypeError):
            self.nma_ring.solve_exact(solver=12345)
        with self.assertRaises(TypeError):
            self.nma_ring.solve_exact(solver=None)

    def test_component_decomposition(self):
        spec_decomp = self.nma_net.solve_exact(decompose_components=True)
        spec_dense = self.nma_net.solve_exact(
            solver=NormalModeSolver.DENSE_LAPACK, decompose_components=False
        )
        np.testing.assert_allclose(
            spec_decomp.eigenvalues, spec_dense.eigenvalues, atol=1e-5
        )
        self.assertEqual(spec_decomp.num_zero_modes, 2)

    def test_spectrum_direct_construction_and_observables(self):
        spec = self.nma_ring.solve_exact()
        omega = np.array([0.1, 1.0, 10.0])
        g_prime = spec.storage_modulus(omega)
        g_double_prime = spec.loss_modulus(omega)
        gp_pair, gpp_pair = spec.moduli(omega)
        np.testing.assert_allclose(gp_pair, g_prime)
        np.testing.assert_allclose(gpp_pair, g_double_prime)
        self.assertEqual(len(g_prime), 3)
        self.assertEqual(len(g_double_prime), 3)
        self.assertTrue(g_prime[0] <= g_prime[1] <= g_prime[2])

        t = np.array([0.0, 1.0, 5.0])
        gt = spec.stress_autocorrelation(t)
        self.assertEqual(len(gt), 3)
        self.assertTrue(gt[0] >= gt[1] >= gt[2])

        # Direct construction from external solver eigenvalues
        ext_spec = NormalModeSpectrum(spec.eigenvalues)
        np.testing.assert_allclose(ext_spec.storage_modulus(omega), g_prime)
        np.testing.assert_allclose(ext_spec.weights, np.ones(3))

    def test_kpm_methods(self):
        moments = self.nma_net.chebyshev_moments(
            KPMConfig(num_moments=64, num_random_vectors=10)
        )
        self.assertEqual(len(moments), 64)
        self.assertAlmostEqual(moments[0], 1.0, delta=0.2)

        grid = np.linspace(0.5, 4.0, 10)
        dos = self.nma_net.density_of_states(
            grid, KPMConfig(num_moments=64, num_random_vectors=10)
        )
        self.assertEqual(len(dos), 10)
        self.assertTrue(np.all(dos >= 0.0))

        kpm_spec = self.nma_net.solve_approx(
            KPMConfig(
                num_moments=128,
                num_random_vectors=15,
                num_quadrature_points=256,
                backend=ComputeBackend.CPU,
            )
        )
        self.assertFalse(kpm_spec.is_complete)

        omega = np.logspace(-1, 1, 5)
        g_prime_kpm, g_double_prime_kpm = kpm_spec.moduli(omega)
        self.assertEqual(len(g_prime_kpm), 5)
        self.assertEqual(len(g_double_prime_kpm), 5)
        self.assertTrue(np.all(np.isfinite(g_prime_kpm)))
        self.assertTrue(np.all(np.isfinite(g_double_prime_kpm)))

        t = np.array([0.0, 0.5, 2.0])
        gt_kpm = kpm_spec.stress_autocorrelation(t)
        self.assertEqual(len(gt_kpm), 3)
        self.assertTrue(gt_kpm[0] >= gt_kpm[1] >= gt_kpm[2])

    def test_approx_spectral_methods(self):
        n = 120
        from_nodes = list(range(n - 1)) + list(range(0, 100, 10))
        to_nodes = list(range(1, n)) + list(range(20, 120, 10))
        nma = NormalModeAnalyzer(from_nodes, to_nodes, n + 2)
        exact_spec = nma.solve_exact()

        omega = np.logspace(-4, 2, 10)
        gp_ex, gpp_ex = exact_spec.moduli(omega)
        t = np.array([0.0, 0.5, 2.0])
        gt_ex = exact_spec.stress_autocorrelation(t)

        for cfg in (
            SLQConfig(num_steps=120, num_random_vectors=30),
            HybridConfig(
                num_steps=120, num_random_vectors=30, num_exact_low_modes=30
            ),
            None,
        ):
            spec = nma.solve_approx(cfg)
            gp, gpp = spec.moduli(omega)
            gt = spec.stress_autocorrelation(t)
            self.assertEqual(len(gp), len(omega))
            self.assertEqual(len(gpp), len(omega))
            self.assertTrue(np.all(np.isfinite(gp)))
            self.assertTrue(np.all(np.isfinite(gpp)))
            self.assertTrue(gt[0] >= gt[1] >= gt[2])
            np.testing.assert_allclose(gt[0], gt_ex[0], rtol=0.15)
            np.testing.assert_allclose(gp, gp_ex, rtol=0.25)
            np.testing.assert_allclose(gpp, gpp_ex, rtol=0.25)

        kpm_spec = nma.solve_approx(
            KPMConfig(num_moments=256, num_random_vectors=20)
        )
        gp_kpm = kpm_spec.storage_modulus(omega)
        self.assertEqual(len(gp_kpm), len(omega))
        self.assertTrue(np.all(np.isfinite(gp_kpm)))

        self.assertIn("backend=ComputeBackend.AUTO", repr(SLQConfig()))
        self.assertIn("backend=ComputeBackend.AUTO", repr(HybridConfig()))
        self.assertIn("backend=ComputeBackend.AUTO", repr(KPMConfig()))
        with self.assertRaises(ValueError):
            nma.solve_approx(SLQConfig(backend=ComputeBackend.METAL))
        with self.assertRaises(ValueError):
            nma.chebyshev_moments(KPMConfig(backend=ComputeBackend.CUDA))
        with self.assertRaises(ValueError):
            nma.chebyshev_moments(KPMConfig(batch_size=0))

    def test_serialization(self):
        spec = self.nma_ring.solve_exact(compute_eigenvectors=True)
        unpickled_nma = pickle.loads(pickle.dumps(self.nma_ring))
        self.assertEqual(unpickled_nma.size, self.nma_ring.size)

        unpickled_spec = pickle.loads(pickle.dumps(spec))
        np.testing.assert_array_equal(
            spec.eigenvalues, unpickled_spec.eigenvalues
        )
        np.testing.assert_array_equal(spec.weights, unpickled_spec.weights)
        np.testing.assert_array_equal(
            spec.eigenvectors, unpickled_spec.eigenvectors
        )
        self.assertEqual(spec.is_complete, unpickled_spec.is_complete)


if __name__ == "__main__":
    unittest.main()
