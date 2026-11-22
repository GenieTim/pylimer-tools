import gc
import math
import os
import pickle
import sys
import tempfile

import numpy as np

from pylimer_tools_cpp import (
    Box,
    CrosslinkerUniverse,
    EntanglementFilteringMode,
    EntanglementMode,
    ExitReason,
    HookeanForceEvaluator,
    LinearSystemMEHPSolver,
    LinkSwappingMode,
    MCUniverseGenerator,
    MEHPForceEvaluator,
    MEHPLinearSolverConfig,
    MEHPNetwork,
    MEHPNonlinearSolverConfig,
    MEHPSimulator,
    MEHPSlipLinkSolverConfig,
    MEHPSolveResult,
    MEHPSolver,
    NLoptAlgorithm,
    NonGaussianSpringForceEvaluator,
    NonlinearOptimizationSolver,
    SLESolver,
    SlipLinkRelaxationSolver,
    SlipMode,
    StrandPairEntanglements,
    StructureSimplificationMode,
    inverse_langevin,
    peel_topological_leaves,
    randomly_sample_entanglements,
    randomly_sample_entanglements_on_demand,
)

if __name__ == "__main__":
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../"))

from tests.pylimer_tools.universeUsingTestCase import UniverseUsingTestCase


class TestMEHPSimulatorBindings(UniverseUsingTestCase):
    def test_enums_and_simulator_instantiations(self) -> None:
        """Test unified enums, stateless evaluators, and MEHPSimulator constructors."""
        self.assertNotEqual(EntanglementMode.NONE, EntanglementMode.MERGED_LINKS)
        self.assertNotEqual(EntanglementMode.MERGED_LINKS, EntanglementMode.SPRING)
        self.assertNotEqual(SlipMode.FIXED, SlipMode.SLIDING)
        self.assertNotEqual(NLoptAlgorithm.LD_MMA, NLoptAlgorithm.LD_LBFGS)

        sim = MEHPSimulator(self.testUniverse)
        gaussian_eval = HookeanForceEvaluator(1.0)
        self.assertIsInstance(gaussian_eval, HookeanForceEvaluator)
        self.assertIsInstance(gaussian_eval, MEHPForceEvaluator)
        self.assertTrue(gaussian_eval.is_linear())
        self.assertAlmostEqual(gaussian_eval.get_kappa(), 1.0)
        sim.set_force_evaluator(gaussian_eval)
        self.assertAlmostEqual(sim.get_force_evaluator().get_kappa(), 1.0)

        non_gaussian_eval = NonGaussianSpringForceEvaluator(1.0, 0.98)
        self.assertIsInstance(non_gaussian_eval, MEHPForceEvaluator)
        self.assertFalse(non_gaussian_eval.is_linear())
        self.assertAlmostEqual(non_gaussian_eval.get_kappa(), 1.0)
        self.assertAlmostEqual(non_gaussian_eval.l, 0.98)

        net = MEHPNetwork.from_universe(self.testUniverse)
        sim_net = MEHPSimulator(net, non_gaussian_eval, False)
        self.assertEqual(sim_net.get_nr_of_iterations(), 0)
        self.assertFalse(sim_net.get_simulation_has_run())
        self.assertEqual(sim_net.get_exit_reason(), ExitReason.UNSET)

        cu = self.testUniverse.to_crosslinker_universe(2)
        sim_cu = MEHPSimulator(self.testUniverse.get_box(), cu)
        self.assertEqual(sim_cu.get_nr_of_nodes(), 4)
        self.assertEqual(sim_cu.get_nr_of_strands(), 3)
        self.assertEqual(sim.get_nr_of_nodes(), 2)

        sim_rand = MEHPSimulator.construct_with_random_entanglements(
            self.saturatedTestUniverse,
            nr_of_entanglements_to_sample=3,
            upper_sampling_cutoff=3.0,
            seed="rand_ctor_seed",
            entanglement_mode=EntanglementMode.MERGED_LINKS,
        )
        self.assertGreater(sim_rand.get_nr_of_springs(), 0)

    def test_mehp_network_builders_and_properties(self) -> None:
        """Test MEHPNetwork factory methods and exposed topology properties."""
        ents = randomly_sample_entanglements(
            self.saturatedTestUniverse,
            nr_of_samples=2,
            upper_cutoff=4.0,
            seed="net_ent_seed",
            crosslinker_type=2,
            filtering_mode=EntanglementFilteringMode.NONE,
        )
        net_none = MEHPNetwork.from_universe(
            self.saturatedTestUniverse,
            ents,
            crosslinker_type=2,
            entanglement_mode=EntanglementMode.NONE,
        )
        net_merged = MEHPNetwork.from_universe(
            self.saturatedTestUniverse,
            ents,
            crosslinker_type=2,
            entanglement_mode=EntanglementMode.MERGED_LINKS,
        )
        net_spring = MEHPNetwork.from_universe(
            self.saturatedTestUniverse,
            ents,
            crosslinker_type=2,
            entanglement_mode=EntanglementMode.SPRING,
        )
        net_rand = MEHPNetwork.from_universe_with_random_entanglements(
            self.saturatedTestUniverse,
            nr_of_entanglements_to_sample=2,
            upper_sampling_cutoff=4.0,
            seed="net_rand_seed",
            entanglement_mode=EntanglementMode.MERGED_LINKS,
        )

        for net in (net_none, net_merged, net_spring, net_rand):
            net.validate()
            self.assertEqual(len(net.L), 3)
            self.assertGreater(net.vol, 0.0)
            self.assertEqual(len(net.coordinates), 3 * net.nr_of_links)
            self.assertEqual(len(net.spring_contour_length), net.nr_of_springs)
            self.assertEqual(len(net.strand_contour_length), net.nr_of_strands)
            self.assertEqual(len(net.spring_partitions), net.nr_of_springs)
            self.assertEqual(len(net.spring_index_a), net.nr_of_springs)
            self.assertEqual(len(net.spring_index_b), net.nr_of_springs)
            self.assertEqual(len(net.spring_box_offset), 3 * net.nr_of_springs)
            self.assertEqual(len(net.strand_index_of_spring), net.nr_of_springs)
            self.assertEqual(len(net.link_is_entanglement), net.nr_of_links)
            self.assertEqual(len(net.spring_is_entanglement), net.nr_of_springs)
            self.assertEqual(len(net.old_atom_ids), net.nr_of_links)
            self.assertEqual(len(net.old_atom_types), net.nr_of_links)
            self.assertIsInstance(net.assume_box_large_enough, bool)
            self.assertIsInstance(net.assume_complete, bool)
            self.assertIsInstance(net.remove_primary_loops, bool)
            self.assertGreaterEqual(net.mean_spring_contour_length, 0.0)
            self.assertGreaterEqual(net.mean_strand_contour_length, 0.0)
            self.assertGreaterEqual(net.default_bond_length, 0.0)
            self.assertGreater(net.default_r02, 0.0)

        self.assertEqual(net_none.nr_of_links, net_none.nr_of_nodes)
        self.assertFalse(np.any(net_none.link_is_entanglement))
        self.assertFalse(np.any(net_none.spring_is_entanglement))
        if len(ents.pairs_of_atoms) > 0:
            self.assertGreater(net_merged.nr_of_links, net_merged.nr_of_nodes)
            self.assertTrue(np.any(net_merged.link_is_entanglement))
            self.assertFalse(np.any(net_merged.spring_is_entanglement))
            self.assertTrue(np.any(net_spring.spring_is_entanglement))

    def test_inverse_langevin_symmetry(self) -> None:
        """Test langevin_inv symmetry: langevin_inv(-x) == -langevin_inv(x)."""
        test_values = [0.0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.999, 1.5]
        for val in test_values:
            pos = inverse_langevin(val)
            neg = inverse_langevin(-val)
            self.assertAlmostEqual(neg, -pos, places=7)
            if val > 0:
                self.assertGreater(pos, 0.0)
                self.assertLess(neg, 0.0)

    def test_box_is_valid_offset_negative(self) -> None:
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

    def test_keep_alive_evaluator_lifetime(self) -> None:
        """Test shared_ptr lifetime keeps evaluator alive with MEHPSimulator."""

        def create_sim_from_network() -> MEHPSimulator:
            net = MEHPNetwork.from_universe(self.testUniverse)
            evaluator = HookeanForceEvaluator(2.5)
            return MEHPSimulator(net, evaluator, False)

        def create_sim_via_setter() -> MEHPSimulator:
            sim = MEHPSimulator(self.testUniverse)
            sim.set_force_evaluator(HookeanForceEvaluator(2.5))
            return sim

        for builder in (create_sim_from_network, create_sim_via_setter):
            sim = builder()
            gc.collect()
            sim.run()
            self.assertTrue(sim.get_simulation_has_run())
            self.assertGreater(sim.get_nr_of_iterations(), 0)

    def test_non_gaussian_stress_contribution_at_zero(self) -> None:
        """Test stress contribution at r=0 returns 0.0 without NaN and validates contour_length."""
        evaluator = NonGaussianSpringForceEvaluator(1.0, 1.0)
        self.assertEqual(evaluator.l, 1.0)
        stress = evaluator.evaluate_stress_contribution([0.0, 0.0, 0.0], 0, 0, 10.0)
        self.assertEqual(stress, 0.0)
        self.assertFalse(math.isnan(stress))
        with self.assertRaises(ValueError):
            evaluator.evaluate_stress_contribution([0.0, 0.0, 0.0], 0, 0, 0.0)

    def test_metal_cg_force_balance_relaxation(self) -> None:
        """Test that MEHPSimulator with METAL_CG completes force relaxation."""
        try:
            sim = MEHPSimulator(self.testUniverseSmall)
            cfg = MEHPLinearSolverConfig()
            cfg.solver = SLESolver.METAL_CG
            sim.run(cfg)
            self.assertGreater(sim.get_nr_of_iterations(), 0)
        except RuntimeError as e:
            if "Metal acceleration is not enabled" in str(
                e
            ) or "Metal GPU backend is not available" in str(e):
                self.skipTest(f"Metal not available: {e}")
            else:
                raise

    def test_metal_bb_force_balance_relaxation(self) -> None:
        """Test that MEHPSimulator with METAL_BARZILAI_BORWEIN completes force relaxation."""
        try:
            sim = MEHPSimulator(self.testUniverseSmall)
            cfg = MEHPLinearSolverConfig()
            cfg.solver = SLESolver.METAL_BARZILAI_BORWEIN
            sim.run(cfg)
            self.assertGreater(sim.get_nr_of_iterations(), 0)
        except RuntimeError as e:
            if "Metal acceleration is not enabled" in str(
                e
            ) or "Metal GPU backend is not available" in str(e):
                self.skipTest(f"Metal not available: {e}")
            else:
                raise

    def test_auto_sle_solver_relaxation(self) -> None:
        """Test that MEHPSimulator with SLESolver.AUTO completes and rejects strings."""
        sim1 = MEHPSimulator(self.testUniverseSmall)
        cfg1 = MEHPLinearSolverConfig()
        cfg1.solver = SLESolver.AUTO
        res1 = sim1.run(cfg1)
        self.assertIsInstance(res1, MEHPSolveResult)
        self.assertGreater(res1.nr_of_iterations, 0)
        self.assertGreater(sim1.get_nr_of_iterations(), 0)

        cfg2 = MEHPLinearSolverConfig()
        with self.assertRaises(TypeError):
            cfg2.solver = "auto"  # type: ignore[assignment]
        with self.assertRaises(TypeError):
            cfg2.solver = "default"  # type: ignore[assignment]

    def test_solver_configs_and_simulator_methods(self) -> None:
        """Test all three solver configs, solver instances, and simulator observables."""
        linear_solver = LinearSystemMEHPSolver()
        slip_solver = SlipLinkRelaxationSolver()
        nonlinear_solver = NonlinearOptimizationSolver()
        self.assertIsInstance(linear_solver, MEHPSolver)
        self.assertIsInstance(slip_solver, MEHPSolver)
        self.assertIsInstance(nonlinear_solver, MEHPSolver)

        linear_cfg = MEHPLinearSolverConfig()
        self.assertEqual(linear_cfg.cg_fallback_max_iterations, 20000)
        linear_cfg.cg_fallback_max_iterations = 10000
        self.assertEqual(linear_cfg.cg_fallback_max_iterations, 10000)
        linear_cfg.in_place_rebuild_threshold = 0.25
        self.assertAlmostEqual(linear_cfg.in_place_rebuild_threshold, 0.25)
        linear_cfg.residual_reduction = 1e-12
        linear_cfg.max_iterations = 10000
        linear_cfg.gradient_descent_learning_rate = 0.01
        linear_cfg.slip_mode = SlipMode.FIXED

        gen = MCUniverseGenerator(20.0, 20.0, 20.0)
        gen.set_seed(321)
        gen.set_mean_squared_bead_distance(1.0)
        gen.config_nr_of_mc_steps(0)
        gen.add_crosslinkers(40, 4, 2, True)
        gen.add_strands(75, [25] * 75, 1)
        gen.link_strands_to_conversion(0.9)
        sim = MEHPSimulator(gen.get_universe())
        sim.config_mean_bond_length(1.0)
        self.assertAlmostEqual(sim.get_default_bond_length(), 1.0)
        sim.config_mean_squared_bond_length(1.44)
        self.assertAlmostEqual(sim.get_default_r0_square(), 1.44)
        sim.config_assume_box_large_enough(False)
        sim.config_assume_network_is_complete(False)
        sim.config_remove_primary_loops(False)
        self.assertFalse(sim.get_is_2d())

        sim.deform_to(Box(22.0, 20.0, 20.0))
        self.assertAlmostEqual(sim.get_volume(), 8800.0)

        res = sim.run(LinearSystemMEHPSolver(linear_cfg))
        self.assertEqual(sim.get_last_solve_result().exit_reason, res.exit_reason)
        self.assertEqual(sim.get_exit_reason(), res.exit_reason)
        self.assertTrue(sim.get_simulation_has_run())
        self.assertGreater(sim.get_nr_of_links(), 0)
        self.assertGreater(sim.get_nr_of_nodes(), 0)
        self.assertGreater(sim.get_nr_of_strands(), 0)
        self.assertGreater(sim.get_nr_of_springs(), 0)
        self.assertGreaterEqual(sim.get_nr_of_active_nodes(), 0)
        self.assertGreaterEqual(sim.get_nr_of_active_springs(), 0)
        self.assertGreaterEqual(sim.get_nr_of_active_strands(), 0)
        self.assertGreater(sim.get_energy(), 0.0)
        self.assertEqual(len(sim.get_residuals()), 3 * sim.get_nr_of_links())
        self.assertLess(sim.get_residual_norm(), 1e-6)
        self.assertEqual(sim.get_stress_tensor().shape, (3, 3))
        self.assertAlmostEqual(
            sim.get_pressure(), float(np.trace(sim.get_stress_tensor()) / 3.0)
        )
        self.assertGreater(sim.get_gamma(), 0.0)
        self.assertGreater(sim.get_gamma_factor(), 0.0)
        self.assertEqual(len(sim.get_gamma_factors()), sim.get_nr_of_springs())
        self.assertEqual(len(sim.get_gamma_factors_in_dir(0)), sim.get_nr_of_springs())
        self.assertEqual(len(sim.get_coordinates()), 3 * sim.get_nr_of_links())
        self.assertEqual(len(sim.get_spring_vectors()), 3 * sim.get_nr_of_springs())
        self.assertEqual(len(sim.get_spring_lengths()), sim.get_nr_of_springs())
        self.assertGreater(sim.get_average_spring_length(), 0.0)
        self.assertEqual(len(sim.get_spring_contour_length()), sim.get_nr_of_springs())
        self.assertEqual(len(sim.get_spring_partitions()), sim.get_nr_of_springs())
        self.assertEqual(len(sim.get_strand_contour_length()), sim.get_nr_of_strands())
        self.assertGreater(sim.get_average_contour_length(), 0.0)
        self.assertIsInstance(sim.get_effective_functionality_of_atoms(), dict)
        self.assertGreaterEqual(sim.count_active_clustered_atoms(), 0.0)
        self.assertIsInstance(sim.get_indices_of_active_nodes(), list)
        self.assertGreaterEqual(sim.get_num_intra_chain_slip_links(), 0)
        xv = sim.get_crosslinker_verse()
        self.assertEqual(xv.get_nr_of_atoms(), sim.get_nr_of_nodes())

        disps = np.array(sim.get_current_displacements(), copy=True)
        sim.set_current_displacements(disps)
        np.testing.assert_allclose(sim.get_current_displacements(), disps)

        # Test 2D toggle zeros z-components
        sim_2d = MEHPSimulator(self.saturatedTestUniverse)
        sim_2d.set_is_2d(True)
        self.assertTrue(sim_2d.get_is_2d())
        sim_2d.run()
        self.assertEqual(sim_2d.get_stress_tensor()[2, 2], 0.0)

        # Test dynamic slip-link insertion on simulator
        sim_sl = MEHPSimulator(self.saturatedTestUniverse)
        added = sim_sl.randomly_add_slip_links(2, cutoff=4.0, seed="add_sl_seed")
        self.assertGreaterEqual(added, 0)
        sim_sl.validate_network()

    def test_solver_precondition_value_errors(self) -> None:
        """Test that violated solver preconditions raise ValueError."""
        non_gaussian_eval = NonGaussianSpringForceEvaluator(1.0, 1.0)
        sim_ng = MEHPSimulator(
            MEHPNetwork.from_universe(self.testUniverse), non_gaussian_eval
        )

        # LinearSystemMEHPSolver rejects nonlinear evaluator
        with self.assertRaises(ValueError):
            sim_ng.run(MEHPLinearSolverConfig())

        # SlipLinkRelaxationSolver rejects nonlinear evaluator
        with self.assertRaises(ValueError):
            sim_ng.run(MEHPSlipLinkSolverConfig())

        # LinearSystemMEHPSolver rejects SlipMode.SLIDING
        sim_hookean = MEHPSimulator(self.testUniverse)
        lin_sliding = MEHPLinearSolverConfig()
        lin_sliding.slip_mode = SlipMode.SLIDING
        with self.assertRaises(ValueError):
            sim_hookean.run(lin_sliding)

        # LinearSystemMEHPSolver rejects deprecated X1F_X2F_THEN_INACTIVE and invalid max_iterations
        lin_bad_simp = MEHPLinearSolverConfig()
        lin_bad_simp.simplification_mode = (
            StructureSimplificationMode.X1F_X2F_THEN_INACTIVE
        )
        with self.assertRaises(ValueError):
            sim_hookean.run(lin_bad_simp)

        lin_bad_iter = MEHPLinearSolverConfig()
        lin_bad_iter.max_iterations = 0
        with self.assertRaises(ValueError):
            sim_hookean.run(lin_bad_iter)

        # NonlinearOptimizationSolver rejects SlipMode.SLIDING
        nl_sliding = MEHPNonlinearSolverConfig()
        nl_sliding.slip_mode = SlipMode.SLIDING
        with self.assertRaises(ValueError):
            sim_hookean.run(nl_sliding)

        # SlipLinkRelaxationSolver with SlipMode.SLIDING rejects SPRING-mode networks
        gen = MCUniverseGenerator(18.0, 18.0, 18.0)
        gen.set_seed(77)
        gen.add_strands(60, [20] * 60)
        gen.add_crosslinkers(30, 4)
        gen.link_strands_to_conversion(0.9)
        es = gen.sample_entanglements(
            20,
            2.5,
            0.0,
            0,
            3.0,
            "spring_err_seed",
            EntanglementFilteringMode.OMIT_INACTIVE,
        )
        if len(es) > 0:
            sim_spring = gen.get_mehp_simulator(es, EntanglementMode.SPRING)
            sl_sliding = MEHPSlipLinkSolverConfig()
            sl_sliding.slip_mode = SlipMode.SLIDING
            with self.assertRaises(ValueError):
                sim_spring.run(sl_sliding)

    def test_cross_method_parity(self) -> None:
        """Verify Linear, Fixed SlipLink, and Nonlinear solvers agree on Hookean networks."""
        gen = MCUniverseGenerator(20.0, 20.0, 20.0)
        gen.set_seed(999)
        gen.set_mean_squared_bead_distance(1.0)
        gen.config_nr_of_mc_steps(0)
        gen.add_crosslinkers(40, 4, 2, True)
        gen.add_strands(75, [25] * 75, 1)
        gen.link_strands_to_conversion(0.9)

        es = gen.sample_entanglements(
            25,
            2.5,
            0.0,
            0,
            3.0,
            "parity_seed",
            filtering_mode=EntanglementFilteringMode.OMIT_INACTIVE,
        )
        base_net = gen.get_mehp_network(es, EntanglementMode.MERGED_LINKS)
        deformed_box = Box(22.0, 19.0, 19.5)

        sim_linear = MEHPSimulator(base_net, HookeanForceEvaluator(1.0), False)
        sim_linear.deform_to(deformed_box)
        lin_cfg = MEHPLinearSolverConfig()
        lin_cfg.residual_reduction = 1e-14
        sim_linear.run(lin_cfg)

        sim_sliplink_fixed = MEHPSimulator(base_net, HookeanForceEvaluator(1.0), False)
        sim_sliplink_fixed.deform_to(deformed_box)
        sl_fixed_cfg = MEHPSlipLinkSolverConfig()
        sl_fixed_cfg.slip_mode = SlipMode.FIXED
        sl_fixed_cfg.xtol = 1e-12
        sl_fixed_cfg.max_nr_of_steps = 50000
        sim_sliplink_fixed.run(sl_fixed_cfg)

        sim_nonlinear = MEHPSimulator(base_net, HookeanForceEvaluator(1.0), False)
        sim_nonlinear.deform_to(deformed_box)
        nl_cfg = MEHPNonlinearSolverConfig()
        nl_cfg.algorithm = NLoptAlgorithm.LD_LBFGS
        nl_cfg.xtol = 1e-12
        nl_cfg.ftol = 1e-14
        nl_cfg.max_nr_of_steps = 50000
        sim_nonlinear.run(nl_cfg)

        for sim_other in (sim_sliplink_fixed, sim_nonlinear):
            self.assertAlmostEqual(
                sim_other.get_energy(), sim_linear.get_energy(), places=4
            )
            self.assertAlmostEqual(
                sim_other.get_gamma_factor(1.0),
                sim_linear.get_gamma_factor(1.0),
                places=4,
            )
            np.testing.assert_allclose(
                sim_other.get_stress_tensor(),
                sim_linear.get_stress_tensor(),
                rtol=1e-3,
                atol=1e-5,
            )
            self.assertAlmostEqual(
                sim_other.get_active_weight_fraction(),
                sim_linear.get_active_weight_fraction(),
                places=6,
            )

        # Sliding slip-links relax further than fixed links
        sim_sliding = MEHPSimulator(base_net, HookeanForceEvaluator(1.0), False)
        sim_sliding.deform_to(deformed_box)
        sl_sliding_cfg = MEHPSlipLinkSolverConfig()
        sl_sliding_cfg.slip_mode = SlipMode.SLIDING
        sl_sliding_cfg.xtol = 1e-10
        sl_sliding_cfg.max_nr_of_steps = 10000
        sim_sliding.run(sl_sliding_cfg)
        self.assertLessEqual(
            sim_sliding.get_energy(), sim_sliplink_fixed.get_energy() + 1e-6
        )

    def test_non_gaussian_with_entanglement_links(self) -> None:
        """Test NonGaussianSpringForceEvaluator combined with EntanglementMode.MERGED_LINKS."""
        gen = MCUniverseGenerator(20.0, 20.0, 20.0)
        gen.set_seed(321)
        gen.set_mean_squared_bead_distance(1.0)
        gen.config_nr_of_mc_steps(0)
        gen.add_crosslinkers(35, 4, 2, True)
        gen.add_strands(65, [25] * 65, 1)
        gen.link_strands_to_conversion(0.9)

        es = gen.sample_entanglements(
            20,
            2.5,
            0.0,
            0,
            3.0,
            "ng_ent_seed",
            filtering_mode=EntanglementFilteringMode.OMIT_INACTIVE,
        )
        net = gen.get_mehp_network(es, EntanglementMode.MERGED_LINKS)
        sim = MEHPSimulator(net, NonGaussianSpringForceEvaluator(1.0, 1.0), False)
        sim.deform_to(Box(22.0, 19.0, 19.0))

        # Auto-dispatches to NonlinearOptimizationSolver since evaluator.is_linear() is False
        res = sim.run()
        self.assertTrue(sim.get_simulation_has_run())
        self.assertIn(
            res.exit_reason,
            (ExitReason.X_TOLERANCE, ExitReason.F_TOLERANCE, ExitReason.MAX_STEPS),
        )
        self.assertGreater(sim.get_energy(), 0.0)
        self.assertGreater(sim.get_gamma_factor(1.0), 0.0)
        stress = sim.get_stress_tensor()
        self.assertTrue(np.all(np.isfinite(stress)))
        np.testing.assert_allclose(stress, stress.T, atol=1e-10)

    def test_universe_free_force_balance(self) -> None:
        """Test Universe-free on-demand entanglement sampling and unified MEHPSimulator."""
        gen = MCUniverseGenerator(22.0, 22.0, 22.0)
        gen.set_seed(123456)
        gen.set_mean_squared_bead_distance(1.0)
        gen.config_nr_of_mc_steps(0)

        gen.add_crosslinkers(60, 4, 2, True)
        gen.add_strands(110, [25] * 110, 1)
        gen.add_monofunctional_strands(10, [15] * 10, 4)
        gen.link_strands_to_conversion(0.88)

        cu = gen.get_crosslinker_universe()
        self.assertIsInstance(cu, CrosslinkerUniverse)
        self.assertEqual(len(cu.xlink_x), 60)
        self.assertEqual(len(cu.beads_in_strand), 120)

        box = gen.get_box()
        self.assertAlmostEqual(box.get_lx(), 22.0)

        es_direct = randomly_sample_entanglements_on_demand(
            cu,
            box,
            nr_of_samples=40,
            upper_cutoff=2.5,
            lower_cutoff=0.2,
            minimum_nr_of_samples=10,
            same_strand_cutoff=3.0,
            seed="py_seed_42",
            filtering_mode=EntanglementFilteringMode.OMIT_INACTIVE,
        )
        self.assertIsInstance(es_direct, StrandPairEntanglements)
        self.assertGreaterEqual(len(es_direct), 10)
        self.assertLessEqual(len(es_direct), 40)

        es_gen = gen.sample_entanglements(
            nr_of_entanglements_to_sample=40,
            upper_cutoff=2.5,
            lower_cutoff=0.2,
            minimum_nr_of_entanglements=10,
            same_strand_cutoff=3.0,
            seed="py_seed_42",
            filtering_mode=EntanglementFilteringMode.OMIT_INACTIVE,
        )
        self.assertEqual(len(es_gen), len(es_direct))
        self.assertEqual(es_gen.strand_idx_1, es_direct.strand_idx_1)
        self.assertEqual(es_gen.strand_idx_2, es_direct.strand_idx_2)

        # Sliding slip-link simulator without Universe
        sim_slip = gen.get_mehp_simulator(
            es_gen, entanglement_mode=EntanglementMode.MERGED_LINKS
        )
        sim_slip.validate_network()
        slip_cfg = MEHPSlipLinkSolverConfig()
        slip_cfg.slip_mode = SlipMode.SLIDING
        slip_cfg.link_swapping_mode = LinkSwappingMode.NO_SWAPPING
        sim_slip.run(slip_cfg)
        self.assertTrue(sim_slip.get_simulation_has_run())

        # Fixed merged entanglement links without Universe
        sim_merged = gen.get_mehp_simulator(
            es_gen, entanglement_mode=EntanglementMode.MERGED_LINKS
        )
        sim_merged.validate_network()
        linear_cfg = MEHPLinearSolverConfig()
        linear_cfg.simplification_mode = StructureSimplificationMode.INACTIVE_THEN_X2F
        sim_merged.run(linear_cfg)
        self.assertGreater(sim_merged.get_active_weight_fraction(), 0.0)
        self.assertGreater(sim_merged.get_soluble_weight_fraction(), 0.0)
        self.assertAlmostEqual(
            sim_merged.get_active_weight_fraction()
            + sim_merged.get_dangling_weight_fraction()
            + sim_merged.get_soluble_weight_fraction(),
            1.0,
            places=6,
        )
        self.assertGreater(sim_merged.get_gamma_factor(1.0), 0.0)

        # Entanglements as compliant springs without Universe
        sim_springs = gen.get_mehp_simulator(
            es_gen, entanglement_mode=EntanglementMode.SPRING
        )
        sim_springs.validate_network()
        self.assertEqual(
            sim_springs.get_nr_of_links() - sim_springs.get_nr_of_nodes(),
            2 * len(es_gen),
        )
        self.assertEqual(
            int(np.sum(sim_springs.get_network().spring_is_entanglement)),
            len(es_gen),
        )
        sim_springs.run(linear_cfg)
        self.assertAlmostEqual(
            sim_springs.get_soluble_weight_fraction(),
            sim_merged.get_soluble_weight_fraction(),
            places=6,
        )
        self.assertAlmostEqual(
            sim_springs.get_active_weight_fraction(),
            sim_merged.get_active_weight_fraction(),
            places=6,
        )

    def test_mehp_sampling_filtering_and_simplification_modes(self) -> None:
        """Test interaction between sampling-time topological filtering and simplification."""
        gen = MCUniverseGenerator(25.0, 25.0, 25.0)
        gen.set_seed(42)
        gen.add_strands(100, [20] * 100)
        gen.add_crosslinkers(50, 4)
        gen.link_strands_to_conversion(0.85)

        es_filt = gen.sample_entanglements(
            30,
            2.5,
            0.0,
            0,
            3.0,
            "seed_f",
            filtering_mode=EntanglementFilteringMode.OMIT_INACTIVE,
        )
        es_unfilt = gen.sample_entanglements(
            30,
            2.5,
            0.0,
            0,
            3.0,
            "seed_f",
            filtering_mode=EntanglementFilteringMode.NONE,
        )
        self.assertLessEqual(len(es_filt), len(es_unfilt))

        sim_filt = gen.get_mehp_simulator(
            es_filt, entanglement_mode=EntanglementMode.MERGED_LINKS
        )
        sim_unfilt = gen.get_mehp_simulator(
            es_unfilt, entanglement_mode=EntanglementMode.MERGED_LINKS
        )

        # Unfiltered network initially contains more entanglement links from dangling chains
        self.assertGreater(
            sim_unfilt.get_nr_of_links() - sim_unfilt.get_nr_of_nodes(),
            sim_filt.get_nr_of_links() - sim_filt.get_nr_of_nodes(),
        )

        linear_cfg = MEHPLinearSolverConfig()
        linear_cfg.simplification_mode = StructureSimplificationMode.INACTIVE_THEN_X2F
        sim_filt.run(linear_cfg)
        sim_unfilt.run(linear_cfg)

        self.assertGreater(sim_filt.get_gamma_factor(1.0), 0.0)
        self.assertGreater(sim_unfilt.get_gamma_factor(1.0), 0.0)
        self.assertGreater(sim_filt.get_active_weight_fraction(), 0.0)
        self.assertGreater(sim_unfilt.get_active_weight_fraction(), 0.0)

    def test_on_demand_entanglement_sampling(self) -> None:
        """Test on-demand Brownian bridge entanglement sampling with bounding box culling."""
        gen = MCUniverseGenerator(22.0, 22.0, 22.0)
        gen.set_seed(123456)
        gen.set_mean_squared_bead_distance(1.0)
        gen.config_nr_of_mc_steps(0)

        gen.add_crosslinkers(60, 4, 2, True)
        gen.add_strands(110, [25] * 110, 1)
        gen.add_monofunctional_strands(10, [15] * 10, 4)
        gen.link_strands_to_conversion(0.88)

        cu = gen.get_crosslinker_universe()
        box = gen.get_box()

        es_direct = randomly_sample_entanglements_on_demand(
            cu,
            box,
            nr_of_samples=40,
            upper_cutoff=2.5,
            lower_cutoff=0.2,
            minimum_nr_of_samples=10,
            same_strand_cutoff=3.0,
            seed="od_seed_42",
            filtering_mode=EntanglementFilteringMode.OMIT_INACTIVE,
        )
        self.assertIsInstance(es_direct, StrandPairEntanglements)
        self.assertGreaterEqual(len(es_direct), 10)
        self.assertLessEqual(len(es_direct), 40)

        es_unified = randomly_sample_entanglements(
            cu,
            box,
            nr_of_samples=40,
            upper_cutoff=2.5,
            lower_cutoff=0.2,
            minimum_nr_of_samples=10,
            same_strand_cutoff=3.0,
            seed="od_seed_42",
            filtering_mode=EntanglementFilteringMode.OMIT_INACTIVE,
        )
        self.assertEqual(len(es_unified), len(es_direct))
        self.assertEqual(es_unified.strand_idx_1, es_direct.strand_idx_1)
        self.assertEqual(es_unified.strand_idx_2, es_direct.strand_idx_2)

        es_gen_method = gen.sample_entanglements(
            nr_of_entanglements_to_sample=40,
            upper_cutoff=2.5,
            lower_cutoff=0.2,
            minimum_nr_of_entanglements=10,
            same_strand_cutoff=3.0,
            seed="od_seed_42",
            filtering_mode=EntanglementFilteringMode.OMIT_INACTIVE,
        )
        self.assertEqual(len(es_gen_method), len(es_direct))
        self.assertEqual(es_gen_method.strand_idx_1, es_direct.strand_idx_1)
        self.assertEqual(es_gen_method.strand_idx_2, es_direct.strand_idx_2)

        sim_merged = gen.get_mehp_simulator(
            es_gen_method, entanglement_mode=EntanglementMode.MERGED_LINKS
        )
        sim_merged.validate_network()
        linear_cfg = MEHPLinearSolverConfig()
        linear_cfg.simplification_mode = StructureSimplificationMode.INACTIVE_THEN_X2F
        sim_merged.run(linear_cfg)
        self.assertGreater(sim_merged.get_active_weight_fraction(), 0.0)
        self.assertGreater(sim_merged.get_gamma_factor(1.0), 0.0)

        # Sampling and building directly via MCUniverseGenerator.get_mehp_simulator
        sim_direct = gen.get_mehp_simulator(
            40,
            2.5,
            0.2,
            10,
            3.0,
            "od_seed_direct",
        )
        sim_direct.validate_network()
        self.assertGreater(sim_direct.get_nr_of_links(), sim_direct.get_nr_of_nodes())

    def test_universe_to_crosslinker_universe_and_leaf_peeling(self) -> None:
        """Test Universe.to_crosslinker_universe, peel_topological_leaves, identify_active_atoms."""
        u = self.testUniverse
        cu = u.to_crosslinker_universe(2)
        self.assertIsInstance(cu, CrosslinkerUniverse)
        self.assertGreater(len(cu.xlink_types), 0)
        self.assertGreater(len(cu.beads_in_strand), 0)

        ignored_strands = u.peel_topological_leaves(2)
        self.assertEqual(len(ignored_strands), len(cu.beads_in_strand))

        active_atoms = u.identify_active_atoms(2)
        self.assertEqual(len(active_atoms), u.get_nr_of_atoms())

        # Test peel_topological_leaves free function
        ignored_free = peel_topological_leaves(cu)
        self.assertEqual(list(ignored_free), list(ignored_strands))

    def test_entanglement_filtering_modes(self) -> None:
        """Test EntanglementFilteringMode enum and discrete sampling with filtering modes."""
        self.assertEqual(int(EntanglementFilteringMode.NONE), 0)
        self.assertEqual(int(EntanglementFilteringMode.OMIT_INACTIVE), 1)
        self.assertEqual(int(EntanglementFilteringMode.PRUNE_UPFRONT), 2)

        u = self.testUniverse
        # Sample with EntanglementFilteringMode.PRUNE_UPFRONT
        es_prune = randomly_sample_entanglements(
            u,
            nr_of_samples=2,
            upper_cutoff=4.0,
            lower_cutoff=0.0,
            minimum_nr_of_samples=0,
            same_strand_cutoff=2.0,
            seed="py_test_seed",
            crosslinker_type=2,
            ignore_crosslinks=True,
            filtering_mode=EntanglementFilteringMode.PRUNE_UPFRONT,
        )
        self.assertIsInstance(es_prune.pairs_of_atoms, list)

        # Sample with EntanglementFilteringMode.OMIT_INACTIVE
        es_omit = randomly_sample_entanglements(
            u,
            nr_of_samples=2,
            upper_cutoff=4.0,
            lower_cutoff=0.0,
            minimum_nr_of_samples=0,
            same_strand_cutoff=2.0,
            seed="py_test_seed",
            crosslinker_type=2,
            ignore_crosslinks=True,
            filtering_mode=EntanglementFilteringMode.OMIT_INACTIVE,
        )
        self.assertIsInstance(es_omit.pairs_of_atoms, list)

        # Sample with EntanglementFilteringMode.NONE
        es_none = randomly_sample_entanglements(
            u,
            nr_of_samples=2,
            upper_cutoff=4.0,
            lower_cutoff=0.0,
            minimum_nr_of_samples=0,
            same_strand_cutoff=2.0,
            seed="py_test_seed",
            crosslinker_type=2,
            ignore_crosslinks=True,
            filtering_mode=EntanglementFilteringMode.NONE,
        )
        self.assertIsInstance(es_none.pairs_of_atoms, list)

        # Passing string to filtering_mode should raise TypeError
        with self.assertRaises(TypeError):
            randomly_sample_entanglements(
                u,
                nr_of_samples=2,
                upper_cutoff=4.0,
                filtering_mode="omit_inactive",  # type: ignore[arg-type]
            )
        with self.assertRaises(TypeError):
            randomly_sample_entanglements(
                u,
                nr_of_samples=2,
                upper_cutoff=4.0,
                filtering_mode=True,  # type: ignore[arg-type]
            )

    def test_cereal_serialization_and_pickle(self) -> None:
        """Test Cereal string/file serialization and pickle round-tripping."""
        net = MEHPNetwork.from_universe(self.saturatedTestUniverse)
        sim = MEHPSimulator(net, NonGaussianSpringForceEvaluator(1.5, 0.9), False)
        sim.deform_to(Box(11.0, 10.0, 10.0))
        sim.run()

        if hasattr(sim, "serialize_to_string"):
            serialized = sim.serialize_to_string()
            self.assertIsInstance(serialized, bytes)
            self.assertGreater(len(serialized), 0)
            restored_str = MEHPSimulator.construct_from_string(serialized)
            self.assertEqual(restored_str.get_nr_of_links(), sim.get_nr_of_links())
            self.assertEqual(restored_str.get_nr_of_springs(), sim.get_nr_of_springs())
            self.assertEqual(
                restored_str.get_nr_of_iterations(), sim.get_nr_of_iterations()
            )
            self.assertTrue(restored_str.get_simulation_has_run())
            self.assertAlmostEqual(
                restored_str.get_energy(), sim.get_energy(), places=10
            )
            self.assertAlmostEqual(
                restored_str.get_gamma_factor(), sim.get_gamma_factor(), places=10
            )
            self.assertAlmostEqual(
                restored_str.get_force_evaluator().get_kappa(), 1.5, places=10
            )
            np.testing.assert_allclose(
                restored_str.get_current_displacements(),
                sim.get_current_displacements(),
            )

        if hasattr(sim, "write_restart_file"):
            with tempfile.TemporaryDirectory() as tmpdir:
                restart_path = os.path.join(tmpdir, "mehp_restart.cereal")
                sim.write_restart_file(restart_path)
                restored_file = MEHPSimulator.read_restart_file(restart_path)
                self.assertEqual(restored_file.get_nr_of_links(), sim.get_nr_of_links())
                self.assertAlmostEqual(
                    restored_file.get_energy(), sim.get_energy(), places=10
                )

        if hasattr(MEHPSimulator, "__getstate__"):
            restored_pkl = pickle.loads(pickle.dumps(sim))
            self.assertEqual(restored_pkl.get_nr_of_links(), sim.get_nr_of_links())
            self.assertEqual(restored_pkl.get_nr_of_springs(), sim.get_nr_of_springs())
            self.assertAlmostEqual(
                restored_pkl.get_energy(), sim.get_energy(), places=10
            )
            self.assertAlmostEqual(
                restored_pkl.get_gamma_factor(), sim.get_gamma_factor(), places=10
            )
            np.testing.assert_allclose(
                restored_pkl.get_current_displacements(),
                sim.get_current_displacements(),
            )
