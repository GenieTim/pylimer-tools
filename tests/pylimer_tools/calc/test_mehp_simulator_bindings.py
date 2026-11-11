import gc
import math
import os
import sys

import numpy as np

from pylimer_tools_cpp import (
    Box,
    CrosslinkerUniverse,
    EntanglementSamplingMethod,
    MCUniverseGenerator,
    MEHPForceRelaxation,
    MEHPForceBalance2,
    SLESolver,
    StrandPairEntanglements,
    StructureSimplificationMode,
    NonGaussianSpringForceEvaluator,
    SimpleSpringMEHPForceEvaluator,
    inverse_langevin,
    randomly_sample_entanglements_analytical,
    randomly_sample_entanglements_on_demand,
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
            if "Metal acceleration is not enabled" in str(
                e
            ) or "Metal GPU backend is not available" in str(e):
                self.skipTest(f"Metal not available: {e}")
            else:
                raise

    def test_metal_bb_force_balance_relaxation(self):
        """Test that MEHPForceBalance2 with METAL_BARZILAI_BORWEIN completes force relaxation."""
        try:
            fb = MEHPForceBalance2(self.testUniverseSmall, 2, False)
            fb.run_force_relaxation(
                sle_solver=SLESolver.METAL_BARZILAI_BORWEIN)
            self.assertGreater(fb.get_nr_of_iterations(), 0)
        except RuntimeError as e:
            if "Metal acceleration is not enabled" in str(
                e
            ) or "Metal GPU backend is not available" in str(e):
                self.skipTest(f"Metal not available: {e}")
            else:
                raise

    def test_auto_sle_solver_relaxation(self):
        """Test that MEHPForceBalance2 with SLESolver.AUTO and string 'auto' completes."""
        fb1 = MEHPForceBalance2(self.testUniverseSmall, 2, False)
        fb1.run_force_relaxation(sle_solver=SLESolver.AUTO)
        self.assertGreater(fb1.get_nr_of_iterations(), 0)

        fb2 = MEHPForceBalance2(self.testUniverseSmall, 2, False)
        fb2.run_force_relaxation(sle_solver="auto")
        self.assertGreater(fb2.get_nr_of_iterations(), 0)

        fb3 = MEHPForceBalance2(self.testUniverseSmall, 2, False)
        fb3.run_force_relaxation(sle_solver="default")
        self.assertGreater(fb3.get_nr_of_iterations(), 0)

    def test_cg_fallback_solver_relaxation(self):
        """Test configurable CG fallback getters/setters and adaptive fallback behavior in Python."""
        fb = MEHPForceBalance2(self.testUniverseSmall, 2, False)
        self.assertEqual(fb.get_cg_fallback_max_iterations(), 20000)
        fb.config_cg_fallback_max_iterations(10000)
        self.assertEqual(fb.get_cg_fallback_max_iterations(), 10000)
        fb.config_cg_fallback_max_iterations(0)
        self.assertEqual(fb.get_cg_fallback_max_iterations(), 0)

    def test_analytical_entanglement_sampling_and_universe_free_force_balance(
            self):
        """Test Universe-free analytical entanglement sampling and MEHPForceBalance / MEHPForceBalance2."""
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

        es_direct = randomly_sample_entanglements_analytical(
            cu,
            box,
            nr_of_samples=40,
            upper_cutoff=2.5,
            lower_cutoff=0.2,
            minimum_nr_of_samples=10,
            same_strand_cutoff=3.0,
            seed="py_seed_42",
            filter_dangling_and_soluble=True,
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
            filter_dangling_and_soluble=True,
        )
        self.assertEqual(len(es_gen), len(es_direct))
        self.assertEqual(es_gen.strand_idx_1, es_direct.strand_idx_1)
        self.assertEqual(es_gen.strand_idx_2, es_direct.strand_idx_2)

        # MEHPForceBalance (v1) without Universe
        fb1 = gen.get_force_balance(es_gen)
        self.assertTrue(fb1.validate_network())

        # MEHPForceBalance2 (merged entanglement links) without Universe
        fb2 = gen.get_force_balance2(
            es_gen, crosslinker_type=2, entanglements_as_springs=False
        )
        fb2.validate_network()
        fb2.run_force_relaxation(
            simplification_mode=StructureSimplificationMode.INACTIVE_THEN_X2F
        )
        self.assertGreater(fb2.get_active_weight_fraction(), 0.0)
        self.assertGreater(fb2.get_soluble_weight_fraction(), 0.0)
        self.assertAlmostEqual(
            fb2.get_active_weight_fraction()
            + fb2.get_dangling_weight_fraction()
            + fb2.get_soluble_weight_fraction(),
            1.0,
            places=6,
        )
        self.assertGreater(fb2.get_gamma_factor(1.0), 0.0)

        # MEHPForceBalance2 (entanglements as springs) without Universe
        fb2_springs = gen.get_force_balance2(
            es_gen, crosslinker_type=2, entanglements_as_springs=True
        )
        fb2_springs.validate_network()
        self.assertEqual(fb2_springs.get_nr_of_extra_atoms(), 2 * len(es_gen))
        self.assertEqual(fb2_springs.get_nr_of_extra_bonds(), len(es_gen))
        fb2_springs.run_force_relaxation(
            simplification_mode=StructureSimplificationMode.INACTIVE_THEN_X2F
        )
        self.assertAlmostEqual(
            fb2_springs.get_soluble_weight_fraction(),
            fb2.get_soluble_weight_fraction(),
            places=6,
        )
        self.assertAlmostEqual(
            fb2_springs.get_active_weight_fraction(),
            fb2.get_active_weight_fraction(),
            places=6,
        )

    def test_mehp2_sampling_filtering_and_simplification_modes(self):
        """Test interaction between sampling-time topological filtering and MEHPForceBalance2 simplification."""
        gen = MCUniverseGenerator(25.0, 25.0, 25.0)
        gen.set_seed(42)
        gen.add_strands(100, [20] * 100)
        gen.add_crosslinkers(50, 4)
        gen.link_strands_to_conversion(0.85)

        # 1. Filtered sampling discards dangling/soluble candidate pairs
        # upfront
        es_filt = gen.sample_entanglements(
            30, 2.5, 0.0, 0, 3.0, "seed_f", filter_dangling_and_soluble=True
        )
        es_unfilt = gen.sample_entanglements(
            30, 2.5, 0.0, 0, 3.0, "seed_f", filter_dangling_and_soluble=False
        )
        self.assertLessEqual(len(es_filt), len(es_unfilt))

        # 2. Compare MEHPForceBalance2 behavior
        fb_filt = gen.get_force_balance2(
            es_filt, crosslinker_type=2, entanglements_as_springs=False
        )
        fb_unfilt = gen.get_force_balance2(
            es_unfilt, crosslinker_type=2, entanglements_as_springs=False
        )

        # Unfiltered network initially contains more extra atoms/links from
        # dangling chains
        self.assertGreater(
            fb_unfilt.get_nr_of_extra_atoms(), fb_filt.get_nr_of_extra_atoms()
        )

        # Running simplification prunes inactive dangling elements
        fb_filt.run_force_relaxation(
            simplification_mode=StructureSimplificationMode.INACTIVE_THEN_X2F
        )
        fb_unfilt.run_force_relaxation(
            simplification_mode=StructureSimplificationMode.INACTIVE_THEN_X2F
        )

        self.assertGreater(fb_filt.get_gamma_factor(1.0), 0.0)
        self.assertGreater(fb_unfilt.get_gamma_factor(1.0), 0.0)
        self.assertGreater(fb_filt.get_active_weight_fraction(), 0.0)
        self.assertGreater(fb_unfilt.get_active_weight_fraction(), 0.0)

    def test_on_demand_entanglement_sampling(self):
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

        # Direct topology function
        es_direct = randomly_sample_entanglements_on_demand(
            cu,
            box,
            nr_of_samples=40,
            upper_cutoff=2.5,
            lower_cutoff=0.2,
            minimum_nr_of_samples=10,
            same_strand_cutoff=3.0,
            seed="od_seed_42",
            filter_dangling_and_soluble=True,
        )
        self.assertIsInstance(es_direct, StrandPairEntanglements)
        self.assertGreaterEqual(len(es_direct), 10)
        self.assertLessEqual(len(es_direct), 40)

        # Via MCUniverseGenerator.sample_entanglements_on_demand
        es_gen_direct = gen.sample_entanglements_on_demand(
            nr_of_entanglements_to_sample=40,
            upper_cutoff=2.5,
            lower_cutoff=0.2,
            minimum_nr_of_entanglements=10,
            same_strand_cutoff=3.0,
            seed="od_seed_42",
            filter_dangling_and_soluble=True,
        )
        self.assertEqual(len(es_gen_direct), len(es_direct))
        self.assertEqual(es_gen_direct.strand_idx_1, es_direct.strand_idx_1)
        self.assertEqual(es_gen_direct.strand_idx_2, es_direct.strand_idx_2)

        # Via MCUniverseGenerator.sample_entanglements(..., method=EntanglementSamplingMethod.ON_DEMAND)
        es_gen_method = gen.sample_entanglements(
            nr_of_entanglements_to_sample=40,
            upper_cutoff=2.5,
            lower_cutoff=0.2,
            minimum_nr_of_entanglements=10,
            same_strand_cutoff=3.0,
            seed="od_seed_42",
            filter_dangling_and_soluble=True,
            method=EntanglementSamplingMethod.ON_DEMAND,
        )
        self.assertEqual(len(es_gen_method), len(es_direct))
        self.assertEqual(es_gen_method.strand_idx_1, es_direct.strand_idx_1)
        self.assertEqual(es_gen_method.strand_idx_2, es_direct.strand_idx_2)

        # Also verify string method continues to work
        es_gen_str = gen.sample_entanglements(
            nr_of_entanglements_to_sample=40,
            upper_cutoff=2.5,
            lower_cutoff=0.2,
            minimum_nr_of_entanglements=10,
            same_strand_cutoff=3.0,
            seed="od_seed_42",
            filter_dangling_and_soluble=True,
            method="on_demand",
        )
        self.assertEqual(len(es_gen_str), len(es_direct))

        # MEHPForceBalance (v1) with pre-sampled entanglements
        fb1 = gen.get_force_balance(es_gen_method)
        self.assertTrue(fb1.validate_network())

        # MEHPForceBalance2 with pre-sampled entanglements
        fb2 = gen.get_force_balance2(
            es_gen_method, crosslinker_type=2, entanglements_as_springs=False
        )
        fb2.validate_network()
        fb2.run_force_relaxation(
            simplification_mode=StructureSimplificationMode.INACTIVE_THEN_X2F
        )
        self.assertGreater(fb2.get_active_weight_fraction(), 0.0)
        self.assertGreater(fb2.get_gamma_factor(1.0), 0.0)

        # Sampling and building directly via get_force_balance2(..., method=EntanglementSamplingMethod.ON_DEMAND)
        fb2_direct = gen.get_force_balance2(
            nr_of_entanglements_to_sample=40,
            upper_cutoff=2.5,
            lower_cutoff=0.2,
            minimum_nr_of_entanglements=10,
            same_strand_cutoff=3.0,
            seed="od_seed_direct",
            filter_dangling_and_soluble=True,
            crosslinker_type=2,
            entanglements_as_springs=False,
            method=EntanglementSamplingMethod.ON_DEMAND,
        )
        fb2_direct.validate_network()
        self.assertGreater(fb2_direct.get_nr_of_extra_atoms(), 0)

