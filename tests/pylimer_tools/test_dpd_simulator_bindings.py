"""Tests for DPDSimulator Python bindings."""

from __future__ import annotations

import os
import tempfile
import unittest

import numpy as np

from pylimer_tools_cpp import DPDSimulator, Universe


class TestDPDSimulatorBindings(unittest.TestCase):
    """Verify DPDSimulator pybind11 methods, PBC flags, slip-spring MC, and restart I/O."""

    def test_dpd_bindings_pbc_mc_and_restart(self) -> None:
        u = Universe(10.0, 10.0, 10.0)
        u.add_atoms(
            [1, 2, 3, 4],
            [1, 1, 1, 1],
            [9.6, 0.4, 9.6, 0.4],
            [1.0, 1.0, 2.0, 2.0],
            [1.0, 1.0, 1.0, 1.0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
        )
        u.add_bonds([1, 3], [2, 4])

        sim = DPDSimulator(
            u,
            crosslinker_type=2,
            slipspring_bond_type=7,
            seed="py_dpd_seed",
            assume_box_large_enough=False,
        )
        self.assertFalse(sim.assumes_box_large_enough())
        sim.config_assume_box_large_enough(True)
        self.assertTrue(sim.assumes_box_large_enough())
        sim.assume_box_large_enough(False)
        self.assertFalse(sim.assumes_box_large_enough())

        sim.config_sigma(3.0)
        sim.config_gamma(2.25)
        self.assertAlmostEqual(sim.get_sigma(), 3.0)
        self.assertAlmostEqual(sim.get_dpd_gamma(), 2.25)
        self.assertAlmostEqual(sim.get_target_kbt(), 2.0)

        sim.config_slipspring_low_cutoff(0.5)
        sim.config_slipspring_high_cutoff(2.5)
        created = sim.create_slip_springs(2)
        self.assertEqual(created, 2)
        self.assertEqual(sim.get_nr_of_slip_springs(), 2)

        u_out = sim.get_universe(with_slipsprings=True)
        edges = u_out.get_edges()
        self.assertEqual(edges["edge_type"][-1], 7)

        shifts = sim.shift_slip_springs(1.5)
        relocs = sim.relocate_slip_springs(1.5)
        self.assertGreaterEqual(shifts, 0)
        self.assertGreaterEqual(relocs, 0)
        sim.validate_state()

        sim.config_num_steps_dpd(2)
        sim.config_num_steps_mc(2)
        sim.run_simulation(n_steps=4, dt=0.05, with_MC=True)
        self.assertEqual(len(sim.get_velocities()), 12)
        self.assertGreaterEqual(sim.get_num_shifts(), 0)
        self.assertGreaterEqual(sim.get_num_relocations(), 0)
        sim.validate_state()

        if hasattr(sim, "write_restart_file"):
            with tempfile.TemporaryDirectory() as tmpdir:
                restart_path = os.path.join(tmpdir, "dpd_restart.bin")
                sim.write_restart_file(restart_path)
                self.assertTrue(os.path.exists(restart_path))
                sim_restored = DPDSimulator.read_restart_file(restart_path)
                np.testing.assert_allclose(
                    sim_restored.get_coordinates(), sim.get_coordinates()
                )
                sim_restored.validate_state()


if __name__ == "__main__":
    unittest.main()
