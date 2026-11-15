import datetime
import os
import pathlib as pl
import tempfile
import time
import unittest

import numpy as np
import pandas as pd
import pandas.testing as pd_testing

from pylimer_tools.utils.cache_utility import do_cache, get_cache_file_name, load_cache
from pylimer_tools.utils.data_utility import get_tail, unify_data_stepsizes
from tests.pylimer_tools.pdComparingTestCase import PandasComparingTestCase


class TestUtilFunctions(PandasComparingTestCase):
    def test_get_tail(self):
        test_dataframe = pd.DataFrame(
            [
                {"a": 1, "b": 2, "c": 3},
                {"a": 4, "b": 5, "c": 6},
                {"a": 1, "b": 2, "c": 3},
                {"a": 1, "b": 2, "c": 3},
            ]
        )
        half_df = test_dataframe.tail(2)
        self.assertEqual(get_tail(test_dataframe), half_df)
        self.assertEqual(
            get_tail(
                test_dataframe,
                max_percentage=1),
            test_dataframe)
        self.assertEqual(
            get_tail(test_dataframe, percentage=0.1,
                     max_percentage=0.5), half_df
        )
        self.assertEqual(
            get_tail(test_dataframe, percentage=0.5,
                     max_percentage=1, min_n=1), half_df
        )
        self.assertEqual(
            get_tail(test_dataframe, percentage=0.1,
                     max_percentage=1, min_n=2), half_df
        )

        self.assertEqual(get_tail([0, 1, 2, 3]), [2, 3])

    def test_unify_step_sizes(self):
        test_dataframe = pd.DataFrame(
            [{"a": 1, "b": 0}, {"a": 2, "b": 0}, {
                "a": 2.5, "b": 0}, {"a": 3, "b": 0}]
        )
        expected_result = pd.DataFrame(
            [{"a": 1.0, "b": 0}, {"a": 2.0, "b": 0}, {"a": 3.0, "b": 0}]
        ).reset_index(drop=True)

        res = unify_data_stepsizes(
            test_dataframe,
            key="a").reset_index(
            drop=True)
        self.assertEqual(expected_result, res)

        with self.assertWarns(Warning):
            unify_data_stepsizes(
                test_dataframe,
                key="a",
                max_expected_step_size=0.5)

    def test_cache_utility(self):
        test_dataframe = pd.DataFrame(
            [{"a": 1, "b": 0}, {"a": 2, "b": 0}, {
                "a": 2.5, "b": 0}, {"a": 3, "b": 0}]
        )
        suffix = "test.out"
        with tempfile.TemporaryDirectory() as tmp_dir:
            file = os.path.join(tmp_dir, "any_file.txt")
            with open(file, "w") as test_file:
                test_file.write("TEST\n")

            past_time = time.time() - 10.0
            os.utime(file, (past_time, past_time))

            # now, test cache capabilities
            do_cache(test_dataframe, file, suffix, tmp_dir=tmp_dir)
            path = pl.Path(get_cache_file_name(file, suffix, tmp_dir=tmp_dir))
            self.assertTrue(path.is_file())
            cached_df = load_cache(file, suffix, tmp_dir=tmp_dir)
            self.assertEqual(test_dataframe, cached_df)

            # and the receival of a warning for non-files
            do_cache(test_dataframe, "non-ex-file", suffix, tmp_dir=tmp_dir)
            path = pl.Path(
                get_cache_file_name(
                    "non-ex-file",
                    suffix,
                    tmp_dir=tmp_dir))
            self.assertTrue(path.is_file())
            with self.assertWarns(Warning):
                loaded_data = load_cache(
                    "non-ex-file", suffix, tmp_dir=tmp_dir)
                self.assertEqual(test_dataframe, loaded_data)

            # test that the cache returns empty if not written yet...
            self.assertIsNone(
                load_cache(
                    file,
                    "non-ex-suffix",
                    tmp_dir=tmp_dir))
            # ...or if the file has been modified
            future_time = time.time() + 10.0
            os.utime(file, (future_time, future_time))
            self.assertIsNone(load_cache(file, suffix, tmp_dir=tmp_dir))

    def test_mc_universe_generator_streaming_stats(self):
        from pylimer_tools_cpp import (
            MEHPForceRelaxation,
            MCUniverseGenerator,
            StrandCategory,
            sample_strand_beads_on_demand,
        )

        gen = MCUniverseGenerator(20.0, 20.0, 20.0)
        gen.set_seed(12345)
        gen.config_nr_of_mc_steps(3)
        gen.config_primary_loop_probability(0.0)
        gen.add_crosslinkers(25, 4, 2)
        gen.add_strands(50, [8] * 50, 1)
        gen.link_strands_to_conversion(0.9)
        gen.add_solvent_chains(4, 6, 1)

        u = gen.get_universe()
        cu = gen.get_crosslinker_universe()
        box = gen.get_box()

        bond_stats = gen.compute_bond_length_stats()
        cu_bond_stats = cu.compute_bond_length_stats(box, 3)
        u_bonds = u.compute_bond_lengths()
        self.assertEqual(bond_stats.nr_of_bonds, len(u_bonds))
        self.assertEqual(cu_bond_stats.nr_of_bonds, len(u_bonds))
        self.assertAlmostEqual(
            bond_stats.mean_bond_length,
            sum(u_bonds) / len(u_bonds),
            places=10,
        )
        self.assertAlmostEqual(
            bond_stats.mean_squared_bond_length,
            sum(b * b for b in u_bonds) / len(u_bonds),
            places=10,
        )

        self.assertEqual(gen.count_atom_types(), u.count_atom_types())
        self.assertEqual(cu.count_atom_types(), u.count_atom_types())
        for deg, cnt in gen.count_atoms_by_degree().items():
            self.assertEqual(cnt, len(u.get_atoms_by_degree(deg)))

        eff_func = gen.determine_effective_functionality_per_type()
        u_eff_func = u.determine_effective_functionality_per_type()
        self.assertEqual(set(eff_func.keys()), set(u_eff_func.keys()))
        for k, v in u_eff_func.items():
            self.assertAlmostEqual(eff_func[k], v, places=10)

        ee_net = gen.compute_end_to_end_stats(StrandCategory.NETWORK_STRAND)
        self.assertAlmostEqual(
            ee_net.mean_squared_end_to_end_distance,
            u.compute_mean_squared_end_to_end_distance(2, True, True),
            places=10,
        )
        ee_all = gen.compute_end_to_end_stats(StrandCategory.ALL)
        self.assertAlmostEqual(
            ee_all.mean_squared_end_to_end_distance,
            u.compute_mean_squared_end_to_end_distance(2, False, True),
            places=10,
        )

        bond_counts = gen.compute_strand_bond_counts(StrandCategory.ALL)
        self.assertEqual(sum(bond_counts), len(u_bonds))

        coords0 = sample_strand_beads_on_demand(
            box, cu, 0, cu.strand_seeds[0], 3)
        self.assertEqual(len(coords0), 3 * cu.beads_in_strand[0])

        fr_cu = MEHPForceRelaxation(box, cu)
        fr_gen = gen.get_force_relaxation()
        self.assertEqual(fr_cu.get_nr_of_springs(), fr_gen.get_nr_of_springs())


if __name__ == "__main__":
    unittest.main()
