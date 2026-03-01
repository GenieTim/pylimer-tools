import datetime
import os
import pathlib as pl
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
        self.assertEqual(get_tail(test_dataframe, max_percentage=1), test_dataframe)
        self.assertEqual(
            get_tail(test_dataframe, percentage=0.1, max_percentage=0.5), half_df
        )
        self.assertEqual(
            get_tail(test_dataframe, percentage=0.5, max_percentage=1, min_n=1), half_df
        )
        self.assertEqual(
            get_tail(test_dataframe, percentage=0.1, max_percentage=1, min_n=2), half_df
        )

        self.assertEqual(get_tail([0, 1, 2, 3]), [2, 3])

    def test_unify_step_sizes(self):
        test_dataframe = pd.DataFrame(
            [{"a": 1, "b": 0}, {"a": 2, "b": 0}, {"a": 2.5, "b": 0}, {"a": 3, "b": 0}]
        )
        expected_result = pd.DataFrame(
            [{"a": 1.0, "b": 0}, {"a": 2.0, "b": 0}, {"a": 3.0, "b": 0}]
        ).reset_index(drop=True)

        res = unify_data_stepsizes(test_dataframe, key="a").reset_index(drop=True)
        self.assertEqual(expected_result, res)

        with self.assertWarns(Warning):
            unify_data_stepsizes(test_dataframe, key="a", max_expected_step_size=0.5)

    def test_cache_utility(self):
        test_dataframe = pd.DataFrame(
            [{"a": 1, "b": 0}, {"a": 2, "b": 0}, {"a": 2.5, "b": 0}, {"a": 3, "b": 0}]
        )
        # make suffix unique so subsequent test runs are consistent
        suffix = "test.out" + datetime.datetime.now().strftime("%d%m%Y%H%M%S%f")
        file = os.path.dirname(__file__) + "/../fixtures/any_file.txt"
        # now, test cache capabilities
        do_cache(test_dataframe, file, suffix)
        path = pl.Path(get_cache_file_name(file, suffix))
        self.assertTrue(path.is_file())
        cached_df = load_cache(file, suffix)
        self.assertEqual(test_dataframe, cached_df)

        # and the receival of a warning for non-files
        do_cache(test_dataframe, "non-ex-file", suffix)
        path = pl.Path(get_cache_file_name("non-ex-file", suffix))
        self.assertTrue(path.is_file())
        with self.assertWarns(Warning):
            loaded_data = load_cache("non-ex-file", suffix)
            self.assertEqual(test_dataframe, loaded_data)

        # test that the cache returns empty if not written yet...
        self.assertIsNone(load_cache(file, "non-ex-suffix"))
        # ...or if the file has been modified
        test_file = open(file, "w")
        test_file.write("TEST\n")
        test_file.write(datetime.datetime.now().strftime("%c"))
        test_file.close()
        self.assertIsNone(load_cache(file, suffix))


if __name__ == "__main__":
    unittest.main()
