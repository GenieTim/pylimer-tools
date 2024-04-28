
import unittest

import pandas as pd
import pandas.testing as pd_testing


class PandasComparingTestCase(unittest.TestCase):
    def assertDataframeEqual(self, a, b, msg="Dataframes are not equal", **kwargs):
        try:
            pd_testing.assert_frame_equal(
                a, b, check_index_type=False, **kwargs)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def assertSeriesEqual(self, a, b, msg="Series are not equal", **kwargs):
        try:
            pd_testing.assert_series_equal(a, b, **kwargs)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)
        self.addTypeEqualityFunc(pd.Series, self.assertSeriesEqual)
