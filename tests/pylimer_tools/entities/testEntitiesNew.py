
import unittest

import igraph
import numpy as np
import pandas as pd
import pandas.testing as pd_testing
from pylimer_tools import Universe

class TestEntitiesNew(unittest.TestCase):

    def assertSeriesEqual(self, a, b, msg):
        try:
            pd_testing.assert_series_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.Series, self.assertSeriesEqual)
        super().setUp()

    def test_universe(self):
        universe = Universe(10, 10, 10)
        self.assertIsInstance(universe, Universe)
