
import unittest

import numpy as np
import pandas as pd
import pandas.testing as pd_testing

from pylimer_tools_cpp.pylimer_tools_cpp import Universe


class TestEntitiesNew(unittest.TestCase):

    def test_universe(self):
        universe = Universe(10, 10, 10)
        self.assertIsInstance(universe, Universe)


if __name__ == '__main__':
    unittest.main()
