
import unittest

import pandas as pd
import pandas.testing as pd_testing
from pylimer_tools.utils.getTail import getTail
from pylimer_tools.utils.unifyDataStepsizes import unifyDataStepsizes


class TestUtilFunctions(unittest.TestCase):

    def assertDataframeEqual(self, a, b, msg):
        try:
            pd_testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

    def test_getTail(self):
        testDf = pd.DataFrame([
            {"a": 1, "b": 2, "c": 3},
            {"a": 4, "b": 5, "c": 6},
            {"a": 1, "b": 2, "c": 3},
            {"a": 1, "b": 2, "c": 3}])
        halfDf = testDf.tail(2)
        self.assertEqual(getTail(testDf, maxPercentage=1), testDf)
        self.assertEqual(getTail(testDf, percentage=0.1,
                                 maxPercentage=0.5), halfDf)
        self.assertEqual(getTail(testDf, percentage=0.5,
                                 maxPercentage=1, minN=1), halfDf)
        self.assertEqual(getTail(testDf, percentage=0.1,
                                 maxPercentage=1, minN=2), halfDf)

    def test_unifyStepSizes(self):
        testDf = pd.DataFrame([
            {"a": 1, "b": 0}, {"a": 2, "b": 0},
            {"a": 2.5, "b": 0}, {"a": 3, "b": 0}])
        expectedResult = pd.DataFrame([
            {"a": 1.0, "b": 0}, {"a": 2.0, "b": 0}, {"a": 3.0, "b": 0}]).reset_index(drop=True)

        res = unifyDataStepsizes(testDf, key="a").reset_index(drop=True)
        self.assertEqual(expectedResult, res)


if __name__ == '__main__':
    unittest.main()
