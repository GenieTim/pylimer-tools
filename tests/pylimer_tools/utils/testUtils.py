
import datetime
import os
import pathlib as pl
import unittest

import numpy as np
import pandas as pd
import pandas.testing as pd_testing
from tests.pylimer_tools.pdComparingTestCase import PandasComparingTestCase
from pylimer_tools.utils.cacheUtility import (doCache, getCacheFileName,
                                              loadCache)
from pylimer_tools.utils.getTail import getTail
from pylimer_tools.utils.unifyDataStepsizes import unifyDataStepsizes


class TestUtilFunctions(PandasComparingTestCase):

    def test_getTail(self):
        testDf = pd.DataFrame([
            {"a": 1, "b": 2, "c": 3},
            {"a": 4, "b": 5, "c": 6},
            {"a": 1, "b": 2, "c": 3},
            {"a": 1, "b": 2, "c": 3}])
        halfDf = testDf.tail(2)
        self.assertEqual(getTail(testDf), halfDf)
        self.assertEqual(getTail(testDf, maxPercentage=1), testDf)
        self.assertEqual(getTail(testDf, percentage=0.1,
                                 maxPercentage=0.5), halfDf)
        self.assertEqual(getTail(testDf, percentage=0.5,
                                 maxPercentage=1, minN=1), halfDf)
        self.assertEqual(getTail(testDf, percentage=0.1,
                                 maxPercentage=1, minN=2), halfDf)

        self.assertEqual(getTail([0, 1, 2, 3]), [2, 3])

    def test_unifyStepSizes(self):
        testDf = pd.DataFrame([
            {"a": 1, "b": 0}, {"a": 2, "b": 0},
            {"a": 2.5, "b": 0}, {"a": 3, "b": 0}])
        expectedResult = pd.DataFrame([
            {"a": 1.0, "b": 0}, {"a": 2.0, "b": 0}, {"a": 3.0, "b": 0}]).reset_index(drop=True)

        res = unifyDataStepsizes(testDf, key="a").reset_index(drop=True)
        self.assertEqual(expectedResult, res)

        with self.assertWarns(Warning):
            unifyDataStepsizes(testDf, key="a", maxExpectedStepSize=0.5)

    def test_cacheUtility(self):
        testDf = pd.DataFrame([
            {"a": 1, "b": 0}, {"a": 2, "b": 0},
            {"a": 2.5, "b": 0}, {"a": 3, "b": 0}])
        # make suffix unique so subsequent test runs are consistent
        suffix = "test.out" + datetime.datetime.now().strftime("%d%m%Y%H%M%S%f")
        file = os.path.dirname(__file__) + "/../fixtures/any_file.txt"
        # now, test cache capabilities
        doCache(testDf, file, suffix)
        path = pl.Path(getCacheFileName(file, suffix))
        self.assertTrue(path.is_file())
        cachedDf = loadCache(file, suffix)
        self.assertEqual(testDf, cachedDf)

        # and the receival of a warning for non-files
        doCache(testDf, "non-ex-file", suffix)
        path = pl.Path(getCacheFileName("non-ex-file", suffix))
        self.assertTrue(path.is_file())
        with self.assertWarns(Warning):
            loadedData = loadCache("non-ex-file", suffix)
            self.assertEqual(testDf, loadedData)

        # test that the cache returns empty if not written yet...
        self.assertIsNone(loadCache(file, "non-ex-suffix"))
        # ...or if the file has been modified
        testFile = open(file, "w")
        testFile.write("TEST\n")
        testFile.write(datetime.datetime.now().strftime("%c"))
        testFile.close()
        self.assertIsNone(loadCache(file, suffix))


if __name__ == '__main__':
    unittest.main()
