
import datetime
import os
import pathlib as pl
import unittest

import numpy as np
import pandas as pd
import pandas.testing as pd_testing
from pylimer_tools.utils.cacheUtility import (doCache, getCacheFileName,
                                              loadCache)
from pylimer_tools.utils.getMolecules import getMolecules
from pylimer_tools.utils.getTail import getTail
from pylimer_tools.utils.unifyDataStepsizes import unifyDataStepsizes


class TestUtilFunctions(unittest.TestCase):

    def assertDataframeEqual(self, a, b, msg):
        try:
            pd_testing.assert_frame_equal(a, b, check_index_type=False)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def assertSeriesEqual(self, a, b, msg):
        try:
            pd_testing.assert_series_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)
        self.addTypeEqualityFunc(pd.Series, self.assertSeriesEqual)

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

    def test_getMolecules(self):
        testAtoms = pd.DataFrame([
            {"id": 1, "nx": 1, "ny": 1, "nz": 1,
                "type": 1, "x": 1, "y": 1, "z": 1},
            {"id": 2, "nx": 1, "ny": 1, "nz": 1,
                "type": 1, "x": 2, "y": 1, "z": 1},
            {"id": 3, "nx": 1, "ny": 1, "nz": 1,
                "type": 1, "x": 3, "y": 1, "z": 1},
            {"id": 5, "nx": 1, "ny": 1, "nz": 1,
                "type": 1, "x": 1, "y": 3, "z": 1},
            {"id": 6, "nx": 1, "ny": 1, "nz": 1,
                "type": 2, "x": 1, "y": 1, "z": 2},
            {"id": 7, "nx": 1, "ny": 1, "nz": 1,
                "type": 2, "x": 1, "y": 1, "z": 3},
        ])
        testBonds = pd.DataFrame([
            {"to": 1, "bondFrom": 2},
            {"to": 3, "bondFrom": 2},
            {"to": 5, "bondFrom": 6},
            {"to": 6, "bondFrom": 7}
        ])

        expectedResult = [[
            {'atom': testAtoms.loc[0, :],
             'bonds_from': pd.DataFrame(columns=['to', 'bondFrom'], dtype=np.int64),
             'bonds_to': pd.DataFrame([{'to': 1, 'bondFrom': 2}])},
            {'atom': testAtoms.loc[1, :], 'bonds_from': pd.DataFrame([
                {'to': 1, 'bondFrom': 2},
                {'to': 3, 'bondFrom': 2}
            ]), 'bonds_to': pd.DataFrame(columns=['to', 'bondFrom'], dtype=np.int64)},
            {'atom': testAtoms.loc[2, :],
             'bonds_from': pd.DataFrame(columns=['to', 'bondFrom'], dtype=np.int64),
             'bonds_to': pd.DataFrame([{"to": 3, "bondFrom": 2}])}
        ], [
            {'atom': testAtoms.loc[3, :],
             'bonds_from': pd.DataFrame(columns=['to', 'bondFrom'], dtype=np.int64),
             'bonds_to': pd.DataFrame(columns=['to', 'bondFrom'], dtype=np.int64)}
        ]]

        resultingMolecules = getMolecules(
            testAtoms, testBonds, [10, 10, 10], 2)

        self.assertEqual(len(resultingMolecules), 2)
        self.assertEqual(len(resultingMolecules[0]), 3)
        self.assertEqual(len(resultingMolecules[1]), 1)
        for moleculeIdx in range(len(expectedResult)):
            self.assertEqual(len(expectedResult[moleculeIdx]), len(
                resultingMolecules[moleculeIdx]))
            for atomIdx in range(len(expectedResult[moleculeIdx])):
                for atomKey in expectedResult[moleculeIdx][atomIdx]:
                    self.assertTrue(
                        atomKey in resultingMolecules[moleculeIdx][atomIdx])
                    atomValueExpected = expectedResult[moleculeIdx][atomIdx][atomKey]
                    atomValueActual = resultingMolecules[moleculeIdx][atomIdx][atomKey]
                    if (isinstance(atomValueExpected, pd.DataFrame)):
                        self.assertEqual(atomValueExpected.reset_index(
                            drop=True), atomValueActual.reset_index(drop=True))
                    else:
                        self.assertEqual(atomValueExpected, atomValueActual)


if __name__ == '__main__':
    unittest.main()
