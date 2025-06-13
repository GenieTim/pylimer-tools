import numpy as np
import pandas as pd
from pylimer_tools.utils.optimize_dataframe import *
from tests.pylimer_tools.pdComparingTestCase import PandasComparingTestCase


class TestOptimizeDf(PandasComparingTestCase):
    def testSeparateOptimizatios(self):
        class Object(object):
            pass

        testObject = Object()
        df = pd.DataFrame(
            [
                {
                    "testFloat": np.float64(1.01),
                    "testInt": np.int64(1e1),
                    "testObject": testObject,
                    "testUint64": np.uint64(np.iinfo(np.uint64).max - 1),
                    "testString": "strong",
                }
            ]
        )
        optimizedDf = optimize(df)

        expectedOptimizedDf = pd.DataFrame(
            [
                {
                    "testFloat": np.float32(1.01),
                    "testInt": np.int8(1e1),
                    "testObject": testObject,
                    "testUint64": np.uint64(np.iinfo(np.uint64).max - 1),
                    "testString": "strong",
                }
            ]
        )

        # contents are equal
        self.assertEqual(df, optimizedDf)
        self.assertEqual(expectedOptimizedDf, optimizedDf)

    def testReduceMemUseage(self):
        class Object(object):
            pass

        testObject = Object()
        dfExtra = pd.DataFrame(
            [
                {
                    "testFloat": np.float64(1.01),
                    "testInt": np.int64(1e1),
                    "testObject": testObject,
                    "testUint64": np.uint64(np.iinfo(np.uint64).max - 1),
                    "testFloat64": np.float64(np.finfo(np.float32).max + 10**20),
                    "testString": "strong",
                }
            ]
        )

        reducedMemDf = reduce_mem_usage(dfExtra)
        expectedReducedMemDf = pd.DataFrame(
            [
                {
                    "testFloat": np.float16(1.01),
                    "testInt": np.int8(1e1),
                    "testObject": testObject,
                    "testUint64": np.uint64(np.iinfo(np.uint64).max - 1),
                    "testFloat64": np.float64(np.finfo(np.float32).max + 10**20),
                    "testString": "strong",
                }
            ]
        )
        self.assertEqual(expectedReducedMemDf, dfExtra)
        self.assertEqual(expectedReducedMemDf, reducedMemDf)

    def testCategoriesAreInitiated(self):
        class Object(object):
            def __init__(self, value):
                self.value = value

            def __eq__(self, other):
                if not isinstance(other, Object):
                    # don't attempt to compare against unrelated types
                    raise NotImplementedError
                return self.value == other.value

            def __hash__(self):
                return hash(self.value)

        testObject1 = Object(0)
        testObject2 = Object(2)

        df = pd.DataFrame(
            {
                "categoryThingy": [
                    testObject1,
                    testObject1,
                    testObject1,
                    testObject1,
                    testObject2,
                    testObject2,
                    testObject1,
                    testObject1,
                ],
                "value": [1.1, 1.2, 3.2, 1.3, 7.9, 8.9, 9.9, 10.0],
            }
        )
        optimizedDf = optimize(df)
        self.assertEqual(optimizedDf["categoryThingy"].dtype, "category")
        reducedDf = reduce_mem_usage(df, obj_to_category=True)
        self.assertEqual(reducedDf["categoryThingy"].dtype, "category")

    def testDatetimeFeatures(self):
        df = pd.DataFrame({"dateThingy": ["2022-08-12", "2021-07-13"]})
        optimizedDf = optimize(df, datetime_features=["dateThingy"])
        self.assertTrue(
            len(
                optimizedDf.loc[
                    :, [np.issubdtype(t, np.datetime64) for t in optimizedDf.dtypes]
                ]
            )
            > 0
        )

    def testDeepIntOptimization(self):
        data = {}
        for i in range(52):
            data[str(i)] = np.uint64(2**i + 1)
        df = pd.DataFrame([data])
        optimizedDf = reduce_mem_usage(df)
        for i in range(52):
            self.assertEqual(optimizedDf[str(i)][0], 2**i + 1)
