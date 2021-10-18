import numpy as np
import pandas as pd

from test.pdComparingTestCase import PandasComparingTestCase
from pylimer_tools.utils.optimizeDf import *

class TestOptimizeDf(PandasComparingTestCase):
  def testSeparateOptimizatios(self):
    class Object(object):
      pass
      
    testObject = Object()
    df = pd.DataFrame([{
      "testFloat": np.float64(1.01),
      "testInt": np.int64(1e1),
      "testObject": testObject
    }])
    optimizedDf = optimize(df)
    expectedOptimizedDf = pd.DataFrame([{
      "testFloat": np.float32(1.01),
      "testInt": np.int8(1e1),
      "testObject": testObject
    }])
    # contents are equal
    self.assertEqual(df, optimizedDf)
    self.assertEqual(expectedOptimizedDf, optimizedDf)

  def testDeepIntOptimization(self):
    data = {}
    for i in range(52):
      data[str(i)] = np.uint64(2**i+1)
    df = pd.DataFrame([data])
    optimizedDf = reduce_mem_usage(df)
    for i in range(52):
      self.assertEqual(optimizedDf[str(i)][0], 2**i+1)
