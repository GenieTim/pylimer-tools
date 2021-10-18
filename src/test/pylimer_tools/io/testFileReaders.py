import os
import unittest
from test.pdComparingTestCase import PandasComparingTestCase

import pandas as pd
from pylimer_tools.io.extractThermoParams import (extractThermoParams,
                                                  getThermoCacheNameSuffix)
from pylimer_tools.io.readLammpData import readLammpData
from pylimer_tools.io.readLammpDump import readLammpDump
from pylimer_tools.utils.cacheUtility import getCacheFileName
from pylimer_tools.utils.optimizeDf import reduce_mem_usage


class TestFileReader(PandasComparingTestCase):

    def test_thermoDataReader(self):
        thermoFile = os.path.join(os.path.dirname(
            __file__), "../fixtures/thermo_file.dat")
        header = "Step Temp E_pair E_mol TotEng Press"
        readData = extractThermoParams(
            thermoFile, header=header, useCache=False)
        self.assertIsInstance(readData, pd.DataFrame)
        # read again. this time from cache
        readData2 = extractThermoParams(thermoFile, header=header)
        self.assertIsInstance(readData2, pd.DataFrame)
        self.assertCountEqual(readData, readData2)
        # read again. this time with a header list (whatever)
        readData3 = extractThermoParams(
            thermoFile, header=[header], useCache=False)
        self.assertIsInstance(readData3, pd.DataFrame)
        self.assertCountEqual(readData3, readData2)
        reducedDf = reduce_mem_usage(readData)
        self.assertIsInstance(reducedDf, pd.DataFrame)
        self.assertFalse(reducedDf.empty)
        # cleanup: delete cache file
        os.remove(getCacheFileName(
            thermoFile, getThermoCacheNameSuffix(header)))

    def testLammpsDataReader(self):
        dataFile = os.path.join(os.path.dirname(
            __file__), "../fixtures/lammps_data_file.out")
        data = readLammpData(dataFile, useCache=False)
        self.assertIsInstance(data, dict)
        expectedKeys = ["N_atoms", "N_Atypes", "N_Btypes", "masses", "Lx", "Ly",
                        "Lz", "xlo", "xhi", "ylo", "yhi", "zlo", "zhi", "atom_data", "bond_data"]
        for key in expectedKeys:
            self.assertTrue(key in data)
            self.assertIsNotNone(data[key])
        # and with cache
        data2 = readLammpData(dataFile, useCache=True)
        for key in data:
            self.assertEqual(data2[key], data[key])

    def testLammpsDumpReader(self):
        dataFile = os.path.join(os.path.dirname(
            __file__), "../fixtures/lammps_dump.lammpstrj")
        data = readLammpDump(dataFile, useCache=False)
        self.assertIsInstance(data, list)
        data2 = readLammpDump(dataFile, useCache=True)
        self.assertEqual(len(data), len(data2))
