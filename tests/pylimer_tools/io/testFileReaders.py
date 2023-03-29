import os
import unittest

import mock
import pandas as pd

from pylimer_tools.io.extractThermoParams import (detectHeaders,
                                                  extractThermoParams,
                                                  getThermoCacheNameSuffix)
from pylimer_tools.io.readLammpsOutputFile import (readAveragesFile,
                                                   readCorrelationFile,
                                                   readDataFile, readDumpFile,
                                                   readHistogramFile,
                                                   readLogFile)
from pylimer_tools.utils.cacheUtility import getCacheFileName
from pylimer_tools.utils.optimizeDf import optimize, reduce_mem_usage
from pylimer_tools_cpp import (DataFileReader, DumpFileReader, Universe,
                               UniverseSequence)
from tests.pylimer_tools.pdComparingTestCase import PandasComparingTestCase


class TestFileReader(PandasComparingTestCase):

    def test_thermoDataReader(self):
        thermoFile = os.path.join(os.path.dirname(
            __file__), "../fixtures/thermo_file.dat")
        # this header does not include "Volume"
        header = "Step Temp E_pair E_mol TotEng Press"
        readData = extractThermoParams(
            thermoFile, header=header, useCache=False)
        self.assertIsInstance(readData, pd.DataFrame)
        self.assertTrue("Volume" in readData.columns)
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
        # test: empty file
        emptyData = extractThermoParams(os.path.join(
            os.path.dirname(__file__), '../fixtures/empty_file.txt'))
        self.assertTrue(emptyData.empty)

        # test: check if we can deduce the header
        print("Testing header detection")
        detectedHeaders = detectHeaders(
            thermoFile, max_nr_of_lines_to_read=150)
        self.assertListEqual(detectedHeaders, [
            "Step Temp E_pair E_mol TotEng Press"
        ])
        readData5 = extractThermoParams(thermoFile, header=detectedHeaders)
        detectedHeaders = detectHeaders(
            thermoFile, max_nr_of_lines_to_read=1500000)
        self.assertListEqual(detectedHeaders, [
            "Step Temp E_pair E_mol TotEng Press",
            "Step Temp E_pair E_mol TotEng Press Volume",
            "Step Temp E_pair E_mol TotEng Press",
            "Step Temp PotEng 2",
            "Step Temp PotEng 2",
        ])
        readData4 = readLogFile(thermoFile, lines_to_read_to_detect_header=500)
        self.assertListEqual(list(readData4.columns), list(readData5.columns))
        self.assertListEqual(list(readData.columns), list(readData4.columns))
        self.assertDataframeEqual(readData4, readData5, check_dtype=False)
        self.assertDataframeEqual(readData, reduce_mem_usage(readData4), check_dtype=False)

    @mock.patch('pylimer_tools.io.extractThermoParams.os.remove')
    def test_cacheDeleteFail(self, mockOsRemove):
        mockOsRemove.side_effect = Exception
        thermoFile = os.path.join(os.path.dirname(
            __file__), "../fixtures/thermo_file.dat")
        header = "Step Temp E_pair E_mol TotEng Press Volume"
        readData = extractThermoParams(
            thermoFile, header=header, useCache=False)
        self.assertIsInstance(readData, pd.DataFrame)
        self.assertEqual(" ".join(readData.columns), header)
        self.assertTrue(mockOsRemove.called)

    @mock.patch('pylimer_tools.io.extractThermoParams.pd.read_csv')
    def test_pdCsvReadFail(self, mockCsvReader):
        mockCsvReader.side_effect = Exception
        thermoFile = os.path.join(os.path.dirname(
            __file__), "../fixtures/thermo_file.dat")
        header = "Step Temp E_pair E_mol TotEng Press Volume"
        with self.assertWarns(Warning):
            readData = extractThermoParams(
                thermoFile, header=header, useCache=False)
            self.assertIsInstance(readData, pd.DataFrame)
            self.assertTrue(readData.empty)
            self.assertTrue(mockCsvReader.called)

    def test_LammpsDataReader(self):
        dataFile = os.path.join(os.path.dirname(
            __file__), "../fixtures/lammps_data_file.out")
        universe = readDataFile(dataFile)
        self.assertIsInstance(universe, Universe)
        self.assertEqual(universe.getNrOfAtoms(), 3000)

        # expectedKeys = ["N_atoms", "N_Atypes", "N_Btypes", "masses", "Lx", "Ly",
        #                 "Lz", "xlo", "xhi", "ylo", "yhi", "zlo", "zhi", "atom_data", "bond_data"]
        # for key in expectedKeys:
        #     self.assertTrue(key in data)
        #     self.assertIsNotNone(data[key])
        # # and with cache
        # data2 = readLammpData(dataFile, useCache=True)
        # for key in data:
        #     self.assertEqual(data2[key], data[key])

    def test_LammpsDumpReader(self):
        dataFile = os.path.join(os.path.dirname(
            __file__), "../fixtures/lammps_data_file_small.out")
        dumpFile = os.path.join(os.path.dirname(
            __file__), "../fixtures/lammps_dump_small.lammpstrj")
        universeSequence = readDumpFile(dataFile, dumpFile)
        self.assertIsInstance(universeSequence, UniverseSequence)
        universe = universeSequence.atIndex(0)
        self.assertIsInstance(universe, Universe)
        self.assertEqual(universeSequence.getLength(), 1)
        universe = universeSequence.atIndex(0)
        self.assertEqual(universe.getNrOfAtoms(), 12)

    def test_avgReader(self):
        dataFile = os.path.join(os.path.dirname(
            __file__), "../fixtures/example_avg_file.out.avg.txt")
        data = readAveragesFile(dataFile)
        self.assertEqual(len(data), 5)
        self.assertEqual(data["TimeStep"].iloc[0], 100)
        self.assertEqual(data["TimeStep"].iloc[4], 500)
        self.assertEqual(data["else"].iloc[4], 9000)

    def test_vecAvgReader(self):
        dataFile = os.path.join(os.path.dirname(
            __file__), "../fixtures/example_vec_avg_file.out.vec-avg.txt")
        data = readAveragesFile(dataFile)
        self.assertIsInstance(data, pd.DataFrame)
        self.assertEqual(len(data), 9)
        self.assertEqual(data["TimeStep"].iloc[0], 100)
        self.assertEqual(data["TimeStep"].iloc[4], 200)
        self.assertEqual(len(data["TimeStep"].unique()), 3)
        self.assertEqual(data["value2"].iloc[4], 4000)
        self.assertEqual(data["value1"].iloc[4], 2.3)
        data2 = readHistogramFile(dataFile)
        self.assertDataframeEqual(data, data2)

    def test_correlationReader(self):
        dataFile = os.path.join(os.path.dirname(
            __file__), "../fixtures/example_correlation_file.out.corr.txt")
        data = readCorrelationFile(dataFile)
        self.assertIsInstance(data, pd.DataFrame)
        self.assertEqual(len(data["Timestep"].unique()), 2)
        self.assertEqual(len(data), 118)
