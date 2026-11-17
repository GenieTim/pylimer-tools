#!/usr/bin/env python
import os
import unittest

import mock
import pandas as pd

from pylimer_tools.io.extract_thermo_data import (
    detect_headers,
    extract_thermo_params,
    get_thermo_cache_name_suffix,
    read_multi_section_separated_value_file,
)
from pylimer_tools.io.read_lammps_output_file import (
    read_averages_file,
    read_correlation_file,
    read_data_file,
    read_dump_file,
    read_histogram_file,
    read_log_file,
)
from pylimer_tools.utils.cache_utility import get_cache_file_name
from pylimer_tools.utils.optimize_dataframe import reduce_mem_usage
from pylimer_tools_cpp import AtomStyle, Universe, UniverseSequence
from tests.pylimer_tools.pdComparingTestCase import PandasComparingTestCase


class TestFileReader(PandasComparingTestCase):
    def test_thermo_data_reader(self):
        thermo_file = os.path.join(
            os.path.dirname(__file__), "../fixtures/thermo_file.dat"
        )
        # this header does not include "Volume"
        header = "Step Temp E_pair E_mol TotEng Press"
        read_data = extract_thermo_params(thermo_file, header=header, use_cache=False)
        self.assertIsInstance(read_data, pd.DataFrame)
        self.assertTrue("Volume" in read_data.columns)
        # read again. this time from cache
        read_data2 = extract_thermo_params(thermo_file, header=header)
        self.assertIsInstance(read_data2, pd.DataFrame)
        self.assertCountEqual(read_data, read_data2)
        # read again. this time with a header list (whatever)
        read_data3 = extract_thermo_params(
            thermo_file, header=[header], use_cache=False
        )
        self.assertIsInstance(read_data3, pd.DataFrame)
        self.assertCountEqual(read_data3, read_data2)
        reduced_df = reduce_mem_usage(read_data)
        self.assertIsInstance(reduced_df, pd.DataFrame)
        self.assertFalse(reduced_df.empty)
        # cleanup: delete cache file
        os.remove(
            get_cache_file_name(thermo_file, get_thermo_cache_name_suffix(header))
        )
        # test: empty file
        empty_data = extract_thermo_params(
            os.path.join(os.path.dirname(__file__), "../fixtures/empty_file.txt")
        )
        self.assertTrue(empty_data.empty)

        # test: check if we can deduce the header
        detected_headers = detect_headers(thermo_file, max_nr_of_lines_to_read=150)
        self.assertListEqual(detected_headers, ["Step Temp E_pair E_mol TotEng Press"])
        # again, this time from cache
        detected_headers2 = detect_headers(thermo_file, max_nr_of_lines_to_read=150)
        self.assertListEqual(detected_headers2, detected_headers)
        # use them to read the file
        read_data5 = extract_thermo_params(thermo_file, header=detected_headers)
        detected_headers = detect_headers(thermo_file, max_nr_of_lines_to_read=1500000)
        self.assertListEqual(
            detected_headers,
            [
                "Step Temp E_pair E_mol TotEng Press",
                "Step Temp E_pair E_mol TotEng Press Volume",
                "Step Temp E_pair E_mol TotEng Press",
                "Step Temp PotEng 2",
                "Step Temp PotEng 2",
            ],
        )
        read_data4 = read_log_file(thermo_file, lines_to_read_to_detect_header=500)
        self.assertListEqual(list(read_data4.columns), list(read_data5.columns))
        self.assertListEqual(list(read_data.columns), list(read_data4.columns))
        self.assertDataframeEqual(read_data4, read_data5, check_dtype=False)
        self.assertDataframeEqual(
            read_data, reduce_mem_usage(read_data4), check_dtype=False
        )

    @mock.patch("pylimer_tools.io.extract_thermo_data.os.remove")
    def test_cache_delete_fail(self, mock_os_remove):
        mock_os_remove.side_effect = Exception
        thermo_file = os.path.join(
            os.path.dirname(__file__), "../fixtures/thermo_file.dat"
        )
        header = "Step Temp E_pair E_mol TotEng Press Volume"
        read_data = extract_thermo_params(thermo_file, header=header, use_cache=False)
        self.assertIsInstance(read_data, pd.DataFrame)
        self.assertEqual(" ".join(read_data.columns), header)
        self.assertTrue(mock_os_remove.called)

    @mock.patch("pylimer_tools.io.extract_thermo_data.pd.read_csv")
    def test_pd_csv_read_fail(self, mock_csv_reader):
        mock_csv_reader.side_effect = Exception
        thermo_file = os.path.join(
            os.path.dirname(__file__), "../fixtures/thermo_file.dat"
        )
        header = "Step Temp E_pair E_mol TotEng Press Volume"
        with self.assertWarns(Warning):
            read_data = extract_thermo_params(
                thermo_file, header=header, use_cache=False
            )
            self.assertIsInstance(read_data, pd.DataFrame)
            self.assertTrue(read_data.empty)
            self.assertTrue(mock_csv_reader.called)

    def test_lammps_data_reader(self):
        data_file = os.path.join(
            os.path.dirname(__file__), "../fixtures/lammps_data_file.out"
        )
        universe = read_data_file(data_file, [AtomStyle.BOND])
        self.assertIsInstance(universe, Universe)
        self.assertEqual(universe.get_nr_of_atoms(), 3000)

        # expectedKeys = ["N_atoms", "N_Atypes", "N_Btypes", "masses", "Lx", "Ly",
        #                 "Lz", "xlo", "xhi", "ylo", "yhi", "zlo", "zhi", "atom_data", "bond_data"]
        # for key in expectedKeys:
        #     self.assertTrue(key in data)
        #     self.assertIsNotNone(data[key])
        # # and with cache
        # data2 = readLammpData(data_file, use_cache=True)
        # for key in data:
        #     self.assertEqual(data2[key], data[key])

    def test_lammps_dump_reader(self):
        data_file = os.path.join(
            os.path.dirname(__file__), "../fixtures/lammps_data_file_small.out"
        )
        dump_file = os.path.join(
            os.path.dirname(__file__), "../fixtures/lammps_dump_small.lammpstrj"
        )
        universe_sequence = read_dump_file(data_file, dump_file)
        self.assertIsInstance(universe_sequence, UniverseSequence)
        universe = universe_sequence.at_index(0)
        self.assertIsInstance(universe, Universe)
        self.assertEqual(universe_sequence.get_length(), 1)
        self.assertEqual(len(universe_sequence), 1)
        universe = universe_sequence.at_index(0)
        self.assertEqual(universe.get_nr_of_atoms(), 12)
        universe_sequence2 = read_dump_file(data_file, dump_file, [AtomStyle.BOND])
        self.assertIsInstance(universe_sequence2, UniverseSequence)
        universe = universe_sequence2.at_index(0)
        self.assertEqual(universe.get_nr_of_atoms(), 12)

    def test_averages_reader(self):
        data_file = os.path.join(
            os.path.dirname(__file__), "../fixtures/example_avg_file.out.avg.txt"
        )
        data = read_averages_file(data_file)
        self.assertEqual(len(data), 5)
        self.assertEqual(data["TimeStep"].iloc[0], 100)
        self.assertEqual(data["TimeStep"].iloc[4], 500)
        self.assertEqual(data["else"].iloc[4], 9000)

    def test_vector_averages_reader(self):
        data_file = os.path.join(
            os.path.dirname(__file__),
            "../fixtures/example_vec_avg_file.out.vec-avg.txt",
        )
        data = read_averages_file(data_file)
        self.assertIsInstance(data, pd.DataFrame)
        self.assertEqual(len(data), 9)
        self.assertEqual(data["TimeStep"].iloc[0], 100)
        self.assertEqual(data["TimeStep"].iloc[4], 200)
        self.assertEqual(len(data["TimeStep"].unique()), 3)
        self.assertEqual(data["value2"].iloc[4], 4000)
        self.assertEqual(data["value1"].iloc[4], 2.3)
        data2 = read_histogram_file(data_file)
        self.assertDataframeEqual(data, data2)

    def test_correlation_reader(self):
        data_file = os.path.join(
            os.path.dirname(__file__),
            "../fixtures/example_correlation_file.out.corr.txt",
        )
        data = read_correlation_file(data_file)
        self.assertIsInstance(data, pd.DataFrame)
        self.assertEqual(len(data["Timestep"].unique()), 2)
        self.assertEqual(len(data), 118)
        # again, from cache
        data2 = read_correlation_file(data_file)
        self.assertDataframeEqual(data, data2)

    def test_multisection_file(self):
        data_file = os.path.join(
            os.path.dirname(__file__), "../fixtures/example_multisection_value.txt"
        )
        data = read_multi_section_separated_value_file(data_file)
        self.assertIsInstance(data, pd.DataFrame)
        self.assertEqual(len(data["Header3"].unique()), 3)
        self.assertEqual(len(data["Step"].unique()), 4)
        self.assertEqual(len(data), 4)
        # again, from cache
        data2 = read_multi_section_separated_value_file(data_file)
        self.assertDataframeEqual(data, data2)


if __name__ == "__main__":
    unittest.main()
