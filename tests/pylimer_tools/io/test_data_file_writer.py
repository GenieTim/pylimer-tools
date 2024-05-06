#!/usr/bin/env python
import os
import random
import unittest

from pylimer_tools.io.read_lammps_output_file import read_data_file
from pylimer_tools_cpp import DataFileWriter, Universe


class DataFileWriterTest(unittest.TestCase):

    def check_image_flags_fixed(
            self, universe: Universe, writer: DataFileWriter):
        file = os.path.join(os.path.dirname(__file__), "..",
                            "fixtures", "tmp-test-data-file.out")
        writer.write_to_file(file)

        read_universe = read_data_file(file)
        self.assertEqual(read_universe.get_nr_of_atoms(),
                         universe.get_nr_of_atoms())

        # then, assert that after writing with the new image flags, things are
        # correct
        molecule = read_universe.get_molecules(2)[0]
        self.assertAlmostEqual(molecule.compute_end_to_end_distance(
        ), molecule.compute_end_to_end_distance_with_derived_image_flags())

    def test_can_write_correcting_image_flags(self):
        random.seed(1290)
        universe = Universe(10, 10, 10)

        # add atoms that require new image flags
        universe.add_atoms(
            [i + 1 for i in range(10)], [1 for _ in range(10)],
            [i + 7. for i in range(10)], [0 for i in range(10)
                                          ], [0 for i in range(10)],
            [random.randint(0, 10) for i in range(10)], [random.randint(
                0, 10) for i in range(10)], [random.randint(0, 10) for i in range(10)]
        )

        universe.add_bonds(
            [i + 1 for i in range(9)],
            [i + 2 for i in range(9)]
        )

        # first, make sure the random image flags actually
        # lead to what we expect not to be true
        molecule = universe.get_molecules(2)[0]
        self.assertNotAlmostEqual(molecule.compute_end_to_end_distance(
        ), molecule.compute_end_to_end_distance_with_derived_image_flags())

        writer = DataFileWriter(universe)
        writer.config_move_into_box(True)
        writer.config_attempt_image_reset(True)
        self.check_image_flags_fixed(universe, writer)

        writer.config_crosslinker_type(2)
        self.check_image_flags_fixed(universe, writer)

        writer.config_move_into_box(False)
        self.check_image_flags_fixed(universe, writer)


if __name__ == '__main__':
    unittest.main()
