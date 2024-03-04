#!/usr/bin/env python
import os
import random
import unittest

import mock
import pandas as pd

from pylimer_tools.io.readLammpsOutputFile import readDataFile
from pylimer_tools_cpp.pylimer_tools_cpp import DataFileWriter, Universe


class DataFileWriterTest(unittest.TestCase):

    def check_imageFlagsFixed(self, universe: Universe, writer: DataFileWriter):
        file = os.path.join(os.path.dirname(__file__), "..",
                            "fixtures", "tmp-test-data-file.out")
        writer.writeToFile(file)

        read_universe = readDataFile(file)
        self.assertEqual(read_universe.getNrOfAtoms(), universe.getNrOfAtoms())

        # then, assert that after writing with the new image flags, things are correct
        molecule = read_universe.getMolecules(2)[0]
        self.assertAlmostEqual(molecule.computeEndToEndDistance(
        ), molecule.computeEndToEndDistanceWithDerivedImageFlags())

    def test_canWriteCorrectingImageFlags(self):
        random.seed(1290)
        universe = Universe(10, 10, 10)

        # add atoms that require new image flags
        universe.addAtoms(
            [i+1 for i in range(10)], [1 for _ in range(10)],
            [i+7. for i in range(10)], [0 for i in range(10)
                                        ], [0 for i in range(10)],
            [random.randint(0, 10) for i in range(10)], [random.randint(
                0, 10) for i in range(10)], [random.randint(0, 10) for i in range(10)]
        )

        universe.addBonds(
            [i+1 for i in range(9)],
            [i+2 for i in range(9)]
        )

        # first, make sure the random image flags actually
        # lead to what we expect not to be true
        molecule = universe.getMolecules(2)[0]
        self.assertNotAlmostEqual(molecule.computeEndToEndDistance(
        ), molecule.computeEndToEndDistanceWithDerivedImageFlags())

        writer = DataFileWriter(universe)
        writer.configMoveIntoBox(True)
        writer.configAttemptImageReset(True)
        self.check_imageFlagsFixed(universe, writer)

        writer.configCrosslinkerType(2)
        self.check_imageFlagsFixed(universe, writer)

        writer.configMoveIntoBox(False)
        self.check_imageFlagsFixed(universe, writer)


if __name__ == '__main__':
    unittest.main()
