
import unittest
import numpy as np

import pandas as pd
from pylimer_tools_cpp import Box, Universe


class TestEntityCalculations(unittest.TestCase):
    def test_universe(self):
        universe = Universe(10, 10, 10)
        self.assertIsInstance(universe, Universe)
        self.assertEqual(universe.getVolume(), 10*10*10)
        universe.setBox(Box(1, 1, 1))
        self.assertEqual(universe.getVolume(), 1)
        universe.setBoxLengths(100, 1, 1)
        self.assertEqual(universe.getVolume(), 100*1*1)

    def test_calculateMeanBondLen(self):
        baseAtom = {
            "id": 1,
            "x": 0,
            "y": 0,
            "z": 0,
            "nx": 0,
            "ny": 0,
            "nz": 0,
            "type": 1
        }

        bondsDf = pd.DataFrame([{"to": 2, "from": 1}, {"to": 3, "from": 2}])

        for dir in ["x", "y", "z"]:
            secondAtom = baseAtom.copy()
            secondAtom[dir] = 1
            secondAtom["id"] = 2

            thirdAtom = baseAtom.copy()
            thirdAtom[dir] = 2
            thirdAtom["id"] = 3

            coordsDf = pd.DataFrame([thirdAtom, secondAtom, baseAtom])
            universe = Universe(2, 2, 2)
            universe.addAtoms(len(coordsDf), coordsDf["id"].tolist(), coordsDf["type"].tolist(),
                              coordsDf["x"].tolist(), coordsDf["y"].tolist(),
                              coordsDf["z"].tolist(),
                              coordsDf["nx"].tolist(), coordsDf["ny"].tolist(), coordsDf["nz"].tolist())
            universe.addBonds(len(bondsDf),
                              bondsDf["from"].tolist(), bondsDf["to"].tolist())
            self.assertEqual(len(universe.getMolecules(-1)), 1)
            molecule = universe.getMolecules(-1)[0]
            self.assertEqual(molecule.getLength(), 3)
            self.assertEqual(molecule.getNrOfBonds(), 2)
            self.assertEqual(np.mean(molecule.computeBondLengths()), 1)

    def test_calculateEndToEndDistance(self):
        baseAtom = {
            "id": 1,
            "x": 0,
            "y": 0,
            "z": 0,
            "nx": 0,
            "ny": 0,
            "nz": 0,
            "type": 1
        }
        universe = Universe(1, 1, 1)
        coordsDf = pd.DataFrame([baseAtom])
        bondsDf = pd.DataFrame([], columns=["to", "bondFrom"])
        universe.addAtoms(len(coordsDf), coordsDf["id"].tolist(), coordsDf["type"].tolist(),
                          coordsDf["x"].tolist(), coordsDf["y"].tolist(),
                          coordsDf["z"].tolist(),
                          coordsDf["nx"].tolist(), coordsDf["ny"].tolist(), coordsDf["nz"].tolist())
        universe.addBonds(len(bondsDf),
                          bondsDf["bondFrom"].tolist(), bondsDf["to"].tolist())
        molecules = universe.getMolecules(-1)
        self.assertEqual(len(molecules), 1)
        self.assertEqual(molecules[0].computeEndToEndDistance(), 0)
        atoms = []
        for i in range(3):
            newAtom = baseAtom.copy()
            newAtom["id"] = i
            atoms.append(newAtom)
        universe = Universe(1, 1, 1)
        coordsDf = pd.DataFrame(atoms)
        bondsDf = pd.DataFrame(
            [{"to": 1, "bondFrom": 0}, {"to": 2, "bondFrom": 1}, {"to": 0, "bondFrom": 2}])
        universe.addAtoms(len(coordsDf), coordsDf["id"].tolist(), coordsDf["type"].tolist(),
                          coordsDf["x"].tolist(), coordsDf["y"].tolist(),
                          coordsDf["z"].tolist(),
                          coordsDf["nx"].tolist(), coordsDf["ny"].tolist(), coordsDf["nz"].tolist())
        universe.addBonds(len(bondsDf),
                          bondsDf["bondFrom"].tolist(), bondsDf["to"].tolist())
        molecules = universe.getMolecules(-1)
        self.assertEqual(len(molecules), 1)
        self.assertEqual(-1, molecules[0].computeEndToEndDistance())

    def test_calculateDistanceThroughPeriodicImage(self):
        baseAtom = {
            "id": 1,
            "x": 0,
            "y": 0,
            "z": 0,
            "nx": 0,
            "ny": 0,
            "nz": 0,
            "type": 1
        }
        for dir in ["x", "y", "z"]:
            secondAtom = baseAtom.copy()
            secondAtom[dir] = 1
            secondAtom["id"] = 2
            secondAtom["n" + dir] = 1
            universe = Universe(1, 1, 1)
            coordsDf = pd.DataFrame([baseAtom, secondAtom])
            bondsDf = pd.DataFrame([{
                "to": 1, "bondFrom": 2
            }])
            universe.addAtoms(len(coordsDf), coordsDf["id"].tolist(), coordsDf["type"].tolist(),
                              coordsDf["x"].tolist(), coordsDf["y"].tolist(),
                              coordsDf["z"].tolist(),
                              coordsDf["nx"].tolist(), coordsDf["ny"].tolist(), coordsDf["nz"].tolist())
            universe.addBonds(len(bondsDf),
                              bondsDf["bondFrom"].tolist(), bondsDf["to"].tolist())
            self.assertEqual(len(universe.getMolecules(-1)), 1)
            self.assertEqual(0.0, np.mean(universe.getMolecules(-1)[
                0].computeBondLengths()))


if __name__ == '__main__':
    unittest.main()
