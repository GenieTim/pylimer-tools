import unittest

import pandas as pd

from pylimer_tools.entities.universum import Universum


class UniverseUsingTestCase(unittest.TestCase):

    testUniverse: Universum
    testUniverseSmall: Universum
    saturatedTestUniverse: Universum

    # The system looks like this (in terms of bonds, not 3D placement):
    # 1-2-3-*6
    # |      |
    # *7-5---|
    # 8
    #
    # *4
    testAtoms = pd.DataFrame([
        {"id": 1, "nx": 1, "ny": 1, "nz": 1,
         "type": 1, "x": 1, "y": 1, "z": 1},
        {"id": 2, "nx": 1, "ny": 1, "nz": 1,
         "type": 1, "x": 2, "y": 1, "z": 1},
        {"id": 3, "nx": 1, "ny": 1, "nz": 1,
         "type": 1, "x": 3, "y": 1, "z": 1},
        {"id": 4, "nx": 1, "ny": 1, "nz": 1,
         "type": 2, "x": 2, "y": 2, "z": 1},
        {"id": 5, "nx": 1, "ny": 1, "nz": 1,
         "type": 1, "x": 1, "y": 3, "z": 1},
        {"id": 6, "nx": 1, "ny": 1, "nz": 1,
         "type": 2, "x": 1, "y": 1, "z": 2},
        {"id": 7, "nx": 1, "ny": 1, "nz": 1,
         "type": 2, "x": 1, "y": 1, "z": 3},
        {"id": 8, "nx": 1, "ny": 1, "nz": 1,
         "type": 1, "x": 2, "y": 2, "z": 2},
    ])
    testBonds = pd.DataFrame([
        {"to": 1, "bondFrom": 2},
        {"to": 3, "bondFrom": 2},
        {"to": 5, "bondFrom": 6},
        {"to": 1, "bondFrom": 7},
        {"to": 5, "bondFrom": 7},
        {"to": 3, "bondFrom": 6},
        {"to": 7, "bondFrom": 8}
    ])

    """
    This system looks like:

    1-2-3

    *7-*6-5
    """
    testAtomsSmall = pd.DataFrame([
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
    testBondsSmall = pd.DataFrame([
        {"to": 1, "bondFrom": 2},
        {"to": 3, "bondFrom": 2},
        {"to": 5, "bondFrom": 6},
        {"to": 6, "bondFrom": 7}
    ])

    def setUp(self):
        # Universe 1: empy
        self.emptyUniverse = Universum([10, 10, 10])
        # Universe 2: small (5 atoms), barely connected in two chains
        self.testUniverseSmall = Universum([10, 10, 10])
        self.testUniverseSmall.addAtomBondData(
            self.testAtomsSmall, self.testBondsSmall)
        # Universe 3: very unsaturated
        self.testUniverse = Universum([10, 10, 10])
        self.testUniverse.addAtomBondData(self.testAtoms, self.testBonds)
        # an additional larget test universe where the stoichiometric inbalance is < 1
        # even when imposing a crosslinker functionality of 1
        # in essence, it is on loop around 4 plus a connction to 6.
        self.saturatedTestUniverse = Universum([10, 10, 10])
        self.saturatedTestUniverse.addAtomBondData(self.testAtoms.append([
            {"id": 9, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 10, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 11, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 12, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 13, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 14, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 15, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 16, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 17, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 18, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 19, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
            {"id": 20, "type": 1, "nx": 1, "ny": 1,
                "nz": 1, "x": 1, "y": 1, "z": 1},
        ]), self.testBonds.append([
            {"to": 9, "bondFrom": 4},
            {"to": 10, "bondFrom": 9},
            {"to": 11, "bondFrom": 10},
            {"to": 12, "bondFrom": 11},
            {"to": 13, "bondFrom": 12},
            {"to": 4, "bondFrom": 13},
            {"to": 14, "bondFrom": 4},
            {"to": 15, "bondFrom": 14},
            {"to": 16, "bondFrom": 15},
            {"to": 17, "bondFrom": 16},
            {"to": 18, "bondFrom": 17},
            {"to": 19, "bondFrom": 18},
            {"to": 20, "bondFrom": 19},
            {"to": 6, "bondFrom": 20},
        ]))
