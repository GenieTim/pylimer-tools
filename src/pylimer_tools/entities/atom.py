from __future__ import annotations

import numpy as np
import pandas as pd


class Atom:
    def __init__(self, data: pd.Series, boxSizes: list) -> None:
        """
        Instantiate the Atom.

        Arguments:
          - data: the data underlying the Atom
          - boxSizes: the size of the box the atom is in. Used for periodic image computations.
        """
        self.data = data
        self.boxDimensions = {
            "x": boxSizes[0], "y": boxSizes[1], "z": boxSizes[2]
        }

    def _getDeltaDistance(self, direction: str, distanceTo: Atom) -> float:
        """
        Calculate the distance in one dimension between two atoms, 
        accounting for periodic displacements.

        Arguments:
          - direction: the dimension to calculate the distance in
          - distanceTo: the other atom to calculate the distance between

        Returns:
          - delta: the distance in the one direction
        """
        delta = abs(self.data[direction] - distanceTo.data[direction])
        if (self.data["n"+direction] != distanceTo.data["n"+direction]):
            delta -= (self.data["n"+direction] -
                      distanceTo.data["n" + direction])*self.boxDimensions[direction]
        return delta

    def getDeltaX(self, secondAtom: Atom) -> float:
        """
        Calculate the distance in the x dimension between this and another atom, 
        accounting for periodic displacements.

        Arguments:
          - secondAtom: the other atom to calculate the distance between

        Returns:
          - delta: the distance in the x direction
        """
        return self._getDeltaDistance("x", secondAtom)

    def getDeltaY(self, secondAtom: Atom) -> float:
        """
        Calculate the distance in the y dimension between this and another atom, 
        accounting for periodic displacements.

        Arguments:
          - secondAtom: the other atom to calculate the distance between

        Returns:
          - delta: the distance in the y direction
        """
        return self._getDeltaDistance("y", secondAtom)

    def getDeltaZ(self, secondAtom: Atom) -> float:
        """
        Calculate the distance in the z dimension between this and another atom, 
        accounting for periodic displacements.

        Arguments:
          - secondAtom: the other atom to calculate the distance between

        Returns:
          - delta: the distance in the z direction
        """
        return self._getDeltaDistance("z", secondAtom)

    def computeDistanceTo(self, secondAtom: Atom) -> float:
        """
        Calculate the the distance between two atoms. 

        Arguments:
            - secondAtom: the atom to compute the distance to

        Returns:
            - meanDistance: the norm of the connecting vector between the two coordinates
        """
        return np.linalg.norm([self.getDeltaX(secondAtom), self.getDeltaY(
            secondAtom), self.getDeltaZ(secondAtom)
        ])

    def getUnderlyingData(self) -> pd.Series:
        """
        Auxilary method to get the pd.Series data associated with this atom

        Returns:
          - data (pd.Series): the data as given to this atom upon instantiation
        """
        return self.data
