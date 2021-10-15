from __future__ import annotations

import numpy as np
import pandas as pd


class Atom:

    name: str

    def __init__(self, data: pd.Series, boxSizes: list, name: str = None):
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
        if (name is None):
            self.name = "atom{}".format(data["id"])
        else:
            self.name = name

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

    def computeVectorTo(self, secondAtom: Atom) -> float:
        """
        Calculate the the vector between two atoms. 

        Arguments:
            - secondAtom: the atom to compute the distance to

        Returns:
            - difference (np.array): the connecting vector between the two coordinates
        """
        # using a numpy array is a bit too much overhead, possibly
        return np.array([self.getDeltaX(secondAtom), self.getDeltaY(
            secondAtom), self.getDeltaZ(secondAtom)
        ])

    def computeDistanceTo(self, secondAtom: Atom) -> float:
        """
        Calculate the the distance between two atoms. 

        Arguments:
            - secondAtom: the atom to compute the distance to

        Returns:
            - meanDistance: the norm of the connecting vector between the two coordinates
        """
        # again, using numpy here is a bit overkill
        return np.linalg.norm(self.computeVectorTo(secondAtom))

    def getUnderlyingData(self) -> pd.Series:
        """
        Auxilary method to get the pd.Series data associated with this atom

        Returns:
          - data (pd.Series): the data as given to this atom upon instantiation
        """
        return self.data
