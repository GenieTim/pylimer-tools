# Using the generator pattern (an iterable)
from pylimer_tools.entities.universum import Universum
from pylimer_tools.io.readLammpData import readLammpData
from pylimer_tools.io.readLammpDump import readLammpDump

"""
This class represents a sequence of Universes, 
as read from a LAMMPS dump file.
"""


class DumpFileUniverseSequence(object):

    def __init__(self, initialStructureFile, dumpFile):
        """
        Initialize this class.

        Arguments:
          - initialStructureFile: the file containing the initial structure, 
                file of type LAMMPS data 
          - dumpFile: the file containing the LAMMPS dump with the atom coordinates 
                at different time-steps.
        """
        self.initialStructure = readLammpData(initialStructureFile)
        self.parsedFile = readLammpDump(dumpFile)
        self.n = len(self.parsedFile)
        self.idx = 0
        self.universes = {}

    def assembleUniverse(self, idx) -> Universum:
        """
        Create a Universe object for a specified time-step.

        Arguments:
            - idx: the index of the data

        Returns:
            - universe (Universum): the universe for that index
        """
        currentTimestepsData = self.parsedFile[idx]
        boxSizes = []
        if ("BOX BOUNDS" in currentTimestepsData):
            boxSizes = [
                currentTimestepsData["BOX BOUNDS"][0][1] -
                currentTimestepsData["BOX BOUNDS"][0][1],
                currentTimestepsData["BOX BOUNDS"][1][1] -
                currentTimestepsData["BOX BOUNDS"][1][1],
                currentTimestepsData["BOX BOUNDS"][2][1] -
                currentTimestepsData["BOX BOUNDS"][2][1],
            ]
        else:
            boxSizes = [
                self.initialStructure["Lx"],
                self.initialStructure["Ly"],
                self.initialStructure["Lz"],
            ]
        universe = Universum(boxSizes)
        newAtoms = currentTimestepsData["ATOMS"]
        universe.addAtomBondData(newAtoms, self.initialStructure["bond_data"])
        self.universes[idx] = universe
        return universe

    def getUniverseAtIndex(self, idx) -> Universum:
        """
        Get the Universe at the given index (as of in the sequence given by the dump file).

        Arguments:
            - idx: The index to get the universe at

        Returns: 
            - universe (Universum): the Universe at the given index
        """
        if (idx in self.universes):
            return self.universes[idx]
        else:
            return self.parseFile(idx)

    def getUniverseAtTimestep(self, timestep) -> Universum:
        """
        Get the Universe at the given timestep.

        Arguments:
            - timestep: The timestep to get the universe at

        Returns: 
            - universe (Universum): the Universe at the given timestep, 
                `None` if there is no entry with this timestep in the dump file.
        """
        for i in range(self.n):
            if (self.parsedFile[i]["TIMESTEP"] == timestep):
                return self.getUniverseAtIndex(i)
        return None

    def __iter__(self):
        return self

    def __next__(self):
        return self.next()

    def next(self):
        if self.idx < self.n:
            cur, self.idx = self.idx, self.idx+1
            return self.getUniverseAtIndex(cur)
        raise StopIteration()
