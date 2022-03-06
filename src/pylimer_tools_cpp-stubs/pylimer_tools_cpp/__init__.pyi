"""
    PylimerTools Cpp
    -----------------

    A collection of utility python functions for handling LAMMPS output and polymers in Python.

    .. currentmodule:: pylimer_tools_cpp
    .. autosummary::
        :toctree: _generate

    """
import typing

import pylimer_tools_cpp.pylimer_tools_cpp

__all__ = [
    "Atom",
    "Box",
    "DANGLING_CHAIN",
    "DataFileReader",
    "DumpFileReader",
    "FREE_CHAIN",
    "Molecule",
    "MoleculeIterator",
    "MoleculeType",
    "NETWORK_STRAND",
    "PRIMARY_LOOP",
    "UNDEFINED",
    "Universe",
    "UniverseSequence",
    "computeStoichiometricInbalance",
    "predictGelationPoint"
]


class Atom():
    def __getstate__(self) -> tuple: ...
    def __init__(self, arg0: int, arg1: int, arg2: float, arg3: float, arg4: float, arg5: int, arg6: int, arg7: int) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def computeVectorTo(self, arg0: Atom, arg1: Box) -> typing.List[float]: 
        """
        Compute the vector to another atom.
        """
    def distanceTo(self, arg0: Atom, arg1: Box) -> float: 
        """
        Compute the distance to another atom.
        """
    def distanceToUnwrapped(self, arg0: Atom, arg1: Box) -> float: 
        """
        Compute the distance to another atom respecting the periodic image flag.
        """
    def getId(self) -> int: 
        """
        Get the id of the atom.
        """
    def getNX(self) -> int: 
        """
        Get the box image that the atom is in in x direction (also known as `ix` or `nx`).
        """
    def getNY(self) -> int: 
        """
        Get the box image that the atom is in in y direction (also known as `iy` or `ny`).
        """
    def getNZ(self) -> int: 
        """
        Get the box image that the atom is in in z direction (also known as `iz` or `nz`).
        """
    def getType(self) -> int: 
        """
        Get the type of the atom.
        """
    def getX(self) -> float: 
        """
        Get the x coordinate of the atom.
        """
    def getY(self) -> float: 
        """
        Get the y coordinate of the atom.
        """
    def getZ(self) -> float: 
        """
        Get the z coordinate of the atom.
        """
    def vectorToUnwrapped(self, arg0: Atom, arg1: Box, arg2: float) -> None: 
        """
        Compute the vector to another atom respecting the periodic image flag.
        """
    pass
class Box():
    """
            The box that the simulation is run in.

            *NOTE*: currently, only rectangular boxes are supported.
            
    """
    def __getstate__(self) -> tuple: ...
    def __init__(self, arg0: float, arg1: float, arg2: float) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def getLx(self) -> float: 
        """
        Get the length of the box in x direction.
        """
    def getLy(self) -> float: 
        """
        Get the length of the box in y direction.
        """
    def getLz(self) -> float: 
        """
        Get the length of the box in z direction.
        """
    def getVolume(self) -> float: 
        """
                    Compute the volume of the box.

                    :math:`V = L_x \cdot L_y \cdot L_z`
                    
        """
    pass
class DataFileReader():
    def __init__(self) -> None: ...
    def getAtomIds(self) -> typing.List[int]: ...
    def getAtomNx(self) -> typing.List[int]: ...
    def getAtomNy(self) -> typing.List[int]: ...
    def getAtomTypes(self) -> typing.List[int]: ...
    def getAtomX(self) -> typing.List[float]: ...
    def getAtomY(self) -> typing.List[float]: ...
    def getBondFrom(self) -> typing.List[int]: ...
    def getBondTo(self) -> typing.List[int]: ...
    def getBondTypes(self) -> typing.List[int]: ...
    def getLx(self) -> float: ...
    def getLy(self) -> float: ...
    def getLz(self) -> float: ...
    def getMasses(self) -> typing.Dict[int, float]: ...
    def getMoleculeIds(self) -> typing.List[int]: ...
    def getNrOfAtomTypes(self) -> int: ...
    def getNrOfAtoms(self) -> int: ...
    def getNrOfBondTypes(self) -> int: ...
    def getNrOfBonds(self) -> int: ...
    def read(self, arg0: str) -> None: ...
    pass
class DumpFileReader():
    def __init__(self) -> None: ...
    def getLength(self) -> int: 
        """
        Get the number of sections in the file
        """
    def getNumericValuesForAt(self, arg0: int, arg1: str, arg2: str) -> typing.List[float]: 
        """
        Get the values for the section `index`, the main header `headerKey` and the column (in the header) `column`.
        """
    def getStringValuesForAt(self, arg0: int, arg1: str, arg2: str) -> typing.List[str]: 
        """
        Get the values for the section `index`, the main header `headerKey` and the column (in the header) `column`.
        """
    def hasKey(self, arg0: str) -> bool: 
        """
        Check whether the first section has the header specified
        """
    def keyHasColumn(self, arg0: str, arg1: str) -> bool: 
        """
        Check whether the header of the first section has the specified column
        """
    def keyHasDirectionalColumn(self, arg0: str, arg1: str, arg2: str) -> bool: 
        """
        Check whether the header of the first section has all the three columns `{dirPraefix}{x|y|z}{dirSuffix}`.
        """
    def read(self, arg0: str) -> None: 
        """
        Read a file
        """
    pass
class Molecule():
    def __getitem__(self, arg0: int) -> Atom: ...
    def __init__(self, arg0: Box, arg1: igraph_s, arg2: MoleculeType) -> None: ...
    def __iter__(self) -> MoleculeIterator: ...
    def __len__(self) -> int: ...
    def computeBondLengths(self) -> typing.List[float]: 
        """
        Computes the length :math:`b` of each bond in the molecule.
        """
    def computeEndToEndDistance(self) -> float: 
        """
        Compute the end-to-end distance (:math:`R_{ee}`) of this molecule. 

        Caution:
            Returns 0.0 if the molecule does not have two or more atoms.
            Returns -1.0 if not exactly 2 ends were found.
        """
    def computeRadiusOfGyration(self) -> float: 
        """
        Computes the radius of gyration, :math:`R_g^2` of this molecule.

        NOTE: 
            The mass of the atoms is not yet taken into account.    
        """
    def getAtomForVertexId(self, arg0: int) -> Atom: 
        """
        Get an atom for a specific vertex.
        """
    def getAtomTypes(self) -> typing.List[int]: 
        """
        Query all types (each one for each atom) ordered by atom vertex id.
        """
    def getAtoms(self) -> typing.List[Atom]: 
        """
        Returns all atom objects enclosed in this molecule.
        """
    def getAtomsWithDegree(self, arg0: int) -> typing.List[Atom]: 
        """
        Get the atoms that have the specified number of bonds.
        """
    def getAtomsWithType(self, arg0: int) -> typing.List[Atom]: 
        """
        Get the atoms with the specified type.
        """
    def getLength(self) -> int: 
        """
        Counts and returns the number of atoms associated with this molecule.
        """
    def getNrOfAtoms(self) -> int: 
        """
        Counts and returns the number of atoms associated with this molecule.
        """
    def getNrOfBonds(self) -> int: 
        """
        Counts and returns the number of bonds associated with this molecule.
        """
    def getType(self) -> MoleculeType: 
        """
        Get the type of this molecule (see 'MoleculeType' enum).
        """
    pass
class MoleculeIterator():
    def __iter__(self) -> MoleculeIterator: ...
    def __next__(self) -> Atom: ...
    pass
class MoleculeType():
    """
    Members:

      UNDEFINED : This value indicates that either the property was not set or not discovered.

      NETWORK_STRAND

      PRIMARY_LOOP

      DANGLING_CHAIN

      FREE_CHAIN
    """
    def __eq__(self, other: object) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: int) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: int) -> None: ...
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @property
    def value(self) -> int:
        """
        :type: int
        """
    DANGLING_CHAIN: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.DANGLING_CHAIN: 3>
    FREE_CHAIN: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.FREE_CHAIN: 4>
    NETWORK_STRAND: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.NETWORK_STRAND: 1>
    PRIMARY_LOOP: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.PRIMARY_LOOP: 2>
    UNDEFINED: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.UNDEFINED: 0>
    __members__: dict # value = {'UNDEFINED': <MoleculeType.UNDEFINED: 0>, 'NETWORK_STRAND': <MoleculeType.NETWORK_STRAND: 1>, 'PRIMARY_LOOP': <MoleculeType.PRIMARY_LOOP: 2>, 'DANGLING_CHAIN': <MoleculeType.DANGLING_CHAIN: 3>, 'FREE_CHAIN': <MoleculeType.FREE_CHAIN: 4>}
    pass
class Universe():
    """
    Represents a full Polymer Network structure, a collection of molecules.
    """
    def __init__(self, arg0: float, arg1: float, arg2: float) -> None: 
        """
        Instantiate this Universe (Collection of Molecules) providing the box lengths.
        """
    def addAtoms(self, arg0: int, arg1: typing.List[int], arg2: typing.List[int], arg3: typing.List[float], arg4: typing.List[float], arg5: typing.List[float], arg6: typing.List[int], arg7: typing.List[int], arg8: typing.List[int]) -> None: 
        """
        Add atoms to the Universe, vertices to the underlying graph.
        """
    def addBonds(self, arg0: int, arg1: typing.List[int], arg2: typing.List[int]) -> None: 
        """
        Add bonds to the underlying atoms, edges to the underlying graph. If the connected atoms are not found, the bonds are silently skipped.
        """
    def determineFunctionalityPerType(self) -> typing.Dict[int, int]: 
        """
        Find the maximum functionality of each atom type in the network.
        """
    def getAtom(self, arg0: int) -> Atom: 
        """
        Find an atom by its ID.
        """
    def getAtomByVertexIdx(self, arg0: int) -> Atom: 
        """
        Find an atom by the ID of the vertex of the underlying graph.
        """
    def getAtomTypes(self) -> typing.List[int]: 
        """
        Get all types (each one for each atom) ordered by atom vertex id.
        """
    def getAtoms(self) -> typing.List[Atom]: 
        """
        Get al atoms.
        """
    def getAtomsWithType(self, arg0: int) -> typing.List[Atom]: 
        """
        Find many atom by their type.
        """
    def getBonds(self) -> typing.Dict[str, typing.List[int]]: 
        """
        Get all bonds.
        """
    def getBox(self) -> Box: 
        """
        Get the underlying bounding box object.
        """
    def getChainsWithCrosslinker(self, arg0: int) -> typing.List[Molecule]: 
        """
                    Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms, without omitting the crosslinkers.
                    In turn, e.g. for a tetrafunctional crosslinker, it will be 4 times in the resulting molecules.
                    
                    **NOTE**: Crosslinkers without bonds to non-crosslinkers are not returned.
        """
    def getMasses(self) -> typing.Dict[int, float]: 
        """
        Get the mass of one atom per type
        """
    def getMolecules(self, arg0: int) -> typing.List[Molecule]: 
        """
                    Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms.
                    
                    Reduces the Universe to a list of molecules. 
                    Specify the crosslinkerType to an existing type id, 
                    then those atoms will be omitted, and this function returns chains instead.
        """
    def getNrOfAtoms(self) -> int: 
        """
        Query the number of atoms in this universe.
        """
    def getNrOfBonds(self) -> int: 
        """
        Query the number of bonds associated with this universe.
        """
    def getNrOfBondsOfAtom(self, arg0: int) -> int: 
        """
        Count the number of immediate neighbours of an atom, specified by its id.
        """
    def getNrOfBondsOfVertex(self, arg0: int) -> int: 
        """
        Count the number of immediate neighbours of an atom, specified by its vertex id.
        """
    def getTimestep(self) -> int: 
        """
        Query the timestep when this universe was captured.
        """
    def getVolume(self) -> float: 
        """
        Query the volume of the underlying bounding box.
        """
    def setBox(self, arg0: Box) -> None: 
        """
        Override the currently assigned box with the one specified.
        """
    def setBoxLengths(self, arg0: float, arg1: float, arg2: float) -> None: 
        """
        Set the box side lengths.
        """
    def setMasses(self, arg0: typing.Dict[int, float]) -> None: 
        """
        Set the mass per type of atom.
        """
    def setTimestep(self, arg0: int) -> None: 
        """
        Set the timestep when this Universe was captured.
        """
    pass
class UniverseSequence():
    """
    This class represents a sequence of Universes, with the Universe's data files only being read on request. Dump files are read at once in order to know how many timesteps/universes are available in total.
    """
    def __init__(self) -> None: ...
    def atIndex(self, arg0: int) -> Universe: 
        """
        Get the Universe at the given index (as of in the sequence given by the dump file).
        """
    def getAll(self) -> typing.List[Universe]: 
        """
        Get all universes initialized back in a list
        """
    def getLength(self) -> int: 
        """
        Get the number of universes in this sequence.
        """
    def initializeFromDataSequence(self, arg0: typing.List[str]) -> None: 
        """
        Reset and initialize the Universes from an ordered list of Lammps data (:code:`write_data`) files.
        """
    def initializeFromDumpFile(self, arg0: str, arg1: str) -> None: 
        """
        Reset and initialize the Universes from a Lammps :code:`dump` output.
        """
    def next(self) -> Universe: 
        """
        Get the Universe that's next in the sequence.
        """
    def resetIterator(self) -> None: 
        """
        Reset the internal iterator, such that a subsequent call to :code:`next()` returns the first one again.
        """
    pass
def computeStoichiometricInbalance(arg0: Universe, arg1: int, arg2: int, arg3: typing.Dict[int, int]) -> float:
    """
    Compute stoichiometric inbalance
    """
def predictGelationPoint(arg0: float, arg1: int, arg2: int) -> float:
    """
    Predict the gelation point of a Universe
    """
DANGLING_CHAIN: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.DANGLING_CHAIN: 3>
FREE_CHAIN: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.FREE_CHAIN: 4>
NETWORK_STRAND: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.NETWORK_STRAND: 1>
PRIMARY_LOOP: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.PRIMARY_LOOP: 2>
UNDEFINED: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.UNDEFINED: 0>
