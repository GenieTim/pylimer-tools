"""
This is the C++ implementation of pylimer_tools. 
The subpackage (actually compiled code) pylimer_tools_cpp.pylimer_tools_cpp 
should be imported automatically when you import pylimer_tools_cpp.
"""
import pylimer_tools_cpp
import typing
from pylimer_tools_cpp.pylimer_tools_cpp import Atom
from pylimer_tools_cpp.pylimer_tools_cpp import Box
from pylimer_tools_cpp.pylimer_tools_cpp import DataFileReader
from pylimer_tools_cpp.pylimer_tools_cpp import DumpFileReader
from pylimer_tools_cpp.pylimer_tools_cpp import Molecule
from pylimer_tools_cpp.pylimer_tools_cpp import MoleculeIterator
from pylimer_tools_cpp.pylimer_tools_cpp import MoleculeType
from pylimer_tools_cpp.pylimer_tools_cpp import Universe
from pylimer_tools_cpp.pylimer_tools_cpp import UniverseSequence
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
    "predictGelationPoint",
    "pylimer_tools_cpp"
]


def computeStoichiometricInbalance(arg0: pylimer_tools_cpp.Universe, arg1: int, arg2: int, arg3: typing.Dict[int, int]) -> float:
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
