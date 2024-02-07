"""
This is the C++ implementation of pylimer_tools. 

The subpackage (actually compiled code) :module:`pylimer_tools_cpp.pylimer_tools_cpp` 
should be imported automatically when you `import pylimer_tools_cpp`, 
meaning there is no need to use the double `pylimer_tools_cpp.pylimer_tools_cpp`,
a single one is sufficient.
"""
from __future__ import annotations
import pylimer_tools_cpp
import typing
from pylimer_tools_cpp.pylimer_tools_cpp import Atom
from pylimer_tools_cpp.pylimer_tools_cpp import AtomStyle
from pylimer_tools_cpp.pylimer_tools_cpp import BalanceRunMode
from pylimer_tools_cpp.pylimer_tools_cpp import Box
from pylimer_tools_cpp.pylimer_tools_cpp import ComputedDoubleValues
from pylimer_tools_cpp.pylimer_tools_cpp import ComputedIntValues
from pylimer_tools_cpp.pylimer_tools_cpp import DPDSimulator
from pylimer_tools_cpp.pylimer_tools_cpp import DataFileReader
from pylimer_tools_cpp.pylimer_tools_cpp import DataFileWriter
from pylimer_tools_cpp.pylimer_tools_cpp import DumpFileReader
from pylimer_tools_cpp.pylimer_tools_cpp import ExitReason
from pylimer_tools_cpp.pylimer_tools_cpp import LazyUniverseSequenceIterator
from pylimer_tools_cpp.pylimer_tools_cpp import LinkSwappingMode
from pylimer_tools_cpp.pylimer_tools_cpp import MCUniverseGenerator
from pylimer_tools_cpp.pylimer_tools_cpp import MEHPForceBalance
from pylimer_tools_cpp.pylimer_tools_cpp import MEHPForceEvaluator
from pylimer_tools_cpp.pylimer_tools_cpp import MEHPForceRelaxation
from pylimer_tools_cpp.pylimer_tools_cpp import Molecule
from pylimer_tools_cpp.pylimer_tools_cpp import MoleculeIterator
from pylimer_tools_cpp.pylimer_tools_cpp import MoleculeType
from pylimer_tools_cpp.pylimer_tools_cpp import NeighbourList
from pylimer_tools_cpp.pylimer_tools_cpp import NonGaussianSpringForceEvaluator
from pylimer_tools_cpp.pylimer_tools_cpp import OutputConfiguration
from pylimer_tools_cpp.pylimer_tools_cpp import SimpleSpringMEHPForceEvaluator
from pylimer_tools_cpp.pylimer_tools_cpp import SimplifiedBalanceNetwork
from pylimer_tools_cpp.pylimer_tools_cpp import SimplifiedNetwork
from pylimer_tools_cpp.pylimer_tools_cpp import StructureSimplificationMode
from pylimer_tools_cpp.pylimer_tools_cpp import Universe
from pylimer_tools_cpp.pylimer_tools_cpp import UniverseSequence
import pylimer_tools_cpp.pylimer_tools_cpp

__all__ = [
    "Atom",
    "AtomStyle",
    "BalanceRunMode",
    "Box",
    "ComputedDoubleValues",
    "ComputedIntValues",
    "DANGLING_CHAIN",
    "DPDSimulator",
    "DataFileReader",
    "DataFileWriter",
    "DumpFileReader",
    "ExitReason",
    "FREE_CHAIN",
    "LazyUniverseSequenceIterator",
    "LinkSwappingMode",
    "MCUniverseGenerator",
    "MEHPForceBalance",
    "MEHPForceEvaluator",
    "MEHPForceRelaxation",
    "Molecule",
    "MoleculeIterator",
    "MoleculeType",
    "NETWORK_STRAND",
    "NeighbourList",
    "NonGaussianSpringForceEvaluator",
    "OutputConfiguration",
    "PRIMARY_LOOP",
    "SimpleSpringMEHPForceEvaluator",
    "SimplifiedBalanceNetwork",
    "SimplifiedNetwork",
    "StructureSimplificationMode",
    "UNDEFINED",
    "Universe",
    "UniverseSequence",
    "computeStoichiometricInbalance",
    "doRandomWalkChainFromTo",
    "inverse_langevin",
    "predictGelationPoint",
    "pylimer_tools_cpp",
    "splitCSV"
]


def computeStoichiometricInbalance(arg0: pylimer_tools_cpp.Universe, arg1: int, arg2: int, arg3: typing.Dict[int, int]) -> float:
    """
    Compute stoichiometric inbalance
    """
def doRandomWalkChainFromTo(box: pylimer_tools_cpp.Box, from_coordinates: Annotated[typing.List[float], FixedSize(3)], to_coordinates: Annotated[typing.List[float], FixedSize(3)], chainLen: int, beadDistance: float = 1.0, seed: str = '') -> typing.Dict[str, typing.List[float]]:
    """
    Do a random walk from one point to another.
    """
def inverse_langevin(x: float) -> float:
    """
    A somewhat accurate (for :math:`x \in (-1, 1)`) implementation of the inverse Langevin.

    Source: https://scicomp.stackexchange.com/a/30251
    """
def predictGelationPoint(arg0: float, arg1: int, arg2: int) -> float:
    """
    Predict the gelation point of a Universe
    """
def splitCSV(arg0: str, arg1: str) -> typing.List[str]:
    """
    Read a file containing a number of CSVs. Returns them split up.
    """
DANGLING_CHAIN: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.DANGLING_CHAIN: 3>
FREE_CHAIN: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.FREE_CHAIN: 4>
NETWORK_STRAND: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.NETWORK_STRAND: 1>
PRIMARY_LOOP: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.PRIMARY_LOOP: 2>
UNDEFINED: pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType # value = <MoleculeType.UNDEFINED: 0>
