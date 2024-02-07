"""
    PylimerTools Cpp
    -----------------

    A collection of utility python functions for handling LAMMPS output and polymers in Python.

    .. autosummary::
        :toctree: _generate

    """
from __future__ import annotations
import pylimer_tools_cpp.pylimer_tools_cpp
import typing
import numpy
_Shape = typing.Tuple[int, ...]

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
    "splitCSV"
]


class Atom():
    """
    A single bead or atom
    """
    def __eq__(self, arg0: Atom) -> bool: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, id: int, type: int, x: float, y: float, z: float, nx: int, ny: int, nz: int) -> None: 
        """
        Construct this atom
        """
    def __setstate__(self, arg0: tuple) -> None: 
        """
        Provides support for pickling
        """
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
    __hash__ = None
    pass
class AtomStyle():
    """
    Members:

      NONE

      ANGLE

      ATOMIC

      BODY

      BOND

      CHARGE

      DIELECTRIC

      DIPOLE

      DPD

      EDPD

      ELECTRON

      ELLIPSOID

      FULL

      LINE

      MDPD

      MOLECULAR

      PERI

      SMD

      SPH

      SPHERE

      BPM_SPHERE

      SPIN

      TDPD

      TEMPLATE

      TRI

      WAVEPACKET

      HYBRID
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
    ANGLE: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.ANGLE: 1>
    ATOMIC: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.ATOMIC: 2>
    BODY: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.BODY: 3>
    BOND: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.BOND: 4>
    BPM_SPHERE: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.BPM_SPHERE: 20>
    CHARGE: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.CHARGE: 5>
    DIELECTRIC: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.DIELECTRIC: 6>
    DIPOLE: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.DIPOLE: 7>
    DPD: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.DPD: 8>
    EDPD: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.EDPD: 9>
    ELECTRON: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.ELECTRON: 10>
    ELLIPSOID: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.ELLIPSOID: 11>
    FULL: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.FULL: 12>
    HYBRID: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.HYBRID: 26>
    LINE: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.LINE: 13>
    MDPD: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.MDPD: 14>
    MOLECULAR: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.MOLECULAR: 15>
    NONE: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.NONE: 0>
    PERI: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.PERI: 16>
    SMD: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.SMD: 17>
    SPH: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.SPH: 18>
    SPHERE: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.SPHERE: 19>
    SPIN: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.SPIN: 21>
    TDPD: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.TDPD: 22>
    TEMPLATE: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.TEMPLATE: 23>
    TRI: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.TRI: 24>
    WAVEPACKET: pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle # value = <AtomStyle.WAVEPACKET: 25>
    __members__: dict # value = {'NONE': <AtomStyle.NONE: 0>, 'ANGLE': <AtomStyle.ANGLE: 1>, 'ATOMIC': <AtomStyle.ATOMIC: 2>, 'BODY': <AtomStyle.BODY: 3>, 'BOND': <AtomStyle.BOND: 4>, 'CHARGE': <AtomStyle.CHARGE: 5>, 'DIELECTRIC': <AtomStyle.DIELECTRIC: 6>, 'DIPOLE': <AtomStyle.DIPOLE: 7>, 'DPD': <AtomStyle.DPD: 8>, 'EDPD': <AtomStyle.EDPD: 9>, 'ELECTRON': <AtomStyle.ELECTRON: 10>, 'ELLIPSOID': <AtomStyle.ELLIPSOID: 11>, 'FULL': <AtomStyle.FULL: 12>, 'LINE': <AtomStyle.LINE: 13>, 'MDPD': <AtomStyle.MDPD: 14>, 'MOLECULAR': <AtomStyle.MOLECULAR: 15>, 'PERI': <AtomStyle.PERI: 16>, 'SMD': <AtomStyle.SMD: 17>, 'SPH': <AtomStyle.SPH: 18>, 'SPHERE': <AtomStyle.SPHERE: 19>, 'BPM_SPHERE': <AtomStyle.BPM_SPHERE: 20>, 'SPIN': <AtomStyle.SPIN: 21>, 'TDPD': <AtomStyle.TDPD: 22>, 'TEMPLATE': <AtomStyle.TEMPLATE: 23>, 'TRI': <AtomStyle.TRI: 24>, 'WAVEPACKET': <AtomStyle.WAVEPACKET: 25>, 'HYBRID': <AtomStyle.HYBRID: 26>}
    pass
class BalanceRunMode():
    """
    Members:

      EIGEN_ALL

      EIGEN_RANDOM

      EIGEN_HEURISTIC

      EIGEN_STRANDS

      ITERATIVE
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
    EIGEN_ALL: pylimer_tools_cpp.pylimer_tools_cpp.BalanceRunMode # value = <BalanceRunMode.EIGEN_ALL: 3>
    EIGEN_HEURISTIC: pylimer_tools_cpp.pylimer_tools_cpp.BalanceRunMode # value = <BalanceRunMode.EIGEN_HEURISTIC: 1>
    EIGEN_RANDOM: pylimer_tools_cpp.pylimer_tools_cpp.BalanceRunMode # value = <BalanceRunMode.EIGEN_RANDOM: 0>
    EIGEN_STRANDS: pylimer_tools_cpp.pylimer_tools_cpp.BalanceRunMode # value = <BalanceRunMode.EIGEN_STRANDS: 2>
    ITERATIVE: pylimer_tools_cpp.pylimer_tools_cpp.BalanceRunMode # value = <BalanceRunMode.ITERATIVE: 4>
    __members__: dict # value = {'EIGEN_ALL': <BalanceRunMode.EIGEN_ALL: 3>, 'EIGEN_RANDOM': <BalanceRunMode.EIGEN_RANDOM: 0>, 'EIGEN_HEURISTIC': <BalanceRunMode.EIGEN_HEURISTIC: 1>, 'EIGEN_STRANDS': <BalanceRunMode.EIGEN_STRANDS: 2>, 'ITERATIVE': <BalanceRunMode.ITERATIVE: 4>}
    pass
class Box():
    """
    The box that the simulation is run in.

    NOTE: 
      currently, only rectangular boxes are supported.
    """
    def __getstate__(self) -> tuple: ...
    def __init__(self, arg0: float, arg1: float, arg2: float) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: 
        """
        Provides support for pickling.
        """
    def applySimpleShear(self, shearMagnitude: float, shearDirection: int = 0) -> None: 
        """
        Apply a simple shear to the box.

        CAUTION:
          currently, this is not supported for all operations.

        For shear magnitude, you specify the angle :math:`\gamma`.

        For the shearDirection parameter, you can specify `0` for x, `1` for y and `2` for z, respectively.
        Specify another integer to disable the shear.
        """
    def getHighX(self) -> float: ...
    def getHighY(self) -> float: ...
    def getHighZ(self) -> float: ...
    def getLowX(self) -> float: ...
    def getLowY(self) -> float: ...
    def getLowZ(self) -> float: ...
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
class ComputedDoubleValues():
    """
    Members:

      TIMESTEP

      TIME

      VOLUME

      PRESSURE

      TEMPERATURE

      STRESS_XX

      STRESS_YY

      STRESS_ZZ

      STRESS_XY

      STRESS_XZ

      STRESS_YZ

      STRESS_NXY

      STRESS_NXZ

      STRESS_NYZ

      MEAN_B

      MAX_B

      MSD
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
    MAX_B: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.MAX_B: 15>
    MEAN_B: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.MEAN_B: 14>
    MSD: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.MSD: 16>
    PRESSURE: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.PRESSURE: 3>
    STRESS_NXY: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.STRESS_NXY: 11>
    STRESS_NXZ: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.STRESS_NXZ: 13>
    STRESS_NYZ: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.STRESS_NYZ: 12>
    STRESS_XX: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.STRESS_XX: 5>
    STRESS_XY: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.STRESS_XY: 8>
    STRESS_XZ: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.STRESS_XZ: 10>
    STRESS_YY: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.STRESS_YY: 6>
    STRESS_YZ: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.STRESS_YZ: 9>
    STRESS_ZZ: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.STRESS_ZZ: 7>
    TEMPERATURE: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.TEMPERATURE: 4>
    TIME: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.TIME: 1>
    TIMESTEP: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.TIMESTEP: 0>
    VOLUME: pylimer_tools_cpp.pylimer_tools_cpp.ComputedDoubleValues # value = <ComputedDoubleValues.VOLUME: 2>
    __members__: dict # value = {'TIMESTEP': <ComputedDoubleValues.TIMESTEP: 0>, 'TIME': <ComputedDoubleValues.TIME: 1>, 'VOLUME': <ComputedDoubleValues.VOLUME: 2>, 'PRESSURE': <ComputedDoubleValues.PRESSURE: 3>, 'TEMPERATURE': <ComputedDoubleValues.TEMPERATURE: 4>, 'STRESS_XX': <ComputedDoubleValues.STRESS_XX: 5>, 'STRESS_YY': <ComputedDoubleValues.STRESS_YY: 6>, 'STRESS_ZZ': <ComputedDoubleValues.STRESS_ZZ: 7>, 'STRESS_XY': <ComputedDoubleValues.STRESS_XY: 8>, 'STRESS_XZ': <ComputedDoubleValues.STRESS_XZ: 10>, 'STRESS_YZ': <ComputedDoubleValues.STRESS_YZ: 9>, 'STRESS_NXY': <ComputedDoubleValues.STRESS_NXY: 11>, 'STRESS_NXZ': <ComputedDoubleValues.STRESS_NXZ: 13>, 'STRESS_NYZ': <ComputedDoubleValues.STRESS_NYZ: 12>, 'MEAN_B': <ComputedDoubleValues.MEAN_B: 14>, 'MAX_B': <ComputedDoubleValues.MAX_B: 15>, 'MSD': <ComputedDoubleValues.MSD: 16>}
    pass
class ComputedIntValues():
    """
    Members:

      STEP

      NUM_SHIFT

      NUM_RELOC

      NUM_ATOMS

      NUM_EXTRA_ATOMS

      NUM_BONDS

      NUM_EXTRA_BONDS

      NUM_BONDS_TO_FORM
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
    NUM_ATOMS: pylimer_tools_cpp.pylimer_tools_cpp.ComputedIntValues # value = <ComputedIntValues.NUM_ATOMS: 3>
    NUM_BONDS: pylimer_tools_cpp.pylimer_tools_cpp.ComputedIntValues # value = <ComputedIntValues.NUM_BONDS: 5>
    NUM_BONDS_TO_FORM: pylimer_tools_cpp.pylimer_tools_cpp.ComputedIntValues # value = <ComputedIntValues.NUM_BONDS_TO_FORM: 7>
    NUM_EXTRA_ATOMS: pylimer_tools_cpp.pylimer_tools_cpp.ComputedIntValues # value = <ComputedIntValues.NUM_EXTRA_ATOMS: 4>
    NUM_EXTRA_BONDS: pylimer_tools_cpp.pylimer_tools_cpp.ComputedIntValues # value = <ComputedIntValues.NUM_EXTRA_BONDS: 6>
    NUM_RELOC: pylimer_tools_cpp.pylimer_tools_cpp.ComputedIntValues # value = <ComputedIntValues.NUM_RELOC: 2>
    NUM_SHIFT: pylimer_tools_cpp.pylimer_tools_cpp.ComputedIntValues # value = <ComputedIntValues.NUM_SHIFT: 1>
    STEP: pylimer_tools_cpp.pylimer_tools_cpp.ComputedIntValues # value = <ComputedIntValues.STEP: 0>
    __members__: dict # value = {'STEP': <ComputedIntValues.STEP: 0>, 'NUM_SHIFT': <ComputedIntValues.NUM_SHIFT: 1>, 'NUM_RELOC': <ComputedIntValues.NUM_RELOC: 2>, 'NUM_ATOMS': <ComputedIntValues.NUM_ATOMS: 3>, 'NUM_EXTRA_ATOMS': <ComputedIntValues.NUM_EXTRA_ATOMS: 4>, 'NUM_BONDS': <ComputedIntValues.NUM_BONDS: 5>, 'NUM_EXTRA_BONDS': <ComputedIntValues.NUM_EXTRA_BONDS: 6>, 'NUM_BONDS_TO_FORM': <ComputedIntValues.NUM_BONDS_TO_FORM: 7>}
    pass
class DPDSimulator():
    """
    A quick-and-dirty implementation of the DPD simulation
    with slip-springs as presented by Langeloth et al.
    """
    def __init__(self, universe: Universe, crosslinker_type: int = 2, slipspring_bond_type: int = 9, is_2D: bool = False, seed: str = '') -> None: 
        """
        Get an instance of this class
        """
    def assumeBoxLargeEnough(self) -> None: 
        """
        Configure whether to run PBC on the bonds or not.

        If your bonds could get larger than half the box length, this must be kept false (default).
        Otherwise, you can set it to true and therewith get some securities.
        """
    def configA(self, A: float = 25.0) -> None: 
        """
        Configure the force-field (pair-style) parameter `A`.
        """
    def configAllowRelocationInNetwork(self, allow_relocation_in_network: bool = False) -> None: 
        """
        Configure whether a relocation step may happen when a slip-spring has ended at a cross-link.

        Side-effect: if true, the relocations may also happen *to* a slip-spring next to a cross-link.

        Arguments:
        - allow_relocation_in_network (bool): Whether to allow relocation in the network or not.
        """
    def configAutoCorrelatorOutput(self, values: typing.List[OutputConfiguration], numcorrin: int = 32, p: int = 16, m: int = 2) -> None: 
        """
        Set which values to compute multiple-tau autocorrelation for.
        If you use this, you should cite `doi:10.1063/1.3491098 <https://pubs.aip.org/aip/jcp/article-abstract/133/15/154103/190247/Efficient-on-the-fly-calculation-of-time?redirectedFrom=fulltext>`_

        Arguments:
             - values: a list of OutputConfiguration structs
             - p
        """
    def configAverageOutput(self, arg0: typing.List[OutputConfiguration]) -> None: 
        """
        Set which values to compute averages for.

        Arguments:
             - values: a list of OutputConfiguration structs
        """
    def configBondFormation(self, num_bonds_to_form: int, max_bonds_per_atom_type: typing.Dict[int, int], bond_formation_dist: float = 1.0, attempt_bond_formation_every: int = 50, atom_type_form_from: int = 2, atom_type_form_to: int = 1) -> None: 
        """
        Configure how to do bond formation during the run.

        Arguments:
        - num_bonds_to_form (int): the nr of bonds to form in total. Use 0 to stop bond formation.
        - num_bonds_per_atom_type (dict): the nr of bonds each atom type may have at most (e.g., 2 for strand atoms, 4 for a tertiary cross-links)
        - bond_formation_dist (float): the maximum distance allowed to form bonds
        - attempt_bond_formation_every (int): attempt to form bonds every this many steps during the simulation run
        - atom_type_form_from (int): the atom type to start forming bonds from. 
        - atom_type_form_to (int): the atom type to start forming bonds to.
        """
    def configLambda(self, l: float = 0.65) -> None: 
        """
        Configure the modified velocity verlet integration parameter `\lambda`.
        """
    def configNumStepsDPD(self, arg0: int) -> None: 
        """
        Configure the number of steps to do in one DPD sequence.
        """
    def configNumStepsMC(self, arg0: int) -> None: 
        """
        Configure the number of steps to do in one MC sequence.
        """
    def configRestartOutput(self, file: str, outputEvery: int = 50000) -> None: 
        """
        Set when to output a restart where.

        Note:
             The filename determines the type of serialisation: 
             .json, .xml are supported; other file endings will lead to binary serialisation (fastest!).

        Caution:
             This method may not be backwards- nor forward-compatible.
             Use the same version of pylimer-tools if you want to be sure that things work.

        Arguments:
             - file: the file path to the restart file to write
             - outputEvery: how often to write the restart file
        """
    def configShiftOneAtATime(self, arg0: bool) -> None: ...
    def configShiftPossibilityEmpty(self, arg0: bool) -> None: ...
    def configSigma(self, sigma: float = 3.0) -> None: 
        """
        Configure the force-field (pair-style) parameter `\sigma`.
        """
    def configSlipspringHighCutoff(self, cutoff: float = 2.0) -> None: 
        """
        Configure the lower cut-off of how far a pair may be distanced for a slip-spring to be created.
        """
    def configSlipspringLowCutoff(self, cutoff: float = 0.5) -> None: 
        """
        Configure the higher cut-off of how far a pair may be distanced for a slip-spring to be created.
        """
    def configSpringConstant(self, k: float = 2.0) -> None: 
        """
        Configure the force-field (bond-style) parameter `k`, the spring constant.
        """
    def configStepOutput(self, arg0: typing.List[OutputConfiguration]) -> None: 
        """
        Set which values to log.

        Arguments:
             - values: a list of OutputConfiguration structs
        """
    def createSlipSprings(self, num: int, bond_type: int = 9) -> int: 
        """
        Randomly add the specified number of slip-springs to neighbours within the specified cut-offs.
        """
    def getBondLengths(self) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]: ...
    def getCoordinates(self) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]: ...
    def getCurrentTimestep(self) -> int: ...
    def getNrOfBondsToForm(self) -> int: 
        """
        Get the number of bonds that are configured to have to be formed.
        """
    def getNumSlipSprings(self) -> int: ...
    def getNumStepsDPD(self) -> int: ...
    def getNumStepsMC(self) -> int: ...
    def getShiftOneAtATime(self) -> bool: ...
    def getShiftPossibilityEmpty(self) -> bool: ...
    def getSlipSpringBondType(self) -> int: ...
    def getSpringConstant(self) -> float: ...
    def getStressTensor(self) -> numpy.ndarray[numpy.float64, _Shape[3, 3]]: ...
    def getTemperature(self) -> float: ...
    def getTimestep(self) -> float: ...
    def getUniverse(self, with_slipsprings: bool = True) -> Universe: 
        """
        Get a universe instance from the current coordinates (and connectivity).

        Arguments:
             - with_slip_springs (bool): whether to include slip-springs in the returned universe.
        """
    @staticmethod
    def readRestartFile(file: str) -> DPDSimulator: 
        """
        Read a restart file in order to continue a simulation.
        """
    def refreshCurrentState(self) -> None: 
        """
        After re-configuring the force-field parameters, 
        this method should be called to update the current stress tensor etc.
        """
    def runSimulation(self, n_steps: int, dt: float = 0.06, with_MC: bool = False) -> None: ...
    def startMeasuringMSDForAtoms(self, atom_ids: typing.List[int]) -> None: 
        """
        Set a new origin for measuing the mean square displacement for a specified set of atoms
        """
    def validateNeighbourList(self, arg0: float) -> None: ...
    def validateState(self) -> None: ...
    def writeRestartFile(self, file: str) -> None: 
        """
        Explicitily force the writing of a restart file, now!

        Arguments:
        - file (str): the file path and name of the restart file to be written.
             Can end in xml, json or anything else (-> binary).
        """
    pass
class DataFileReader():
    """
    A reader for LAMMPS's `write_data` files.
    """
    def __init__(self) -> None: ...
    def getAtomIds(self) -> typing.List[int]: ...
    def getAtomNx(self) -> typing.List[int]: ...
    def getAtomNy(self) -> typing.List[int]: ...
    def getAtomNz(self) -> typing.List[int]: ...
    def getAtomTypes(self) -> typing.List[int]: ...
    def getAtomX(self) -> typing.List[float]: ...
    def getAtomY(self) -> typing.List[float]: ...
    def getAtomZ(self) -> typing.List[float]: ...
    def getBondFrom(self) -> typing.List[int]: ...
    def getBondTo(self) -> typing.List[int]: ...
    def getBondTypes(self) -> typing.List[int]: ...
    def getHighX(self) -> float: ...
    def getHighY(self) -> float: ...
    def getHighZ(self) -> float: ...
    def getLowX(self) -> float: ...
    def getLowY(self) -> float: ...
    def getLowZ(self) -> float: ...
    def getLx(self) -> float: ...
    def getLy(self) -> float: ...
    def getLz(self) -> float: ...
    def getMasses(self) -> typing.Dict[int, float]: ...
    def getMoleculeIds(self) -> typing.List[int]: ...
    def getNrOfAtomTypes(self) -> int: ...
    def getNrOfAtoms(self) -> int: ...
    def getNrOfBondTypes(self) -> int: ...
    def getNrOfBonds(self) -> int: ...
    def read(self, path_of_file_to_read: str, atom_style: AtomStyle = AtomStyle.ANGLE, atom_style2: AtomStyle = AtomStyle.NONE, atom_style_3: AtomStyle = AtomStyle.NONE) -> None: 
        """
        Actually read a LAMMPS's `write_data` file.

        Arguments:
           - `path_of_file_to_read`: The path to the file to read
           - `atom_style`: The format of the "Atoms" section, see https://docs.lammps.org/read_data.html
           - `atom_style2`: The format of the "Atoms" section if the previous parameter is equal to AtomStyle::HYBRID
           - `atom_style3`: The format of the "Atoms" section if the second to last parameter is equal to AtomStyle::HYBRID
        """
    pass
class DataFileWriter():
    def __init__(self, universe: Universe) -> None: 
        """
        Initialize the writer with the universe to write.
        """
    def configAttemptImageReset(self, attemptImageReset: bool = True) -> None: 
        """
        Set whether to change the outuput coordinates to lie in the box or not.

        Default: false.
        """
    def configCrosslinkerType(self, crosslinkerType: int = 2) -> None: 
        """
        Set which atom type represents cross-linkers. 
        Needed in case the moleculeIdx in the output file should have any meaning.
        (e.g. with :func:`~pylimer_tools_cpp.pylimer_tools_cpp.DataFileWriter.configMoleculeIdxForSwap`).

        Default: 2.
        """
    def configIncludeAngles(self, includeAngles: bool = True) -> None: 
        """
        Set whether to include the angles from the universe in the file or not.

        Default: true.
        """
    def configIncludeDihedralAngles(self, includeDihedralAngles: bool = True) -> None: 
        """
        Set whether to include the dihedral angles from the universe in the file or not.

        Default: true.
        """
    def configMoleculeIdxForSwap(self, enableSwappability: bool = True) -> None: 
        """
        Swappable chains implies that their `moleculeIdx` in the LAMMPS data file is not 
        identical per chain, but identical per position in the chain.
        That's how you can have bond swapping with constant chain length distribution.

        Default: false.
        """
    def configMoveIntoBox(self, moveIntoBox: bool = True) -> None: 
        """
        Set whether to change the outuput coordinates to lie in the box or not.

        Default: false (used to be true).
        """
    def configReindexAtoms(self, reindexAtoms: bool = True) -> None: 
        """
        Set whether to reindex the atoms or not. 
        Re-indexing leads to atom ids being in the range of 1 to the number of atoms.

        Default: false.
        """
    def setCustomAtomFormat(self, atomFormat: str = '\t$atomId\t$moleculeId\t$atomType\t$x\t$y\t$z\t$nx\t$ny\t$nz') -> None: 
        """
         Specify a custom format for the atom section.
         Placeholder options are:
             - $atomId
             - $moleculeId
             - $atomType
             - $x
             - $y
             - $z
             - $nx
             - $ny
             - $nz

        Additionally, you can use the keys used in 
        :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.setPropertyValue`
        as placeholders (as long as they are alphanumeric only; prefix in the format with '$' as well).
        Specifically useful if you need a different (or hybrid) atom style in LAMMPS.
        """
    def setUniverseToWrite(self, universe: Universe) -> None: 
        """
        Re-set the universe to write.
        """
    def writeToFile(self, arg0: str) -> None: ...
    pass
class DumpFileReader():
    """
    A reader for LAMMPS's `dump` files.
    """
    def __init__(self, pathOfFileToRead: str) -> None: ...
    def getLength(self) -> int: 
        """
        Get the number of sections (time-steps) in the file
        """
    def getNumericValuesForAt(self, arg0: int, arg1: str, arg2: str) -> typing.List[float]: 
        """
        Get the values for the section `index`, the main header `headerKey` and the column (in the header) `column`.
        """
    def getStringValuesForAt(self, rowIndex: int, headerKey: str, columnIndex: str) -> typing.List[str]: 
        """
        Get the values for the section `index`, the main header `headerKey` and the column (in the header) `column`.
        """
    def hasKey(self, headerKey: str) -> bool: 
        """
        Check whether the first section has the header specified
        """
    def keyHasColumn(self, headerKey: str, columnName: str) -> bool: 
        """
        Check whether the header of the first section has the specified column
        """
    def keyHasDirectionalColumn(self, headerKey: str, dirPraefix: str = '', dirSuffix: str = '') -> bool: 
        """
        Check whether the header of the first section has all the three columns `{dirPraefix}{x|y|z}{dirSuffix}`.
        """
    def read(self) -> None: 
        """
        Read the whole file
        """
    pass
class ExitReason():
    """
    Members:

      UNSET

      MAX_STEPS

      F_TOLERANCE

      X_TOLERANCE

      FAILURE

      OTHER
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
    FAILURE: pylimer_tools_cpp.pylimer_tools_cpp.ExitReason # value = <ExitReason.FAILURE: 4>
    F_TOLERANCE: pylimer_tools_cpp.pylimer_tools_cpp.ExitReason # value = <ExitReason.F_TOLERANCE: 1>
    MAX_STEPS: pylimer_tools_cpp.pylimer_tools_cpp.ExitReason # value = <ExitReason.MAX_STEPS: 3>
    OTHER: pylimer_tools_cpp.pylimer_tools_cpp.ExitReason # value = <ExitReason.OTHER: 5>
    UNSET: pylimer_tools_cpp.pylimer_tools_cpp.ExitReason # value = <ExitReason.UNSET: 0>
    X_TOLERANCE: pylimer_tools_cpp.pylimer_tools_cpp.ExitReason # value = <ExitReason.X_TOLERANCE: 2>
    __members__: dict # value = {'UNSET': <ExitReason.UNSET: 0>, 'MAX_STEPS': <ExitReason.MAX_STEPS: 3>, 'F_TOLERANCE': <ExitReason.F_TOLERANCE: 1>, 'X_TOLERANCE': <ExitReason.X_TOLERANCE: 2>, 'FAILURE': <ExitReason.FAILURE: 4>, 'OTHER': <ExitReason.OTHER: 5>}
    pass
class LazyUniverseSequenceIterator():
    """
    An iterator to iterate throught the universes in :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.UniverseSequence`.
    """
    def __iter__(self) -> LazyUniverseSequenceIterator: ...
    def __next__(self) -> Universe: ...
    pass
class LinkSwappingMode():
    """
    Members:

      NO_SWAPPING

      SLIPLINKS_ONLY

      ALL

      ALL_CYCLE

      ALL_MC

      ALL_MC_CYCLE
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
    ALL: pylimer_tools_cpp.pylimer_tools_cpp.LinkSwappingMode # value = <LinkSwappingMode.ALL: 2>
    ALL_CYCLE: pylimer_tools_cpp.pylimer_tools_cpp.LinkSwappingMode # value = <LinkSwappingMode.ALL_CYCLE: 3>
    ALL_MC: pylimer_tools_cpp.pylimer_tools_cpp.LinkSwappingMode # value = <LinkSwappingMode.ALL_MC: 4>
    ALL_MC_CYCLE: pylimer_tools_cpp.pylimer_tools_cpp.LinkSwappingMode # value = <LinkSwappingMode.ALL_MC_CYCLE: 5>
    NO_SWAPPING: pylimer_tools_cpp.pylimer_tools_cpp.LinkSwappingMode # value = <LinkSwappingMode.NO_SWAPPING: 0>
    SLIPLINKS_ONLY: pylimer_tools_cpp.pylimer_tools_cpp.LinkSwappingMode # value = <LinkSwappingMode.SLIPLINKS_ONLY: 1>
    __members__: dict # value = {'NO_SWAPPING': <LinkSwappingMode.NO_SWAPPING: 0>, 'SLIPLINKS_ONLY': <LinkSwappingMode.SLIPLINKS_ONLY: 1>, 'ALL': <LinkSwappingMode.ALL: 2>, 'ALL_CYCLE': <LinkSwappingMode.ALL_CYCLE: 3>, 'ALL_MC': <LinkSwappingMode.ALL_MC: 4>, 'ALL_MC_CYCLE': <LinkSwappingMode.ALL_MC_CYCLE: 5>}
    pass
class MCUniverseGenerator():
    """
    A :obj:`pylimer_tools_cpp.pylimer_tools_cpp.Universe` generator using a Monte-Carlo procedure.
    """
    def __init__(self, Lx: float, Ly: float, Lz: float) -> None: ...
    def addAndLinkStrands(self, nrOfStrands: int, strandLengths: typing.List[int], crosslinkerConversion: float, crosslinkerFunctionality: int, strandAtomType: int = 1) -> None: 
        """
        Actually add strands, link them to the previously added cross-linkers.
        """
    def addCrosslinkers(self, nrOfCrosslinkers: int, crosslinkerAtomType: int = 2) -> None: 
        """
        Add the cross-linkers.
        """
    def addSolventChains(self, nrOfSolventChains: int, solventChainLength: int, solventAtomType: int = 3) -> None: 
        """
        Randomly distribute additional, free chains.
        """
    def getUniverse(self) -> Universe: 
        """
        Fetch the current (or final) state of the universe.
        """
    def setBeadDistance(self, distance: float) -> None: 
        """
        Set the optimal distance between beads.
        """
    def setSeed(self, seed: int) -> None: 
        """
        Set the seed for the random generator.
        """
    pass
class MEHPForceBalance():
    """
    A small simulation tool for quickly minimizing the force between the cross-linker beads.
     
    """
    def __copy__(self) -> MEHPForceBalance: ...
    def __init__(self, universe: Universe, crosslinkerType: int = 2, is2D: bool = False, kappa: float = 1.0, remove2functionalCrosslinkers: bool = True) -> None: 
        """
        Instantiate the simulator for a certain universe.

        :param universe: the universe to simulate with
        :param crosslinkerType: The atom type of the cross-linkers. Needed to reduce the network.
        :param is2D: Whether to ignore the z direction.
        :param forceEvaluator: The force evaluator to use
        """
    def addSlipLinks(self, strandIdx1: typing.List[int], strandIdx2: typing.List[int], x: typing.List[float], y: typing.List[float], z: typing.List[float], alpha1: typing.List[float], alpha2: typing.List[float], clampAlpha: bool = False) -> None: 
        """
        Add the slip-links
        """
    def addSliplinksBasedOnCycles(self, maxLoopLength: int = -1) -> int: ...
    def deformTo(self, newBox: Box) -> None: ...
    def evaluateDistanceBetween(self, network: SimplifiedBalanceNetwork, displacements: numpy.ndarray[numpy.float64, _Shape[m, 1]], linkIndexA: int, linkIndexB: int, is2D: bool = False) -> numpy.ndarray[numpy.float64, _Shape[3, 1]]: ...
    def getAverageSpringLength(self) -> float: 
        """
        Get the average length of the springs. Note that in contrast to :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceBalance.getGammaFactor()`,
        this value is normalized by the number of springs rather than the number of chains.
        """
    def getCrosslinkerVerse(self) -> Universe: 
        """
        Returns the universe [of cross-linkers] with the positions of the current state of the simulation.
        """
    def getDanglingWeightFraction(self, tolerance: float = 0.1) -> float: 
        """
        Compute the weight fraction of non-active springs

        Caution: ignores atom masses.
        """
    def getDefaultNrOfChains(self) -> int: 
        """
        Returns the value effectively used in :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceBalance.getGammaFactor()` for normalizing the distances.`.
        """
    def getDefaultR0Square(self) -> float: 
        """
        Returns the value effectively used in :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceBalance.getGammaFactor()` for :math:`\langle R_{0,\eta}^2\rangle`.
        """
    def getDisplacementResidualNorm(self, arg0: float) -> float: 
        """
        Get the current link displacement residual norm.
        """
    def getDisplacements(self) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]: 
        """
        Get the current link displacements.
        """
    def getEffectiveFunctionalityOfAtoms(self, tolerance: float = 0.1) -> typing.Dict[int, int]: 
        """
        Returns the number of active springs connected to each atom, atomId used as index

        :param tolerance: springs under this length are considered inactive
        """
    def getExitReason(self) -> ExitReason: 
        """
        Returns the reason for termination of the simulation
        """
    def getForceOn(self, linkIdx: int, oneOverSpringPartitionUpperLimit: float = 1.0) -> numpy.ndarray[numpy.float64, _Shape[3, 3]]: ...
    def getGammaFactor(self, r0squared: float = -1.0, nrOfChains: int = -1) -> float: 
        """
        Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

        :math:`\Gamma = \langle\gamma_{\eta}\rangle`, with :math:`\gamma_{\eta} = \frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}`,
        which you can use as :math:`G_{\mathrm{ANT}} = \Gamma \nu k_B T`,
        where :math:`\eta` is the index of a particular strand, 
        :math:`R_{0}^2` is the melt mean square end to end distance, in phantom systems :math:`$= N_{\eta}*b^2$`
        :math:`N_{\eta}` is the number of atoms in this strand :math:`\eta`, 
        :math:`b` its mean square bond length,
        :math:`T` the temperature and 
        :math:`k_B` Boltzmann's constant.

        :param r0squared: The denominator in the equation of :math:`\Gamma`. If :math:`-1.0` (default), the network is used for determination (which is not accurate). For phantom systems, the correct value is :math:`Nb^2`.
             For other systems, the value could be determined by `~pylimer_tools_cpp.pylimer_tools_cpp.Universe.computeMeanEndToEndDistance` on the melt system.
        :param nrOfChains: the value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of chains. 
        """
    def getIdsOfActiveNodes(self, tolerance: float = 0.1, minimumNrOfActiveConnections: int = 2, maximumNrOfActiveConnections: int = -1, usePartial: bool = False) -> typing.List[int]: 
        """
        Get the atom ids of the nodes that are considered active.

        :param tolerance: springs under this length are considered inactive. A node is active if it has > 2 active springs.
        :param minimumNrOfActiveConnections:  A node is active if it has equal or more than this number of active springs.
        :param maximumNrOfActiveConnections:  A node is active if it has equal or less than this number of active springs.
             Use a value < 0 to indicate that there is no maximum number of active connections.
        """
    def getNeighbourLinkIndices(self, network: SimplifiedBalanceNetwork, linkIdx: int) -> typing.List[int]: ...
    def getNrOfActiveNodes(self, tolerance: float = 0.1, minimumNrOfActiveConnections: int = 2, maximumNrOfActiveConnections: int = -1, usePartial: bool = False) -> int: 
        """
         Get the number of active nodes remaining after running the simulation.

        :param tolerance: springs under this length are considered inactive.
        :param minimumNrOfActiveConnections:  A node is active if it has equal or more than this number of active springs.
        :param maximumNrOfActiveConnections:  A node is active if it has equal or less than this number of active springs.
             Use a value < 0 to indicate that there is no maximum number of active connections.
        :param usePartial: Whether to use the partial spring distances rather than the total (set to true if you want primary loop contributors)
        """
    def getNrOfActivePartialSprings(self, tolerance: float = 0.1) -> int: 
        """
         Get the number of active partial springs remaining after running the simulation.

        :param tolerance: springs under this length are considered inactive
        """
    def getNrOfActiveSprings(self, tolerance: float = 0.1) -> int: 
        """
         Get the number of active springs remaining after running the simulation.

        :param tolerance: springs under this length are considered inactive
        """
    def getNrOfIterations(self) -> int: 
        """
        Returns the number of iterations used for force relaxation.
        """
    def getNrOfNodes(self) -> int: 
        """
        Get the number of nodes considered in this simulation.
        """
    def getNrOfSprings(self) -> int: 
        """
        Get the number of springs considered in this simulation.

        :param tolerance: springs under this length are considered inactive
        """
    def getPressure(self) -> float: 
        """
        Returns the pressure at the current state of the simulation.
        """
    def getSolubleWeightFraction(self, tolerance: float = 0.1) -> float: 
        """
        Compute the weight fraction of springs connected to active
        springs (any depth). 

        Caution: ignores atom masses.
        """
    def getSpringPartitions(self) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]: 
        """
        Get the current spring partitions.
        """
    def getSpringpartitionIndicesOfSliplink(self, network: SimplifiedBalanceNetwork, linkIdx: int) -> typing.List[int]: ...
    def getStressTensor(self, oneOverSpringPartitionUpperLimit: float = 1.0) -> numpy.ndarray[numpy.float64, _Shape[3, 3]]: 
        """
        Returns the stress tensor at the current state of the simulation.
        """
    def getStressTensorLinkBased(self, oneOverSpringPartitionUpperLimit: float = 1.0, xlinksOnly: bool = False) -> numpy.ndarray[numpy.float64, _Shape[3, 3]]: 
        """
        Returns the stress tensor at the current state of the simulation.
        """
    def inspectDisplacementToMeanPositionUpdate(self, linkIdx: int, oneOverSpringPartitionUpperLimit: float = 1.0) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]: ...
    def inspectLinkDisplacementToMeanPositionUpdate(self, linkIdx: int, damping: float = 1.0) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]: ...
    def inspectParametrisationOptimsationForLink(self, linkIdx: int, displacements: numpy.ndarray[numpy.float64, _Shape[m, 1]], springPartitions: numpy.ndarray[numpy.float64, _Shape[m, 1]], maxNrOfSteps: int = 100, alpha_tol: float = 1e-09, minNrOfSteps: int = 1, oneOverSpringPartitionUpperLimit: float = 1.0) -> typing.Tuple[numpy.ndarray[numpy.float64, _Shape[m, 1]], numpy.ndarray[numpy.float64, _Shape[m, 1]], int, float, float, float, float]: ...
    def inspectSpringPartitionUpdate(self, linkIdx: int) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]: ...
    def moveSlipLinksToTheirBestBranch(self, arg0: SimplifiedBalanceNetwork, arg1: numpy.ndarray[numpy.float64, _Shape[m, 1]], arg2: numpy.ndarray[numpy.float64, _Shape[m, 1]], arg3: float, arg4: int, arg5: bool) -> None: ...
    def randomlyAddSlipLinks(self, nrOfSlipLinksToSample: int, cutoff: float = 2.0, minimumNrOfSliplinks: int = 0, sameStrandCutoff: float = 2, excludeCrosslinks: bool = False, seed: int = -1) -> int: ...
    def runForceRelaxation(self, runMode: BalanceRunMode = BalanceRunMode.ITERATIVE, damping: float = 1.0, maxNrOfSteps: int = 250000, xTolerance: float = 1e-12, initialResidualNorm: float = -1.0, simplificationMode: StructureSimplificationMode = StructureSimplificationMode.NO_SIMPLIFICATION, inactiveRemovalCutoff: float = -1.0, outputFrequency: int = 50, doInnerIterations: bool = False, allowSlipLinksToPassEachOther: LinkSwappingMode = LinkSwappingMode.NO_SWAPPING, swappingFrequency: int = 10, oneOverSpringPartitionUpperLimit: float = 1.0, nrOfCrosslinkSwapsAllowedPerSliplink: int = -1) -> None: 
        """
        Run the simulation.
        Note that the final state of the minimization is persisted and reused if you use this method again.
        This is useful if you want to run a global optimization first and add a local one afterwards.
        As a consequence though, you cannot simply benchmark only this method; you must include the setup.

        :param runMode: Choice of the mode to run the simulation with.
        :param damping: For certain run modes, a damping factor helps to improve performance.
        :param maxNrOfSteps: The maximum number of steps to do during the simulation.
        :param xTolerance: The tolerance of the displacements as an exit condition.
        :param innerMaxNrOfSteps: The maximum number of steps to do per iteration during the slip-link displacements.
        :param innerXTolerance: The tolerance of the displacements of the slip-link as an inner exit condition.
        :param innerAlphaTolerance: The tolerance of the contour-length when slipping the slip-link as an inner exit condition.
        """
    def setDisplacements(self, arg0: numpy.ndarray[numpy.float64, _Shape[m, 1]]) -> None: 
        """
        Set the current link displacements.
        """
    def setSpringContourLengths(self, arg0: numpy.ndarray[numpy.float64, _Shape[m, 1]]) -> None: 
        """
        Set/overwrite the contour lengths.
        """
    def setSpringPartitions(self, arg0: numpy.ndarray[numpy.float64, _Shape[m, 1]]) -> None: 
        """
        Set the current spring partitions.
        """
    def swapSlipLinksInclXlinks(self, arg0: SimplifiedBalanceNetwork, arg1: numpy.ndarray[numpy.float64, _Shape[m, 1]], arg2: numpy.ndarray[numpy.float64, _Shape[m, 1]], arg3: float, arg4: bool) -> None: ...
    @property
    def network(self) -> SimplifiedBalanceNetwork:
        """
        :type: SimplifiedBalanceNetwork
        """
    pass
class MEHPForceEvaluator():
    """
    The base interface to change the way the force is evaluated during a MEHP run.
    """
    def __init__(self) -> None: ...
    def evaluateStressContribution(self, springDistances: float, i: int, j: int, spring_index: int) -> float: 
        """
        An evaluation of the stress-contribution.

        :param springDistances: the three coordinate differences for one spring.
        :param i: the row index of the stress tensor
        :param j: the column index of the stress tensor
        """
    @property
    def is2D(self) -> bool:
        """
        :type: bool
        """
    @is2D.setter
    def is2D(self, arg1: bool) -> None:
        pass
    @property
    def network(self) -> SimplifiedNetwork:
        """
        :type: SimplifiedNetwork
        """
    pass
class MEHPForceRelaxation():
    """
    A small simulation tool for quickly minimizing the force between the cross-linker beads.
     
    """
    def __init__(self, universe: Universe, crosslinkerType: int = 2, is2D: bool = False, forceEvaluator: MEHPForceEvaluator = None, kappa: float = 1.0, remove2functionalCrosslinkers: bool = True, removeDanglingChains: bool = False) -> None: 
        """
        Instantiate the simulator for a certain universe.

        :param universe: the universe to simulate with
        :param crosslinkerType: The atom type of the cross-linkers. Needed to reduce the network.
        :param is2D: Whether to ignore the z direction.
        :param forceEvaluator: The force evaluator to use
        :param kappa: The spring constant
        :param remove2functionalCrosslinkers: Whether to replace two-functional cross-links with a "normal" chain bead
        :param removeDanglingChains: Whether to remove dangling chains before running the simulation. 
             **Caution*: Removing the dangling chains will result in incorrect results fo the computation of 
             :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getSolubleWeightFraction()` and
             :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getDanglingWeightFraction()`
        """
    def configRerunEpsilon(self, epsilon: float = 0.001) -> None: 
        """
        Configure the offset from the lower and upper bounds for the simulation to suggest another run (
             See: :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.requiresAnotherRun()`
        ).
        """
    def getAverageSpringLength(self) -> float: 
        """
        Get the average length of the springs. Note that in contrast to :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getGammaFactor()`,
        this value is normalized by the number of springs rather than the number of chains.
        """
    def getCrosslinkerVerse(self) -> Universe: 
        """
        Returns the universe [of cross-linkers] with the positions of the current state of the simulation.
        """
    def getDanglingWeightFraction(self, tolerance: float = 0.1) -> float: 
        """
        Compute the weight fraction of non-active springs

        Caution: ignores atom masses.
        """
    def getDefaultNrOfChains(self) -> int: 
        """
        Returns the value effectively used in :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getGammaFactor()` for normalizing the distances.`.
        """
    def getDefaultR0Square(self) -> float: 
        """
        Returns the value effectively used in :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.getGammaFactor()` for :math:`\langle R_{0,\eta}^2\rangle`.
        """
    def getEffectiveFunctionalityOfAtoms(self, tolerance: float = 0.1) -> typing.Dict[int, int]: 
        """
        Returns the number of active springs connected to each atom, atomId used as index

        :param tolerance: springs under this length are considered inactive
        """
    def getExitReason(self) -> ExitReason: 
        """
        Returns the reason for termination of the simulation
        """
    def getForce(self) -> float: 
        """
        Returns the force at the current state of the simulation.
        """
    def getGammaFactor(self, r0squared: float = -1.0, nrOfChains: int = -1) -> float: 
        """
        Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

        :math:`\Gamma = \langle\gamma_{\eta}\rangle`, with :math:`\gamma_{\eta} = \frac{\bar{r_{\eta}}^2}{R_{0,\eta}^2}`,
        which you can use as :math:`G_{\mathrm{ANT}} = \Gamma \nu k_B T`,
        where :math:`\eta` is the index of a particular strand, 
        :math:`R_{0}^2` is the melt mean square end to end distance, in phantom systems :math:`$= N_{\eta}*b^2$`
        :math:`N_{\eta}` is the number of atoms in this strand :math:`\eta`, 
        :math:`b` its mean square bond length,
        :math:`T` the temperature and 
        :math:`k_B` Boltzmann's constant.

        :param r0squared: The denominator in the equation of :math:`\Gamma`. If :math:`-1.0` (default), the network is used for determination (which is not accurate). For phantom systems, the correct value is :math:`Nb^2`.
             For other systems, the value could be determined by `~pylimer_tools_cpp.pylimer_tools_cpp.Universe.computeMeanEndToEndDistance` on the melt system.
        :param nrOfChains: the value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of chains. 
        """
    def getIdsOfActiveNodes(self, tolerance: float = 0.1, minimumNrOfActiveConnections: int = 2, maximumNrOfActiveConnections: int = -1) -> typing.List[int]: 
        """
        Get the atom ids of the nodes that are considered active.

        :param tolerance: springs under this length are considered inactive. A node is active if it has > 2 active springs.
        :param minimumNrOfActiveConnections:  A node is active if it has equal or more than this number of active springs.
        :param maximumNrOfActiveConnections:  A node is active if it has equal or less than this number of active springs.
             Use a value < 0 to indicate that there is no maximum number of active connections.
        """
    def getNrOfActiveNodes(self, tolerance: float = 0.1, minimumNrOfActiveConnections: int = 2, maximumNrOfActiveConnections: int = -1) -> int: 
        """
         Get the number of active nodes remaining after running the simulation.

        :param tolerance: springs under this length are considered inactive.
        :param minimumNrOfActiveConnections:  A node is active if it has equal or more than this number of active springs.
        :param maximumNrOfActiveConnections:  A node is active if it has equal or less than this number of active springs.
             Use a value < 0 to indicate that there is no maximum number of active connections.
        """
    def getNrOfActiveSprings(self, tolerance: float = 0.1) -> int: 
        """
         Get the number of active springs remaining after running the simulation.

        :param tolerance: springs under this length are considered inactive
        """
    def getNrOfIterations(self) -> int: 
        """
        Returns the number of iterations used for force relaxation.
        """
    def getNrOfNodes(self) -> int: 
        """
        Get the number of nodes considered in this simulation.
        """
    def getNrOfSprings(self) -> int: 
        """
        Get the number of springs considered in this simulation.

        :param tolerance: springs under this length are considered inactive
        """
    def getPressure(self) -> float: 
        """
        Returns the pressure at the current state of the simulation.
        """
    def getResidualNorm(self) -> float: 
        """
        Returns the residual norm at the current state of the simulation.
        """
    def getResiduals(self) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]: 
        """
        Returns the residuals at the current state of the simulation.
        """
    def getSolubleWeightFraction(self, tolerance: float = 0.1) -> float: 
        """
        Compute the weight fraction of springs connected to active
        springs (any depth). 

        Caution: ignores atom masses.
        """
    def getSpringDistances(self) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]: 
        """
        Get the current coordinate differences for all the springs.

        Returns:
             - distances: a vector of size 3*nrOfSprings, with each x, y, z values of the springs
        """
    def getSpringLengths(self) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]: 
        """
        Get the current lengths for all the springs.

        Returns:
             - distances: a vector of size nrOfSprings, with each the norm of the distances
        """
    def getStressTensor(self) -> numpy.ndarray[numpy.float64, _Shape[3, 3]]: 
        """
        Returns the stress tensor at the current state of the simulation.
        """
    def requiresAnotherRun(self) -> bool: 
        """
        For performance reasons, the objective is only minimised within the distances of one box.
        This means, that there is a possibility, e.g. for a single strand longer than two boxes, 
        that it would not be globally minimised.

        If the final displacement of one of the atoms is close 
        (1e-3, configurable via :func:`~pylimer_tools_cpp.pylimer_tools_cpp.MEHPForceRelaxation.configRerunEpsilon()`) 
        to the imposed min/max, after minimizing,
        this method would return true.
        """
    def runForceRelaxation(self, algorithm: str = 'LD_MMA', maxNrOfSteps: int = 250000, xTolerance: float = 1e-12, fTolerance: float = 1e-09) -> None: 
        """
        Run the simulation.
        Note that the final state of the minimization is persisted and reused if you use this method again.
        This is useful if you want to run a global optimization first and add a local one afterwards.
        As a consequence though, you cannot simply benchmark only this method; you must include the setup.

        :param algorithm: The algorithm to use for the force relaxation. Choices: see `NLopt Algorithms <https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/>`_
        :param maxNrOfSteps: The maximum number of steps to do during the simulation.
        :param xTolerance: The tolerance of the displacements as an exit condition.
        :param fTolerance: The tolerance of the force as an exit condition.
        :param is2d: Specify true if you want to evaluate the force relation only in x and y direction.
        """
    def setForceEvaluator(self, arg0: MEHPForceEvaluator) -> None: 
        """
        Reset the currently used force evaluator.
        """
    pass
class Molecule():
    """
    An (ideally) connected series of atoms/beads.
    """
    def __copy__(self) -> Molecule: ...
    def __getitem__(self, arg0: int) -> Atom: 
        """
        Access an atom by its vertex index.
        """
    def __init__(self, arg0: Box, arg1: igraph_s, arg2: MoleculeType, arg3: typing.Dict[int, float]) -> None: ...
    def __iter__(self) -> MoleculeIterator: 
        """
        Iterate through the atoms in this molecule.
        No specific order is guaranteed.
        """
    def __len__(self) -> int: 
        """
        Get the number of atoms in this molecule.
        """
    def computeBondLengths(self) -> typing.List[float]: 
        """
        Computes the length :math:`b` of each bond in the molecule, respecting periodic boundaries.
        """
    def computeEndToEndDistance(self) -> float: 
        """
        Compute the end-to-end distance (:math:`R_{ee}`) of this molecule. 

        CAUTION:
           Returns 0.0 if the molecule does not have two or more atoms.
           Returns -1.0 if not exactly 2 ends were found.
           Computes the distance between 2 atoms with functionality 1, 
           ignoring whether they are cross-linkers or not.
        """
    def computeEndToEndDistanceWithDerivedImageFlags(self) -> float: 
        """
        Compute the end-to-end distance (:math:`R_{ee}`) of this molecule,
        but ignoring the image flags attached to the atoms. 
        This only works for Molecules that can be lined up with 
        :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule.getAtomsLinedUp()`,
        as it needs the atoms sorted such that the periodic box can still be respected somewhat.

        CAUTION:
           Returns 0.0 if the molecule does not have two or more atoms.
           Requires bonds to be shorter than half the box length.
           Computes the distance between 2 atoms with functionality 1, 
           ignoring whether they are cross-linkers or not.
        """
    def computeRadiusOfGyration(self) -> float: 
        """
        Computes the radius of gyration, :math:`R_g^2` of this molecule.

        :math:`{R_g}^2 = \frac{1}{M} \sum_i m_i (r_i - r_{cm})^2`,
        where :math:`M` is the total mass of the molecule, :math:`r_{cm}`
        are the coordinates of the center of mass of the molecule and the
        sum is over all contained atoms.
        """
    def computeRadiusOfGyrationWithDerivedImageFlags(self) -> float: 
        """
        Computes the radius of gyration, :math:`R_g^2` of this molecule,
        but ignoring the image flags attached to the atoms.
        This only works for Molecules that can be lined up with 
        :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule.getAtomsLinedUp()`,
        as it needs the atoms sorted such that the periodic box can still be respected somewhat.
        In other words, this function computes the radius of gyration 
        assuming the distance between two lined-up beads 
        is smaller than half the periodic box in each direction.

        See also: :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule.computeRadiusOfGyration()`.
        """
    def computeTotalMass(self) -> float: 
        """
        Computes the total mass of this molecule.
        """
    def getAtomById(self, arg0: int) -> Atom: 
        """
        Get an atom by its id.
        """
    def getAtomForVertexId(self, arg0: int) -> Atom: 
        """
        Get an atom for a specific vertex.
        """
    def getAtomIdByIdx(self, arg0: int) -> int: 
        """
        Get the id of the atom by the vertex id of the underlying graph.
        """
    def getAtomTypes(self) -> typing.List[int]: 
        """
        Query all types (each one for each atom) ordered by atom vertex id.
        """
    def getAtoms(self) -> typing.List[Atom]: 
        """
        Returns all atom objects enclosed in this molecule.
        """
    def getAtomsConnectedTo(self, arg0: int) -> typing.List[Atom]: 
        """
        Get the atoms connected to a specified vertex id.
        """
    def getAtomsLinedUp(self, crosslinkType: int = 2, assumed_coordinates: bool = False) -> typing.List[Atom]: 
        """
        Returns all atom objects enclosed in this molecule based on the connectivity.

        This method works only for lone chains, atoms and loops, 
        as it throws an error if the molecule does not allow such a "line-up", 
        for example because of cross-links.

        Use the `crosslinkType` parameter to force the atoms in a primary loop 
        to start with the cross-link.
        """
    def getAtomsOfDegree(self, arg0: int) -> typing.List[Atom]: 
        """
        Get the atoms that have the specified number of bonds.
        """
    def getAtomsOfType(self, arg0: int) -> typing.List[Atom]: 
        """
        Get the atoms with the specified type.
        """
    def getBonds(self) -> typing.Dict[str, typing.List[int]]: 
        """
        Get all bonds. Returns a dict with three properties: 'bond_from', 'bond_to' and 'bond_type'.
        """
    def getConnectedAtoms(self, arg0: Atom) -> typing.List[Atom]: 
        """
        Get the atoms connected to a specified atom.

        Internally uses :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule.getAtomsConnectedTo`
        """
    def getEdges(self) -> typing.Dict[str, typing.List[int]]: 
        """
        Get all bonds. Returns a dict with three properties: 'edge_from', 'edge_to' and 'edge_type'.
        The order is not necessarily related to any structural property.

        NOTE:
           The integer values returned refer to the vertex ids, not the atom ids.
           Use :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule.getAtomIdByIdx` to translate them to atom ids, or 
           :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule.getBonds` to have that done for you.
        """
    def getIdxByAtomId(self, arg0: int) -> int: 
        """
        Get the vertex id of the underlying graph for an atom with a specified id.
        """
    def getKey(self) -> str: 
        """
        Get a unique identifier for this molecule.
        """
    def getLength(self) -> int: 
        """
        Counts and returns the number of atoms associated with this 
        molecule.
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
        Get the type of this molecule (see :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.MoleculeType` enum).
        """
    pass
class MoleculeIterator():
    """
    An iterator to iterate throught the atoms in :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule`.
    """
    def __iter__(self) -> MoleculeIterator: ...
    def __next__(self) -> Atom: ...
    pass
class MoleculeType():
    """
    Members:

      UNDEFINED : This value indicates that either the property was not set or not discovered.

      NETWORK_STRAND : 
               A network strand is a strand in a network.
          

      PRIMARY_LOOP : 
               A primary loop is a network strand looping from and to the same cross-linker.
          

      DANGLING_CHAIN : 
               A dangling chain is a network strand where only one end is attached to a cross-linker.
          

      FREE_CHAIN : 
               A free chain is a strand not connected to any cross-linker.
          
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
class NeighbourList():
    """
    Gives access to somewhat fast queries on the neighbourhood of atoms
    """
    def __init__(self, atoms: typing.List[Atom], box: Box, cutoff: float) -> None: 
        """
        Instantiates a new neighbour list
        """
    def getAtomsCloseTo(self, atom: Atom, upperCutoff: float = 1.0, lowerCutoff: float = 0.0, unwrapped: bool = True) -> typing.List[Atom]: 
        """
        List all atoms that are close to a given one. 

        It is possible to request it within a new cutoff, 
        though the underlying neighbour list will not be regenerated.
        For performance reasons, it is recommended to initialize a 
        new NeighbourList if you require a different cutoff, depending on your use case.

        You can use a negative value for the newCutoff to use the cutoff used for 
        filling the neighbour list buckets.
        """
    pass
class NonGaussianSpringForceEvaluator(MEHPForceEvaluator):
    """
    This is equal to a spring evaluator for Langevin chains.

    The force for a certain spring is given by:
    :math:`f = 0.5 \cdot \frac{1}{l} \scriptL^{-1}(\frac{r}{N\cdot l})`, 
    where :math:`r` is the spring [between cross-linkers] length 
    and :math:`\scriptL^{-1}` the inverse langevin function.

    Please note that the inverse langevin is only approximated.

    Recommended optimization algorithm: "LD_MMA"

    :param kappa: the spring constant :math:`\kappa`
    :param N: The number of links in a spring
    :param l: The  the length of a spring in the chain
    """
    def __init__(self, kappa: float = 1.0, N: float = 1.0, l: float = 1.0) -> None: 
        """
        Initialize this ForceEvaluator
        """
    pass
class OutputConfiguration():
    def __init__(self) -> None: 
        """
        Get an instance of this struct
        """
    @property
    def append(self) -> bool:
        """
             Whether to append to the file or truncate it
            

        :type: bool
        """
    @append.setter
    def append(self, arg0: bool) -> None:
        """
        Whether to append to the file or truncate it
        """
    @property
    def doubleValues(self) -> typing.List[ComputedDoubleValues]:
        """
        :type: typing.List[ComputedDoubleValues]
        """
    @doubleValues.setter
    def doubleValues(self, arg0: typing.List[ComputedDoubleValues]) -> None:
        pass
    @property
    def filename(self) -> str:
        """
        The file to write to. Empty means standard output (console).

        :type: str
        """
    @filename.setter
    def filename(self, arg0: str) -> None:
        """
        The file to write to. Empty means standard output (console).
        """
    @property
    def intValues(self) -> typing.List[ComputedIntValues]:
        """
        :type: typing.List[ComputedIntValues]
        """
    @intValues.setter
    def intValues(self, arg0: typing.List[ComputedIntValues]) -> None:
        pass
    @property
    def outputEvery(self) -> int:
        """
        How often to write the values to the output. 
              For averages, this value also says how many values will be averaged.
             

        :type: int
        """
    @outputEvery.setter
    def outputEvery(self, arg0: int) -> None:
        """
        How often to write the values to the output. 
              For averages, this value also says how many values will be averaged.
             
        """
    @property
    def useEvery(self) -> int:
        """
             for autocorrelation/averaging, how often to include values
            

        :type: int
        """
    @useEvery.setter
    def useEvery(self, arg0: int) -> None:
        """
        for autocorrelation/averaging, how often to include values
        """
    pass
class SimpleSpringMEHPForceEvaluator(MEHPForceEvaluator):
    """
    This is equal to a spring evaluator for Gaussian chains.

    The force for a certain spring is given by:
    :math:`f = 0.5 \cdot \kappa r`, 
    where :math:`r` is the spring [between cross-linkers] length.

    Recommended optimization algorithm: "LD_LBFGS"

    :param kappa: the spring constant :math:`\kappa`
    """
    def __init__(self, kappa: float = 1.0) -> None: ...
    pass
class SimplifiedBalanceNetwork():
    """
    A more efficient structure of the network for use in MEHP force balance.
    Consists usually only of the cross- and slip-links.
    """
    @property
    def boxLengths(self) -> float:
        """
        :type: float
        """
    @property
    def coordinates(self) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.float64, _Shape[m, 1]]
        """
    @property
    def linkIndicesOfSprings(self) -> typing.List[typing.List[int]]:
        """
        :type: typing.List[typing.List[int]]
        """
    @property
    def linkIsSliplink(self) -> numpy.ndarray[bool, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[bool, _Shape[m, 1]]
        """
    @property
    def localToGlobalSpringIndex(self) -> typing.List[typing.List[int]]:
        """
        :type: typing.List[typing.List[int]]
        """
    @property
    def nrOfCrossLinks(self) -> int:
        """
        :type: int
        """
    @property
    def nrOfCrosslinkSwapsEndured(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def nrOfLinks(self) -> int:
        """
        :type: int
        """
    @property
    def nrOfPartialSprings(self) -> int:
        """
        :type: int
        """
    @property
    def nrOfSprings(self) -> int:
        """
        :type: int
        """
    @property
    def oldAtomIds(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def partialToFullSpringIndex(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def springContourLength(self) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.float64, _Shape[m, 1]]
        """
    @property
    def springCoordinateIndexA(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def springCoordinateIndexB(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def springIndexA(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def springIndexB(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def springIndicesOfLinks(self) -> typing.List[typing.List[int]]:
        """
        :type: typing.List[typing.List[int]]
        """
    @property
    def springPartCoordinateIndexA(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def springPartCoordinateIndexB(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def springPartIndexA(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def springPartIndexB(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def volume(self) -> float:
        """
        :type: float
        """
    pass
class SimplifiedNetwork():
    """
    A more efficient structure of the network for use in MEHP.
    Consists usually only of the cross-links.
    """
    @property
    def boxLengths(self) -> float:
        """
        :type: float
        """
    @property
    def coordinates(self) -> numpy.ndarray[numpy.float64, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.float64, _Shape[m, 1]]
        """
    @property
    def nrOfNodes(self) -> int:
        """
        :type: int
        """
    @property
    def nrOfSprings(self) -> int:
        """
        :type: int
        """
    @property
    def oldAtomIds(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def springCoordinateIndexA(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def springCoordinateIndexB(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def springIndexA(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def springIndexB(self) -> numpy.ndarray[numpy.int32, _Shape[m, 1]]:
        """
        :type: numpy.ndarray[numpy.int32, _Shape[m, 1]]
        """
    @property
    def volume(self) -> float:
        """
        :type: float
        """
    pass
class StructureSimplificationMode():
    """
    Members:

      NO_SIMPLIFICATION

      X2F_ONLY

      INACTIVE_ONLY

      ALL_TIM

      ALL_ANDREI
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
    ALL_ANDREI: pylimer_tools_cpp.pylimer_tools_cpp.StructureSimplificationMode # value = <StructureSimplificationMode.ALL_ANDREI: 4>
    ALL_TIM: pylimer_tools_cpp.pylimer_tools_cpp.StructureSimplificationMode # value = <StructureSimplificationMode.ALL_TIM: 3>
    INACTIVE_ONLY: pylimer_tools_cpp.pylimer_tools_cpp.StructureSimplificationMode # value = <StructureSimplificationMode.INACTIVE_ONLY: 2>
    NO_SIMPLIFICATION: pylimer_tools_cpp.pylimer_tools_cpp.StructureSimplificationMode # value = <StructureSimplificationMode.NO_SIMPLIFICATION: 0>
    X2F_ONLY: pylimer_tools_cpp.pylimer_tools_cpp.StructureSimplificationMode # value = <StructureSimplificationMode.X2F_ONLY: 1>
    __members__: dict # value = {'NO_SIMPLIFICATION': <StructureSimplificationMode.NO_SIMPLIFICATION: 0>, 'X2F_ONLY': <StructureSimplificationMode.X2F_ONLY: 1>, 'INACTIVE_ONLY': <StructureSimplificationMode.INACTIVE_ONLY: 2>, 'ALL_TIM': <StructureSimplificationMode.ALL_TIM: 3>, 'ALL_ANDREI': <StructureSimplificationMode.ALL_ANDREI: 4>}
    pass
class Universe():
    """
    Represents a full Polymer Network structure, a collection of molecules.
    """
    def __copy__(self) -> Universe: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, Lx: float, Ly: float, Lz: float) -> None: 
        """
        Instantiate this Universe (Collection of Molecules) providing the box lengths.
        """
    def __setstate__(self, arg0: tuple) -> None: ...
    def addAngles(self, angles_from: typing.List[int], angles_via: typing.List[int], angles_to: typing.List[int], angle_types: typing.List[int]) -> None: 
        """
        Add angles to the Universe. No relation to the underlying graph, just a method to preserve read & write capabilities
        """
    def addAtoms(self, ids: typing.List[int], types: typing.List[int], x: typing.List[float], y: typing.List[float], z: typing.List[float], nx: typing.List[int], ny: typing.List[int], nz: typing.List[int]) -> None: 
        """
        Add atoms to the Universe, vertices to the underlying graph.
        """
    @typing.overload
    def addBonds(self, bonds_from: typing.List[int], bonds_to: typing.List[int]) -> None: 
        """
        Add bonds to the underlying atoms, edges to the underlying graph. If the connected atoms are not found, the bonds are silently skipped.

        Add bonds to the underlying atoms, edges to the underlying graph. 
        """
    @typing.overload
    def addBonds(self, num_bonds: int, bonds_from: typing.List[int], bonds_to: typing.List[int], bond_types: typing.List[int], ignore_non_existent_atoms: bool = False, simplify_universe: bool = True) -> None: ...
    def addBondsWithTypes(self, bonds_from: typing.List[int], bonds_to: typing.List[int], bond_types: typing.List[int]) -> None: 
        """
        Add bonds to the underlying atoms, edges to the underlying graph. If the connected atoms are not found, the bonds are silently skipped.
        """
    def addDihedralAngles(self, angles_from: typing.List[int], angles_via1: typing.List[int], angles_via2: typing.List[int], angles_to: typing.List[int], angle_types: typing.List[int]) -> None: 
        """
        Add dihedral angles to the Universe. No relation to the underlying graph, just a method to preserve read & write capabilities
        """
    def computeBondLengths(self) -> typing.List[float]: 
        """
        Computes the length :math:`b` of each bond in the molecule, respecting periodic boundaries.
        """
    def computeDxs(self, atomIdsTo: typing.List[int], atomIdsFrom: typing.List[int]) -> typing.List[float]: 
        """
        Compute the dx distance for certain bonds (length in x direction).
        """
    def computeDys(self, atomIdsTo: typing.List[int], atomIdsFrom: typing.List[int]) -> typing.List[float]: 
        """
        Compute the dy distance for certain bonds (length in y direction).
        """
    def computeDzs(self, atomIdsTo: typing.List[int], atomIdsFrom: typing.List[int]) -> typing.List[float]: 
        """
        Compute the dz distance for certain bonds (length in z direction).
        """
    def computeEndToEndDistances(self, arg0: int) -> typing.List[float]: 
        """
        Compute the end-to-end distance of each strand in the network.

        NOTE:
             Internally, this uses the :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule.computeEndToEndDistance`.
             All its cautionary facts apply.
        """
    def computeMeanEndToEndDistance(self, arg0: int) -> float: 
        """
        Computes the mean of the end-to-end distances of each strand in the network.

        NOTE:
             Internally, this uses the :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule.computeEndToEndDistance`.
             All its cautionary facts apply.
             Invalid strands (where said function returns 0.0 or -1.0) are ignored.
        """
    def computeMeanSquareEndToEndDistance(self, crossLinkerType: int, onlyThoseWithTwoCrosslinkers: bool = False) -> float: 
        """
        Computes the mean square of the end-to-end distances of each strand in the network.

        NOTE:
             Internally, this uses the :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Molecule.computeEndToEndDistance`.
             All its cautionary facts apply.
             Invalid strands (where said function returns 0.0 or -1.0) are ignored.
        """
    def computeMeanStrandLength(self, crossLinkerType: int) -> float: 
        """
        Compute the mean number of beads per strand.
        """
    def computeNumberAverageMolecularWeight(self, crossLinkerType: int) -> float: 
        """
        Compute the number average molecular weight.

        NOTE: 
              Cross-linkers are ignored completely.
        """
    def computePolydispersityIndex(self, crossLinkerType: int) -> float: 
        """
        Compute the polydispersity indiex: 
        the weight average molecular weight over the number average molecular weight.
        """
    def computeTotalMass(self) -> float: 
        """
        Compute the total mass of this network/universe in whatever mass unit was used when 
        :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.setMasses()` was called.
        """
    def computeWeightAverageMolecularWeight(self, crossLinkerType: int) -> float: 
        """
        Compute the weight average molecular weight.

        NOTE: 
              Cross-linkers are ignored completely.
        """
    def computeWeightFractions(self) -> typing.Dict[int, float]: 
        """
        Compute the weight fractions of each atom type in the network.

        If no masses are stored, 
        """
    def countAtomTypes(self) -> typing.Dict[int, int]: 
        """
        Count how often each atom type is present.
        """
    def countAtomsInSkinDistance(self, distances: typing.List[float], unwrapped: bool = False) -> typing.List[int]: 
        """
        This is a function that may help you to compute the radial distribution function.
        It loops the 

        Parameters:
             - distances: the edges of the bins
             - unwrapped: whether to measure the distance in unwrapped coordinates or as PBC-corrected distance
        """
    def detectAngles(self) -> typing.Dict[str, typing.List[int]]: 
        """
        Returns just as 
                  :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.getAngles`, 
                  but all angles that are detected in the network, rather than the one already set.
                  Note that the angle types are determined by 
                  :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.hashAngleType`,
                  which serves angle types that should be mapped by you back to smaller numbers, 
                  before serving them to :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.addAngles`.
                 
        """
    def detectDihedralAngles(self) -> typing.Dict[str, typing.List[int]]: 
        """
        Returns just as 
                  :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.getDihedralAngles`, 
                  but all dihedral angles that are detected in the network, rather than the one already set.
                  Note that the angle types are determined by 
                  :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.hashDihedralAngleType`,
                  which serves angle types that should be mapped by you back to smaller numbers, 
                  before serving them to :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.addDiheralAngles`.
                 
        """
    def determineEffectiveFunctionalityPerType(self) -> typing.Dict[int, float]: 
        """
        Find the average functionality of each atom type in the network.
        """
    def determineFunctionalityPerType(self) -> typing.Dict[int, int]: 
        """
        Find the maximum functionality of each atom type in the network.
        """
    def findLoops(self, crossLinkerType: int, maxLength: int = -1, skipSelfLoops: bool = False) -> typing.Dict[int, typing.List[typing.List[Atom]]]: 
        """
        Decompose the Universe into loops.
        The primary index specifies the degree of the loop.

        CAUTION:
           There are exponentially many paths between two cross-linkers of a network,
           and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
           You can use the maxLength parameter to restrict the algorithm to only search for loops up to a certain length.
           Use a negative value to find all loops and paths.
        """
    def findMinimalOrderLoopFrom(self, loopStart: int, loopStep1: int, maxLength: int = -1, skipSelfLoops: bool = False) -> typing.List[Atom]: 
        """
        Decompose the Universe into loops.
        The primary index specifies the degree of the loop.

        CAUTION:
           There are exponentially many paths between two cross-linkers of a network,
           and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
           You can use the maxLength parameter to restrict the algorithm to only search for loops up to a certain length.
           Use a negative value to find all loops and paths.
        """
    def getAngles(self) -> typing.Dict[str, typing.List[int]]: 
        """
        Get all angles added to this network.

        Returns a dict with three properties: 'angle_from', 'angle_via' and 'angle_to'.

        NOTE:
            The integer values returned refer to the the atom ids, not the vertex ids.
            Use :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.getIdxByAtomId` to translate them to vertex ids.
        """
    def getAtom(self, atomId: int) -> Atom: 
        """
        Find an atom by its ID.
        """
    def getAtomByVertexIdx(self, vertexId: int) -> Atom: 
        """
        Find an atom by the ID of the vertex of the underlying graph.
        """
    def getAtomForVertexId(self, vertexId: int) -> Atom: 
        """
        Get an atom for a specific vertex.
        """
    def getAtomIdByIdx(self, atomId: int) -> int: 
        """
        Get the id of the atom by the vertex id of the underlying graph.
        """
    def getAtomTypes(self) -> typing.List[int]: 
        """
        Get all types (each one for each atom) ordered by atom vertex id.
        """
    def getAtoms(self) -> typing.List[Atom]: 
        """
        Get all atoms.
        """
    def getAtomsConnectedTo(self, arg0: int) -> typing.List[Atom]: 
        """
        Get the atoms connected to a specified vertex id.
        """
    def getAtomsOfDegree(self, arg0: int) -> typing.List[Atom]: 
        """
        Get the atoms that have the specified number of bonds.
        """
    def getAtomsOfType(self, arg0: int) -> typing.List[Atom]: 
        """
        Find many atom by their type.
        """
    def getBonds(self) -> typing.Dict[str, typing.List[int]]: 
        """
        Get all bonds. Returns a dict with three properties: 'bond_from', 'bond_to' and 'bond_type'.
        The order is not necessarily related to any structural characteristic.
        """
    def getBox(self) -> Box: 
        """
        Get the underlying bounding box object.
        """
    def getChainsWithCrosslinker(self, crossLinkerType: int) -> typing.List[Molecule]: 
        """
        Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms, without omitting the cross-linkers.
        In turn, e.g. for a tetrafunctional cross-linker, it will be 4 times in the resulting molecules.

        NOTE:
           Cross-linkers without bonds to non-cross-linkers are not returned 
           (i.e., cross-linker-cross-linker bonds, or single cross-linkers, are not counted as strands).
        """
    def getClusters(self) -> typing.List[Universe]: 
        """
        Get the components of the universe that are not connected to each other.
        Returns a list of :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe`s.
        Unconnected, free atoms/beads become their own :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe`.
        """
    def getConnectedAtoms(self, arg0: Atom) -> typing.List[Atom]: 
        """
        Get the atoms connected to a specified atom.

        Internally uses :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.getAtomsConnectedTo`
        """
    def getEdges(self) -> typing.Dict[str, typing.List[int]]: 
        """
        Get all edges. Returns a dict with three properties: 'edge_from', 'edge_to' and 'edge_type'.
        The order is not necessarily related to any structural characteristic.

        NOTE:
           The integer values returned refer to the vertex ids, not the atom ids.
           Use :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.getAtomIdByIdx` to translate them to atom ids, or
           :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.getBonds` to have that done for you.
        """
    def getIdxByAtomId(self, atomId: int) -> int: 
        """
        Get the vertex id of the underlying graph for an atom with a specified id.
        """
    def getMasses(self) -> typing.Dict[int, float]: 
        """
        Get the mass of one atom per type
        """
    def getMolecules(self, atomTypeToOmit: int) -> typing.List[Molecule]: 
        """
        Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms.

        Reduces the Universe to a list of molecules. 
        Specify the crossLinkerType to an existing type id, 
        then those atoms will be omitted, and this function returns chains instead.
        """
    def getNetworkOfCrosslinker(self, crossLinkerType: int) -> Universe: 
        """
        Reduce the network to contain only cross-linkers, replacing all the strands with a single bond.
        Useful e.g. to reduce the memory useage and runtime of 
        :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.findLoops()` or 
        :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.hasInfiniteStrand()`.

        Further use :func:`~pylimer_tools_cpp.pylimer_tools_cpp.Universe.simplify()` to remove primary loops.
        """
    def getNrOfAngles(self) -> int: 
        """
        Query the number of angles that have been added to this universe.
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
    def getNrOfDihedralAngles(self) -> int: 
        """
        Query the number of dihedralangles that have been added to this universe.
        """
    def getTimestep(self) -> int: 
        """
        Query the timestep when this universe was captured.
        """
    def getVolume(self) -> float: 
        """
        Query the volume of the underlying bounding box.
        """
    def hasInfiniteStrand(self, arg0: int, arg1: int) -> bool: 
        """
        Checks whether there is a strand (with cross-linker) in the universe that loops through periodic images without coming back.

         CAUTION:
            There are exponentially many paths between two cross-linkers of a network,
            and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
        """
    def hashAngleType(self, angle_from: int, angle_via: int, angle_to: int) -> int: 
        """
        Convert the three integers 
        """
    def hashDihedralAngleType(self, angle_from: int, angle_via1: int, angle_via2: int, angle_to: int) -> int: 
        """
        Convert the four integers 
        """
    def removeAtoms(self, atomIds: typing.List[int]) -> None: 
        """
        Remove atoms and all associated bonds by their atom ids. 
        """
    def removeBonds(self, bonds_from: typing.List[int], bonds_to: typing.List[int]) -> None: 
        """
        Remove bonds by their connected atom ids. 
        """
    def removeBondsOfType(self, bond_type: int) -> None: 
        """
        Remove bonds with a specific type. 
        """
    def replaceAtom(self, atomId: int, replacementAtom: Atom) -> None: 
        """
        Replace the properties of an atom with the properties of another given atom.
        """
    def setBox(self, box: Box, rescaleAtoms: bool = False) -> None: 
        """
        Override the currently assigned box with the one specified.
        """
    def setBoxLengths(self, Lx: float, Ly: float, Lz: float, rescaleAtoms: bool = False) -> None: 
        """
        Set the box side lengths.
        """
    def setMasses(self, massPerType: typing.Dict[int, float]) -> None: 
        """
        Set the mass per type of atom.
        """
    def setTimestep(self, timestep: int) -> None: 
        """
        Set the timestep when this Universe was captured.
        """
    def setVertexProperty(self, vertexId: int, propertyName: str, value: float) -> None: 
        """
        Set a specific property for a specific vertex.
        """
    def simplify(self) -> None: 
        """
        Remove self links and double bonds. This function is called automatically after adding bonds.
        """
    pass
class UniverseSequence():
    """
    This class represents a sequence of Universes, with the Universe's data
    only being read on request. Dump files are read at once in order
    to know how many timesteps/universes are available in total 
    (but the universes' data is not read on first look through the file).
    This, while it can lead to two (or more) reads of the whole file, 
    is a measure in order to enable low memory useage if needed (i.e. for large dump files).
    Use Python's iterator to have this UniverseSequence only ever retain one universe in memory.
    Alternatively, use :func:`~pylimer_tools_cpp.pylimer_tools_cpp.UniverseSequence.forgetAtIndex`
    to have the UniverseSequence forget about already read universes.
    """
    def __getitem__(self, arg0: int) -> Universe: 
        """
        Get a universe by its index.
        """
    def __init__(self) -> None: ...
    def __iter__(self) -> LazyUniverseSequenceIterator: 
        """
        Lazily (memory-efficiently) iterate through all the universes in this sequence.
        This is the standard Python iteration way. 

        Example:

        .. code::

            for (universe in universeSequence):
                 # do something with the universe
                 pass


        Note: 
            this iterator is supposed to be memory-efficient. Therefore, no cache is kept;
            iterating twice will lead to the file(s) being read twice 
            (plus, for dump files, a third time initially to determine the number of universes in the file).
        """
    def __len__(self) -> int: 
        """
        Get the number of universes
        """
    def atIndex(self, index: int) -> Universe: 
        """
        Get the Universe at the given index (as of in the sequence given by the dump file).
        """
    def computeMsdForAtoms(self, atom_ids: typing.List[int], nr_of_origins: int = 25, reduce_memory: bool = False) -> typing.Dict[int, float]: 
        """
        Compute the mean square displacement for atoms with the specified ids
        """
    def forgetAtIndex(self, index: int) -> None: 
        """
        Clear the memory of the Universe at the given index (as of in the 
                   sequence given by the dump file).
        """
    def getAll(self) -> typing.List[Universe]: 
        """
        Get all universes initialized back in a list.
        For big dump files or lots of data files, this might lead to memory issues.
        Use :func:`~pylimer_tools_cpp.pylimer_tools_cpp.UniverseSequence.__iter__`
        to have
        or :func:`~pylimer_tools_cpp.pylimer_tools_cpp.UniverseSequence.atIndex`
        and :func:`~pylimer_tools_cpp.pylimer_tools_cpp.UniverseSequence.forgetAtIndex`
        to craft a more memory-efficient retrieval mechanism.
        """
    def getLength(self) -> int: 
        """
        Get the number of universes in this sequence.
        """
    def initializeFromDataSequence(self, data_files: typing.List[str]) -> None: 
        """
        Reset and initialize the Universes from an ordered list of Lammps data (:code:`write_data`) files.
        """
    def initializeFromDumpFile(self, initial_data_file: str, dump_file: str) -> None: 
        """
        Reset and initialize the Universes from a Lammps :code:`dump` output. 

        NOTE:
             If you have not output the id of the atoms in the dump file, they will be assigned sequentially. 
             If you have not output the type of the atoms in the dump file, they will be set to -1 if they cannot be infered from the data file.
        """
    def next(self) -> Universe: 
        """
        Get the Universe that's next in the sequence.
        """
    def resetIterator(self) -> None: 
        """
        Reset the internal iterator, such that a subsequent call to 
        :func:`~pylimer_tools_cpp.pylimer_tools_cpp.UniverseSequence.next` returns the first one again.
        """
    def setDataFileAtomStyle(self, atom_styles: typing.List[AtomStyle]) -> None: 
        """
        Set the format of the data files to be read. See :obj:`~pylimer_tools_cpp.pylimer_tools_cpp.AtomStyle`.
        """
    pass
def computeStoichiometricInbalance(arg0: Universe, arg1: int, arg2: int, arg3: typing.Dict[int, int]) -> float:
    """
    Compute stoichiometric inbalance
    """
def doRandomWalkChainFromTo(box: Box, from_coordinates: Annotated[typing.List[float], FixedSize(3)], to_coordinates: Annotated[typing.List[float], FixedSize(3)], chainLen: int, beadDistance: float = 1.0, seed: str = '') -> typing.Dict[str, typing.List[float]]:
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
