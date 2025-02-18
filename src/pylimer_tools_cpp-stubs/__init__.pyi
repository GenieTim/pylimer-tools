"""

    PylimerTools Cpp
    -----------------

    A collection of utility python functions for handling LAMMPS output and polymers in Python.

    .. autosummary::
        :toctree: _generate


"""
from __future__ import annotations
import numpy
import scipy.sparse
import typing
__all__ = [
    'Atom',
    'AtomPairEntanglements',
    'AtomStyle',
    'AveFileReader',
    'Box',
    'ComputedDoubleValues',
    'ComputedIntValues',
    'DANGLING_CHAIN',
    'DPDSimulator',
    'DataFileReader',
    'DataFileWriter',
    'DumpFileReader',
    'ExitReason',
    'FREE_CHAIN',
    'LazyUniverseSequenceIterator',
    'LinkSwappingMode',
    'MCUniverseGenerator',
    'MEHPForceBalance',
    'MEHPForceEvaluator',
    'MEHPForceRelaxation',
    'Molecule',
    'MoleculeIterator',
    'MoleculeType',
    'NETWORK_STRAND',
    'NeighbourList',
    'NonGaussianSpringForceEvaluator',
    'NormalModeAnalyzer',
    'OutputConfiguration',
    'PRIMARY_LOOP',
    'SimpleSpringMEHPForceEvaluator',
    'SimplifiedBalanceNetwork',
    'SimplifiedNetwork',
    'StructureSimplificationMode',
    'UNDEFINED',
    'Universe',
    'UniverseSequence',
    'compute_stoichiometric_imbalance',
    'do_linear_walk_chain_from_to',
    'do_random_walk',
    'do_random_walk_chain_from_to',
    'do_random_walk_chain_from_to_mc',
    'inverse_langevin',
    'predict_gelation_point',
    'randomly_sample_entanglements',
    'split_csv',
    'version_information']


class Atom:
    """

           A single bead or atom

    """
    __hash__: typing.ClassVar[None] = None

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __eq__(self, arg0: Atom) -> bool:
        ...

    def __getstate__(self) -> tuple:
        ...

    def __init__(self, id: int, type: int, x: float, y: float,
                 z: float, nx: int, ny: int, nz: int) -> None:
        """
        Construct this atom
        """

    def __setstate__(self, arg0: tuple) -> None:
        """
        Provides support for pickling
        """

    def compute_vector_to(self, to_atom: Atom, pbc_box: Box) -> numpy.ndarray:
        """
                    Compute the vector to another atom.
        """

    def distance_to(self, to_atom: Atom, pbc_box: Box) -> float:
        """
                    Compute the distance to another atom.
        """

    def distance_to_unwrapped(self, arg0: Atom, arg1: Box) -> float:
        """
        Compute the distance to another atom respecting the periodic image flag.
        """

    def get_coordinates(self) -> numpy.ndarray:
        ...

    def get_id(self) -> int:
        """
                    Get the id of the atom.
        """

    def get_nx(self) -> int:
        """
        Get the box image that the atom is in in x direction (also known as `ix` or `nx`).
        """

    def get_ny(self) -> int:
        """
        Get the box image that the atom is in in y direction (also known as `iy` or `ny`).
        """

    def get_nz(self) -> int:
        """
        Get the box image that the atom is in in z direction (also known as `iz` or `nz`).
        """

    def get_type(self) -> int:
        """
                    Get the type of the atom.
        """

    def get_unwrapped_coordinates(self, arg0: Box) -> numpy.ndarray:
        ...

    def get_x(self) -> float:
        """
                    Get the x coordinate of the atom.
        """

    def get_y(self) -> float:
        """
                    Get the y coordinate of the atom.
        """

    def get_z(self) -> float:
        """
                    Get the z coordinate of the atom.
        """

    def vector_to_unwrapped(self, arg0: Atom, arg1: Box) -> numpy.ndarray:
        """
        Compute the vector to another atom respecting the periodic image flags.
        """


class AtomPairEntanglements:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self) -> None:
        """
        Get an instance of this struct
        """
    @property
    def pair_of_atom(self) -> list[int]:
        """
              An index in the pairs_of_atoms if the atom is part of a pair, -1 else.
        """
    @pair_of_atom.setter
    def pair_of_atom(self, arg0: list[int]) -> None:
        ...

    @property
    def pairs_of_atoms(self) -> list[tuple[int, int]]:
        """
              A list of pairs of atom ids that are close together and could be entanglements
        """
    @pairs_of_atoms.setter
    def pairs_of_atoms(self, arg0: list[tuple[int, int]]) -> None:
        ...


class AtomStyle:
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
    ANGLE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.ANGLE: 1>
    ATOMIC: typing.ClassVar[AtomStyle]  # value = <AtomStyle.ATOMIC: 2>
    BODY: typing.ClassVar[AtomStyle]  # value = <AtomStyle.BODY: 3>
    BOND: typing.ClassVar[AtomStyle]  # value = <AtomStyle.BOND: 4>
    # value = <AtomStyle.BPM_SPHERE: 20>
    BPM_SPHERE: typing.ClassVar[AtomStyle]
    CHARGE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.CHARGE: 5>
    DIELECTRIC: typing.ClassVar[AtomStyle]  # value = <AtomStyle.DIELECTRIC: 6>
    DIPOLE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.DIPOLE: 7>
    DPD: typing.ClassVar[AtomStyle]  # value = <AtomStyle.DPD: 8>
    EDPD: typing.ClassVar[AtomStyle]  # value = <AtomStyle.EDPD: 9>
    ELECTRON: typing.ClassVar[AtomStyle]  # value = <AtomStyle.ELECTRON: 10>
    ELLIPSOID: typing.ClassVar[AtomStyle]  # value = <AtomStyle.ELLIPSOID: 11>
    FULL: typing.ClassVar[AtomStyle]  # value = <AtomStyle.FULL: 12>
    HYBRID: typing.ClassVar[AtomStyle]  # value = <AtomStyle.HYBRID: 26>
    LINE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.LINE: 13>
    MDPD: typing.ClassVar[AtomStyle]  # value = <AtomStyle.MDPD: 14>
    MOLECULAR: typing.ClassVar[AtomStyle]  # value = <AtomStyle.MOLECULAR: 15>
    NONE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.NONE: 0>
    PERI: typing.ClassVar[AtomStyle]  # value = <AtomStyle.PERI: 16>
    SMD: typing.ClassVar[AtomStyle]  # value = <AtomStyle.SMD: 17>
    SPH: typing.ClassVar[AtomStyle]  # value = <AtomStyle.SPH: 18>
    SPHERE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.SPHERE: 19>
    SPIN: typing.ClassVar[AtomStyle]  # value = <AtomStyle.SPIN: 21>
    TDPD: typing.ClassVar[AtomStyle]  # value = <AtomStyle.TDPD: 22>
    TEMPLATE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.TEMPLATE: 23>
    TRI: typing.ClassVar[AtomStyle]  # value = <AtomStyle.TRI: 24>
    # value = <AtomStyle.WAVEPACKET: 25>
    WAVEPACKET: typing.ClassVar[AtomStyle]
    __members__: typing.ClassVar[dict[str, AtomStyle]]  # value = {'NONE': <AtomStyle.NONE: 0>, 'ANGLE': <AtomStyle.ANGLE: 1>, 'ATOMIC': <AtomStyle.ATOMIC: 2>, 'BODY': <AtomStyle.BODY: 3>, 'BOND': <AtomStyle.BOND: 4>, 'CHARGE': <AtomStyle.CHARGE: 5>, 'DIELECTRIC': <AtomStyle.DIELECTRIC: 6>, 'DIPOLE': <AtomStyle.DIPOLE: 7>, 'DPD': <AtomStyle.DPD: 8>, 'EDPD': <AtomStyle.EDPD: 9>, 'ELECTRON': <AtomStyle.ELECTRON: 10>, 'ELLIPSOID': <AtomStyle.ELLIPSOID: 11>, 'FULL': <AtomStyle.FULL: 12>, 'LINE': <AtomStyle.LINE: 13>, 'MDPD': <AtomStyle.MDPD: 14>, 'MOLECULAR': <AtomStyle.MOLECULAR: 15>, 'PERI': <AtomStyle.PERI: 16>, 'SMD': <AtomStyle.SMD: 17>, 'SPH': <AtomStyle.SPH: 18>, 'SPHERE': <AtomStyle.SPHERE: 19>, 'BPM_SPHERE': <AtomStyle.BPM_SPHERE: 20>, 'SPIN': <AtomStyle.SPIN: 21>, 'TDPD': <AtomStyle.TDPD: 22>, 'TEMPLATE': <AtomStyle.TEMPLATE: 23>, 'TRI': <AtomStyle.TRI: 24>, 'WAVEPACKET': <AtomStyle.WAVEPACKET: 25>, 'HYBRID': <AtomStyle.HYBRID: 26>}

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...


class AveFileReader:
    """

              Alternative implementation of the data file reader implemented in
              :func:`pylimer_tools.read_lammps_output_file.read_averages_file`.

              This implementation is better for certain use cases, worse for others.
              In the end, only performance and memory usage are different.
              For moderately sized and small files, we recommend to use the Python interface instead.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, file_path: str) -> None:
        ...

    def autocorrelate_column(self, column_index: int,
                             delta_indices: list[int]) -> list[float]:
        """
                  Do autocorrelation on one particular column for a specified set of delta indices.

                  Assumes the data is equally spaced.
        """

    def autocorrelate_column_difference(
            self, column_index1: int, column_index2: int, delta_indices: list[int]) -> list[float]:
        """
                  Do autocorrelation on the difference between two particular columns for a specified set of delta indices.

                  Assumes the data is equally spaced.
        """

    def get_column_names(self) -> list[str]:
        ...

    def get_data(self) -> list[list[float]]:
        ...

    def get_nr_of_columns(self) -> int:
        ...

    def get_nr_of_rows(self) -> int:
        ...


class Box:
    """

            The box that the simulation is run in.

            NOTE:
              currently, only rectangular boxes are supported.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __getstate__(self) -> tuple:
        ...

    @typing.overload
    def __init__(self, arg0: float, arg1: float, arg2: float) -> None:
        ...

    @typing.overload
    def __init__(self, arg0: float, arg1: float, arg2: float,
                 arg3: float, arg4: float, arg5: float) -> None:
        ...

    def __setstate__(self, arg0: tuple) -> None:
        """
        Provides support for pickling.
        """

    def apply_pbc(self, distances: numpy.ndarray) -> numpy.ndarray:
        """
              Apply periodic boundary conditions (PBC): adjust the specified distances to fit into this box.
        """

    def apply_simple_shear(self, shear_magnitude: float,
                           shear_direction: int = 0) -> None:
        """
                  Apply a simple shear to the box.

                  CAUTION:
                    currently, this is not supported for all operations.

                  For shear magnitude, you specify the angle :math:`\\gamma`.

                  For the shearDirection parameter, you can specify `0` for x, `1` for y and `2` for z, respectively.
                  Specify another integer to disable the shear.
        """

    def get_bounding_box(self) -> Box:
        """
             Get an orthogonal box that encloses this box.
             For non-sheared boxes, the resulting box is identical to the current box.
        """

    def get_high_x(self) -> float:
        ...

    def get_high_y(self) -> float:
        ...

    def get_high_z(self) -> float:
        ...

    def get_l(self) -> numpy.ndarray:
        """
                  Get the three lengths of the box in an array/list.
        """

    def get_low_x(self) -> float:
        ...

    def get_low_y(self) -> float:
        ...

    def get_low_z(self) -> float:
        ...

    def get_lx(self) -> float:
        """
                    Get the length of the box in x direction.
        """

    def get_ly(self) -> float:
        """
                    Get the length of the box in y direction.
        """

    def get_lz(self) -> float:
        """
                    Get the length of the box in z direction.
        """

    def get_offset(self, distances: numpy.ndarray) -> numpy.ndarray:
        """
             Compute the offset required to compensate for periodic boundary conditions.

             Useful e.g. if you are using absolute coordinates for distances, but
             still need an infinite network,
             e.g., if the bonds need to be able to get longer than half the box.
        """

    def get_volume(self) -> float:
        """
                    Compute the volume of the box.

                    :math:`V = L_x \\cdot L_y \\cdot L_z`
        """

    def is_valid_offset(self, potential_offset: numpy.ndarray,
                        abs_precision: float = 1e-05) -> bool:
        """
                  Check whether the passed offest is a valid one in this box.
        """


class ComputedDoubleValues:
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
    MAX_B: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.MAX_B: 15>
    # value = <ComputedDoubleValues.MEAN_B: 14>
    MEAN_B: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.MSD: 16>
    MSD: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.PRESSURE: 3>
    PRESSURE: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.STRESS_NXY: 11>
    STRESS_NXY: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.STRESS_NXZ: 13>
    STRESS_NXZ: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.STRESS_NYZ: 12>
    STRESS_NYZ: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.STRESS_XX: 5>
    STRESS_XX: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.STRESS_XY: 8>
    STRESS_XY: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.STRESS_XZ: 10>
    STRESS_XZ: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.STRESS_YY: 6>
    STRESS_YY: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.STRESS_YZ: 9>
    STRESS_YZ: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.STRESS_ZZ: 7>
    STRESS_ZZ: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.TEMPERATURE: 4>
    TEMPERATURE: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.TIME: 1>
    TIME: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.TIMESTEP: 0>
    TIMESTEP: typing.ClassVar[ComputedDoubleValues]
    # value = <ComputedDoubleValues.VOLUME: 2>
    VOLUME: typing.ClassVar[ComputedDoubleValues]
    __members__: typing.ClassVar[dict[str, ComputedDoubleValues]]  # value = {'TIMESTEP': <ComputedDoubleValues.TIMESTEP: 0>, 'TIME': <ComputedDoubleValues.TIME: 1>, 'VOLUME': <ComputedDoubleValues.VOLUME: 2>, 'PRESSURE': <ComputedDoubleValues.PRESSURE: 3>, 'TEMPERATURE': <ComputedDoubleValues.TEMPERATURE: 4>, 'STRESS_XX': <ComputedDoubleValues.STRESS_XX: 5>, 'STRESS_YY': <ComputedDoubleValues.STRESS_YY: 6>, 'STRESS_ZZ': <ComputedDoubleValues.STRESS_ZZ: 7>, 'STRESS_XY': <ComputedDoubleValues.STRESS_XY: 8>, 'STRESS_XZ': <ComputedDoubleValues.STRESS_XZ: 10>, 'STRESS_YZ': <ComputedDoubleValues.STRESS_YZ: 9>, 'STRESS_NXY': <ComputedDoubleValues.STRESS_NXY: 11>, 'STRESS_NXZ': <ComputedDoubleValues.STRESS_NXZ: 13>, 'STRESS_NYZ': <ComputedDoubleValues.STRESS_NYZ: 12>, 'MEAN_B': <ComputedDoubleValues.MEAN_B: 14>, 'MAX_B': <ComputedDoubleValues.MAX_B: 15>, 'MSD': <ComputedDoubleValues.MSD: 16>}

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...


class ComputedIntValues:
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
    NUM_ATOMS: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_ATOMS: 3>
    # value = <ComputedIntValues.NUM_BONDS: 5>
    NUM_BONDS: typing.ClassVar[ComputedIntValues]
    # value = <ComputedIntValues.NUM_BONDS_TO_FORM: 7>
    NUM_BONDS_TO_FORM: typing.ClassVar[ComputedIntValues]
    # value = <ComputedIntValues.NUM_EXTRA_ATOMS: 4>
    NUM_EXTRA_ATOMS: typing.ClassVar[ComputedIntValues]
    # value = <ComputedIntValues.NUM_EXTRA_BONDS: 6>
    NUM_EXTRA_BONDS: typing.ClassVar[ComputedIntValues]
    # value = <ComputedIntValues.NUM_RELOC: 2>
    NUM_RELOC: typing.ClassVar[ComputedIntValues]
    # value = <ComputedIntValues.NUM_SHIFT: 1>
    NUM_SHIFT: typing.ClassVar[ComputedIntValues]
    # value = <ComputedIntValues.STEP: 0>
    STEP: typing.ClassVar[ComputedIntValues]
    # value = {'STEP': <ComputedIntValues.STEP: 0>, 'NUM_SHIFT':
    # <ComputedIntValues.NUM_SHIFT: 1>, 'NUM_RELOC':
    # <ComputedIntValues.NUM_RELOC: 2>, 'NUM_ATOMS':
    # <ComputedIntValues.NUM_ATOMS: 3>, 'NUM_EXTRA_ATOMS':
    # <ComputedIntValues.NUM_EXTRA_ATOMS: 4>, 'NUM_BONDS':
    # <ComputedIntValues.NUM_BONDS: 5>, 'NUM_EXTRA_BONDS':
    # <ComputedIntValues.NUM_EXTRA_BONDS: 6>, 'NUM_BONDS_TO_FORM':
    # <ComputedIntValues.NUM_BONDS_TO_FORM: 7>}
    __members__: typing.ClassVar[dict[str, ComputedIntValues]]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...


class DPDSimulator:
    """

              A quick-and-dirty implementation of the DPD simulation
              with slip-springs as presented by Langeloth et al.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @staticmethod
    def read_restart_file(file: str) -> DPDSimulator:
        """
                  Read a restart file in order to continue a simulation.
        """

    def __init__(self, universe: Universe, crosslinker_type: int = 2,
                 slipspring_bond_type: int = 9, is_2d: bool = False, seed: str = '') -> None:
        """
        Get an instance of this class
        """

    def assume_box_large_enough(self) -> None:
        """
                  Configure whether to run PBC on the bonds or not.

                  If your bonds could get larger than half the box length, this must be kept false (default).
                  Otherwise, you can set it to true and therewith get some securities.
        """

    def configAutoCorrelatorOutput(
            self, values: list[OutputConfiguration], num_corr_in: int = 32, p: int = 16, m: int = 2) -> None:
        """
                  Set which values to compute multiple-tau autocorrelation for.
                  If you use this, you should cite `doi:10.1063/1.3491098 <https://pubs.aip.org/aip/jcp/article-abstract/133/15/154103/190247/Efficient-on-the-fly-calculation-of-time?redirectedFrom=fulltext>`_

                  Arguments:
                       - values: a list of OutputConfiguration structs
                       - ...
        """

    def configAverageOutput(self, arg0: list[OutputConfiguration]) -> None:
        """
                  Set which values to compute averages for.

                  Arguments:
                       - values: a list of OutputConfiguration structs
        """

    def config_a(self, A: float = 25.0) -> None:
        """
                  Configure the force-field (pair-style) parameter `A`.
        """

    def config_allow_relocation_in_network(
            self, allow_relocation_in_network: bool = False) -> None:
        """
                  Configure whether a relocation step may happen when a slip-spring has ended at a cross-link.

                  Side-effect: if true, the relocations may also happen *to* a slip-spring next to a cross-link.

                  Arguments:
                  - allow_relocation_in_network (bool): Whether to allow relocation in the network or not.
        """

    def config_bond_formation(self, num_bonds_to_form: int, max_bonds_per_atom_type: dict[int, int], bond_formation_dist: float = 1.0,
                              attempt_bond_formation_every: int = 50, atom_type_form_from: int = 2, atom_type_form_to: int = 1) -> None:
        """
                  Configure how to do bond formation during the run.

                  Arguments:
                  - num_bonds_to_form (int): the nr of bonds to form in total. Use 0 to stop bond formation.
                  - num_bonds_per_atom_type (dict): the nr of bonds each atom type may have at most (e.g., 2 for strand atoms, 4 for a tertiary crosslinkers)
                  - bond_formation_dist (float): the maximum distance allowed to form bonds
                  - attempt_bond_formation_every (int): attempt to form bonds every this many steps during the simulation run
                  - atom_type_form_from (int): the atom type to start forming bonds from.
                  - atom_type_form_to (int): the atom type to start forming bonds to.
        """

    def config_box_deformation(self, target_box: Box) -> None:
        """
                  Configure where to (incrementally) deform the box to during the next simulation run.
        """

    def config_lambda(self, l: float = 0.65) -> None:
        """
                  Configure the modified velocity verlet integration parameter `\\lambda`.
        """

    def config_num_steps_dpd(self, arg0: int) -> None:
        """
                  Configure the number of steps to do in one DPD sequence.
        """

    def config_num_steps_mc(self, arg0: int) -> None:
        """
                  Configure the number of steps to do in one MC sequence.
        """

    def config_restart_output(self, file: str,
                              output_every: int = 50000) -> None:
        """
                  Set when to output a restart where.

                  Note:
                       The filename determines the type of serialization:
                       .json, .xml are supported; other file endings will lead to binary serialization (fastest!).

                  Caution:
                       This method may not be backwards- nor forward-compatible.
                       Use the same version of pylimer-tools if you want to be sure that things work.

                  Arguments:
                       - file: the file path to the restart file to write
                       - outputEvery: how often to write the restart file
        """

    def config_shift_one_at_a_time(self, arg0: bool) -> None:
        ...

    def config_shift_possibility_empty(self, arg0: bool) -> None:
        ...

    def config_sigma(self, sigma: float = 3.0) -> None:
        """
                  Configure the force-field (pair-style) parameter `\\sigma`.
        """

    def config_slipspring_high_cutoff(self, cutoff: float = 2.0) -> None:
        """
                  Configure the lower cut-off of how far a pair may be distanced for a slip-spring to be created.
        """

    def config_slipspring_low_cutoff(self, cutoff: float = 0.5) -> None:
        """
                  Configure the higher cut-off of how far a pair may be distanced for a slip-spring to be created.
        """

    def config_spring_constant(self, k: float = 2.0) -> None:
        """
                  Configure the force-field (bond-style) parameter `k`, the spring constant.
        """

    def config_step_output(self, arg0: list[OutputConfiguration]) -> None:
        """
                  Set which values to log.

                  Arguments:
                       - values: a list of OutputConfiguration structs
        """

    def create_slip_springs(self, num: int, bond_type: int = 9) -> int:
        """
                  Randomly add the specified number of slip-springs to neighbours within the specified cut-offs.
        """

    def get_bond_lengths(self) -> numpy.ndarray:
        ...

    def get_coordinates(self) -> numpy.ndarray:
        ...

    def get_current_timestep(self) -> int:
        ...

    def get_nr_of_atoms(self) -> int:
        ...

    def get_nr_of_bonds(self) -> int:
        ...

    def get_nr_of_bonds_to_form(self) -> int:
        """
                  Get the number of bonds that are configured to have to be formed.
        """

    def get_nr_of_extra_atoms(self) -> int:
        ...

    def get_nr_of_extra_bonds(self) -> int:
        ...

    def get_nr_of_slip_springs(self) -> int:
        ...

    def get_nr_of_steps_dpd(self) -> int:
        ...

    def get_nr_of_steps_mc(self) -> int:
        ...

    def get_shift_one_at_a_time(self) -> bool:
        ...

    def get_shift_possibility_empty(self) -> bool:
        ...

    def get_slip_spring_bond_type(self) -> int:
        ...

    def get_spring_constant(self) -> float:
        ...

    def get_stress_tensor(self) -> numpy.ndarray:
        ...

    def get_temperature(self) -> float:
        ...

    def get_timestep(self) -> float:
        ...

    def get_universe(self, with_slipsprings: bool = True) -> Universe:
        """
             Get a universe instance from the current coordinates (and connectivity).

             Arguments:
                  - with_slip_springs (bool): whether to include slip-springs in the returned universe.
        """

    def get_volume(self) -> float:
        ...

    def refresh_current_state(self) -> None:
        """
                  After re-configuring the force-field parameters,
                  this method should be called to update the current stress tensor etc.
        """

    def run_simulation(self, n_steps: int, dt: float = 0.06,
                       with_MC: bool = False) -> None:
        ...

    def start_measuring_msd_for_atoms(self, atom_ids: list[int]) -> None:
        """
                  Set a new origin for measuing the mean square displacement for a specified set of atoms
        """

    def validate_neighbour_list(self, arg0: float) -> None:
        ...

    def validate_state(self) -> None:
        ...

    def write_restart_file(self, file: str) -> None:
        """
                  Explicitily force the writing of a restart file, now!

                  Arguments:
                  - file (str): the file path and name of the restart file to be written.
                       Can end in xml, json or anything else (-> binary).
        """


class DataFileReader:
    """

           A reader for LAMMPS's `write_data` files.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self) -> None:
        ...

    def get_atom_ids(self) -> list[int]:
        ...

    def get_atom_nx(self) -> list[int]:
        ...

    def get_atom_ny(self) -> list[int]:
        ...

    def get_atom_nz(self) -> list[int]:
        ...

    def get_atom_types(self) -> list[int]:
        ...

    def get_atom_x(self) -> list[float]:
        ...

    def get_atom_y(self) -> list[float]:
        ...

    def get_atom_z(self) -> list[float]:
        ...

    def get_bond_from(self) -> list[int]:
        ...

    def get_bond_to(self) -> list[int]:
        ...

    def get_bond_types(self) -> list[int]:
        ...

    def get_high_x(self) -> float:
        ...

    def get_high_y(self) -> float:
        ...

    def get_high_z(self) -> float:
        ...

    def get_low_x(self) -> float:
        ...

    def get_low_y(self) -> float:
        ...

    def get_low_z(self) -> float:
        ...

    def get_lx(self) -> float:
        ...

    def get_ly(self) -> float:
        ...

    def get_lz(self) -> float:
        ...

    def get_masses(self) -> dict[int, float]:
        ...

    def get_molecule_ids(self) -> list[int]:
        ...

    def get_nr_of_atom_types(self) -> int:
        ...

    def get_nr_of_atoms(self) -> int:
        ...

    def get_nr_of_bond_types(self) -> int:
        ...

    def get_nr_of_bonds(self) -> int:
        ...

    def read(self, path_of_file_to_read: str, atom_style: AtomStyle = ...,
             atom_style2: AtomStyle = ..., atom_style_3: AtomStyle = ...) -> None:
        """
               Actually read a LAMMPS's `write_data` file.

               Arguments:
                  - `path_of_file_to_read`: The path to the file to read
                  - `atom_style`: The format of the "Atoms" section, see https://docs.lammps.org/read_data.html
                  - `atom_style2`: The format of the "Atoms" section if the previous parameter is equal to AtomStyle::HYBRID
                  - `atom_style3`: The format of the "Atoms" section if the second to last parameter is equal to AtomStyle::HYBRID
        """


class DataFileWriter:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, universe: Universe) -> None:
        """
                   Initialize the writer with the universe to write.
        """

    def config_atom_style(self, arg0: AtomStyle) -> None:
        """
                  Set the (LAMMPS) atom style to use for writing the atoms.
        """

    def config_attempt_image_reset(
            self, attempt_image_reset: bool = True) -> None:
        """
                   Set whether to change the outuput coordinates to lie in the box or not.

                   Default: false.
        """

    def config_crosslinker_type(self, crosslinker_type: int = 2) -> None:
        """
                   Set which atom type represents cross-linkers.
                   Needed in case the moleculeIdx in the output file should have any meaning.
                   (e.g. with :func:`~pylimer_tools_cpp.DataFileWriter.configMoleculeIdxForSwap`).

                   Default: 2.
        """

    def config_include_angles(self, include_angles: bool = True) -> None:
        """
                   Set whether to include the angles from the universe in the file or not.

                   Default: true.
        """

    def config_include_dihedral_angles(
            self, include_dihedral_angles: bool = True) -> None:
        """
                   Set whether to include the dihedral angles from the universe in the file or not.

                   Default: true.
        """

    def config_include_velocities(
            self, include_velocities: bool = True) -> None:
        """
                   Set whether to include the velocities from the universe (if any) in the file or not.

                   Default: true.
        """

    def config_molecule_idx_for_swap(
            self, enableSwappability: bool = True) -> None:
        """
                        Swappable chains implies that their `moleculeIdx` in the LAMMPS data file is not
                        identical per chain, but identical per position in the chain.
                        That's how you can have bond swapping with constant chain length distribution.

                        Default: false.
        """

    def config_move_into_box(self, move_into_box: bool = True) -> None:
        """
                   Set whether to change the outuput coordinates to lie in the box or not.

                   Default: false (used to be true).
        """

    def config_reindex_atoms(self, reindex_atoms: bool = True) -> None:
        """
                   Set whether to reindex the atoms or not.
                   Re-indexing leads to atom ids being in the range of 1 to the number of atoms.

                   Default: false.
        """

    def set_custom_atom_format(
            self, atom_format: str = '\t$atomId\t$moleculeId\t$atomType\t$x\t$y\t$z\t$nx\t$ny\t$nz') -> None:
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
                  :func:`~pylimer_tools_cpp.Universe.setPropertyValue`
                  as placeholders (as long as they are alphanumeric only; prefix in the format with '$' as well).
                  Specifically useful if you need a different (or hybrid) atom style in LAMMPS.

                  Be sure to still call :func:`~pylimer_tools_cpp.DataFileWriter.configAtomStyle`,
                  so that the file can be read correctly again.
        """

    def set_universe_to_write(self, universe: Universe) -> None:
        """
                   Re-set the universe to write.
        """

    def write_to_file(self, file: str) -> None:
        """
                  Actually do the writing to the disk.

                  Arguments:
                       - file (str): the path and file name to write to
        """


class DumpFileReader:
    """

           A reader for LAMMPS's `dump` files.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, path_of_file_to_read: str) -> None:
        ...

    def get_length(self) -> int:
        """
        Get the number of sections (time-steps) in the file
        """

    def get_numeric_values_for_at(
            self, arg0: int, arg1: str, arg2: str) -> list[float]:
        """
        Get the values for the section `index`, the main header `headerKey` and the column (in the header) `column`.
        """

    def get_string_values_for_at(
            self, rowIndex: int, headerKey: str, columnIndex: str) -> list[str]:
        """
        Get the values for the section `index`, the main header `headerKey` and the column (in the header) `column`.
        """

    def has_key(self, headerKey: str) -> bool:
        """
        Check whether the first section has the header specified
        """

    def key_has_column(self, headerKey: str, columnName: str) -> bool:
        """
        Check whether the header of the first section has the specified column
        """

    def key_has_directional_column(
            self, headerKey: str, dirPraefix: str = '', dirSuffix: str = '') -> bool:
        """
        Check whether the header of the first section has all the three columns `{dirPraefix}{x|y|z}{dirSuffix}`.
        """

    def read(self) -> None:
        """
        Read the whole file
        """


class ExitReason:
    """
    Members:

      UNSET

      MAX_STEPS

      F_TOLERANCE

      X_TOLERANCE

      FAILURE

      OTHER
    """
    FAILURE: typing.ClassVar[ExitReason]  # value = <ExitReason.FAILURE: 5>
    # value = <ExitReason.F_TOLERANCE: 1>
    F_TOLERANCE: typing.ClassVar[ExitReason]
    MAX_STEPS: typing.ClassVar[ExitReason]  # value = <ExitReason.MAX_STEPS: 3>
    OTHER: typing.ClassVar[ExitReason]  # value = <ExitReason.OTHER: 7>
    UNSET: typing.ClassVar[ExitReason]  # value = <ExitReason.UNSET: 0>
    # value = <ExitReason.X_TOLERANCE: 2>
    X_TOLERANCE: typing.ClassVar[ExitReason]
    # value = {'UNSET': <ExitReason.UNSET: 0>, 'MAX_STEPS':
    # <ExitReason.MAX_STEPS: 3>, 'F_TOLERANCE': <ExitReason.F_TOLERANCE: 1>,
    # 'X_TOLERANCE': <ExitReason.X_TOLERANCE: 2>, 'FAILURE':
    # <ExitReason.FAILURE: 5>, 'OTHER': <ExitReason.OTHER: 7>}
    __members__: typing.ClassVar[dict[str, ExitReason]]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...


class LazyUniverseSequenceIterator:
    """

           An iterator to iterate throught the universes in :obj:`~pylimer_tools_cpp.UniverseSequence`.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __iter__(self) -> LazyUniverseSequenceIterator:
        ...

    def __next__(self) -> Universe:
        ...


class LinkSwappingMode:
    """
    Members:

      NO_SWAPPING

      SLIPLINKS_ONLY

      ALL

      ALL_CYCLE

      ALL_MC

      ALL_MC_CYCLE

      ALL_MC_TRY

      ALL_MC_TRY_CYCLE
    """
    ALL: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL: 2>
    # value = <LinkSwappingMode.ALL_CYCLE: 3>
    ALL_CYCLE: typing.ClassVar[LinkSwappingMode]
    # value = <LinkSwappingMode.ALL_MC: 4>
    ALL_MC: typing.ClassVar[LinkSwappingMode]
    # value = <LinkSwappingMode.ALL_MC_CYCLE: 5>
    ALL_MC_CYCLE: typing.ClassVar[LinkSwappingMode]
    # value = <LinkSwappingMode.ALL_MC_TRY: 6>
    ALL_MC_TRY: typing.ClassVar[LinkSwappingMode]
    # value = <LinkSwappingMode.ALL_MC_TRY_CYCLE: 7>
    ALL_MC_TRY_CYCLE: typing.ClassVar[LinkSwappingMode]
    # value = <LinkSwappingMode.NO_SWAPPING: 0>
    NO_SWAPPING: typing.ClassVar[LinkSwappingMode]
    # value = <LinkSwappingMode.SLIPLINKS_ONLY: 1>
    SLIPLINKS_ONLY: typing.ClassVar[LinkSwappingMode]
    # value = {'NO_SWAPPING': <LinkSwappingMode.NO_SWAPPING: 0>,
    # 'SLIPLINKS_ONLY': <LinkSwappingMode.SLIPLINKS_ONLY: 1>, 'ALL':
    # <LinkSwappingMode.ALL: 2>, 'ALL_CYCLE': <LinkSwappingMode.ALL_CYCLE: 3>,
    # 'ALL_MC': <LinkSwappingMode.ALL_MC: 4>, 'ALL_MC_CYCLE':
    # <LinkSwappingMode.ALL_MC_CYCLE: 5>, 'ALL_MC_TRY':
    # <LinkSwappingMode.ALL_MC_TRY: 6>, 'ALL_MC_TRY_CYCLE':
    # <LinkSwappingMode.ALL_MC_TRY_CYCLE: 7>}
    __members__: typing.ClassVar[dict[str, LinkSwappingMode]]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...


class MCUniverseGenerator:
    """

           A :obj:`pylimer_tools_cpp.Universe` generator using a Monte-Carlo procedure.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, lx: float, ly: float, lz: float) -> None:
        ...

    def add_crosslinkers(self, nr_of_crosslinkers: int, crosslinker_functionality: int = 4,
                         crosslinker_type: int = 2, white_noise: bool = True) -> None:
        """
                    Add cross-linkers at random positions.

                    :param nr_of_crosslinkers: Number of cross-linkers to add.
                    :param crosslinker_functionality: Functionality of the cross-linkers (default: 4).
                    :param crosslinker_type: Atom type of the cross-linkers (default: 2).
                    :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the cross-linkers (default: true).
        """

    def add_crosslinkers_at(self, coordinates: numpy.ndarray,
                            crosslinker_functionality: int = 4, crosslinker_type: int = 2) -> None:
        """
                  Add cross-linkers at specific coordinates.

                  :param coordinates: Coordinates of the cross-linkers.
                  :param crosslinker_functionality: Functionality of the cross-linkers (default: 4).
                  :param crosslinker_type: Atom type of the cross-linkers (default: 2).
        """

    def add_end_functionalized_strands(
            self, nr_of_strands: int, strand_length: list[int], crosslinker_functionality: int = 4, crosslinker_type: int = 2, strand_atom_type: int = 1, white_noise: bool = True) -> None:
        """
                    Add strands with cross-linkers at the ends.

                    :param nr_of_strands: Number of strands to add.
                    :param strand_length: Length of each strand.
                    :param crosslinker_functionality: Functionality of the cross-linkers (default: 4).
                    :param crosslinker_type: Atom type of the cross-linkers (default: 2).
                    :param strand_atom_type: Atom type of the beads that are not at the ends (default: 1).
                    :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the cross-linkers (default: true).
        """

    def add_monofunctional_strands(self, nr_of_monofunctional_strands: int,
                                   monofunctional_strand_length: list[int], monofunctional_strand_atom_type: int = 4) -> None:
        """
                 Add multiple monofunctional strands with specified bead types.
        """

    def add_randomly_functionalized_strands(self, nr_of_strands: int, strand_length: list[int], functionalization_probability: float,
                                            crosslinker_functionality: int = 4, crosslinker_type: int = 2, strand_atom_type: int = 1, white_noise: bool = True) -> None:
        """
                  Add strands with randomly distributed cross-linkers in between.

                  :param nr_of_strands: Number of strands to add.
                  :param strand_length: Length of each strand.
                  :param functionalization_probability: Proportion of beads that are made cross-link.
                  :param crosslinker_functionality: Functionality of the cross-linkers (default: 4).
                  :param crosslinker_type: Atom type of the cross-linkers (default: 2).
                  :param strand_atom_type: Atom type of the beads that stay (default: 1).
                  :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the cross-linkers (default: true).
        """

    def add_solvent_chains(self, nr_of_solvent_chains: int, solvent_chain_length: int,
                           solvent_atom_type: int = 3, white_noise: bool = True) -> None:
        """
                    Randomly distribute additional, free chains.
        """

    def add_strands(self, nr_of_strands: int,
                    strand_lengths: list[int], strand_atom_type: int = 1) -> None:
        """
                    Add strands.
                    Adds them unconnected at first.
                    To link them to cross-linkers, use some of the :code:`link_strand_` methods.

                    :param nr_of_strands: Number of strands to add.
                    :param strand_lengths: A list of integers representing the number of beads of each of the strands.
                    :param strand_atom_type: Type of atoms for the strands (default: 1).
        """

    def config_nr_of_mc_steps(self, n_steps: int = 2000) -> None:
        """
        Set the number of Monte-Carlo steps during bond length equilibration.
        """

    def config_primary_loop_probability(
            self, probability: float = 1.0) -> None:
        """
                 Configure an additional weight reduction for the primary loop formation.

                 Defaults to 1., which means the general :math:`P(\\vec{R})` is used without any bias.
                 This results in more primary loops for shorter chains than for longer ones.

                 Set to 0. to disable the formation of primary loops.
        """

    def config_secondary_loop_probability(
            self, probability: float = 1.0) -> None:
        """
                 Configure an additional weight reduction for the secondary loop formation.

                 Defaults to 1., which means the general :math:`P(\\vec{R})` is used without any bias.
                 This results in more secondary loops for shorter chains than for longer ones.

                 Set to 0. to disable the formation of secondary loops.
        """

    def get_mean_bead_distance(self) -> float:
        """
        Get the currently configured mean bead distance.
        """

    def get_mean_squared_bead_distance(self) -> float:
        """
        Get the currently configured mean squared bead distance.
        """

    def get_nr_of_atoms(self) -> int:
        """
        Get the current number of atoms that the universe would/will have.
        """

    def get_nr_of_bonds(self) -> int:
        """
        Get the current number of bonds that the universe would/will have.
        """

    def get_universe(self) -> Universe:
        """
                    Fetch the current (or final) state of the universe.

                    Use this method to actually (MC) place beads between the cross-links and retrieve the generated structure.
        """

    def link_strand(self, strand_idx: int, c_infinity: float = 1.0) -> None:
        """
                  Link one particular strand to a cross-linker.
                  This strand will have one bond made, to an appropriate cross-linker,
                  as chosen by the parameters associated with the strand.

                  :param strand_idx: Index of the strand to be linked.
                  :param c_infinity: As needed for the end-to-end distribution, given by :math:`\\langle R^2\\rangle_0 = C_{\\infty} N b^2`.
        """

    def link_strand_to(self, strand_idx: int, link_idx: int) -> None:
        """
                  Link a strand to a specific cross-linker.
                  This assumes that you keep track of the order in which you added the cross-linkers and strands.

                  :param strand_idx: Index of the strand to be linked.
                  :param link_idx: Index of the cross-linker to be linked.
        """

    def link_strands_to_conversion(
            self, crosslinker_conversion: float, c_infinity: float = 1.0) -> None:
        """
                    Actually link the previously added strands to the previously added cross-linkers,
                    until a certain cross-link conversion is reached.

                    :param crosslinker_conversion: Target conversion of cross-linkers (0: no connections to cross-links; 1: all cross-linkers fully connected).
                    :param c_infinity: As needed for the end-to-end distribution, given by :math:`\\langle R^2\\rangle_0 = C_{\\infty} N b^2`.
        """

    def link_strands_to_soluble_fraction(
            self, soluble_fraction: float, c_infinity: float = 1.0) -> None:
        """
                    Actually link the previously added strands to the previously added cross-linkers,
                    until a certain soluble fraction is reached.

                    :param soluble_fraction: Target soluble_fraction (0: no connections to cross-links; 1: all cross-linkers fully connected).
                    :param c_infinity: As needed for the end-to-end distribution, given by :math:`\\langle R^2\\rangle_0 = C_{\\infty} N b^2`.
        """

    def relax_crosslinks(self) -> None:
        """
                 Run force relaxation with the cross-linkers and their strands,
                 to have the cross-links in their statistically most probable position.
        """

    def remove_soluble_fraction(self, rescale_box: bool = True) -> None:
        """
                    Remove soluble fraction (as determined by phantom force relaxation) of the strands and cross-links.
        """

    def set_bead_distance(self, distance: float,
                          update_mean_squared: bool = True) -> None:
        """
                 Set the mean distance between beads when doing MC stepping.
                 Also used for the target cross-linker partner sampling.

                 NOTE: Mainly the mean squared bead distance is effectively used in the Monte-Carlo simulation.

                 :param distance: Mean distance between beads.
                 :param update_mean_squared: Whether to update the mean squared distance as well, deduced from the assumed gaussian distribution in 3D (default: true).
        """

    def set_mean_squared_bead_distance(
            self, mean_squared_distance: float, update_mean: bool = True) -> None:
        """
        Set the mean squared distance between beads.
                 :param mean_squared_distance: Mean squared distance between beads.
                 :param update_mean: Whether to update the mean bead distance as well, deduced from the assumed gaussian distribution in 3D (default: true).
        """

    def set_seed(self, seed: int) -> None:
        """
        Set the seed for the random generator.
        """

    def validate(self) -> None:
        """
                 Check whether the internal state of the generator is valid.
                 Throws errors if not.
                 Should in principle always be valid when called from Python; if not, there is a bug in the code.
        """


class MEHPForceBalance:
    """

        A small simulation tool for quickly minimizing the force between the cross-linker beads.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @staticmethod
    def construct_with_random_sliplinks(universe: Universe, nr_of_sliplinks_to_sample: int, upper_sampling_cutoff: float = 1.2, lower_sampling_cutoff: float = 0.0,
                                        minimum_nr_of_sliplinks: int = 0, same_strand_cutoff: float = 3, seed: str = '', crosslinker_type: int = 2, is_2d: bool = False) -> MEHPForceBalance:
        """
                  Instantiate this simulator with randomly chosen slip-links.
        """

    def __copy__(self) -> MEHPForceBalance:
        ...

    def __getstate__(self) -> tuple:
        ...

    def __init__(self, universe: Universe, crosslinker_type: int = 2, is_2d: bool = False,
                 remove_2functional_crosslinkers: bool = False, remove_dangling_chains: bool = False) -> None:
        """
                  Instantiate the simulator for a certain universe.

                  :param universe: the universe to simulate with
                  :param crosslinker_type: The atom type of the cross-linkers. Needed to reduce the network.
                  :param is2D: Whether to ignore the z direction.
                  :param kappa: the spring constant
                  :param remove_2functionalCrosslinkers: whether to keep or remove the 2-functional crosslinkers when setting up the network
                  :param remove_dangling_chains: whether to keep or remove obviously dangling chains when setting up the network
        """

    def __setstate__(self, arg0: tuple) -> None:
        ...

    def add_sliplinks(self, strand_idx_1: list[int], strand_idx_2: list[int], x: list[float], y: list[float],
                      z: list[float], alpha_1: list[float], alpha_2: list[float], clamp_alpha: bool = False) -> None:
        """
                  Add new slip-links
        """

    def add_sliplinks_based_on_cycles(self, maxLoopLength: int = -1) -> int:
        """
                  Detect and add slip-links based on detected entanglements.

                  WARNING:
                       Does not work yet.
        """

    def config_assume_box_large_enough(
            self, box_large_enough: bool = False) -> None:
        """
                  Configure whether to run PBC on the bonds or not.

                  If your bonds could get larger than half the box length, this must be kept false (default).
                  Otherwise, you can set it to true and therewith get some securities.
        """

    def config_entanglement_type(self, type: int = -1) -> None:
        """
                 To have certain cross-links behave as entanglements in the removal process,
                 you can specify the here a type, that you have used in the universe to specify:
                 - the type of entanglement atoms (expected with functionality f = 3),
                 - and the entanglement-bonds between the entanglement atoms.

                 I.e., say you want to model some entanglements as non-slipping,
                 bonds between two strand beads resulting in f = 3 beads, for example,
                 you can call this method to have the "StructureSimplificationMode" also remove these atoms,
                 if they have a functionality of 2 or less while still being connected to its partner bead.
        """

    def config_mean_bond_length(self, b: float = 1.0) -> None:
        """
             Configure the :math:`b` used e.g. for the topological Gamma-factor.
        """

    def config_simplification_frequency(self, frequency: int = 10) -> None:
        """
                 Config every how many steps to simplify the structure.
                 Default: 10.
        """

    def config_spring_breaking_distance(
            self, distance_over_contour_length: float = -1) -> None:
        """
                  Configure the "force" (distance over contour length) at which the bonds break.
                  Can be used to model the effect of fracture, to reduce the stiffening happening upon deformation.
                  Springs breaking will happen before the simplification procedure is run.
                  Negative values will disable spring breaking.
                  Default: -1..
        """

    def config_spring_constant(self, kappa: float = 1.0) -> None:
        ...

    def config_step_output(self, arg0: list[OutputConfiguration]) -> None:
        """
                  Set which values to log.

                  Arguments:
                       - values: a list of OutputConfiguration structs
        """

    def deform_to(self, new_box: Box) -> None:
        """
                  Perform a deformation of the system box to a different box.
                  All coordinates etc. will be scaled as needed.
        """

    def evaluate_partial_spring_distance(
            self, network: SimplifiedBalanceNetwork, displacements: numpy.ndarray, spring_idx: int) -> numpy.ndarray:
        ...

    def evaluate_partial_spring_distance_from(
            self, network: SimplifiedBalanceNetwork, displacements: numpy.ndarray, spring_idx: int, link_idx: int) -> numpy.ndarray:
        ...

    def evaluate_partial_spring_distance_to(self, network: SimplifiedBalanceNetwork,
                                            displacements: numpy.ndarray, spring_idx: int, link_idx: int) -> numpy.ndarray:
        ...

    def get_average_spring_length(self) -> float:
        """
                   Get the average length of the springs. Note that in contrast to :func:`~pylimer_tools_cpp.MEHPForceBalance.getGammaFactor()`,
                   this value is normalized by the number of springs rather than the number of chains.
        """

    def get_crosslinker_universe(self) -> Universe:
        """
                  Returns the universe [of cross-linkers] with the positions of the current state of the simulation.
        """

    def get_current_partial_spring_lengths(self) -> list[float]:
        """
                  Get the partial spring distances.
        """

    def get_current_partial_spring_vectors(self) -> numpy.ndarray:
        """
                  Get the partial spring vectors.
        """

    def get_current_spring_vectors(self) -> numpy.ndarray:
        ...

    def get_dangling_weight_fraction(self, tolerance: float = 0.001) -> float:
        """
                  Compute the weight fraction of non-active springs

                  Caution: ignores atom masses.
        """

    def get_default_mean_bond_length(self) -> float:
        """
                   Returns the value effectively used in :func:`~pylimer_tools_cpp.MEHPForceBalance.getGammaFactor()` for
                   :math:`b` in :math:`\\langle R_{0,\\eta}^2 = N_{\\eta} b^2\\rangle`.
        """

    def get_displacement_residual_norm(
            self, one_over_spring_partition_upper_limit: float = 1.0) -> float:
        """
                  Get the current link displacement residual norm.
        """

    def get_displacements(self) -> numpy.ndarray:
        """
                  Get the current link displacements.
        """

    def get_effective_functionality_of_atoms(
            self, tolerance: float = 0.001) -> dict[int, int]:
        """
                  Returns the number of active springs connected to each atom, atomId used as index

                  :param tolerance: springs under this length are considered inactive
        """

    def get_exit_reason(self) -> ExitReason:
        """
                   Returns the reason for termination of the simulation
        """

    def get_force_magnitude_vector(self, arg0: float) -> numpy.ndarray:
        """
                  Evaluate the norm of the force on each (slip- or cross-) link.
        """

    def get_force_on(self, link_idx: int,
                     one_over_spring_partition_upper_limit: float = 1.0) -> numpy.ndarray:
        """
                  Evaluate the force on a particular (slip- or cross-) link.
        """

    def get_gamma_factor(self, b02: float = -1.0, nr_of_chains: int = -1,
                         one_over_spring_partition_upper_limit: float = 1.0) -> float:
        """
                  Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

                  :math:`\\Gamma = \\langle\\gamma_{\\eta}\\rangle`, with :math:`\\gamma_{\\eta} = \\frac{\\bar{r_{\\eta}}^2}{R_{0,\\eta}^2}`,
                  which you can use as :math:`G_{\\mathrm{ANT}} = \\Gamma \\nu k_B T`,
                  where :math:`\\eta` is the index of a particular strand,
                  :math:`R_{0}^2` is the melt mean square end to end distance, in phantom systems :math:`$= N_{\\eta}*b^2$`
                  :math:`N_{\\eta}` is the number of atoms in this strand :math:`\\eta`,
                  :math:`b` its mean square bond length,
                  :math:`T` the temperature and
                  :math:`k_B` Boltzmann's constant.

                  :param b02: the melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\\infinity b^2`.
                  :param nr_of_chains: the value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of springs.
        """

    def get_gamma_factors(
            self, b02: float, one_over_spring_partition_upper_limit: float = 1.0) -> numpy.ndarray:
        """
                  Evaluates the gamma factor for each strand (i.e., the squared distance divided by the contour length multiplied by b02)
        """

    def get_gamma_factors_in_dir(self, b02: float, direction: int,
                                 one_over_spring_partition_upper_limit: float = 1.0) -> numpy.ndarray:
        """
                       Evaluates the gamma factor for each strand in the specified direction (i.e., the squared distance divided by the contour length multiplied by b02)

                       :param b02: the melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\\infinity b^2`.
                       :param direction: the direction in which to compute the gamma factors (0: x, 1: y, 2: z)
        """

    def get_ids_of_active_nodes(self, tolerance: float = 0.001) -> list[int]:
        """
                  Get the atom ids of the nodes that are considered active.
                  Only cross-link ids are returned (not e.g. entanglement links).

                  :param tolerance: springs under this length are considered inactive. A node is active if it has > 1 active springs.
        """

    def get_neighbour_link_indices(
            self, network: SimplifiedBalanceNetwork, link_idx: int) -> list[int]:
        ...

    def get_nr_of_active_nodes(self, tolerance: float = 0.001) -> int:
        """
                  Get the number of active nodes (incl. entanglement nodes [atoms with type = entanglementType, present in the universe when creating this simulator],
                  excl. entanglement links [the slip-links created internally when e.g. constructing the simulator with random slip-links]).

                  :param tolerance: springs under this length are considered inactive. A node is active if it has > 1 active springs.
        """

    def get_nr_of_active_partial_springs(
            self, tolerance: float = 0.001) -> int:
        """
                   Get the number of active partial springs remaining after running the simulation.

                  :param tolerance: springs under this length are considered inactive
        """

    def get_nr_of_active_springs(self, tolerance: float = 0.001) -> int:
        """
                   Get the number of active springs remaining after running the simulation.

                  :param tolerance: springs under this length are considered inactive
        """

    def get_nr_of_active_springs_in_dir(
            self, direction: int, tolerance: float = 0.001) -> int:
        """
                        Get the number of active springs remaining after running the simulation.

                       :param direction: the direction in which to compute the active springs (0: x, 1: y, 2: z)
                       :param tolerance: springs under this length are considered inactive
        """

    def get_nr_of_atoms(self) -> int:
        ...

    def get_nr_of_bonds(self) -> int:
        ...

    def get_nr_of_extra_atoms(self) -> int:
        ...

    def get_nr_of_extra_bonds(self) -> int:
        ...

    def get_nr_of_intra_chain_sliplinks(self) -> int:
        ...

    def get_nr_of_iterations(self) -> int:
        """
                  Returns the number of iterations used for force relaxation so far.
        """

    def get_nr_of_nodes(self) -> int:
        """
                   Get the number of nodes (crosslinkers) considered in this simulation.
        """

    def get_nr_of_springs(self) -> int:
        """
                  Get the number of springs considered in this simulation.

                  :param tolerance: springs under this length are considered inactive
        """

    def get_overall_spring_lengths(self) -> list[float]:
        """
                  Get the sum of the lengths of the partial springs of each spring.
        """

    def get_pressure(self) -> float:
        """
                  Returns the pressure at the current state of the simulation.
        """

    def get_soluble_weight_fraction(self, tolerance: float = 0.001) -> float:
        """
                  Compute the weight fraction of springs connected to active
                  springs (any depth).

                  Caution: ignores atom masses.
        """

    def get_spring_partitions(self) -> numpy.ndarray:
        """
                  Get the current spring partitions (the fraction of the contour length associated with each partial spring).
        """

    def get_springpartition_indices_of_sliplink(
            self, network: SimplifiedBalanceNetwork, link_idx: int) -> list[int]:
        ...

    def get_stress_on(self, link_idx: int,
                      one_over_spring_partition_upper_limit: float = 1.0) -> numpy.ndarray:
        """
                  Evaluate the stress on a particular (slip- or cross-) link.
        """

    def get_stress_tensor(
            self, one_over_spring_partition_upper_limit: float = 1.0) -> numpy.ndarray:
        """
                  Returns the stress tensor at the current state of the simulation.

                  The units are in :math:`[\\text{units of }\\kappa]/[\\text{distance units}]`,
                  where the units of :math:`\\kappa` should be :math:`[\\text{force}]/[\\text{distance units}]^2`.
                  Make sure to multiply by :math:`\\kappa` or configure it appropriately.
        """

    def get_stress_tensor_link_based(
            self, one_over_spring_partition_upper_limit: float = 1.0, xlinks_only: bool = False) -> numpy.ndarray:
        """
                  Returns the stress tensor at the current state of the simulation.
        """

    def get_weighted_partial_spring_lengths(
            self, one_over_spring_partition_upper_limit: float = 1.0) -> numpy.ndarray:
        """
                  Get the current partial spring lengths (norm of vector) divided by the spring partition times the contour length.
        """

    def inspect_displacement_to_mean_position_update(
            self, link_idx: int, one_over_spring_partition_upper_limit: float = 1.0) -> numpy.ndarray:
        """
                  Helper method to debug and/or understand what happens to certain links when being displaced.
        """

    def inspect_parametrisation_optimsation_for_link(self, link_idx: int, displacements: numpy.ndarray, spring_partitions: numpy.ndarray, max_nr_of_steps: int = 100, alpha_tol: float = 1e-09,
                                                     min_nr_of_steps: int = 1, one_over_spring_partition_upper_limit: float = 1.0) -> tuple[numpy.ndarray, numpy.ndarray, int, float, float, float, float]:
        """
                  Helper method to debug and/or understand what happens to certain links
                  when being displaced and their partition updated.
        """

    def inspect_spring_partition_update(self, link_idx: int) -> numpy.ndarray:
        """
                  Helper method to debug and/or understand what happens to certain links
                  when the spring partition is being updated.
        """

    def move_sliplinks_to_their_best_branch(self, arg0: SimplifiedBalanceNetwork, arg1: numpy.ndarray,
                                            arg2: numpy.ndarray, arg3: float, arg4: int, arg5: bool, arg6: bool) -> None:
        ...

    def randomly_add_sliplinks(self, nr_of_sliplinks_to_sample: int, cutoff: float = 2.0, minimum_nr_of_sliplinks: int = 0,
                               same_strand_cutoff: float = 2, exclude_crosslinks: bool = False, seed: int = -1) -> int:
        """
                  Randomly sample and add slip-links based on certain criteria.
        """

    def run_force_relaxation(self, max_nr_of_steps: int = 250000, x_tolerance: float = 1e-12, initial_residual_norm: float = -1.0, simplification_mode: StructureSimplificationMode = ..., inactive_removal_cutoff: float = 0.001, do_inner_iterations: bool = False,
                             allow_sliplinks_to_pass_each_other: LinkSwappingMode = ..., swapping_frequency: int = 10, one_over_spring_partition_upper_limit: float = 1.0, nr_of_crosslink_swaps_allowed_per_sliplink: int = -1, disable_slipping: bool = False) -> None:
        """
                  Run the simulation.
                  Note that the final state of the minimization is persisted and reused if you use this method again.
                  This is useful if you want to run a global optimization first and add a local one afterwards.
                  As a consequence though, you cannot simply benchmark only this method; you must include the setup.

                  :param max_nr_of_steps: The maximum number of steps to do during the simulation.
                  :param x_tolerance: The tolerance of the displacements as an exit condition.
                  :param initial_residual_norm: The residual norm relative to which the relative tolerance is specified. Negative values mean, it will be replaced with the current one.
                  :param simplification_mode: How to simplify the structure during the minimization.
                  :param inactive_removal_cutoff: The tolerance in distance units of the partial spring length to count as active, relevant if simplification mode is specified to be something other than NO_SIMPLIFICATION.
                  :param do_inner_iterations: Whether to do inner iterations; usually, they are not helpful.
                  :param allow_sliplinks_to_pass_each_other: Whether slip-links can pass each other.
                  :param swapping_frequency: How often slip-links attempt to swap.
                  :param one_over_spring_partition_upper_limit: Super-secret parameter. Use 1, gradually increase (and then -1) if you want to publish.
                  :param nr_of_crosslink_swaps_allowed_per_sliplink: Use to steer whether slip-links can cross cross-links when swapping is enabled.
                  :param disable_slipping: Whether slip-links should be prohibited from slipping.
        """

    def set_displacements(self, arg0: numpy.ndarray) -> None:
        """
                  Set the current link displacements.
        """

    def set_spring_contour_lengths(self, arg0: numpy.ndarray) -> None:
        """
                  Set/overwrite the contour lengths.
        """

    def set_spring_partitions(self, arg0: numpy.ndarray) -> None:
        """
                  Set the current spring partitions.
        """

    def swap_sliplinks_incl_xlinks(self, arg0: SimplifiedBalanceNetwork,
                                   arg1: numpy.ndarray, arg2: numpy.ndarray, arg3: float, arg4: bool) -> None:
        ...

    def validate_network(self) -> bool:
        """
                  Validates the internal structures.

                  Throws an error if something is not ok.
                  Otherwise, it returns true.

                  Can be used e.g. as :code:`assert fb.validate_network()`.
        """
    @property
    def network(self) -> SimplifiedBalanceNetwork:
        ...


class MEHPForceEvaluator:
    """

         The base interface to change the way the force is evaluated during a MEHP run.

    """
    is_2d: bool

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self) -> None:
        ...

    def evaluate_stress_contribution(
            self, spring_distances: float, i: int, j: int, spring_index: int) -> float:
        """
                  An evaluation of the stress-contribution.

                  :param springDistances: the three coordinate differences for one spring.
                  :param i: the row index of the stress tensor
                  :param j: the column index of the stress tensor
        """
    @property
    def network(self) -> SimplifiedNetwork:
        ...


class MEHPForceRelaxation:
    """

        A small simulation tool for quickly minimizing the force between the cross-linker beads.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __getstate__(self) -> tuple:
        ...

    def __init__(self, universe: Universe, crosslinker_type: int = 2, is_2d: bool = False, force_evaluator: MEHPForceEvaluator = None,
                 kappa: float = 1.0, remove_2functional_crosslinkers: bool = False, remove_dangling_chains: bool = False) -> None:
        """
                  Instantiate the simulator for a certain universe.

                  :param universe: the universe to simulate with
                  :param crosslinker_type: The atom type of the cross-linkers. Needed to reduce the network.
                  :param is2d: Whether to ignore the z direction.
                  :param force_evaluator: The force evaluator to use
                  :param kappa: The spring constant
                  :param remove_2functional_crosslinkers: Whether to replace two-functional crosslinkers with a "normal" chain bead
                  :param remove_dangling_chains: Whether to remove dangling chains before running the simulation.
                       **Caution*: Removing the dangling chains will result in incorrect results fo the computation of
                       :func:`~pylimer_tools_cpp.MEHPForceRelaxation.getSolubleWeightFraction()` and
                       :func:`~pylimer_tools_cpp.MEHPForceRelaxation.getDanglingWeightFraction()`
        """

    def __setstate__(self, arg0: tuple) -> None:
        ...

    def assume_box_large_enough(self, box_large_enough: bool = False) -> None:
        """
                  Configure whether to run PBC on the bonds or not.

                  If your bonds could get larger than half the box length, this must be kept false (default).
                  Otherwise, you can set it to true and therewith get some securities.
        """

    def config_rerun_epsilon(self, epsilon: float = 0.001) -> None:
        """
                  Configure the offset from the lower and upper bounds for the simulation to suggest another run (
                       See: :func:`~pylimer_tools_cpp.MEHPForceRelaxation.requiresAnotherRun()`
                  ).
        """

    def config_step_output(
            self, output_configuration: list[OutputConfiguration]) -> None:
        """
                  Set which values to log.

                  Arguments:
                       - values: a list of OutputConfiguration structs
        """

    def count_active_clustered_atoms(self, tolerance: float = 0.001) -> float:
        """
                  Counts the active clustered atoms in the system.

                  :param tolerance: springs under this length are considered inactive.
        """

    def get_active_chains(self, tolerance: float = 0.001) -> list[Molecule]:
        """
                  Get the cross-linker chains that are active.
        """

    def get_average_spring_length(self) -> float:
        """
                   Get the average length of the springs. Note that in contrast to :func:`~pylimer_tools_cpp.MEHPForceRelaxation.getGammaFactor()`,
                   this value is normalized by the number of springs rather than the number of chains.
        """

    def get_crosslinker_universe(self) -> Universe:
        """
                  Returns the universe [of cross-linkers] with the positions of the current state of the simulation.
        """

    def get_dangling_weight_fraction(self, tolerance: float = 0.001) -> float:
        """
                  Compute the weight fraction of non-active springs

                  Caution: ignores atom masses.
        """

    def get_default_r0_square(self) -> float:
        """
                   Returns the value effectively used in :func:`~pylimer_tools_cpp.MEHPForceRelaxation.getGammaFactor()` for :math:`\\langle R_{0,\\eta}^2\\rangle`.
        """

    def get_effective_functionality_of_atoms(
            self, tolerance: float = 0.001) -> dict[int, int]:
        """
                  Returns the number of active springs connected to each atom, atomId used as index

                  :param tolerance: springs under this length are considered inactive
        """

    def get_exit_reason(self) -> ExitReason:
        """
                   Returns the reason for termination of the simulation
        """

    def get_force(self) -> float:
        """
                  Returns the force at the current state of the simulation.
        """

    def get_gamma_factor(self, b0_squared: float = -1.0,
                         nr_of_chains: int = -1) -> float:
        """
                  Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:

                  :math:`\\Gamma = \\langle\\gamma_{\\eta}\\rangle`, with :math:`\\gamma_{\\eta} = \\frac{\\bar{r_{\\eta}}^2}{R_{0,\\eta}^2}`,
                  which you can use as :math:`G_{\\mathrm{ANT}} = \\Gamma \\nu k_B T`,
                  where :math:`\\eta` is the index of a particular strand,
                  :math:`R_{0}^2` is the melt mean square end to end distance, in phantom systems :math:`= N_{\\eta} b^2$`,
                  :math:`N_{\\eta}` is the number of atoms in this strand :math:`\\eta`,
                  :math:`b` its mean square bond length,
                  :math:`\\nu` the volume fraction,
                  :math:`T` the temperature and
                  :math:`k_B` Boltzmann's constant.

                  :param b0_squared: Part of the denominator in the equation of :math:`\\Gamma`.
                       If :math:`-1.0` (default), the network is used for determination (which is not accurate), the system is assumed to be phantom.
                       For real systems, the value could be determined by :func:`~pylimer_tools_cpp.Universe.compute_mean_squared_end_to_end_distance()`
                       on the melt system, with subsequent division by the nr of bonds in the chain.
                  :param nr_of_chains: the value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of chains.
        """

    def get_gamma_factors(self, b0_squared: float = -1.0) -> numpy.ndarray:
        """
                  Computes the gamma factor for each spring as part of the ANT/MEHP formulism.

                  :math:`\\gamma_{\\eta} = \\frac{\\bar{r_{\\eta}}^2}{R_{0,\\eta}^2}`, with (here)
                  :math:`R_{0,\\eta}^2 = N_\\eta \\cdot ` the parameter `b0_squared`.
                  You can obtain this parameter e.g. by doing melt simulations at different lengths,
                  it's the slope you obtain.

                  :param b0_squared: Part of the denominator in the equation of :math:`\\Gamma`.
                       If :math:`-1.0` (default), the network is used for determination (which is not accurate), the system is assumed to be phantom.
                       For real systems, the value could be determined by :func:`~pylimer_tools_cpp.Universe.compute_mean_squared_end_to_end_distance()`
                       on the melt system, with subsequent division by the nr of bonds in the chain.

                  See also :func:`~pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor` for the mean of these.
        """

    def get_ids_of_active_nodes(self, tolerance: float = 0.001, minimum_nr_of_active_connections: int = 2,
                                maximum_nr_of_active_connections: int = -1) -> list[int]:
        """
                  Get the atom ids of the nodes that are considered active.

                  Arguments:
                   - :param tolerance: springs under this length are considered inactive. A node is active if it has > 2 active springs.
        """

    def get_nr_of_active_nodes(self, tolerance: float = 0.001,
                               minimumNrOfActiveConnections: int = 2, maximumNrOfActiveConnections: int = -1) -> int:
        """
                   Get the number of active nodes remaining after running the simulation.

                  :param tolerance: springs under this length are considered inactive.
                  :param minimumNrOfActiveConnections:  A node is active if it has equal or more than this number of active springs.
                  :param maximumNrOfActiveConnections:  A node is active if it has equal or less than this number of active springs.
                       Use a value < 0 to indicate that there is no maximum number of active connections.
                  :param usePartial: Whether to use the partial spring distances rather than the total (set to true if you want primary loop contributors)
        """

    def get_nr_of_active_springs(self, tolerance: float = 0.001) -> int:
        """
                   Get the number of active springs remaining after running the simulation.

                  :param tolerance: springs under this length are considered inactive
        """

    def get_nr_of_iterations(self) -> int:
        """
                  Returns the number of iterations used for force relaxation.
        """

    def get_nr_of_nodes(self) -> int:
        """
                   Get the number of nodes considered in this simulation.
        """

    def get_nr_of_springs(self) -> int:
        """
                  Get the number of springs considered in this simulation.

                  :param tolerance: springs under this length are considered inactive
        """

    def get_pressure(self) -> float:
        """
                  Returns the pressure at the current state of the simulation.
        """

    def get_residual_norm(self) -> float:
        """
                  Returns the residual norm at the current state of the simulation.
        """

    def get_residuals(self) -> numpy.ndarray:
        """
                  Returns the residuals at the current state of the simulation.
        """

    def get_soluble_weight_fraction(self, tolerance: float = 0.001) -> float:
        """
                  Compute the weight fraction of springs connected to active
                  springs (any depth).

                  Caution: ignores atom masses.
        """

    def get_spring_distances(self) -> numpy.ndarray:
        """
                  Get the current coordinate differences for all the springs.

                  Returns:
                       - distances: a vector of size 3*nrOfSprings, with each x, y, z values of the springs
        """

    def get_spring_lengths(self) -> numpy.ndarray:
        """
                  Get the current lengths for all the springs.

                  Returns:
                       - distances: a vector of size nrOfSprings, with each the norm of the distances
        """

    def get_stress_tensor(self) -> numpy.ndarray:
        """
                  Returns the stress tensor at the current state of the simulation.
        """

    def requires_another_run(self) -> bool:
        """
                  For performance reasons, the objective is only minimised within the distances of one box.
                  This means, that there is a possibility, e.g. for a single strand longer than two boxes,
                  that it would not be globally minimised.

                  If the final displacement of one of the atoms is close
                  (1e-3, configurable via :func:`~pylimer_tools_cpp.MEHPForceRelaxation.configRerunEpsilon()`)
                  to the imposed min/max, after minimizing,
                  this method would return true.
        """

    def run_force_relaxation(self, algorithm: str = 'LD_MMA', max_nr_of_steps: int = 250000,
                             x_tolerance: float = 1e-12, f_tolerance: float = 1e-09) -> None:
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

    def set_force_evaluator(self, force_evaluator: MEHPForceEvaluator) -> None:
        """
                  Reset the currently used force evaluator.
        """


class Molecule:
    """

           An (ideally) connected series of atoms/beads.

    """
    __hash__: typing.ClassVar[None] = None

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __contains__(self, arg0: Atom) -> bool:
        """
                  Check whether a particular atom is contained in this molecule.
        """

    def __copy__(self) -> Molecule:
        ...

    def __eq__(self, arg0: Molecule) -> bool:
        ...

    def __getitem__(self, arg0: int) -> Atom:
        """
               Access an atom by its vertex index.
        """

    def __init__(self, arg0: Box, arg1: igraph_s,
                 arg2: MoleculeType, arg3: dict[int, float]) -> None:
        ...

    def __iter__(self) -> MoleculeIterator:
        """
               Iterate through the atoms in this molecule.
               No specific order is guaranteed.
        """

    def __len__(self) -> int:
        """
               Get the number of atoms in this molecule.
        """

    def compute_bond_lengths(self) -> list[float]:
        """
        Computes the length :math:`b` of each bond in the molecule, respecting periodic boundaries.
        """

    def compute_end_to_end_distance(self) -> float:
        """
                    Compute the end-to-end distance (:math:`R_{ee}`) of this molecule.

                    CAUTION:
                       Returns 0.0 if the molecule does not have two or more atoms.
        """

    def compute_end_to_end_distance_with_derived_image_flags(self) -> float:
        """
                    Compute the end-to-end distance (:math:`R_{ee}`) of this molecule,
                    but ignoring the image flags attached to the atoms.
                    This only works for Molecules that can be lined up with
                    :func:`~pylimer_tools_cpp.Molecule.getAtomsLinedUp()`,
                    as it needs the atoms sorted such that the periodic box can still be respected somewhat.

                    CAUTION:
                       Returns 0.0 if the molecule does not have two or more atoms.
                       Requires bonds to be shorter than half the box length.
        """

    def compute_end_to_end_vector(self) -> numpy.ndarray:
        """
                    Compute the end-to-end vector (:math:`\\overrightarrow{R}_{ee}`) of this molecule.

                    CAUTION:
                       Returns 0.0 if the molecule does not have two or more atoms.
        """

    def compute_end_to_end_vector_with_derived_image_flags(
            self) -> numpy.ndarray:
        """
                    Compute the end-to-end vector (:math:`\\overrightarrow{R}_{ee}`) of this molecule,
                    but ignoring the image flags attached to the atoms.
                    This only works for Molecules that can be lined up with
                    :func:`~pylimer_tools_cpp.Molecule.getAtomsLinedUp()`,
                    as it needs the atoms sorted such that the periodic box can still be respected somewhat.

                    CAUTION:
                       Returns 0.0 if the molecule does not have two or more atoms.
                       Requires bonds to be shorter than half the box length.
        """

    def compute_radius_of_gyration(self) -> float:
        """
                    Computes the radius of gyration, :math:`R_g^2` of this molecule.

                    :math:`{R_g}^2 = \\frac{1}{M} \\sum_i m_i (r_i - r_{cm})^2`,
                    where :math:`M` is the total mass of the molecule, :math:`r_{cm}`
                    are the coordinates of the center of mass of the molecule and the
                    sum is over all contained atoms.
        """

    def compute_radius_of_gyration_with_derived_image_flags(self) -> float:
        """
                    Computes the radius of gyration, :math:`R_g^2` of this molecule,
                    but ignoring the image flags attached to the atoms.
                    This only works for Molecules that can be lined up with
                    :func:`~pylimer_tools_cpp.Molecule.getAtomsLinedUp()`,
                    as it needs the atoms sorted such that the periodic box can still be respected somewhat.
                    In other words, this function computes the radius of gyration
                    assuming the distance between two lined-up beads
                    is smaller than half the periodic box in each direction.

                    See also: :func:`~pylimer_tools_cpp.Molecule.computeRadiusOfGyration()`.
        """

    def compute_total_length(self) -> float:
        """
             Computes the sum of the lengths of all bonds.
             In most cases, this is equal to the contour length.
        """

    def compute_total_mass(self) -> float:
        """
                    Computes the total mass of this molecule.
        """

    def compute_total_vector(self, crosslinker_type: int = 2,
                             close_loop: bool = True) -> numpy.ndarray:
        """
                       Computes the sum of all bond vectors.
        """

    def compute_vector_from_to(self, atom_id_from: int, atom_id_to: int,
                               crosslinker_type: int = 2, require_order: bool = True) -> numpy.ndarray:
        """
                       Computes the sum of all bond vectors between two specified atoms.
        """

    def get_atom_by_id(self, atom_id: int) -> Atom:
        """
                    Get an atom by its id.
        """

    def get_atom_by_vertex_idx(self, vertex_idx: int) -> Atom:
        """
                    Get an atom for a specific vertex.
        """

    def get_atom_id_by_vertex_idx(self, vertex_idx: int) -> int:
        """
                    Get the id of the atom by the vertex id of the underlying graph.
        """

    def get_atom_types(self) -> list[int]:
        """
        Query all types (each one for each atom) ordered by atom vertex id.
        """

    def get_atoms(self) -> list[Atom]:
        """
                    Returns all atom objects enclosed in this molecule, ordered by vertex id.
        """

    def get_atoms_by_degree(self, degree: int) -> list[Atom]:
        """
                    Get the atoms that have the specified number of bonds.
        """

    def get_atoms_by_type(self, type: int) -> list[Atom]:
        """
                    Get the atoms with the specified type.
        """

    def get_atoms_connected_to(self, atom: Atom) -> list[Atom]:
        """
                    Get the atoms connected to a specified atom.

                    Internally uses :func:`~pylimer_tools_cpp.Molecule.getAtomsConnectedTo`
        """

    def get_atoms_connected_to_vertex(self, vertex_idx: int) -> list[Atom]:
        """
                    Get the atoms connected to a specified vertex id.
        """

    def get_atoms_lined_up(self, crosslinker_type: int = 2,
                           assumed_coordinates: bool = False, close_loop: bool = False) -> list[Atom]:
        """
                    Returns all atom objects enclosed in this molecule based on the connectivity.

                    This method works only for lone chains, atoms and loops,
                    as it throws an error if the molecule does not allow such a "line-up",
                    for example because of crosslinkers.

                    Use the `crosslinker_type` parameter to force the atoms in a primary loop
                    to start with the cross-link.
        """

    def get_bonds(self) -> dict[str, list[int]]:
        """
                    Get all bonds. Returns a dict with three properties: 'bond_from', 'bond_to' and 'bond_type'.
        """

    def get_edges(self) -> dict[str, list[int]]:
        """
                    Get all bonds. Returns a dict with three properties: 'edge_from', 'edge_to' and 'edge_type'.
                    The order is not necessarily related to any structural property.

                    NOTE:
                       The integer values returned refer to the vertex ids, not the atom ids.
                       Use :func:`~pylimer_tools_cpp.Molecule.getAtomIdByIdx` to translate them to atom ids, or
                       :func:`~pylimer_tools_cpp.Molecule.getBonds` to have that done for you.
        """

    def get_key(self) -> str:
        """
                    Get a unique identifier for this molecule.
        """

    def get_nr_of_atoms(self) -> int:
        """
        Counts and returns the number of atoms associated with this molecule.
        """

    def get_nr_of_bonds(self) -> int:
        """
        Counts and returns the number of bonds associated with this molecule.
        """

    def get_strand_ends(self, crosslinker_type: int = 2,
                        close_loop: bool = False) -> list[Atom]:
        """
                  Get the ends of the given strand (= molecule).
                  In case of a primary loop, the cross-link is returned, if there is one.
                  Use the argument `close_loop` to decide, whether this should be returned once or twice.

                  NOTE:
                       Currently only works for linear strands.
        """

    def get_strand_type(self) -> MoleculeType:
        """
                   Get the type of this molecule (see :obj:`~pylimer_tools_cpp.MoleculeType` enum).

                   Note that this type might be unset; currently, only
                   :func:`~pylimer_tools_cpp.Universe.get_chains_with_crosslinker` assigns them automatically.
        """

    def get_vertex_idx_by_atom_id(self, atom_id: int) -> int:
        """
        Get the vertex index of the underlying graph for an atom with a specified id.
        """


class MoleculeIterator:
    """

           An iterator to iterate throught the atoms in :obj:`~pylimer_tools_cpp.Molecule`.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __iter__(self) -> MoleculeIterator:
        ...

    def __next__(self) -> Atom:
        ...


class MoleculeType:
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
    DANGLING_CHAIN: typing.ClassVar[MoleculeType]  # value = <MoleculeType.DANGLING_CHAIN: 3>
    # value = <MoleculeType.FREE_CHAIN: 4>
    FREE_CHAIN: typing.ClassVar[MoleculeType]
    # value = <MoleculeType.NETWORK_STRAND: 1>
    NETWORK_STRAND: typing.ClassVar[MoleculeType]
    # value = <MoleculeType.PRIMARY_LOOP: 2>
    PRIMARY_LOOP: typing.ClassVar[MoleculeType]
    # value = <MoleculeType.UNDEFINED: 0>
    UNDEFINED: typing.ClassVar[MoleculeType]
    # value = {'UNDEFINED': <MoleculeType.UNDEFINED: 0>, 'NETWORK_STRAND':
    # <MoleculeType.NETWORK_STRAND: 1>, 'PRIMARY_LOOP':
    # <MoleculeType.PRIMARY_LOOP: 2>, 'DANGLING_CHAIN':
    # <MoleculeType.DANGLING_CHAIN: 3>, 'FREE_CHAIN':
    # <MoleculeType.FREE_CHAIN: 4>}
    __members__: typing.ClassVar[dict[str, MoleculeType]]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...


class NeighbourList:
    """
    Gives access to somewhat fast queries on the neighbourhood of atoms
    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, atoms: list[Atom], box: Box, cutoff: float) -> None:
        """
        Instantiates a new neighbour list
        """

    def get_atoms_close_to(self, atom: Atom, upper_cutoff: float = 1.0, lower_cutoff: float = 0.0,
                           unwrapped: bool = False, expect_self: bool = False) -> list[Atom]:
        """
                  List all atoms that are close to a given one.

                  It is possible to request it within a new cutoff,
                  though the underlying neighbour list will not be regenerated.
                  For performance reasons, it is recommended to initialize a
                  new NeighbourList if you require a different cutoff, depending on your use case.

                  You can use a negative value for the upper_cutoff to use the cutoff used for
                  filling the neighbour list buckets.
        """

    def remove_atom(self, atom: Atom, debug_hint: str = '') -> None:
        """
                  Remove an atom from this neighbour list.
                  It will not show up when querying for neighbours,
                  but its neighbours cannot be queried either.
        """


class NonGaussianSpringForceEvaluator(MEHPForceEvaluator):
    """

         This is equal to a spring evaluator for Langevin chains.

         The force for a certain spring is given by:
         :math:`f = 0.5 \\cdot \\frac{1}{l} \\scriptL^{-1}(\\frac{r}{N\\cdot l})`,
         where :math:`r` is the spring [between cross-linkers] length
         and :math:`\\scriptL^{-1}` the inverse langevin function.

         Please note that the inverse langevin is only approximated.

         Recommended optimization algorithm: "LD_MMA"

         :param kappa: the spring constant :math:`\\kappa`
         :param N: The number of links in a spring
         :param l: The  the length of a spring in the chain

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, kappa: float = 1.0, N: float = 1.0,
                 l: float = 1.0) -> None:
        """
        Initialize this ForceEvaluator
        """


class NormalModeAnalyzer:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __getstate__(self) -> tuple:
        ...

    def __init__(self, spring_from: list[int], spring_to: list[int]) -> None:
        """
        Initialize NormalModeAnalyzer
        """

    def __setstate__(self, arg0: tuple) -> None:
        ...

    def evaluate_loss_modulus(self, omega: numpy.ndarray) -> numpy.ndarray:
        """
        Evaluate the loss modulus :math:`G''`. Yet misses the conversion factor.
        """

    def evaluate_storage_modulus(self, omega: numpy.ndarray) -> numpy.ndarray:
        """
        Evaluate the storage modulus :math:`G'`. Yet misses the conversion factor.
        """

    def evaluate_stress_autocorrelation(
            self, t: numpy.ndarray) -> numpy.ndarray:
        """
        Evaluate stress autocorrelation
        """

    def find_all_eigenvalues(self, compute_eigenvectors: bool = False) -> None:
        """
        Find all eigenvalues using a dense solver
        """

    def find_sparse_eigenvalues(
            self, nr_of_eigenvalues: int, compute_eigenvectors: bool = False) -> None:
        """
        Find the `k` smallest eigenvalues using a sparse solver
        """

    def get_eigenvalues(self) -> numpy.ndarray:
        """
        Get the eigenvalues
        """

    def get_eigenvectors(self) -> numpy.ndarray:
        """
        Get eigenvectors
        """

    def get_matrix(self) -> scipy.sparse.csc_matrix:
        """
        Get the assembled connectivity matrix.
        """

    def get_matrix_size(self) -> int:
        """
        Get the size of the matrix (the maximum of eigenvalues that could be queried)
        """

    def get_nr_of_soluble_clusters(self) -> int:
        """
        Get the number of soluble clusters (Eigenvalues = 0)
        """

    def set_eigenvalues(self, eigenvalues: numpy.ndarray) -> None:
        """
        Set the eigenvalues, e.g. if you use an external solver
        """

    def set_eigenvectors(self, eigenvectors: numpy.ndarray) -> None:
        """
        Set eigenvectors, e.g. if you use an external solver
        """


class OutputConfiguration:
    double_values: list[ComputedDoubleValues]
    int_values: list[ComputedIntValues]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self) -> None:
        """
        Get an instance of this struct
        """
    @property
    def append(self) -> bool:
        """
             Whether to append to the file or truncate it
        """
    @append.setter
    def append(self, arg0: bool) -> None:
        ...

    @property
    def filename(self) -> str:
        """
             The file to write to. Empty means standard output (console).
        """
    @filename.setter
    def filename(self, arg0: str) -> None:
        ...

    @property
    def output_every(self) -> int:
        """
             How often to write the values to the output.
             For averages, this value also says how many values will be averaged.
        """
    @output_every.setter
    def output_every(self, arg0: int) -> None:
        ...

    @property
    def use_every(self) -> int:
        """
             For autocorrelation/averaging, how often to include values
        """
    @use_every.setter
    def use_every(self, arg0: int) -> None:
        ...


class SimpleSpringMEHPForceEvaluator(MEHPForceEvaluator):
    """

         This is equal to a spring evaluator for Gaussian chains.

         The force for a certain spring is given by:
         :math:`f = 0.5 \\cdot \\kappa r`,
         where :math:`r` is the spring [between cross-linkers] length.

         Recommended optimization algorithm: "LD_LBFGS"

         :param kappa: the spring constant :math:`\\kappa`

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __init__(self, kappa: float = 1.0) -> None:
        ...


class SimplifiedBalanceNetwork:
    """

         A more efficient structure of the network for use in MEHP force balance.
         Consists usually only of the cross- and slip-links.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def box_lengths(self) -> float:
        ...

    @property
    def coordinates(self) -> numpy.ndarray:
        ...

    @property
    def link_indices_of_springs(self) -> list[list[int]]:
        ...

    @property
    def link_is_sliplink(self) -> numpy.ndarray:
        ...

    @property
    def local_to_global_spring_index(self) -> list[list[int]]:
        ...

    @property
    def nr_of_crosslink_swaps_endured(self) -> numpy.ndarray:
        ...

    @property
    def nr_of_crosslinks(self) -> int:
        ...

    @property
    def nr_of_links(self) -> int:
        ...

    @property
    def nr_of_partial_springs(self) -> int:
        ...

    @property
    def nr_of_springs(self) -> int:
        ...

    @property
    def old_atom_ids(self) -> numpy.ndarray:
        ...

    @property
    def partial_to_full_spring_index(self) -> numpy.ndarray:
        ...

    @property
    def spring_contour_length(self) -> numpy.ndarray:
        ...

    @property
    def spring_coordinate_index_a(self) -> numpy.ndarray:
        ...

    @property
    def spring_coordinate_index_b(self) -> numpy.ndarray:
        ...

    @property
    def spring_index_a(self) -> numpy.ndarray:
        ...

    @property
    def spring_index_b(self) -> numpy.ndarray:
        ...

    @property
    def spring_indices_of_links(self) -> list[list[int]]:
        ...

    @property
    def spring_part_coordinate_index_a(self) -> numpy.ndarray:
        ...

    @property
    def spring_part_coordinate_index_b(self) -> numpy.ndarray:
        ...

    @property
    def spring_part_index_a(self) -> numpy.ndarray:
        ...

    @property
    def spring_part_index_b(self) -> numpy.ndarray:
        ...

    @property
    def volume(self) -> float:
        ...


class SimplifiedNetwork:
    """

         A more efficient structure of the network for use in MEHP.
         Consists usually only of the crosslinkers.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    @property
    def box_lengths(self) -> float:
        ...

    @property
    def coordinates(self) -> numpy.ndarray:
        ...

    @property
    def nr_of_crosslinks(self) -> int:
        ...

    @property
    def nr_of_nodes(self) -> int:
        ...

    @property
    def nr_of_springs(self) -> int:
        ...

    @property
    def old_atom_ids(self) -> numpy.ndarray:
        ...

    @property
    def spring_coordinate_index_a(self) -> numpy.ndarray:
        ...

    @property
    def spring_coordinate_index_b(self) -> numpy.ndarray:
        ...

    @property
    def spring_index_a(self) -> numpy.ndarray:
        ...

    @property
    def spring_index_b(self) -> numpy.ndarray:
        ...

    @property
    def volume(self) -> float:
        ...


class StructureSimplificationMode:
    """
    Members:

      NO_SIMPLIFICATION

      X2F_ONLY

      INACTIVE_ONLY

      ALL_TIM

      ALL_ANDREI
    """
    ALL_ANDREI: typing.ClassVar[StructureSimplificationMode]  # value = <StructureSimplificationMode.ALL_ANDREI: 4>
    # value = <StructureSimplificationMode.ALL_TIM: 3>
    ALL_TIM: typing.ClassVar[StructureSimplificationMode]
    # value = <StructureSimplificationMode.INACTIVE_ONLY: 2>
    INACTIVE_ONLY: typing.ClassVar[StructureSimplificationMode]
    # value = <StructureSimplificationMode.NO_SIMPLIFICATION: 0>
    NO_SIMPLIFICATION: typing.ClassVar[StructureSimplificationMode]
    # value = <StructureSimplificationMode.X2F_ONLY: 1>
    X2F_ONLY: typing.ClassVar[StructureSimplificationMode]
    # value = {'NO_SIMPLIFICATION':
    # <StructureSimplificationMode.NO_SIMPLIFICATION: 0>, 'X2F_ONLY':
    # <StructureSimplificationMode.X2F_ONLY: 1>, 'INACTIVE_ONLY':
    # <StructureSimplificationMode.INACTIVE_ONLY: 2>, 'ALL_TIM':
    # <StructureSimplificationMode.ALL_TIM: 3>, 'ALL_ANDREI':
    # <StructureSimplificationMode.ALL_ANDREI: 4>}
    __members__: typing.ClassVar[dict[str, StructureSimplificationMode]]

    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __eq__(self, other: typing.Any) -> bool:
        ...

    def __getstate__(self) -> int:
        ...

    def __hash__(self) -> int:
        ...

    def __index__(self) -> int:
        ...

    def __init__(self, value: int) -> None:
        ...

    def __int__(self) -> int:
        ...

    def __ne__(self, other: typing.Any) -> bool:
        ...

    def __repr__(self) -> str:
        ...

    def __setstate__(self, state: int) -> None:
        ...

    def __str__(self) -> str:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def value(self) -> int:
        ...


class Universe:
    """
    Represents a full Polymer Network structure, a collection of molecules.
    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __contains__(self, arg0: Atom) -> bool:
        """
                  Check whether a particular atom is contained in this universe.
        """

    def __copy__(self) -> Universe:
        ...

    def __getitem__(self, arg0: int) -> Atom:
        """
               Access an atom by its vertex index.
        """

    def __getstate__(self) -> tuple:
        ...

    def __init__(self, Lx: float, Ly: float, Lz: float) -> None:
        """
        Instantiate this Universe (Collection of Molecules) providing the box lengths.
        """

    def __len__(self) -> int:
        """
               Get the number of atoms in this universe.
        """

    def __setstate__(self, arg0: tuple) -> None:
        ...

    def add_angles(self, angles_from: list[int], angles_via: list[int],
                   angles_to: list[int], angle_types: list[int]) -> None:
        """
        Add angles to the Universe. No relation to the underlying graph, just a method to preserve read & write capabilities
        """

    def add_atoms(self, ids: list[int], types: list[int], x: list[float], y: list[float],
                  z: list[float], nx: list[int], ny: list[int], nz: list[int]) -> None:
        """
        Add atoms to the Universe, vertices to the underlying graph.
        """
    @typing.overload
    def add_bonds(self, bonds_from: list[int], bonds_to: list[int]) -> None:
        """
        Add bonds to the underlying atoms, edges to the underlying graph. If the connected atoms are not found, the bonds are silently skipped.
        """
    @typing.overload
    def add_bonds(self, nr_of_bonds: int, bonds_from: list[int], bonds_to: list[int], bond_types: list[int],
                  ignore_non_existent_atoms: bool = False, simplify_universe: bool = True) -> None:
        """
        Add bonds to the underlying atoms, edges to the underlying graph.
        """

    def add_bonds_with_types(
            self, bonds_from: list[int], bonds_to: list[int], bond_types: list[int]) -> None:
        """
        Add bonds to the underlying atoms, edges to the underlying graph. If the connected atoms are not found, the bonds are silently skipped.
        """

    def add_dihedral_angles(self, angles_from: list[int], angles_via1: list[int],
                            angles_via2: list[int], angles_to: list[int], angle_types: list[int]) -> None:
        """
        Add dihedral angles to the Universe. No relation to the underlying graph, just a method to preserve read & write capabilities
        """

    def compute_angles(self) -> list[float]:
        """
        Computes the angle :math:`\\theta` of each angle in the molecule, respecting periodic boundaries.
        """

    def compute_bond_lengths(self) -> list[float]:
        """
        Computes the length :math:`b` of each bond in the molecule, respecting periodic boundaries.
        """

    def compute_dxs(self, atomIdsTo: list[int],
                    atomIdsFrom: list[int]) -> list[float]:
        """
        Compute the dx distance for certain bonds (length in x direction).
        """

    def compute_dys(self, atomIdsTo: list[int],
                    atomIdsFrom: list[int]) -> list[float]:
        """
        Compute the dy distance for certain bonds (length in y direction).
        """

    def compute_dzs(self, atomIdsTo: list[int],
                    atomIdsFrom: list[int]) -> list[float]:
        """
        Compute the dz distance for certain bonds (length in z direction).
        """

    def compute_end_to_end_distances(
            self, crosslinker_type: int, derive_image_flags: bool = False) -> list[float]:
        """
                  Compute the end-to-end distance of each strand in the network.

                  NOTE:
                       Internally, this uses either :func:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance`
                       or :func:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags`,
                       depending on `derive_image_flags`.
                       Invalid strands (where said function returns 0.0 or -1.0) are ignored.
        """

    def compute_mean_end_to_end_distance(
            self, crosslinker_type: int, derive_image_flags: bool = False) -> float:
        """
                  Computes the mean of the end-to-end distances of each strand in the network.

                  NOTE:
                       Internally, this uses either :func:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance`
                       or :func:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags`,
                       depending on `derive_image_flags`.
                       Invalid strands (where said function returns 0.0 or -1.0) are ignored.
        """

    def compute_mean_squared_end_to_end_distance(
            self, crosslinker_type: int, only_those_with_two_crosslinkers: bool = False, derive_image_flags: bool = False) -> float:
        """
                  Computes the mean square of the end-to-end distances of each strand (incl. cross-links) in the network.

                  NOTE:
                       Internally, this uses either :func:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance`
                       or :func:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags`,
                       depending on `derive_image_flags`.
                       Invalid strands (where said function returns 0.0 or -1.0) are ignored.
        """

    def compute_mean_strand_length(self, crosslinker_type: int) -> float:
        """
                      Compute the mean number of beads per strand.
        """

    def compute_number_average_molecular_weight(
            self, crosslinker_type: int) -> float:
        """
                      Compute the number average molecular weight.

                      NOTE:
                            Cross-linkers are ignored completely.
        """

    def compute_polydispersity_index(self, crosslinker_type: int) -> float:
        """
                      Compute the polydispersity index:
                      the weight average molecular weight over the number average molecular weight.
        """

    def compute_temperature(self, dimensions: int = 3,
                            k_b: float = 1.0) -> float:
        """
        Use the velocities per atom to compute the temperature from the kinetic energy of the system.
        """

    def compute_total_mass(self) -> float:
        """
                  Compute the total mass of this network/universe in whatever mass unit was used when
                  :func:`~pylimer_tools_cpp.Universe.setMasses()` was called.
        """

    def compute_weight_average_molecular_weight(
            self, crosslinker_type: int) -> float:
        """
                      Compute the weight average molecular weight.

                      NOTE:
                            Cross-linkers are ignored completely.
        """

    def compute_weight_fractions(self) -> dict[int, float]:
        """
                    Compute the weight fractions of each atom type in the network.

                    If no masses are stored, assumes a mass of 1 for each atom.

                    If the total mass is 0., returns the total mass per atom type.
        """

    def contract_vertices_along_bond_type(self, bond_type: int) -> Universe:
        """
                  Merge vertices along a specific bond type.

                  May result in new self-loops; use :func:`~pylimer_tools_cpp.Universe.simplify()` to remove them.
        """

    def count_atom_types(self) -> dict[int, int]:
        """
                  Count how often each atom type is present.
        """

    def count_atoms_in_skin_distance(
            self, distances: list[float], unwrapped: bool = False) -> list[int]:
        """
                  This is a function that may help you to compute the radial distribution function.
                  It loops the

                  Parameters:
                       - distances: the edges of the bins
                       - unwrapped: whether to measure the distance in unwrapped coordinates or as PBC-corrected distance
        """

    def count_loop_lengths(self, max_length: int = -1) -> dict[int, int]:
        """
                  Find all loops (below a specific length) and count the number of atoms involved in them.
                  Returns the count, how many loops per length are found.
        """

    def detect_angles(self) -> dict[str, list[int]]:
        """
        Returns just as
                  :func:`~pylimer_tools_cpp.Universe.getAngles`,
                  but all angles that are detected in the network, rather than the one already set.
                  Note that the angle types are determined by
                  :func:`~pylimer_tools_cpp.Universe.hashAngleType`,
                  which serves angle types that should be mapped by you back to smaller numbers,
                  before serving them to :func:`~pylimer_tools_cpp.Universe.addAngles`.
        """

    def detect_dihedral_angles(self) -> dict[str, list[int]]:
        """
        Returns just as
                  :func:`~pylimer_tools_cpp.Universe.getDihedralAngles`,
                  but all dihedral angles that are detected in the network, rather than the one already set.
                  Note that the angle types are determined by
                  :func:`~pylimer_tools_cpp.Universe.hashDihedralAngleType`,
                  which serves angle types that should be mapped by you back to smaller numbers,
                  before serving them to :func:`~pylimer_tools_cpp.Universe.addDiheralAngles`.
        """

    def determine_effective_functionality_per_type(self) -> dict[int, float]:
        """
                    Find the average functionality of each atom type in the network.
        """

    def determine_functionality_per_type(self) -> dict[int, int]:
        """
                    Find the maximum functionality of each atom type in the network.
        """

    def find_loops(self, crosslinker_type: int, max_length: int = -1,
                   skip_self_loops: bool = False) -> dict[int, list[list[Atom]]]:
        """
                    Decompose the Universe into loops.
                    The primary index specifies the degree of the loop.

                    CAUTION:
                       There are exponentially many paths between two cross-linkers of a network,
                       and you may run out of memory when using this function, if your Universe/Network is lattice-like.
                       You can use the `max_length` parameter to restrict the algorithm to only search for loops up to a certain length.
                       Use a negative value to find all loops and paths.
        """

    def find_minimal_order_loop_from(self, loop_start: int, loop_step1: int,
                                     max_length: int = -1, skip_self_loops: bool = False) -> list[Atom]:
        """
                    Decompose the Universe into loops.
                    The primary index specifies the degree of the loop.

                    CAUTION:
                       There are exponentially many paths between two cross-linkers of a network,
                       and you may run out of memory when using this function, if your Universe/Network is lattice-like.
                       You can use the `max_length` parameter to restrict the algorithm to only search for loops up to a certain length.
                       Use a negative value to find all loops and paths.
        """

    def get_angles(self) -> dict[str, list[int]]:
        """
                   Get all angles added to this network.

                   Returns a dict with three properties: 'angle_from', 'angle_via' and 'angle_to'.

                   NOTE:
                       The integer values returned refer to the the atom ids, not the vertex ids.
                       Use :func:`~pylimer_tools_cpp.Universe.get_idx_by_atom_id` to translate them to vertex ids.
        """

    def get_atom(self, atom_id: int) -> Atom:
        """
        Find an atom by its ID.
        """

    def get_atom_by_vertex_id(self, vertex_id: int) -> Atom:
        """
        Find an atom by the ID of the vertex of the underlying graph.
        """

    def get_atom_id_by_vertex_idx(self, vertex_id: int) -> int:
        """
        Get the id of the atom by the vertex id of the underlying graph.
        """

    def get_atom_types(self) -> list[int]:
        """
                  Get all types (each one for each atom) ordered by atom vertex id.
        """

    def get_atoms(self) -> list[Atom]:
        """
                    Get all atoms.
        """

    def get_atoms_by_degree(self, functionality: int) -> list[Atom]:
        """
                    Get the atoms that have the specified number of bonds.
        """

    def get_atoms_by_type(self, atom_type: int) -> list[Atom]:
        """
                    Query all atoms by their type.
        """

    def get_atoms_connected_to(self, atom: Atom) -> list[Atom]:
        """
                    Get the atoms connected to a specified atom.

                    Internally uses :func:`~pylimer_tools_cpp.Universe.getAtomsConnectedTo`
        """

    def get_atoms_connected_to_vertex(self, vertex_idx: int) -> list[Atom]:
        """
                    Get the atoms connected to a specified vertex id.
        """

    def get_bonds(self) -> dict[str, list[int]]:
        """
                    Get all bonds. Returns a dict with three properties: 'bond_from', 'bond_to' and 'bond_type'.
                    The order is not necessarily related to any structural characteristic.
        """

    def get_box(self) -> Box:
        """
                    Get the underlying bounding box object.
        """

    def get_chains_with_crosslinker(
            self, crosslinker_type: int) -> list[Molecule]:
        """
                    Decompose the Universe into strands (molecules, which could be either chains, or even lonely atoms), without omitting the cross-linkers
                    (as in :func:`~pylimer_tools_cpp.Universe.getMolecules(crosslinker_type)`).
                    In turn, e.g. for a tetrafunctional cross-linker, it will be 4 times in the resulting molecules.

                    NOTE:
                       Cross-linkers without bonds to non-cross-linkers are not returned
                       (i.e., single cross-linkers, are not counted as strands).
        """

    def get_clusters(self) -> list[Universe]:
        """
                  Get the components of the universe that are not connected to each other.
                  Returns a list of :obj:`~pylimer_tools_cpp.Universe`s.
                  Unconnected, free atoms/beads become their own :obj:`~pylimer_tools_cpp.Universe`.
        """

    def get_edge_ids_from(self, vertex_id: int) -> list[int]:
        ...

    def get_edge_ids_from_to(self, vertex_id_from: int,
                             vertex_id_to: int) -> list[int]:
        """
              Get the edge ids of the edges between two specific vertices
        """

    def get_edges(self) -> dict[str, list[int]]:
        """
                    Get all edges. Returns a dict with three properties: 'edge_from', 'edge_to' and 'edge_type'.
                    The order is not necessarily related to any structural characteristic.

                    NOTE:
                       The integer values returned refer to the vertex ids, not the atom ids.
                       Use :func:`~pylimer_tools_cpp.Universe.get_atom_id_by_idx` to translate them to atom ids, or
                       :func:`~pylimer_tools_cpp.Universe.get_bonds` to have that done for you.
        """

    def get_masses(self) -> dict[int, float]:
        """
                    Get the mass of one atom per type
        """

    def get_molecules(self, atom_type_to_omit: int) -> list[Molecule]:
        """
                  Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms.

                  Reduces the Universe to a list of molecules.
                  Specify the crosslinker_type to an existing type id,
                  then those atoms will be omitted, and this function returns chains instead.
        """

    def get_network_of_crosslinker(self, crosslinker_type: int) -> Universe:
        """
                    Reduce the network to contain only cross-linkers, replacing all the strands with a single bond.
                    Useful e.g. to reduce the memory useage and runtime of
                    :func:`~pylimer_tools_cpp.Universe.findLoops()` or
                    :func:`~pylimer_tools_cpp.Universe.hasInfiniteStrand()`.

                    Further use :func:`~pylimer_tools_cpp.Universe.simplify()` to remove primary loops.
        """

    def get_nr_of_angles(self) -> int:
        """
                    Query the number of angles that have been added to this universe.
        """

    def get_nr_of_atoms(self) -> int:
        """
                    Query the number of atoms in this universe.
        """

    def get_nr_of_bonds(self) -> int:
        """
                    Query the number of bonds associated with this universe.
        """

    def get_nr_of_bonds_of_atom(self, arg0: int) -> int:
        """
        Count the number of immediate neighbors of an atom, specified by its id.
        """

    def get_nr_of_bonds_of_vertex(self, arg0: int) -> int:
        """
        Count the number of immediate neighbors of an atom, specified by its vertex id.
        """

    def get_nr_of_dihedral_angles(self) -> int:
        """
                    Query the number of dihedral angles that have been added to this universe.
        """

    def get_timestep(self) -> int:
        """
                    Query the timestep when this universe was captured.
        """

    def get_vertex_degrees(self) -> list[int]:
        """
                  Get the degree (functionality) of each vertex.
        """

    def get_vertex_idx_by_atom_id(self, atom_id: int) -> int:
        """
        Get the vertex id of the underlying graph for an atom with a specified id.
        """

    def get_volume(self) -> float:
        """
                    Query the volume of the underlying bounding box.
        """

    def has_infinite_strand(self, arg0: int, arg1: int) -> bool:
        """
                   Checks whether there is a strand (with cross-linker) in the universe that loops through periodic images without coming back.

                    CAUTION:
                       There are exponentially many paths between two cross-linkers of a network,
                       and you may run out of memory when using this function, if your Universe/Network is lattice-like.
        """

    def hash_angle_type(self, angle_from: int,
                        angle_via: int, angle_to: int) -> int:
        """
                  Convert the three integers to one long number/hash.
                  Used internally for duplicate detection.
        """

    def hash_dihedral_angle_type(
            self, angle_from: int, angle_via1: int, angle_via2: int, angle_to: int) -> int:
        """
                  Convert the four integers to one long number/hash.
                  Used internally for duplicate detection.
        """

    def interpolate_edges(self, crosslinker_type: int,
                          interpolation_factor: float) -> list[tuple[int, int]]:
        """
                  Get more or less edges than currently present, interpolating between junctions.
        """

    def remove_all_angles(self) -> None:
        ...

    def remove_all_dihedral_angles(self) -> None:
        ...

    def remove_atoms(self, atom_ids: list[int]) -> None:
        """
                  Remove atoms and all associated bonds by their atom ids.
        """

    def remove_bonds(self, bonds_from: list[int], bonds_to: list[int]) -> None:
        """
                  Remove bonds by their connected atom ids.
        """

    def remove_bonds_by_type(self, bond_type: int) -> None:
        """
                  Remove bonds with a specific type.
        """

    def replace_atom(self, atom_id: int, replacement_atom: Atom) -> None:
        """
                  Replace the properties of an atom with the properties of another given atom.
        """

    def replace_atom_type(self, atom_id: int, new_type: int) -> None:
        """
                  Replace the type of an atom with another type.
        """

    def resample_velocities(self, mean: float, variance: float,
                            seed: str = '', is_2d: bool = False) -> None:
        ...

    def set_box(self, box: Box, rescale_atoms: bool = False) -> None:
        """
                  Override the currently assigned box with the one specified.
        """

    def set_box_lengths(self, lx: float, ly: float, lz: float,
                        rescale_atoms: bool = False) -> None:
        """
                  Override the currently assigned box with one with the side lengths specified.
        """

    def set_mass(self, atom_type: int, mass: float) -> None:
        """
        Set the mass for a specific atom type.
        """

    def set_masses(self, mass_per_type: dict[int, float]) -> None:
        """
        Set the mass per type of atom.
        """

    def set_timestep(self, timestep: int) -> None:
        """
        Set the time-step when this Universe was captured.
        """

    def set_vertex_property(self, vertex_id: int,
                            property_name: str, value: float) -> None:
        """
                  Set a specific property for a specific vertex.
        """

    def simplify(self) -> None:
        """
        Remove self links and double bonds. This function is called automatically after adding bonds.
        """


class UniverseSequence:
    """

         This class represents a sequence of Universes, with the Universe's data
         only being read on request. Dump files are read at once in order
         to know how many timesteps/universes are available in total
         (but the universes' data is not read on first look through the file).
         This, while it can lead to two (or more) reads of the whole file,
         is a measure in order to enable low memory useage if needed (i.e. for large dump files).
         Use Python's iterator to have this UniverseSequence only ever retain one universe in memory.
         Alternatively, use :func:`~pylimer_tools_cpp.UniverseSequence.forgetAtIndex`
         to have the UniverseSequence forget about already read universes.

    """
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...

    def __getitem__(self, arg0: int) -> Universe:
        """
        Get a universe by its index.
        """

    def __init__(self) -> None:
        ...

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

    def at_index(self, index: int) -> Universe:
        """
        Get the Universe at the given index (as of in the sequence given by the dump file).
        """

    def compute_distance_autocorrelation_from_to(
            self, atom_ids_from: list[int], atom_ids_to: list[int], nr_of_origins: int = 25, reduce_memory: bool = False) -> dict[int, float]:
        """
                  Compute the autocorrelation of the dot product of the distance vector from certain to other atoms.

                  For example, this can be used to compute Eq. 4.51 from Masao Doi, Introduction to Polymer Physics, p. 74.
        """

    def compute_distance_from_to_atoms(
            self, atom_ids_from: list[int], atom_ids_to: list[int], reduce_memory: bool = False) -> list[float]:
        """
                  Compute the root square norm of all the (unwrapped!) distances for the given pair of atoms.

                  Can be used to somewhat faster compute e.g. all the end-to-end or bond distances.
                  Pay attention that the image flags are correct, otherwise, this data may not be useable.
        """

    def compute_msd_for_atom_properties(self, atom_ids: list[int], x_property: str, y_property: str,
                                        z_property: str, nr_of_origins: int = 25, reduce_memory: bool = False) -> dict[int, float]:
        """
                  Compute the mean square displacement for atoms with the specified ids
        """

    def compute_msd_for_atoms(
            self, atom_ids: list[int], nr_of_origins: int = 25, reduce_memory: bool = False) -> dict[int, float]:
        """
                  Compute the mean square displacement for atoms with the specified ids
        """

    def compute_vector_from_to_atoms(
            self, atom_ids_from: list[int], atom_ids_to: list[int], reduce_memory: bool = False) -> list[numpy.ndarray]:
        """
                  Compute the (unwrapped!) distances for the given pair of atoms.

                  Can be used to somewhat faster compute e.g. all the end-to-end or bond vectors.
                  Pay attention that the image flags are correct, otherwise, this data may not be useable.
        """

    def forget_at_index(self, index: int) -> None:
        """
        Clear the memory of the Universe at the given index (as of in the
                   sequence given by the dump file).
        """

    def get_all(self) -> list[Universe]:
        """
                    Get all universes initialized back in a list.
                    For big dump files or lots of data files, this might lead to memory issues.
                    Use :func:`~pylimer_tools_cpp.UniverseSequence.__iter__`
                    to have
                    or :func:`~pylimer_tools_cpp.UniverseSequence.atIndex`
                    and :func:`~pylimer_tools_cpp.UniverseSequence.forgetAtIndex`
                    to craft a more memory-efficient retrieval mechanism.
        """

    def get_length(self) -> int:
        """
                    Get the number of universes in this sequence.
        """

    def initialize_from_data_sequence(self, data_files: list[str]) -> None:
        """
        Reset and initialize the Universes from an ordered list of Lammps data (:code:`write_data`) files.
        """

    def initialize_from_dump_file(
            self, initial_data_file: str, dump_file: str) -> None:
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

    def reset_iterator(self) -> None:
        """
                  Reset the internal iterator, such that a subsequent call to
                  :func:`~pylimer_tools_cpp.UniverseSequence.next` returns the first one again.
        """

    def set_data_file_atom_style(self, atom_styles: list[AtomStyle]) -> None:
        """
                  Set the format of the data files to be read. See :obj:`~pylimer_tools_cpp.AtomStyle`.
        """


def compute_stoichiometric_imbalance(
        arg0: Universe, arg1: int, arg2: int, arg3: dict[int, int]) -> float:
    """
    Compute stoichiometric imbalance
    """


def do_linear_walk_chain_from_to(box: Box, from_coordinates: numpy.ndarray,
                                 to_coordinates: numpy.ndarray, chain_len: int, include_ends: bool = False) -> numpy.ndarray:
    """
                Get coordinates linearly interpolated from one point to another (both exclusive).

                :param box: The box for doing PBC correction on the from/to.
                :param from_coordinates: Coordinates of the start point.
                :param to_coordinates: Coordinates of the end point.
                :param chain_len: Number of coordinates to generate between the start and end-point.
                :param include_ends: Whether to include the start and end points in the output (default: false).
                   If yes, chain_len + 2 coordinates will be returned,
                   where the first will be from_coordinates and the last will be to_coordinates.
    """


def do_random_walk(chain_len: int, bead_distance: float = 1.0,
                   mean_squared_bead_distance: float = 1.0, seed: str = '') -> numpy.ndarray:
    """
                Do a random walk, return the coordinates of each point visited.
    """


def do_random_walk_chain_from_to(box: Box, from_coordinates: numpy.ndarray, to_coordinates: numpy.ndarray, chain_len: int,
                                 bead_distance: float = 1.0, mean_squared_bead_distance: float = 1.0, seed: str = '') -> numpy.ndarray:
    """
                Do a random walk from one point to another.
    """


def do_random_walk_chain_from_to_mc(box: Box, from_coordinates: numpy.ndarray, to_coordinates: numpy.ndarray, chain_len: int,
                                    bead_distance: float = 1.0, mean_squared_bead_distance: float = 1.0, seed: str = '', n_iterations: int = 10000) -> numpy.ndarray:
    """
                Do a random walk from one point to another.
                Then, relax the points in between using a Metropolis-Monte Carlo simulation.
    """


def inverse_langevin(x: float) -> float:
    """
         A somewhat accurate (for :math:`x \\in (-1, 1)`) implementation of the inverse Langevin.

         Source: https://scicomp.stackexchange.com/a/30251
    """


def predict_gelation_point(arg0: float, arg1: int, arg2: int) -> float:
    """
    Predict the gelation point of a Universe
    """


def randomly_sample_entanglements(universe: Universe, nr_of_samples: int, upper_cutoff: float, lower_cutoff: float = 0, minimum_nr_of_samples: int = 0,
                                  same_strand_cutoff: float = 3.0, seed: str = '', crosslinker_type: int = 2, ignore_crosslinks: bool = True) -> AtomPairEntanglements:
    """
        Randomly find pairs of atoms that are close together and could be
        entanglements

        Arguments:
        :param universe: The universe of atoms from which to sample entanglements from.
        :param nr_of_samples: The number of pairs of atoms to randomly sample.
        :param upper_cutoff: The maximum distance between atoms for a pair to be considered a potential entanglement.
        :param lower_cutoff: The minimum distance between atoms for a pair to be considered a potential entanglement.
        :param minimum_nr_of_samples: The minimum number of entanglements to be found.
        :param same_strand_cutoff: The maximum distance between atoms on the same strand for a pair to be considered a potential entanglement.
        :param seed: A seed for the random number generator.
        :param crosslinker_type: The type of crosslinker to consider when finding entanglements. Used for the splitting into strands.
        :param ignore_crosslinks: Whether to ignore crosslinks when finding entanglements.
          Careful: if you don't ignore them, the same-strand policy might not work correctly,
          since each cross-link should actually be associated with more than one strand.
    """


def split_csv(arg0: str, arg1: str) -> list[str]:
    """
    Read a file containing a number of CSVs. Returns them split up.
    """


def version_information() -> str:
    """
        Returns  a string of the the current version, incl. git hash and date of compilation.
    """


DANGLING_CHAIN: MoleculeType  # value = <MoleculeType.DANGLING_CHAIN: 3>
FREE_CHAIN: MoleculeType  # value = <MoleculeType.FREE_CHAIN: 4>
NETWORK_STRAND: MoleculeType  # value = <MoleculeType.NETWORK_STRAND: 1>
PRIMARY_LOOP: MoleculeType  # value = <MoleculeType.PRIMARY_LOOP: 2>
UNDEFINED: MoleculeType  # value = <MoleculeType.UNDEFINED: 0>
