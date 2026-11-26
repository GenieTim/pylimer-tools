"""

    pylimer_tools_cpp
    -----------------

    A collection of utility python functions for handling LAMMPS output and polymers in Python.

    .. autosummary::
        :toctree: _generate

    
"""
from __future__ import annotations
import collections.abc
import enum
import numpy
import numpy.typing
import scipy.sparse
import typing
__all__: list[str] = ['Atom', 'AtomPairEntanglements', 'AtomStyle', 'AveFileReader', 'BackTrackStatus', 'BondLengthStats', 'Box', 'ComputeBackend', 'ComputedDoubleValues', 'ComputedIntValues', 'Correlator', 'CrosslinkerUniverse', 'DPDSimulator', 'DataFileReader', 'DataFileWriter', 'DumpFileReader', 'EndToEndStats', 'EntanglementFilteringMode', 'EntanglementMode', 'ExitReason', 'HookeanForceEvaluator', 'HybridConfig', 'KPMConfig', 'LazyUniverseSequenceIterator', 'LinearMaxDistanceProvider', 'LinearSystemMEHPSolver', 'LinkSwappingMode', 'MCUniverseGenerator', 'MEHPForceEvaluator', 'MEHPLinearSolverConfig', 'MEHPNetwork', 'MEHPNonlinearSolverConfig', 'MEHPSimulator', 'MEHPSlipLinkSolverConfig', 'MEHPSolveResult', 'MEHPSolver', 'MaxDistanceProvider', 'Molecule', 'MoleculeIterator', 'MoleculeType', 'NLoptAlgorithm', 'NoMaxDistanceProvider', 'NonGaussianSpringForceEvaluator', 'NonlinearOptimizationSolver', 'NormalModeAnalyzer', 'NormalModeSolver', 'NormalModeSpectrum', 'OutputConfiguration', 'ResistanceDistanceAnalysis', 'RouseBrownianSimulator', 'RouseSimulationResult', 'RouseSolverType', 'SLESolver', 'SLQConfig', 'STOP', 'SlipLinkRelaxationSolver', 'SlipMode', 'StrandCategory', 'StrandPairEntanglements', 'StructureSimplificationMode', 'TRACK_BACKWARD', 'TRACK_FORWARD', 'Universe', 'UniverseSequence', 'ZScoreMaxDistanceProvider', 'do_linear_walk_chain_from_to', 'do_random_walk', 'do_random_walk_chain_from_to', 'do_random_walk_chain_from_to_mc', 'inverse_langevin', 'peel_topological_leaves', 'randomly_sample_entanglements', 'randomly_sample_entanglements_on_demand', 'randomly_sample_entanglements_v2', 'sample_strand_beads_on_demand', 'split_csv', 'version_information']
class Atom:
    """
    
           A single bead or atom.
      
    """
    __hash__: typing.ClassVar[None] = None
    def __eq__(self, arg0: Atom) -> bool:
        ...
    def __getstate__(self) -> tuple[int, int, float, float, float, int, int, int, dict[str, float]]:
        ...
    @typing.overload
    def __init__(self, id: typing.SupportsInt | typing.SupportsIndex, type: typing.SupportsInt | typing.SupportsIndex, x: typing.SupportsFloat | typing.SupportsIndex, y: typing.SupportsFloat | typing.SupportsIndex, z: typing.SupportsFloat | typing.SupportsIndex, nx: typing.SupportsInt | typing.SupportsIndex, ny: typing.SupportsInt | typing.SupportsIndex, nz: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
                 Construct this atom.
        
                 :param id: Unique identifier for the atom
                 :param type: Type classification of the atom
                 :param x: X coordinate position
                 :param y: Y coordinate position  
                 :param z: Z coordinate position
                 :param nx: Periodic image flag in x direction
                 :param ny: Periodic image flag in y direction
                 :param nz: Periodic image flag in z direction
        """
    @typing.overload
    def __init__(self, properties: collections.abc.Mapping[str, typing.SupportsFloat | typing.SupportsIndex]) -> None:
        """
                 Construct this atom from a properties dictionary.
        
                 The dictionary should contain at least the following keys:
                 - "id": Unique identifier for the atom
                 - "type": Type classification of the atom  
                 - "x", "y", "z": Coordinate positions
                 - "nx", "ny", "nz": Periodic image flags
                 
                 Any additional properties will be stored as extra data.
        
                 :param properties: Dictionary containing atom properties
        """
    def __repr__(self) -> str:
        """
                 Return a string representation of the Atom.
                 
                 :return: String representation showing atom properties
        """
    def __setstate__(self, arg0: tuple) -> None:
        """
        Provides support for pickling
        """
    def compute_vector_to(self, to_atom: Atom, pbc_box: Box) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                    Compute the vector to another atom.
                    
                    :param to_atom: The target atom
                    :param pbc_box: The periodic boundary conditions box
                    :return: Vector pointing from this atom to the target atom
        """
    def distance_to(self, to_atom: Atom, pbc_box: Box) -> float:
        """
                    Compute the distance to another atom.
                    
                    :param to_atom: The target atom
                    :param pbc_box: The periodic boundary conditions box
                    :return: Euclidean distance between the atoms
        """
    def distance_to_unwrapped(self, arg0: Atom, arg1: Box) -> float:
        """
                 Compute the distance to another atom respecting the periodic image flags.
                 
                 :param to_atom: The target atom
                 :return: Unwrapped distance to the target atom
        """
    def get_coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                 Get the coordinates of this atom as a vector.
                 
                 :return: A vector containing the x, y, z coordinates
        """
    def get_extra_data(self) -> dict[str, float]:
        """
                 Get all extra data properties stored with this atom (e.g., charge, dipole, etc.).
                 
                 :return: Dictionary containing all extra properties
        """
    def get_id(self) -> int:
        """
                    Get the ID of the atom.
                    
                    :return: The atom's unique identifier
        """
    def get_nx(self) -> int:
        """
                 Get the box image that the atom is in in x direction (also known as `ix` or `nx`).
                 
                 :return: The periodic image flag in x direction
        """
    def get_ny(self) -> int:
        """
                 Get the box image that the atom is in in y direction (also known as `iy` or `ny`).
                 
                 :return: The periodic image flag in y direction
        """
    def get_nz(self) -> int:
        """
                 Get the box image that the atom is in in z direction (also known as `iz` or `nz`).
                 
                 :return: The periodic image flag in z direction
        """
    def get_property(self, property: str) -> float:
        """
                 Get a specific property value from the extra data.
                 
                 :param property: The name of the property to retrieve
                 :return: The value of the specified property
                 :raises: std::out_of_range if the property doesn't exist
        """
    def get_type(self) -> int:
        """
                    Get the type of the atom.
                    
                    :return: The atom's type classification
        """
    def get_unwrapped_coordinates(self, arg0: Box) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                 Get the unwrapped coordinates of this atom.
                 
                 :param box: The simulation box to use for unwrapping
                 :return: A vector containing the unwrapped x, y, z coordinates
        """
    def get_unwrapped_x(self, box: Box) -> float:
        """
                 Get the unwrapped x coordinate of the atom.
                 
                 :param box: The simulation box to use for unwrapping
                 :return: The unwrapped x coordinate
        """
    def get_unwrapped_y(self, box: Box) -> float:
        """
                 Get the unwrapped y coordinate of the atom.
                 
                 :param box: The simulation box to use for unwrapping
                 :return: The unwrapped y coordinate
        """
    def get_unwrapped_z(self, box: Box) -> float:
        """
                 Get the unwrapped z coordinate of the atom.
                 
                 :param box: The simulation box to use for unwrapping
                 :return: The unwrapped z coordinate
        """
    def get_x(self) -> float:
        """
                    Get the x coordinate of the atom.
                    
                    :return: The x coordinate
        """
    def get_y(self) -> float:
        """
                    Get the y coordinate of the atom.
                    
                    :return: The y coordinate
        """
    def get_z(self) -> float:
        """
                    Get the z coordinate of the atom.
                    
                    :return: The z coordinate
        """
    def mean_position_with(self, other_atom: Atom, pbc_box: Box) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                 Compute the mean position between this atom and another atom, considering periodic boundaries.
                 
                 :param other_atom: The other atom
                 :param pbc_box: The periodic boundary conditions box
                 :return: Vector representing the mean position
        """
    def mean_position_with_unwrapped(self, other_atom: Atom, pbc_box: Box) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                 Compute the mean position between this atom and another atom using unwrapped coordinates.
                 
                 :param other_atom: The other atom
                 :param pbc_box: The periodic boundary conditions box  
                 :return: Vector representing the mean position (unwrapped)
        """
    def vector_to_unwrapped(self, arg0: Atom, arg1: Box) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                 Compute the vector to another atom respecting the periodic image flags.
                 
                 :param to_atom: The target atom
                 :return: Unwrapped vector to the target atom
        """
class AtomPairEntanglements:
    """
    
          A struct to store pairs of atoms that are close together and could be entanglements.
        
    """
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
    def pair_of_atom(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def pairs_of_atoms(self) -> list[tuple[int, int]]:
        """
              A list of pairs of atom ids that are close together and could be entanglements
        """
    @pairs_of_atoms.setter
    def pairs_of_atoms(self, arg0: collections.abc.Sequence[tuple[typing.SupportsInt | typing.SupportsIndex, typing.SupportsInt | typing.SupportsIndex]]) -> None:
        ...
class AtomStyle(enum.Enum):
    """
    An enumeration of the LAMMPS atom styles.
    """
    ANGLE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.ANGLE: 1>
    ATOMIC: typing.ClassVar[AtomStyle]  # value = <AtomStyle.ATOMIC: 2>
    BODY: typing.ClassVar[AtomStyle]  # value = <AtomStyle.BODY: 3>
    BOND: typing.ClassVar[AtomStyle]  # value = <AtomStyle.BOND: 4>
    BPM_SPHERE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.BPM_SPHERE: 5>
    CHARGE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.CHARGE: 6>
    DIELECTRIC: typing.ClassVar[AtomStyle]  # value = <AtomStyle.DIELECTRIC: 7>
    DIPOLE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.DIPOLE: 8>
    DPD: typing.ClassVar[AtomStyle]  # value = <AtomStyle.DPD: 9>
    EDPD: typing.ClassVar[AtomStyle]  # value = <AtomStyle.EDPD: 10>
    ELECTRON: typing.ClassVar[AtomStyle]  # value = <AtomStyle.ELECTRON: 11>
    ELLIPSOID: typing.ClassVar[AtomStyle]  # value = <AtomStyle.ELLIPSOID: 12>
    FULL: typing.ClassVar[AtomStyle]  # value = <AtomStyle.FULL: 13>
    HYBRID: typing.ClassVar[AtomStyle]  # value = <AtomStyle.HYBRID: 28>
    LINE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.LINE: 14>
    MDPD: typing.ClassVar[AtomStyle]  # value = <AtomStyle.MDPD: 15>
    MOLECULAR: typing.ClassVar[AtomStyle]  # value = <AtomStyle.MOLECULAR: 16>
    NONE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.NONE: 0>
    PERI: typing.ClassVar[AtomStyle]  # value = <AtomStyle.PERI: 17>
    RHEO: typing.ClassVar[AtomStyle]  # value = <AtomStyle.RHEO: 18>
    RHEO_THERMAL: typing.ClassVar[AtomStyle]  # value = <AtomStyle.RHEO_THERMAL: 19>
    SMD: typing.ClassVar[AtomStyle]  # value = <AtomStyle.SMD: 20>
    SPH: typing.ClassVar[AtomStyle]  # value = <AtomStyle.SPH: 21>
    SPHERE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.SPHERE: 22>
    SPIN: typing.ClassVar[AtomStyle]  # value = <AtomStyle.SPIN: 23>
    TDPD: typing.ClassVar[AtomStyle]  # value = <AtomStyle.TDPD: 24>
    TEMPLATE: typing.ClassVar[AtomStyle]  # value = <AtomStyle.TEMPLATE: 25>
    TRI: typing.ClassVar[AtomStyle]  # value = <AtomStyle.TRI: 26>
    WAVEPACKET: typing.ClassVar[AtomStyle]  # value = <AtomStyle.WAVEPACKET: 27>
class AveFileReader:
    """
    
              Alternative implementation of the data file reader implemented in 
              :func:`pylimer_tools.io.read_lammps_output_file.read_averages_file`.
    
              This implementation is better for certain use cases, worse for others.
              In the end, only performance and memory usage are different.
              For moderately sized and small files, we recommend to use the Python interface instead.
         
    """
    def __init__(self, file_path: str) -> None:
        """
                 Initialize the AveFileReader with a file path.
                 
                 :param file_path: Path to the averages file to read
        """
    def __repr__(self) -> str:
        """
                 Return a string representation of the AveFileReader.
                 
                 :return: String representation showing file dimensions
        """
    def autocorrelate_column(self, column_index: typing.SupportsInt | typing.SupportsIndex, delta_indices: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> list[float]:
        """
                  Do autocorrelation on one particular column for a specified set of delta indices.
        
                  Assumes the data is equally spaced.
                  
                  :param column_index: Index of the column to autocorrelate
                  :param delta_indices: List of delta indices for the autocorrelation
                  :return: Autocorrelation values for the specified deltas
        """
    def autocorrelate_column_difference(self, column_index1: typing.SupportsInt | typing.SupportsIndex, column_index2: typing.SupportsInt | typing.SupportsIndex, delta_indices: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> list[float]:
        """
                  Do autocorrelation on the difference between two particular columns for a specified set of delta indices.
        
                  Assumes the data is equally spaced.
                  
                  :param column_index1: Index of the first column
                  :param column_index2: Index of the second column  
                  :param delta_indices: List of delta indices for the autocorrelation
                  :return: Autocorrelation values for the column differences at specified deltas
        """
    def get_column_names(self) -> list[str]:
        """
                 Get the names of all columns in the file.
                 
                 :return: List of column names
        """
    def get_data(self) -> list[list[float]]:
        """
                 Get all data from the file.
                 
                 :return: 2D array containing all the numerical data
        """
    def get_nr_of_columns(self) -> int:
        """
                 Get the number of columns in the file.
                 
                 :return: Number of columns
        """
    def get_nr_of_rows(self) -> int:
        """
                 Get the number of data rows in the file.
                 
                 :return: Number of rows
        """
class BackTrackStatus(enum.Enum):
    """
    
         Enum for controlling the strand linking process in linkStrandsCallback.
         
    """
    STOP: typing.ClassVar[BackTrackStatus]  # value = <BackTrackStatus.STOP: 0>
    TRACK_BACKWARD: typing.ClassVar[BackTrackStatus]  # value = <BackTrackStatus.TRACK_BACKWARD: 2>
    TRACK_FORWARD: typing.ClassVar[BackTrackStatus]  # value = <BackTrackStatus.TRACK_FORWARD: 1>
class BondLengthStats:
    """
    
        Summary statistics for realized bond lengths.
      
    """
    def __init__(self) -> None:
        ...
    @property
    def max_bond_length(self) -> float:
        ...
    @property
    def mean_bond_length(self) -> float:
        ...
    @property
    def mean_squared_bond_length(self) -> float:
        ...
    @property
    def min_bond_length(self) -> float:
        ...
    @property
    def nr_of_bonds(self) -> int:
        ...
class Box:
    """
    
            The box that the simulation is run in.
    
            .. note:: 
              Currently, only rectangular boxes are supported.
            
    """
    def __getstate__(self) -> tuple[float, float, float, float, float, float, float, int]:
        ...
    @typing.overload
    def __init__(self, arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex, arg2: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex, arg2: typing.SupportsFloat | typing.SupportsIndex, arg3: typing.SupportsFloat | typing.SupportsIndex, arg4: typing.SupportsFloat | typing.SupportsIndex, arg5: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def __repr__(self) -> str:
        """
                 Return a string representation of the Box.
                 
                 :return: String representation showing box dimensions and shear parameters
        """
    def __setstate__(self, arg0: tuple) -> None:
        """
        Provides support for pickling.
        """
    def apply_pbc(self, distances: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
              Apply periodic boundary conditions (PBC): adjust the specified distances to fit into this box.
              
              :param distances: The distances to adjust
              :return: The adjusted distances
        """
    def apply_simple_shear(self, shear_magnitude: typing.SupportsFloat | typing.SupportsIndex, shear_direction: typing.SupportsInt | typing.SupportsIndex = 0) -> None:
        """
                  Apply a simple shear to the box.
        
                  .. warning::
                    Currently, this is not supported for all operations.
        
                  For shear magnitude, you specify the angle :math:`\\gamma`.
        
                  :param shear_magnitude: The shear magnitude (angle :math:`\\gamma`)
                  :param shear_direction: Direction of shear: 0 for x, 1 for y, 2 for z. 
                                         Use any other integer to disable shear.
        """
    def get_bounding_box(self) -> Box:
        """
             Get an orthogonal box that encloses this box.
             
             For non-sheared boxes, the resulting box is identical to the current box.
        
             :return: A new Box object representing the bounding box
        """
    def get_high_x(self) -> float:
        """
                    Get the upper bound of the box in x direction.
        
                    :return: The upper x-coordinate boundary
        """
    def get_high_y(self) -> float:
        """
                    Get the upper bound of the box in y direction.
        
                    :return: The upper y-coordinate boundary
        """
    def get_high_z(self) -> float:
        """
                    Get the upper bound of the box in z direction.
        
                    :return: The upper z-coordinate boundary
        """
    def get_l(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                  Get the three lengths of the box in an array/list.
        
                  :return: Array containing [Lx, Ly, Lz] box dimensions
        """
    def get_low_x(self) -> float:
        """
                    Get the lower bound of the box in x direction.
        
                    :return: The lower x-coordinate boundary
        """
    def get_low_y(self) -> float:
        """
                    Get the lower bound of the box in y direction.
        
                    :return: The lower y-coordinate boundary
        """
    def get_low_z(self) -> float:
        """
                    Get the lower bound of the box in z direction.
        
                    :return: The lower z-coordinate boundary
        """
    def get_lx(self) -> float:
        """
                    Get the length of the box in x direction.
        
                    :return: The x-dimension length of the box
        """
    def get_ly(self) -> float:
        """
                    Get the length of the box in y direction.
        
                    :return: The y-dimension length of the box
        """
    def get_lz(self) -> float:
        """
                    Get the length of the box in z direction.
        
                    :return: The z-dimension length of the box
        """
    def get_offset(self, distances: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
             Compute the offset required to compensate for periodic boundary conditions.
        
             Useful e.g. if you are using absolute coordinates for distances, but 
             still need an infinite network, 
             e.g., if the bonds need to be able to get longer than half the box.
             
             :param distances: The distances to compute offset for
             :return: The computed offset
        """
    def get_volume(self) -> float:
        """
                    Compute the volume of the box.
        
                    :math:`V = L_x \\cdot L_y \\cdot L_z`
        
                    :return: The volume of the box
        """
    def is_valid_offset(self, potential_offset: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], abs_precision: typing.SupportsFloat | typing.SupportsIndex = 1e-05) -> bool:
        """
                  Check whether the passed offset is a valid one in this box.
                  
                  :param potential_offset: The offset to validate
                  :param abs_precision: Absolute precision for the validation
                  :return: True if the offset is valid, False otherwise
        """
class ComputeBackend(enum.Enum):
    """
    Hardware compute backend for approximate spectral methods
    """
    AUTO: typing.ClassVar[ComputeBackend]  # value = <ComputeBackend.AUTO: 0>
    CPU: typing.ClassVar[ComputeBackend]  # value = <ComputeBackend.CPU: 1>
    CUDA: typing.ClassVar[ComputeBackend]  # value = <ComputeBackend.CUDA: 3>
    METAL: typing.ClassVar[ComputeBackend]  # value = <ComputeBackend.METAL: 2>
class ComputedDoubleValues(enum.IntEnum):
    """
    Floating point output quantities
    """
    GAMMA: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.GAMMA: 14>
    MAX_B: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.MAX_B: 17>
    MEAN_B: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.MEAN_B: 16>
    MSD: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.MSD: 18>
    PRESSURE: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.PRESSURE: 3>
    RESIDUAL: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.RESIDUAL: 15>
    STRESS_NXY: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.STRESS_NXY: 11>
    STRESS_NXZ: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.STRESS_NXZ: 13>
    STRESS_NYZ: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.STRESS_NYZ: 12>
    STRESS_XX: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.STRESS_XX: 5>
    STRESS_XY: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.STRESS_XY: 8>
    STRESS_XZ: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.STRESS_XZ: 10>
    STRESS_YY: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.STRESS_YY: 6>
    STRESS_YZ: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.STRESS_YZ: 9>
    STRESS_ZZ: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.STRESS_ZZ: 7>
    TEMPERATURE: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.TEMPERATURE: 4>
    TIME: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.TIME: 1>
    TIMESTEP: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.TIMESTEP: 0>
    VOLUME: typing.ClassVar[ComputedDoubleValues]  # value = <ComputedDoubleValues.VOLUME: 2>
    @classmethod
    def __new__(cls, value):
        ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """
class ComputedIntValues(enum.IntEnum):
    """
    Integer output quantities
    """
    NUM_ATOMS: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_ATOMS: 3>
    NUM_BONDS: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_BONDS: 5>
    NUM_BONDS_TO_FORM: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_BONDS_TO_FORM: 7>
    NUM_EXTRA_ATOMS: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_EXTRA_ATOMS: 4>
    NUM_EXTRA_BONDS: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_EXTRA_BONDS: 6>
    NUM_RELOC: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_RELOC: 2>
    NUM_SHIFT: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_SHIFT: 1>
    STEP: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.STEP: 0>
    @classmethod
    def __new__(cls, value):
        ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """
class Correlator:
    """
    
        Multi-tau scalar autocorrelator computing :math:`f(\\tau) = \\langle A(t) A(t + \\tau) \\rangle`.
    
        Uses a hierarchy of coarse-grained shift registers (Likhtman/Ramirez multi-tau algorithm)
        to evaluate time correlations across exponentially spanning lag times in O(1) amortized time per sample.
        
    """
    def __getstate__(self) -> tuple[bytes]:
        ...
    def __init__(self, num_correlators: typing.SupportsInt | typing.SupportsIndex = 32, p: typing.SupportsInt | typing.SupportsIndex = 16, m: typing.SupportsInt | typing.SupportsIndex = 2) -> None:
        """
                 Construct a multi-tau scalar Correlator.
        
                 :param num_correlators: Number of hierarchical correlator levels (default: 32)
                 :param p: Number of lag points per correlator level (default: 16)
                 :param m: Averaging factor between consecutive levels; recommended ``p % m == 0`` (default: 2)
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def add(self, value: typing.SupportsFloat | typing.SupportsIndex, level: typing.SupportsInt | typing.SupportsIndex = 0) -> None:
        """
                 Add a scalar observation to correlator level ``level``.
        
                 :param value: Finite scalar observation
                 :param level: Correlator level index (default: 0)
        """
    def add_batch(self, values: typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]) -> None:
        """
                 Add a 1D array of scalar observations at level 0.
        
                 :param values: 1D array of finite scalar observations
        """
    def evaluate(self, norm: bool = False) -> None:
        """
                 Evaluate the autocorrelation function at all populated lag times.
        
                 :param norm: If True, subtract the squared sample mean :math:`\\langle A \\rangle^2` (default: False)
        """
    def reset(self) -> None:
        """
        Reset accumulated state while keeping the current correlator dimensions.
        """
    def set_size(self, num_correlators: typing.SupportsInt | typing.SupportsIndex = 32, p: typing.SupportsInt | typing.SupportsIndex = 16, m: typing.SupportsInt | typing.SupportsIndex = 2) -> None:
        """
                 Reconfigure correlator dimensions and reset accumulated state.
        
                 :param num_correlators: Number of hierarchical correlator levels (default: 32)
                 :param p: Number of lag points per correlator level (default: 16)
                 :param m: Averaging factor between consecutive levels (default: 2)
        """
    @property
    def accval(self) -> float:
        """
        Accumulated sum of level-0 samples.
        """
    @property
    def correlation(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Evaluated correlation values (call evaluate() first).
        """
    @property
    def m(self) -> int:
        """
        Averaging factor between consecutive levels.
        """
    @property
    def mean(self) -> float:
        """
        Running sample mean at level 0.
        """
    @property
    def npcorr(self) -> int:
        """
        Number of populated correlation points after evaluate().
        """
    @property
    def num_correlators(self) -> int:
        """
        Number of correlator levels configured.
        """
    @property
    def num_samples(self) -> int:
        """
        Total number of samples added at level 0.
        """
    @property
    def p(self) -> int:
        """
        Number of lag points per correlator level.
        """
    @property
    def raw_f(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Full pre-allocated correlation buffer.
        """
    @property
    def raw_t(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Full pre-allocated lag time buffer.
        """
    @property
    def times(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Evaluated lag times (call evaluate() first).
        """
class CrosslinkerUniverse:
    """
    
          A struct representing the topological and spatial state of crosslinkers and strands
          without explicitly instantiating all individual strand beads.
        
    """
    def __init__(self) -> None:
        """
        Get an instance of this struct
        """
    def classify_strand(self, strand_idx: typing.SupportsInt | typing.SupportsIndex) -> StrandCategory:
        """
              Classify a strand into its topological StrandCategory.
        
              :param strand_idx: Index of the strand.
              :return: StrandCategory of the strand.
        """
    def compute_bond_length_stats(self, box: Box, n_mc_steps: typing.SupportsInt | typing.SupportsIndex = 0) -> BondLengthStats:
        """
              Stream over all strands in O(1) strand memory to compute realized bond length statistics.
        
              :param box: Simulation box.
              :param n_mc_steps: Number of MC equilibration steps per strand (default: 0).
              :return: BondLengthStats object.
        """
    def compute_end_to_end_stats(self, box: Box, category: StrandCategory = ..., n_mc_steps: typing.SupportsInt | typing.SupportsIndex = 0) -> EndToEndStats:
        """
              Stream over strands matching `category` to compute realized end-to-end distance statistics.
        
              :param box: Simulation box.
              :param category: StrandCategory filter (default: NETWORK_STRAND).
              :param n_mc_steps: Number of MC equilibration steps per strand (default: 0).
              :return: EndToEndStats object.
        """
    def compute_strand_bond_counts(self, category: StrandCategory = ...) -> list[int]:
        """
              Return the number of bonds for each strand matching `category`.
        
              :param category: StrandCategory filter (default: ALL).
              :return: List of bond counts per matching strand.
        """
    def count_atom_types(self) -> dict[int, int]:
        """
              Count atoms per atom type in O(N_xlinks + N_strands).
        
              :return: Dict mapping atom type to count.
        """
    def count_atoms_by_degree(self) -> dict[int, int]:
        """
              Count atoms per bond degree in O(N_xlinks + N_strands).
        
              :return: Dict mapping degree to atom count.
        """
    def determine_effective_functionality_per_type(self) -> dict[int, float]:
        """
              Compute effective functionality (mean degree) per atom type in O(N_xlinks + N_strands).
        
              :return: Dict mapping atom type to mean degree.
        """
    @property
    def bead_distance_in_strand(self) -> list[float]:
        ...
    @bead_distance_in_strand.setter
    def bead_distance_in_strand(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def beads_in_strand(self) -> list[int]:
        ...
    @beads_in_strand.setter
    def beads_in_strand(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def mean_squared_bead_distance_in_strand(self) -> list[float]:
        ...
    @mean_squared_bead_distance_in_strand.setter
    def mean_squared_bead_distance_in_strand(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def strand_bead_type(self) -> list[int]:
        ...
    @strand_bead_type.setter
    def strand_bead_type(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def strand_from(self) -> list[int]:
        ...
    @strand_from.setter
    def strand_from(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def strand_seeds(self) -> list[int]:
        ...
    @strand_seeds.setter
    def strand_seeds(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def strand_to(self) -> list[int]:
        ...
    @strand_to.setter
    def strand_to(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def strands_of_xlink(self) -> list[list[int]]:
        ...
    @strands_of_xlink.setter
    def strands_of_xlink(self, arg0: collections.abc.Sequence[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]]) -> None:
        ...
    @property
    def xlink_chain_id(self) -> list[int]:
        ...
    @xlink_chain_id.setter
    def xlink_chain_id(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def xlink_types(self) -> list[int]:
        ...
    @xlink_types.setter
    def xlink_types(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def xlink_x(self) -> list[float]:
        ...
    @xlink_x.setter
    def xlink_x(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def xlink_y(self) -> list[float]:
        ...
    @xlink_y.setter
    def xlink_y(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def xlink_z(self) -> list[float]:
        ...
    @xlink_z.setter
    def xlink_z(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
class DPDSimulator:
    """
    
         A quick-and-dirty implementation of the dissipative particle dynamics (DPD) simulation
         with slip-springs as presented by :cite:t:`langeloth_recovering_2013` 
         and :cite:t:`schneider_simulation_2021`.
         
    """
    @staticmethod
    def read_restart_file(file: str) -> DPDSimulator:
        """
                  Read a restart file in order to continue a simulation.
        
                  :param file: The file path to the restart file to read
        """
    def __init__(self, universe: Universe, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, slipspring_bond_type: typing.SupportsInt | typing.SupportsIndex = 9, is_2d: bool = False, seed: str = '', assume_box_large_enough: bool = False) -> None:
        """
        Get an instance of this class
        """
    def assume_box_large_enough(self, assume: bool = True) -> None:
        """
                  Configure whether to assume minimum-image PBC is sufficient for all bonds.
        
                  If your bonds could exceed half the box length, keep this ``False`` (default) so per-bond periodic box offsets are tracked.
                  Set to ``True`` to apply minimum-image PBC on every bond evaluation.
        """
    def assumes_box_large_enough(self) -> bool:
        """
                  Return ``True`` if minimum-image PBC is assumed sufficient for bonds, ``False`` if per-bond periodic box offsets are tracked.
        """
    def config_a(self, A: typing.SupportsFloat | typing.SupportsIndex = 25.0) -> None:
        """
                  Configure the force-field (pair-style) parameter `A`.
        """
    def config_allow_relocation_in_network(self, allow_relocation_in_network: bool = False) -> None:
        """
                  Configure whether a relocation step may happen when a slip-spring has ended at a crosslink.
        
                  Side-effect: if true, the relocations may also happen *to* a slip-spring next to a crosslink.
        
                  :param allow_relocation_in_network (bool): Whether to allow relocation in the network or not.
        """
    def config_assume_box_large_enough(self, assume: bool = True) -> None:
        """
                  Configure whether to assume minimum-image PBC is sufficient for all bonds.
        
                  If your bonds could exceed half the box length, keep this ``False`` (default) so per-bond periodic box offsets are tracked.
                  Set to ``True`` to apply minimum-image PBC on every bond evaluation.
        """
    def config_auto_correlator_output(self, values: collections.abc.Sequence[OutputConfiguration], num_corr_in: typing.SupportsInt | typing.SupportsIndex = 32, p: typing.SupportsInt | typing.SupportsIndex = 16, m: typing.SupportsInt | typing.SupportsIndex = 2) -> None:
        """
                  Set which values to compute multiple-tau autocorrelation for.
                  If you use this, you should cite :cite:t:`ramirez_efficient_2010`.
        
                  :param values: a list of OutputConfiguration structs
                  :param num_corr_in: Number of correlations in
                  :param p: Parameter p for the autocorrelator
                  :param m: Parameter m for the autocorrelator
        """
    def config_average_output(self, values: collections.abc.Sequence[OutputConfiguration]) -> None:
        """
                  Set which values to compute averages for.
        
                  :param values: A list of OutputConfiguration structs specifying what to average
        """
    def config_bond_formation(self, num_bonds_to_form: typing.SupportsInt | typing.SupportsIndex, max_bonds_per_atom_type: collections.abc.Mapping[typing.SupportsInt | typing.SupportsIndex, typing.SupportsInt | typing.SupportsIndex], bond_formation_dist: typing.SupportsFloat | typing.SupportsIndex = 1.0, attempt_bond_formation_every: typing.SupportsInt | typing.SupportsIndex = 50, atom_type_form_from: typing.SupportsInt | typing.SupportsIndex = 2, atom_type_form_to: typing.SupportsInt | typing.SupportsIndex = 1) -> None:
        """
                  Configure how to do bond formation during the run.
        
                  :param num_bonds_to_form (int): The nr of bonds to form in total. Use 0 to stop bond formation.
                  :param num_bonds_per_atom_type (dict): The nr of bonds each atom type may have at most (e.g., 2 for strand atoms, 4 for a tertiary crosslinkers)
                  :param bond_formation_dist (float): The maximum distance allowed to form bonds
                  :param attempt_bond_formation_every (int): attempt to form bonds every this many steps during the simulation run
                  :param atom_type_form_from (int): The atom type to start forming bonds from.
                  :param atom_type_form_to (int): The atom type to start forming bonds to.
        """
    def config_box_deformation(self, target_box: Box) -> None:
        """
                  Configure where to (incrementally) deform the box to during the next simulation run.
        """
    def config_gamma(self, gamma: typing.SupportsFloat | typing.SupportsIndex = 4.5) -> None:
        """
                  Configure the DPD dissipative friction parameter `\\gamma` independently of `\\sigma`.
        """
    def config_lambda(self, l: typing.SupportsFloat | typing.SupportsIndex = 0.65) -> None:
        """
                  Configure the modified velocity verlet integration parameter `\\lambda`.
        """
    def config_num_steps_dpd(self, num_steps: typing.SupportsInt | typing.SupportsIndex = 500) -> None:
        """
                  Configure the number of steps to do in one DPD sequence.
        
                  :param num_steps: Number of DPD steps per sequence
        """
    def config_num_steps_mc(self, num_steps: typing.SupportsInt | typing.SupportsIndex = 500) -> None:
        """
                  Configure the number of steps to do in one MC sequence.
        
                  :param num_steps: Number of Monte Carlo steps per sequence
        """
    def config_restart_output(self, file: str, output_every: typing.SupportsInt | typing.SupportsIndex = 50000) -> None:
        """
                  Set when to output a restart file.
        
                  Note:
                       The filename determines the type of serialization:
                       .json, .xml are supported; other file endings will lead to binary serialization (fastest!).
        
                  Caution:
                       This method may not be backwards- nor forward-compatible.
                       Use the same version of pylimer-tools if you want to be sure that things work.
        
                  :param file: The file path to the restart file to write
                  :param output_every: How often to write the restart file (default: 50000)
        """
    def config_shift_one_at_a_time(self, shift_one_at_a_time: bool = False) -> None:
        """
                  Configure whether to shift atoms one at a time.
        
                  This setting affects Monte Carlo move behavior in the simulation.
        """
    def config_shift_possibility_empty(self, shift_possibility_empty: bool = True) -> None:
        """
                  Configure the possibility of shifting to empty positions.
        
                  This setting affects Monte Carlo moves in the simulation.
        """
    def config_sigma(self, sigma: typing.SupportsFloat | typing.SupportsIndex = 3.0) -> None:
        """
                  Configure the force-field (pair-style) parameter `\\sigma` (and set `\\gamma = 0.5 \\sigma^2`).
        """
    def config_slipspring_high_cutoff(self, cutoff: typing.SupportsFloat | typing.SupportsIndex = 2.0) -> None:
        """
                  Configure the higher cut-off of how far a pair may be distanced for a slip-spring to be created.
        """
    def config_slipspring_low_cutoff(self, cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.5) -> None:
        """
                  Configure the lower cut-off of how far a pair may be distanced for a slip-spring to be created.
        """
    def config_spring_constant(self, k: typing.SupportsFloat | typing.SupportsIndex = 2.0) -> None:
        """
                  Configure the force-field (bond-style) parameter `k`, the spring constant.
        """
    def config_step_output(self, values: collections.abc.Sequence[OutputConfiguration]) -> None:
        """
                  Set which values to log.
        
                  :param values: a list of OutputConfiguration structs
        """
    def create_slip_springs(self, num: typing.SupportsInt | typing.SupportsIndex, bond_type: typing.SupportsInt | typing.SupportsIndex = 0) -> int:
        """
                  Randomly add up to the specified number of slip-springs to eligible neighbours within the configured cut-offs.
        
                  :param num: Maximum number of slip-springs to create
                  :param bond_type: Bond type for created slip-springs (``0`` inherits ``slipspring_bond_type``)
                  :return: Number of slip-springs actually created
        """
    def get_bond_lengths(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Get the lengths of all bonds in the system.
        
                  :return: Vector containing the length of each bond
        """
    def get_coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Get the current particle coordinates.
        
                  :return: Vector of particle coordinates (x1,y1,z1,x2,y2,z2,...)
        """
    def get_current_timestep(self) -> int:
        """
                  Get the current timestep number.
        
                  :return: The current timestep index
        """
    def get_dpd_gamma(self) -> float:
        """
                  Get the configured DPD dissipative friction coefficient `\\gamma`.
        """
    def get_nr_of_atoms(self) -> int:
        """
                  Get the total number of atoms in the system.
        
                  :return: Total number of atoms
        """
    def get_nr_of_bonds(self) -> int:
        """
                  Get the number of regular bonds (excluding slip-springs).
        
                  :return: Number of regular bonds
        """
    def get_nr_of_bonds_to_form(self) -> int:
        """
                  Get the number of bonds that are configured to have to be formed.
        """
    def get_nr_of_extra_atoms(self) -> int:
        """
                  Get the number of extra atoms (always 0 for DPD simulations).
        
                  :return: Number of extra atoms
        """
    def get_nr_of_extra_bonds(self) -> int:
        """
                  Get the number of extra bonds (slip-springs).
        
                  :return: Number of slip-springs
        """
    def get_nr_of_slip_springs(self) -> int:
        """
                  Get the current number of slip-springs in the system.
        
                  :return: Number of slip-springs
        """
    def get_nr_of_steps_dpd(self) -> int:
        """
                  Get the configured number of DPD steps per sequence.
        
                  :return: Number of DPD steps
        """
    def get_nr_of_steps_mc(self) -> int:
        """
                  Get the configured number of Monte Carlo steps per sequence.
        
                  :return: Number of MC steps
        """
    def get_num_relocations(self) -> int:
        """
                  Get the number of accepted slip-spring relocation moves in the most recent MC sequence.
        """
    def get_num_shifts(self) -> int:
        """
                  Get the number of accepted slip-spring shift moves in the most recent MC sequence.
        """
    def get_shift_one_at_a_time(self) -> bool:
        """
                  Get whether slip-springs are shifted one at a time.
        
                  :return: True if shifting one at a time, False otherwise
        """
    def get_shift_possibility_empty(self) -> bool:
        """
                  Get whether shifting to empty positions is allowed.
        
                  :return: True if shifting to empty positions is allowed
        """
    def get_sigma(self) -> float:
        """
                  Get the configured DPD random force noise amplitude `\\sigma`.
        """
    def get_slip_spring_bond_type(self) -> int:
        """
                  Get the bond type identifier used for slip-springs.
        
                  :return: Bond type for slip-springs
        """
    def get_spring_constant(self) -> float:
        """
                  Get the current spring constant value.
        
                  :return: The current spring constant for bond interactions
        """
    def get_stress_tensor(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 3]"]:
        """
                  Get the current stress tensor.
        
                  :return: 3x3 stress tensor matrix
        """
    def get_target_kbt(self) -> float:
        """
                  Get the target thermal energy :math:`k_B T = \\sigma^2 / (2 \\gamma)` configured on the DPD thermostat.
        
                  :return: Target thermal energy :math:`k_B T`
        """
    def get_temperature(self) -> float:
        """
                  Get the current system temperature.
        
                  :return: The current temperature of the system
        """
    def get_timestep(self) -> float:
        """
                  Get the timestep used in the simulation.
        
                  :return: The simulation timestep value
        """
    def get_universe(self, with_slipsprings: bool = True) -> Universe:
        """
             Get a universe instance from the current coordinates (and connectivity).
        
             :param with_slipsprings: Whether to include slip-springs in the returned universe (default: True)
        """
    def get_velocities(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Get the current particle velocities.
        
                  :return: Vector of particle velocities (vx1,vy1,vz1,vx2,vy2,vz2,...)
        """
    def get_volume(self) -> float:
        """
                  Get the current system volume.
        
                  :return: Current simulation box volume
        """
    def refresh_current_state(self) -> None:
        """
                  After re-configuring the force-field parameters,
                  this method should be called to update the current stress tensor etc.
        """
    def relocate_slip_springs(self, kbt: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> int:
        """
                  Perform one sweep of Monte Carlo slip-spring chain-end relocation moves.
        
                  :param kbt: Thermal energy scale :math:`k_B T` for Metropolis acceptance (must be ``>= 0``)
                  :return: Number of accepted relocation moves
        """
    def run_simulation(self, n_steps: typing.SupportsInt | typing.SupportsIndex, dt: typing.SupportsFloat | typing.SupportsIndex = 0.06, with_MC: bool = False) -> None:
        ...
    def shift_slip_springs(self, kbt: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> int:
        """
                  Perform one sweep of Monte Carlo slip-spring shift moves along polymer chains.
        
                  :param kbt: Thermal energy scale :math:`k_B T` for Metropolis acceptance (must be ``>= 0``)
                  :return: Number of accepted shift moves
        """
    def start_measuring_msd_for_atoms(self, atom_ids: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        """
                  Set a new origin for measuing the mean square displacement for a specified set of atoms
        """
    def validate_neighbour_list(self, cutoff: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
                  Validate the neighbor list consistency for debugging purposes.
        
                  :param cutoff: Cutoff distance for validation
        """
    def validate_state(self) -> None:
        """
                  Validate the current simulation state for debugging purposes.
                  
                  Checks internal data structure consistency and throws exceptions if issues are found.
        """
    def write_restart_file(self, file: str) -> None:
        """
                  Explicitly force the writing of a restart file, now!
        
                  :param file: The file path and name of the restart file to be written.
                               Can end in .xml, .json or anything else (-> binary)
        """
class DataFileReader:
    """
    
           A reader for LAMMPS's `write_data` files.
      
    """
    def __init__(self) -> None:
        """
                 Initialize a new data file parser.
        """
    def get_atom_ids(self) -> list[int]:
        """
                 Get all atom IDs from the data file.
                 
                 :return: Vector of atom IDs
        """
    def get_atom_nx(self) -> list[int]:
        """
                 Get periodic image flags in x direction for all atoms.
                 
                 :return: Vector of nx values (image flags)
        """
    def get_atom_ny(self) -> list[int]:
        """
                 Get periodic image flags in y direction for all atoms.
                 
                 :return: Vector of ny values (image flags)
        """
    def get_atom_nz(self) -> list[int]:
        """
                 Get periodic image flags in z direction for all atoms.
                 
                 :return: Vector of nz values (image flags)
        """
    def get_atom_types(self) -> list[int]:
        """
                 Get all atom types from the data file.
                 
                 :return: Vector of atom types
        """
    def get_atom_x(self) -> list[float]:
        """
                 Get x coordinates of all atoms.
                 
                 :return: Vector of x coordinates
        """
    def get_atom_y(self) -> list[float]:
        """
                 Get y coordinates of all atoms.
                 
                 :return: Vector of y coordinates
        """
    def get_atom_z(self) -> list[float]:
        """
                 Get z coordinates of all atoms.
                 
                 :return: Vector of z coordinates
        """
    def get_bond_from(self) -> list[int]:
        """
                 Get starting atom IDs for all bonds.
                 
                 :return: Vector of starting atom IDs for bonds
        """
    def get_bond_to(self) -> list[int]:
        """
                 Get ending atom IDs for all bonds.
                 
                 :return: Vector of ending atom IDs for bonds
        """
    def get_bond_types(self) -> list[int]:
        """
                 Get all bond types from the data file.
                 
                 :return: Vector of bond types
        """
    def get_high_x(self) -> float:
        """
                 Get the upper bound of the box in x direction.
                 
                 :return: Upper x boundary
        """
    def get_high_y(self) -> float:
        """
                 Get the upper bound of the box in y direction.
                 
                 :return: Upper y boundary
        """
    def get_high_z(self) -> float:
        """
                 Get the upper bound of the box in z direction.
                 
                 :return: Upper z boundary
        """
    def get_low_x(self) -> float:
        """
                 Get the lower bound of the box in x direction.
                 
                 :return: Lower x boundary
        """
    def get_low_y(self) -> float:
        """
                 Get the lower bound of the box in y direction.
                 
                 :return: Lower y boundary
        """
    def get_low_z(self) -> float:
        """
                 Get the lower bound of the box in z direction.
                 
                 :return: Lower z boundary
        """
    def get_lx(self) -> float:
        """
                 Get the box length in x direction.
                 
                 :return: Box length in x direction
        """
    def get_ly(self) -> float:
        """
                 Get the box length in y direction.
                 
                 :return: Box length in y direction
        """
    def get_lz(self) -> float:
        """
                 Get the box length in z direction.
                 
                 :return: Box length in z direction
        """
    def get_masses(self) -> dict[int, float]:
        """
                 Get the mass values for each atom type.
                 
                 :return: Map of atom types to their masses
        """
    def get_molecule_ids(self) -> list[int]:
        """
                 Get all molecule IDs from the data file.
                 
                 :return: Vector of molecule IDs
        """
    def get_nr_of_atom_types(self) -> int:
        """
                 Get the number of atom types in the data file.
                 
                 :return: Number of atom types
        """
    def get_nr_of_atoms(self) -> int:
        """
                 Get the number of atoms in the data file.
                 
                 :return: Number of atoms
        """
    def get_nr_of_bond_types(self) -> int:
        """
                 Get the number of bond types in the data file.
                 
                 :return: Number of bond types
        """
    def get_nr_of_bonds(self) -> int:
        """
                 Get the number of bonds in the data file.
                 
                 :return: Number of bonds
        """
    def read(self, path_of_file_to_read: str, atom_style: AtomStyle = ..., atom_style_2: AtomStyle = ..., atom_style_3: AtomStyle = ...) -> None:
        """
               Actually read a LAMMPS's `write_data` file.
        
               :param path_of_file_to_read: The path to the file to read
               :param atom_style: The format of the "Atoms" section, see https://docs.lammps.org/read_data.html
               :param atom_style2: The format of the "Atoms" section if the previous parameter is equal to AtomStyle::HYBRID
               :param atom_style_3: The format of the "Atoms" section if the second to last parameter is equal to AtomStyle::HYBRID
        """
class DataFileWriter:
    """
    
        A class to write a LAMMPS data file from a universe.
        
        .. attention::
            The resulting file is not guaranteed to be a completely valid LAMMPS file.
            In particular, this writer does not force you to set masses for all atom types,
            and it does not have limits on the image flags.
            You can use :meth:`~pylimer_tools_cpp.Universe.set_masses` to set the masses for all atom types,
            and either :meth:`~pylimer_tools_cpp.Universe.set_vertex_property` or 
            :meth:`~pylimer_tools_cpp.DataFileWriter.config_attempt_image_reset` and
            :meth:`~pylimer_tools_cpp.DataFileWriter.config_move_into_box` to ensure that the image flags are correct.
        
    """
    def __init__(self, universe: Universe) -> None:
        """
                Initialize the writer with the universe to write.
                
                :param universe: The universe to write to the data file
        """
    def config_atom_style(self, atom_style: AtomStyle = ...) -> None:
        """
                Set the (LAMMPS) atom style to use for writing the atoms.
        
                :param atom_style: The LAMMPS atom style to use (default: AtomStyle.ANGLE)
        """
    def config_attempt_image_reset(self, attempt_image_reset: bool = False) -> None:
        """
                Set whether to attempt to reset image flags so that output coordinates lie in the box.
        
                :param attempt_image_reset: Whether to attempt image flag reset (default: False)
        """
    def config_crosslinker_type(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2) -> None:
        """
                Set which atom type represents crosslinkers.
                Needed in case the moleculeIdx in the output file should have any meaning.
                (e.g. with :meth:`~pylimer_tools_cpp.DataFileWriter.config_molecule_idx_for_swap`).
        
                :param crosslinker_type: The atom type representing crosslinkers (default: 2)
        """
    def config_include_angles(self, include_angles: bool = True) -> None:
        """
                Set whether to include the angles from the universe in the file or not.
        
                :param include_angles: Whether to include angles (default: True)
        """
    def config_include_dihedral_angles(self, include_dihedral_angles: bool = True) -> None:
        """
                Set whether to include the dihedral angles from the universe in the file or not.
        
                :param include_dihedral_angles: Whether to include dihedral angles (default: True)
        """
    def config_include_velocities(self, include_velocities: bool = True) -> None:
        """
                Set whether to include the velocities from the universe (if any) in the file or not.
        
                :param include_velocities: Whether to include velocities (default: True)
        """
    def config_molecule_idx_for_swap(self, enable_swappability: bool = False) -> None:
        """
                Swappable chains implies that their `moleculeIdx` in the LAMMPS data file is not
                identical per chain, but identical per position in the chain.
                That's how you can have bond swapping with constant chain length distribution.
        
                :param enable_swappability: Whether to enable molecule index swappability (default: False)
        """
    def config_move_into_box(self, move_into_box: bool = False) -> None:
        """
                Set whether to change the output coordinates to lie in the box or not.
        
                :param move_into_box: Whether to move coordinates into box (default: False)
        """
    def config_reindex_atoms(self, reindex_atoms: bool = False) -> None:
        """
                Set whether to reindex the atoms or not.
                Re-indexing leads to atom IDs being in the range of 1 to the number of atoms.
        
                :param reindex_atoms: Whether to reindex atoms (default: False)
        """
    def set_custom_atom_format(self, atom_format: str = '\t$atomId\t$moleculeId\t$atomType\t$x\t$y\t$z\t$nx\t$ny\t$nz') -> None:
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
                :meth:`~pylimer_tools_cpp.Universe.set_property_value`
                as placeholders (as long as they are alphanumeric only; prefix in the format with '$' as well).
                Other placeholders are available if the universe was read from a LAMMPS data file with an
                atom style with additional data.
        
                This method is specifically useful if you need a different (or hybrid) atom style in LAMMPS.
        
                Be sure to still call :meth:`~pylimer_tools_cpp.DataFileWriter.config_atom_style`,
                so that the file can be read correctly again.
        
                :param atom_format: Custom format string for atoms (default: tab-separated standard format)
        """
    def set_universe_to_write(self, universe: Universe) -> None:
        """
                Re-set the universe to write.
                
                :param universe: The new universe to write to the data file
        """
    def write_to_file(self, file: str) -> None:
        """
                Actually do the writing to the disk.
        
                :param file: The path and file name to write to
        """
class DumpFileReader:
    """
    
           A reader for LAMMPS's `dump` files.
      
    """
    def __init__(self, path_of_file_to_read: str) -> None:
        """
                 Initialize the dump file reader.
                 
                 :param path_of_file_to_read: Path to the dump file to read
        """
    def __repr__(self) -> str:
        """
                 Return a string representation of the DumpFileReader.
                 
                 :return: String representation showing number of sections
        """
    def get_length(self) -> int:
        """
                 Get the number of sections (time-steps) in the file.
                 
                 :return: Number of timesteps/sections in the dump file
        """
    def get_numeric_values_for_at(self, rowIndex: typing.SupportsInt | typing.SupportsIndex, headerKey: str, columnIndex: str) -> list[float]:
        """
                 Get numeric values for the section at index, the main header headerKey and the column columnIndex.
                 
                 :param rowIndex: Index of the section/timestep to query
                 :param headerKey: Name of the header section to query  
                 :param columnIndex: Index of the column within the header
                 :return: Vector of numeric values for the specified location
        """
    def get_string_values_for_at(self, rowIndex: typing.SupportsInt | typing.SupportsIndex, headerKey: str, columnIndex: str) -> list[str]:
        """
                 Get string values for the section at index, the main header headerKey and the column columnIndex.
                 
                 :param rowIndex: Index of the section/timestep to query
                 :param headerKey: Name of the header section to query
                 :param columnIndex: Index of the column within the header
                 :return: Vector of string values for the specified location
        """
    def has_key(self, headerKey: str) -> bool:
        """
                 Check whether the first section has the header specified.
                 
                 :param headerKey: Name of the header to check for
                 :return: True if the header exists in the first section
        """
    def key_has_column(self, headerKey: str, columnName: str) -> bool:
        """
                 Check whether the header of the first section has the specified column.
                 
                 :param headerKey: Name of the header section to check
                 :param columnName: Name of the column to look for
                 :return: True if the column exists in the specified header
        """
    def key_has_directional_column(self, header_key: str, dir_prefix: str = '', dir_suffix: str = '') -> bool:
        """
                 Check whether the header of the first section has all the three columns {dir_prefix}{x|y|z}{dir_suffix}.
                 
                 :param header_key: Name of the header section to check
                 :param dir_prefix: Prefix for the directional columns (default: "")
                 :param dir_suffix: Suffix for the directional columns (default: "")
                 :return: True if all three directional columns (x, y, z) exist
        """
    def read(self) -> None:
        """
                 Read the whole file.
        """
class EndToEndStats:
    """
    
        Summary statistics for realized end-to-end distances.
      
    """
    def __init__(self) -> None:
        ...
    @property
    def mean_end_to_end_distance(self) -> float:
        ...
    @property
    def mean_squared_end_to_end_distance(self) -> float:
        ...
    @property
    def nr_of_excluded_zero_length(self) -> int:
        ...
    @property
    def nr_of_strands(self) -> int:
        ...
class EntanglementFilteringMode(enum.IntEnum):
    """
    Enum representing the filtering mode for entanglement sampling.
    """
    NONE: typing.ClassVar[EntanglementFilteringMode]  # value = <EntanglementFilteringMode.NONE: 0>
    OMIT_INACTIVE: typing.ClassVar[EntanglementFilteringMode]  # value = <EntanglementFilteringMode.OMIT_INACTIVE: 1>
    PRUNE_UPFRONT: typing.ClassVar[EntanglementFilteringMode]  # value = <EntanglementFilteringMode.PRUNE_UPFRONT: 2>
    @classmethod
    def __new__(cls, value):
        ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """
class EntanglementMode(enum.Enum):
    """
    Representation of topological entanglements in an MEHP network.
    """
    MERGED_LINKS: typing.ClassVar[EntanglementMode]  # value = <EntanglementMode.MERGED_LINKS: 1>
    NONE: typing.ClassVar[EntanglementMode]  # value = <EntanglementMode.NONE: 0>
    SPRING: typing.ClassVar[EntanglementMode]  # value = <EntanglementMode.SPRING: 2>
class ExitReason(enum.IntEnum):
    """
    
    An enum representing the reason for exiting
    the simulation or optimization procedure.
    """
    FAILURE: typing.ClassVar[ExitReason]  # value = <ExitReason.FAILURE: 5>
    F_TOLERANCE: typing.ClassVar[ExitReason]  # value = <ExitReason.F_TOLERANCE: 1>
    INTERRUPT: typing.ClassVar[ExitReason]  # value = <ExitReason.INTERRUPT: 6>
    MAX_STEPS: typing.ClassVar[ExitReason]  # value = <ExitReason.MAX_STEPS: 3>
    NO_STEPS_POSSIBLE: typing.ClassVar[ExitReason]  # value = <ExitReason.NO_STEPS_POSSIBLE: 4>
    OTHER: typing.ClassVar[ExitReason]  # value = <ExitReason.OTHER: 7>
    UNSET: typing.ClassVar[ExitReason]  # value = <ExitReason.UNSET: 0>
    X_TOLERANCE: typing.ClassVar[ExitReason]  # value = <ExitReason.X_TOLERANCE: 2>
    @classmethod
    def __new__(cls, value):
        ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """
class HookeanForceEvaluator(MEHPForceEvaluator):
    """
    
          Linear Gaussian (Hookean) spring constitutive evaluator:
          :math:`E_s = \\frac{\\kappa}{2 N_s} \\|\\Delta\\mathbf{r}_s\\|^2`.
    
          :param kappa: Spring stiffness prefactor (>= 0).
        
    """
    def __init__(self, kappa: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> None:
        ...
class HybridConfig:
    """
    Configuration for Hybrid exact low-modes + deflated SLQ.
    """
    backend: ComputeBackend
    def __init__(self, num_steps: typing.SupportsInt | typing.SupportsIndex = 200, num_random_vectors: typing.SupportsInt | typing.SupportsIndex = 20, num_exact_low_modes: typing.SupportsInt | typing.SupportsIndex = 100, backend: ComputeBackend = ...) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def num_exact_low_modes(self) -> int:
        ...
    @num_exact_low_modes.setter
    def num_exact_low_modes(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_random_vectors(self) -> int:
        ...
    @num_random_vectors.setter
    def num_random_vectors(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_steps(self) -> int:
        ...
    @num_steps.setter
    def num_steps(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
class KPMConfig:
    """
    Configuration for Kernel Polynomial Method (KPM) spectral approximation.
    """
    backend: ComputeBackend
    def __init__(self, num_moments: typing.SupportsInt | typing.SupportsIndex = 1024, num_random_vectors: typing.SupportsInt | typing.SupportsIndex = 20, num_quadrature_points: typing.SupportsInt | typing.SupportsIndex = 2048, batch_size: typing.SupportsInt | typing.SupportsIndex = 4, backend: ComputeBackend = ...) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def batch_size(self) -> int:
        ...
    @batch_size.setter
    def batch_size(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_moments(self) -> int:
        ...
    @num_moments.setter
    def num_moments(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_quadrature_points(self) -> int:
        ...
    @num_quadrature_points.setter
    def num_quadrature_points(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_random_vectors(self) -> int:
        ...
    @num_random_vectors.setter
    def num_random_vectors(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
class LazyUniverseSequenceIterator:
    """
    
           An iterator to iterate throught the universes in :obj:`~pylimer_tools_cpp.UniverseSequence`.
      
    """
    def __iter__(self) -> LazyUniverseSequenceIterator:
        ...
    def __next__(self) -> Universe:
        ...
class LinearMaxDistanceProvider(MaxDistanceProvider):
    """
    
        For MC generation, converts the :math:`N` to a maximum distance within which to sample.
        The distance will be calculated as :math:`N \\times \\text{max_distance_multiplier}`.
        Useful only for performance improvements in large systems.
    
        :param max_distance_multiplier: Multiplier for the maximum distance.
        
    """
    def __init__(self, max_distance_multiplier: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def __repr__(self) -> str:
        """
                 Return a string representation of the LinearMaxDistanceProvider.
                 
                 :return: String representation
        """
    def get_max_distance(self, N: typing.SupportsFloat | typing.SupportsIndex) -> float:
        """
                 Get the maximum distance for a given N.
        
                 :param N: Number of segments.
        """
class LinearSystemMEHPSolver(MEHPSolver):
    """
    Equilibrium solver for linear (Hookean) MEHP networks via sparse linear system assembly.
    """
    config: MEHPLinearSolverConfig
    def __init__(self, config: MEHPLinearSolverConfig = ...) -> None:
        ...
    def get_config(self) -> MEHPLinearSolverConfig:
        ...
    def set_config(self, config: MEHPLinearSolverConfig) -> None:
        ...
class LinkSwappingMode(enum.Enum):
    """
    How slip-links may act when they reach each-other or even a crosslink.
    """
    ALL: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL: 2>
    ALL_CYCLE: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL_CYCLE: 3>
    ALL_MC: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL_MC: 4>
    ALL_MC_CYCLE: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL_MC_CYCLE: 5>
    ALL_MC_TRY: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL_MC_TRY: 6>
    ALL_MC_TRY_CYCLE: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL_MC_TRY_CYCLE: 7>
    NO_SWAPPING: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.NO_SWAPPING: 0>
    SLIPLINKS_ONLY: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.SLIPLINKS_ONLY: 1>
class MCUniverseGenerator:
    """
    
           A :obj:`pylimer_tools_cpp.Universe` generator using a Monte-Carlo procedure.
           This generator creates phantom bead-spring polymer networks with crosslinkers and Gaussian chain and bond statistics.
    
           .. note::
    
              A note on the characteristic ratio :math:`C_{\\infty}` and the bead distance :math:`b`:
              If you use parameters from :mod:`pylimer_tools.io.bead_spring_parameter_provider` with 
              type :obj:`pylimer_tools.io.bead_spring_parameter_provider.ParameterType.GAUSSIAN`,
              there is no need to set :math:`C_{\\infty}` in the MCUniverseGenerator methods,
              as the bead distance :math:`b` you set via :meth:`~pylimer_tools_cpp.MCUniverseGenerator.set_bead_distance`
              from :meth:`pylimer_tools.io.bead_spring_parameter_provider.Parameters.get("<b>")`
              already incorporates the characteristic ratio.
    
              This is the recommended way to use the MCUniverseGenerator, since using e.g.
              Kuhn segment parameters conflicts with the generator producing Gaussian bond statistics.
    
           Please cite :cite:t:`gusev_molecular_2024` and/or :cite:t:`bernhard_phantom_2025` if you use this method in your work.
      
    """
    def __copy__(self) -> MCUniverseGenerator:
        ...
    def __getstate__(self) -> tuple[str]:
        ...
    def __init__(self, lx: typing.SupportsFloat | typing.SupportsIndex = 10.0, ly: typing.SupportsFloat | typing.SupportsIndex = 10.0, lz: typing.SupportsFloat | typing.SupportsIndex = 10.0) -> None:
        """
                 Initialize a new MCUniverseGenerator with specified box dimensions.
        
                 :param lx: Box length in x-direction (default: 10.0).
                 :param ly: Box length in y-direction (default: 10.0).
                 :param lz: Box length in z-direction (default: 10.0).
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def add_crosslinkers(self, nr_of_crosslinkers: typing.SupportsInt | typing.SupportsIndex, crosslinker_functionality: typing.SupportsInt | typing.SupportsIndex = 4, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, white_noise: bool = True) -> None:
        """
                    Add crosslinkers at random positions.
        
                    :param nr_of_crosslinkers: Number of crosslinkers to add.
                    :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
                    :param crosslinker_type: Atom type of the crosslinkers (default: 2).
                    :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).
        """
    def add_crosslinkers_at(self, coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], crosslinker_functionality: typing.SupportsInt | typing.SupportsIndex = 4, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2) -> None:
        """
                  Add crosslinkers at specific coordinates.
        
                  :param coordinates: Coordinates of the crosslinkers as a flat array [x1, y1, z1, x2, y2, z2, ...].
                  :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
                  :param crosslinker_type: Atom type of the crosslinkers (default: 2).
        """
    def add_end_functionalized_strands(self, nr_of_strands: typing.SupportsInt | typing.SupportsIndex, strand_length: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], crosslinker_functionality: typing.SupportsInt | typing.SupportsIndex = 4, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, strand_atom_type: typing.SupportsInt | typing.SupportsIndex = 1, white_noise: bool = True) -> None:
        """
                    Add strands with crosslinkers at the ends.
        
                    :param nr_of_strands: Number of strands to add.
                    :param strand_length: Length of each strand.
                    :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
                    :param crosslinker_type: Atom type of the crosslinkers (default: 2).
                    :param strand_atom_type: Atom type of the beads that are not at the ends (default: 1).
                    :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).
        """
    def add_functionalized_strands(self, nr_of_strands: typing.SupportsInt | typing.SupportsIndex, strand_length: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], crosslink_selector: collections.abc.Callable[[int, int, int], tuple[bool, typing.SupportsInt | typing.SupportsIndex]], default_crosslinker_functionality: typing.SupportsInt | typing.SupportsIndex, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, strand_atom_type: typing.SupportsInt | typing.SupportsIndex = 1, white_noise: bool = True) -> None:
        """
                  Add strands with custom crosslink placement using a selector function.
        
                  This is a general implementation that allows for flexible crosslink placement patterns.
                  The crosslink_selector function is called for each bead position and should return
                  a tuple (should_convert, functionality) indicating whether that bead should be
                  converted to a crosslink and what functionality it should have.
        
                  .. tip::
        
                     This method provides the most flexibility for creating custom crosslink patterns.
                     For common patterns, consider using the more specific methods like
                     :meth:`~pylimer_tools_cpp.MCUniverseGenerator.add_randomly_functionalized_strands` or
                     :meth:`~pylimer_tools_cpp.MCUniverseGenerator.add_regularly_spaced_functionalized_strands`.
        
                  Example usage::
        
                      # Create a custom pattern with crosslinks every 5 beads starting at position 2
                      def my_selector(strand_index, bead_index, total_beads):
                          if bead_index >= 2 and (bead_index - 2) % 5 == 0:
                              return (True, 4)  # Convert to crosslink with functionality 4
                          return (False, 0)     # Don't convert
        
                      generator.add_functionalized_strands(
                          nr_of_strands=10,
                          strand_length=[20] * 10,
                          crosslink_selector=my_selector,
                          default_crosslinker_functionality=4
                      )
        
                  :param nr_of_strands: Number of strands to add.
                  :param strand_length: Length of each strand (list of integers).
                  :param crosslink_selector: Function that takes (strand_index, bead_index, total_beads) and returns (should_convert, functionality).
                  :param default_crosslinker_functionality: Default functionality for crosslinks (used for positioning).
                  :param crosslinker_type: Atom type of the crosslinkers (default: 2).
                  :param strand_atom_type: Atom type of the beads that stay (default: 1).
                  :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).     
        """
    def add_monofunctional_strands(self, nr_of_monofunctional_strands: typing.SupportsInt | typing.SupportsIndex, monofunctional_strand_length: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], monofunctional_strand_atom_type: typing.SupportsInt | typing.SupportsIndex = 4) -> None:
        """
                 Add multiple monofunctional strands with specified bead types.
        
                 :param nr_of_monofunctional_strands: Number of monofunctional strands to add.
                 :param monofunctional_strand_length: Vector specifying the length of each strand in beads.
                 :param monofunctional_strand_atom_type: Atom type for the strand beads (default: 4).
        """
    def add_randomly_functionalized_strands(self, nr_of_strands: typing.SupportsInt | typing.SupportsIndex, strand_length: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], functionalization_probability: typing.SupportsFloat | typing.SupportsIndex, crosslinker_functionality: typing.SupportsInt | typing.SupportsIndex = 4, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, strand_atom_type: typing.SupportsInt | typing.SupportsIndex = 1, white_noise: bool = True) -> None:
        """
                  Add strands with randomly distributed crosslinkers in between.
        
                  :param nr_of_strands: Number of strands to add.
                  :param strand_length: Length of each strand.
                  :param functionalization_probability: Proportion of beads that are made crosslink.
                  :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
                  :param crosslinker_type: Atom type of the crosslinkers (default: 2).
                  :param strand_atom_type: Atom type of the beads that stay (default: 1).
                  :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).     
        """
    def add_regularly_spaced_functionalized_strands(self, nr_of_strands: typing.SupportsInt | typing.SupportsIndex, strand_length: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], spacing_between_crosslinks: typing.SupportsInt | typing.SupportsIndex, offset_to_first_crosslink: typing.SupportsInt | typing.SupportsIndex = 0, crosslinker_functionality: typing.SupportsInt | typing.SupportsIndex = 4, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, strand_atom_type: typing.SupportsInt | typing.SupportsIndex = 1, white_noise: bool = True) -> None:
        """
                  Add strands with regularly spaced crosslinkers in between.
        
                  :param nr_of_strands: Number of strands to add.
                  :param strand_length: Length of each strand.
                  :param spacing_between_crosslinks: Number of beads between crosslinks.
                  :param offset_to_first_crosslink: Offset from start of strand to first crosslink (0-based, default: 0).
                  :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
                  :param crosslinker_type: Atom type of the crosslinkers (default: 2).
                  :param strand_atom_type: Atom type of the beads that stay (default: 1).
                  :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).     
        """
    def add_solvent_chains(self, nr_of_solvent_chains: typing.SupportsInt | typing.SupportsIndex, solvent_chain_length: typing.SupportsInt | typing.SupportsIndex, solvent_atom_type: typing.SupportsInt | typing.SupportsIndex = 3, white_noise: bool = True) -> None:
        """
                    Randomly distribute additional, free chains.
        
                    :param nr_of_solvent_chains: Number of solvent chains to add.
                    :param solvent_chain_length: Length of each solvent chain in beads.
                    :param solvent_atom_type: Atom type for solvent chain beads (default: 3).
                    :param white_noise: Whether to use white noise (true) or blue noise (false) for positioning (default: true).
        """
    def add_star_crosslinkers(self, nr_of_stars: typing.SupportsInt | typing.SupportsIndex, functionality: typing.SupportsInt | typing.SupportsIndex, beads_per_strand: typing.SupportsInt | typing.SupportsIndex, crosslinker_atom_type: typing.SupportsInt | typing.SupportsIndex = 2, strand_atom_type: typing.SupportsInt | typing.SupportsIndex = 1, white_noise: bool = True) -> None:
        """
                 Add star-like crosslinkers with pre-connected strands (useful e.g. for tetra-PEG networks).
                 Each star crosslinker will have the specified functionality with strands already attached.
        
                 .. tip::
        
                    To have a certain polydispersity in the arms of __one__ star crosslinker,
                    the stars can alternatively be created using :meth:`~pylimer_tools_cpp.MCUniverseGenerator.add_crosslinkers`,
                    :meth:`~pylimer_tools_cpp.MCUniverseGenerator.add_strands` and :meth:`~pylimer_tools_cpp.MCUniverseGenerator.link_strand_to`.
        
                 :param nr_of_stars: Number of star crosslinkers to add
                 :param functionality: Functionality of each star crosslinker (number of strands)
                 :param beads_per_strand: Number of beads in each strand
                 :param crosslinker_atom_type: Atom type for the crosslinker
                 :param strand_atom_type: Atom type for the strand beads
                 :param white_noise: Whether to use white noise positioning
        """
    def add_strands(self, nr_of_strands: typing.SupportsInt | typing.SupportsIndex, strand_lengths: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], strand_atom_type: typing.SupportsInt | typing.SupportsIndex = 1) -> None:
        """
                    Add strands.
                    Adds them unconnected at first.
                    To link them to crosslinkers, use some of the :code:`link_strand_` methods.
        
                    :param nr_of_strands: Number of strands to add.
                    :param strand_lengths: A list of integers representing the number of beads of each of the strands.
                    :param strand_atom_type: Type of atoms for the strands (default: 1).
        """
    def compute_bond_length_stats(self) -> BondLengthStats:
        """
                    Stream over strands in O(1) strand memory to compute realized bond length statistics
                    without constructing a Universe.
        
                    :return: BondLengthStats object.
        """
    def compute_end_to_end_stats(self, category: StrandCategory = ...) -> EndToEndStats:
        """
                    Stream over strands matching `category` to compute end-to-end distance statistics
                    without constructing a Universe.
        
                    :param category: StrandCategory filter (default: NETWORK_STRAND).
                    :return: EndToEndStats object.
        """
    def compute_strand_bond_counts(self, category: StrandCategory = ...) -> list[int]:
        """
                    Return the number of bonds for each strand matching `category`.
        
                    :param category: StrandCategory filter (default: ALL).
                    :return: List of bond counts per matching strand.
        """
    def config_nr_of_mc_steps(self, n_steps: typing.SupportsInt | typing.SupportsIndex = 0) -> None:
        """
               Set the number of Monte-Carlo steps during bond length equilibration.
        
               :param n_steps: Number of MC steps to perform (default: 0)
        """
    def config_primary_loop_probability(self, probability: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> None:
        """
                 Configure an additional weight reduction for the primary loop formation.
        
                 Defaults to 1., which means the general :math:`P(\\vec{R})` is used without any bias.
                 This results in more primary loops for shorter chains than for longer ones.
        
                 Set to 0. to disable the formation of primary loops.
        """
    def config_secondary_loop_probability(self, probability: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> None:
        """
                 Configure an additional weight reduction for the secondary loop formation.
        
                 Defaults to 1., which means the general :math:`P(\\vec{R})` is used without any bias.
                 This results in more secondary loops for shorter chains than for longer ones.
        
                 Set to 0. to disable the formation of secondary loops.
        """
    def count_atom_types(self) -> dict[int, int]:
        """
                    Count atoms per atom type without constructing a Universe.
        
                    :return: Dict mapping atom type to count.
        """
    def count_atoms_by_degree(self) -> dict[int, int]:
        """
                    Count atoms per bond degree without constructing a Universe.
        
                    :return: Dict mapping degree to atom count.
        """
    def determine_effective_functionality_per_type(self) -> dict[int, float]:
        """
                    Compute effective functionality (mean degree) per atom type without constructing a Universe.
        
                    :return: Dict mapping atom type to mean degree.
        """
    def disable_max_distance(self) -> None:
        """
                 Disable the maximum distance provider to allow unlimited distance sampling.
                 This may slow down performance for large systems but ensures complete sampling.
        """
    def get_box(self) -> Box:
        """
                    Get the simulation box of the generator.
        
                    :return: The Box object.
        """
    def get_crosslinker_universe(self) -> CrosslinkerUniverse:
        """
                    Get a copy of the internal CrosslinkerUniverse (crosslinker coordinates and strand topology)
                    without instantiating individual strand beads or constructing a Universe.
        
                    :return: The CrosslinkerUniverse struct.
        """
    def get_current_crosslinker_conversion(self) -> float:
        """
                 Get the current conversion of crosslinkers, i.e., the fraction of crosslinkers
                 that have been linked to strands.
        
                 :return: Crosslinker conversion as a fraction between 0 and 1.
        """
    def get_current_nr_of_atoms(self) -> int:
        """
                 Get the current number of atoms that the universe would/will have.
        
                 :return: Number of atoms in the generated universe.
        """
    def get_current_nr_of_bonds(self) -> int:
        """
                 Get the current number of bonds that the universe would/will have.
        
                 :return: Number of bonds in the generated universe.
        """
    def get_current_strand_conversion(self) -> float:
        """
                 Get the current conversion of strands, i.e., the fraction of strand ends that have been consumed.
        
                 :return: Strand conversion as a fraction between 0 and 1.
        """
    def get_mean_bead_distance(self) -> float:
        """
                  Get the currently configured mean bead distance.
        
                  :return: The currently configured mean distance between beads
        """
    def get_mean_squared_bead_distance(self) -> float:
        """
                  Get the currently configured mean squared bead distance.
        
                  :return: The currently configured mean squared distance between beads
        """
    def get_mehp_network(self, entanglements: StrandPairEntanglements = ..., entanglement_mode: EntanglementMode = ..., crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, remove_dangling_chains: bool = False, include_free_chains: bool = True) -> MEHPNetwork:
        """
                 Build a canonical MEHPNetwork from the current simplified universe and optional sampled entanglements.
        
                 :param entanglements: Optional pre-sampled strand-pair entanglements.
                 :param entanglement_mode: Entanglement representation mode (default: EntanglementMode.MERGED_LINKS).
                 :param crosslinker_type: Atom type identifying crosslinkers (default: 2).
                 :param remove_dangling_chains: Omit dangling chains when True (default: False).
                 :param include_free_chains: Include free/soluble chains when True (default: True).
                 :return: Configured MEHPNetwork instance.
        """
    @typing.overload
    def get_mehp_simulator(self, entanglements: StrandPairEntanglements = ..., entanglement_mode: EntanglementMode = ..., force_evaluator: MEHPForceEvaluator = None, is_2d: bool = False, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, remove_dangling_chains: bool = False, include_free_chains: bool = True) -> MEHPSimulator:
        """
              Construct a unified MEHPSimulator from the current simplified universe and optional pre-sampled entanglements,
              bypassing the creation of a Universe.
        
              :param entanglements: Optional pre-sampled strand-pair entanglements.
              :param entanglement_mode: Entanglement representation mode (default: EntanglementMode.MERGED_LINKS).
              :param force_evaluator: Optional constitutive force evaluator (defaults to HookeanForceEvaluator(1.0) when None).
              :param is_2d: Treat system as 2D when True (default: False).
              :param crosslinker_type: Atom type identifying crosslinkers (default: 2).
              :param remove_dangling_chains: Omit dangling chains when True (default: False).
              :param include_free_chains: Include free/soluble chains when True (default: True).
              :return: Configured MEHPSimulator instance.
        """
    @typing.overload
    def get_mehp_simulator(self, nr_of_entanglements_to_sample: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_entanglements: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = '', filter_dangling_and_soluble: bool = True, entanglement_mode: EntanglementMode = ..., force_evaluator: MEHPForceEvaluator = None, is_2d: bool = False, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2) -> MEHPSimulator:
        """
              Sample random entanglements from the simplified universe and construct a unified MEHPSimulator,
              bypassing the creation of a Universe.
        
              :param nr_of_entanglements_to_sample: Target number of entanglements to sample.
              :param upper_cutoff: Upper pair sampling cutoff.
              :param lower_cutoff: Lower pair sampling cutoff (default: 0.0).
              :param minimum_nr_of_entanglements: Minimum entanglements required if possible (default: 0).
              :param same_strand_cutoff: Minimum chemical distance for same-strand pairs (default: 3.0).
              :param seed: PRNG seed string (empty uses generator seed).
              :param filter_dangling_and_soluble: Omit entanglements on inactive strands (default: True).
              :param entanglement_mode: Entanglement representation mode (default: EntanglementMode.MERGED_LINKS).
              :param force_evaluator: Optional constitutive force evaluator (defaults to HookeanForceEvaluator(1.0) when None).
              :param is_2d: Treat system as 2D when True (default: False).
              :param crosslinker_type: Atom type identifying crosslinkers (default: 2).
              :return: Configured MEHPSimulator instance.
        """
    def get_nr_of_atoms(self) -> int:
        """
                  Get the current number of atoms that the universe would/will have.
        
                  :return: Number of atoms in the generated universe
        """
    def get_nr_of_bonds(self) -> int:
        """
                  Get the current number of bonds that the universe would/will have.
        
                  :return: Number of bonds in the generated universe
        """
    def get_universe(self) -> Universe:
        """
                    Fetch the current (or final) state of the universe.
        
                    Use this method to actually (MC) place beads between the crosslinks and retrieve the generated structure.
        
                    :return: The generated Universe object containing all atoms, bonds, and their coordinates.
        """
    def link_strand(self, strand_idx: typing.SupportsInt | typing.SupportsIndex, c_infinity: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> None:
        """
                  Link one particular strand to a crosslinker.
                  This strand will have one bond made, to an appropriate crosslinker, 
                  as chosen by the parameters associated with the strand.
        
                  :param strand_idx: Index of the strand to be linked.
                  :param c_infinity: As needed for the end-to-end distribution, given by :math:`\\langle R^2\\rangle_0 = C_{\\infty} N b^2` (default: 1.0).
        """
    def link_strand_to(self, strand_idx: typing.SupportsInt | typing.SupportsIndex, link_idx: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
                  Link a strand to a specific crosslinker.
                  This assumes that you keep track of the order in which you added the crosslinkers and strands,
                  as those will determine the indices.
        
                  .. caution::
        
                       Be aware that some few methods may change the cross-link or strand indices.
        
                  :param strand_idx: Index of the strand to be linked.
                  :param link_idx: Index of the crosslinker to be linked.
        """
    def link_strand_to_strand(self, c_infinity: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> bool:
        """
                 Link two free strand ends together directly.
                 This method finds free strand ends and combines them into a single strand based on
                 end-to-end distance probability distributions.
        
                 .. caution::
        
                       - The strand indices will change after this method is called.
                       - The probability of linking to the end of a free strand is the same as 
                         linking to an end with a distance 0.
        
                 :param c_infinity: Statistical parameter for end-to-end distance calculation
                 :return: True if a successful link was made, False if no suitable pair was found
        """
    def link_strands_callback(self, linking_controller: collections.abc.Callable[[MCUniverseGenerator, int], BackTrackStatus], c_infinity: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> None:
        """
                    Link strands to crosslinkers using a custom callback function to control when to stop.
        
                    The callback function receives the current MCUniverseGenerator state and the current step number,
                    and should return a BackTrackStatus value:
                    - STOP: Stop the linking process
                    - TRACK_FORWARD: Continue linking, make more bonds
                    - TRACK_BACKWARD: Track backward in the linking process, i.e., remove bonds again
        
                    :param linking_controller: Callback function that controls the linking process. 
                                              Function signature: (MCUniverseGenerator, int) -> BackTrackStatus
                    :param c_infinity: As needed for the end-to-end distribution, given by :math:`\\langle R^2\\rangle_0 = C_{\\infty} N b^2`.
        """
    def link_strands_to_conversion(self, crosslinker_conversion: typing.SupportsFloat | typing.SupportsIndex, c_infinity: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> None:
        """
                    Actually link the previously added strands to the previously added crosslinkers,
                    until a certain crosslink conversion is reached.
        
                    :param crosslinker_conversion: Target conversion of crosslinkers (0: no connections to crosslinks; 1: all crosslinkers fully connected).
                    :param c_infinity: As needed for the end-to-end distribution, given by :math:`\\langle R^2\\rangle_0 = C_{\\infty} N b^2`.
        """
    def link_strands_to_soluble_fraction(self, soluble_fraction: typing.SupportsFloat | typing.SupportsIndex, c_infinity: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> None:
        """
                    Actually link the previously added strands to the previously added crosslinkers,
                    until a certain soluble fraction is reached.
        
                    :param soluble_fraction: Target soluble fraction (0: all material is soluble; 1: no material is soluble).
                    :param c_infinity: As needed for the end-to-end distribution, given by :math:`\\langle R^2\\rangle_0 = C_{\\infty} N b^2` (default: 1.0).
        """
    def link_strands_to_strands_to_conversion(self, target_strand_conversion: typing.SupportsFloat | typing.SupportsIndex, c_infinity: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> None:
        """
                 Link free strand ends to each other until target conversion is reached.
                 This method repeatedly calls :meth:`~pylimer_tools_cpp.MCUniverseGenerator.link_strand_to_strand` 
                 in a slightly more efficient manner than you could from Python,
                 until the specified conversion of free strand ends is achieved.
        
                 .. warning::
        
                    The strands are merged during this process, i.e., from two strands results one strand.
                    Consequently, calling :meth:`~pylimer_tools_cpp.MCUniverseGenerator.get_current_strands_conversion` after this method will not return
                    the requested `target_strand_conversion`, since that one is taken relative to the number of strands when calling this method.
        
                 .. caution::
        
                    Beware of the notes on :meth:`~pylimer_tools_cpp.MCUniverseGenerator.link_strand_to_strand`.
        
                 :param target_strand_conversion: Target conversion of free strand ends (0.0 to 1.0)
                 :param c_infinity: Statistical parameter for end-to-end distance calculation
        """
    def relax_crosslinks(self) -> None:
        """
                 Run force relaxation with the crosslinkers and their strands,
                 to have the crosslinks in their statistically most probable position.
        """
    def remove_soluble_fraction(self, rescale_box: bool = True) -> None:
        """
                    Remove soluble fraction (as determined by phantom force relaxation) of the strands and crosslinks.
        
                    :param rescale_box: Whether to rescale the box dimensions to maintain constant density (default: true).
        """
    def sample_entanglements(self, nr_of_entanglements_to_sample: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_entanglements: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = '', filtering_mode: EntanglementFilteringMode = ...) -> StrandPairEntanglements:
        """
                 Sample entanglements directly from the generator's CrosslinkerUniverse
                 using per-strand bounding boxes and on-demand Brownian bridge bead
                 resampling without constructing a full Universe.
        
                 :param nr_of_entanglements_to_sample: Target number of entanglements to sample.
                 :param upper_cutoff: Upper cutoff distance between entangled beads.
                 :param lower_cutoff: Lower cutoff distance between entangled beads (default: 0.0).
                 :param minimum_nr_of_entanglements: Minimum number of entanglements required (default: 0).
                 :param same_strand_cutoff: Minimum contour separation on the same strand (default: 3.0).
                 :param seed: Optional random seed (uses generator RNG if empty).
                 :param filtering_mode: Filtering mode for entanglements (default: EntanglementFilteringMode.OMIT_INACTIVE).
                 :return: Sampled StrandPairEntanglements.
        """
    def set_bead_distance(self, distance: typing.SupportsFloat | typing.SupportsIndex, update_mean_squared: bool = True) -> None:
        """
                 Set the mean distance between beads when doing MC stepping.
                 Also used for the target crosslinker partner sampling.
        
                 NOTE: Mainly the mean squared bead distance is effectively used in the Monte-Carlo simulation.
        
                 :param distance: Mean distance between beads.
                 :param update_mean_squared: Whether to update the mean squared distance as well, deduced from the assumed gaussian distribution in 3D (default: true).
        """
    def set_mean_squared_bead_distance(self, mean_squared_distance: typing.SupportsFloat | typing.SupportsIndex, update_mean: bool = True) -> None:
        """
                  Set the mean squared distance between beads.
        
                  :param mean_squared_distance: Mean squared distance between beads
                  :param update_mean: Whether to update the mean bead distance as well, deduced from the assumed gaussian distribution in 3D (default: true)
        """
    def set_seed(self, seed: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
                 Set the seed for the random generator.
        
                 :param seed: Random seed value for reproducible results.
        """
    def use_linear_max_distance(self, multiplier: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
                  Converts the :math:`N` to a maximum distance within which to sample.
                  The distance will be calculated as :math:`N \\times \\text{multiplier}`.
                  For example, commonly, a maximum distance is :math:`N \\times <b>`, 
                  in which case the multiplier is :math:`<b>`.
        
                  Useful only for performance improvements in large systems.
        """
    def use_zscore_max_distance(self, std_multiplier: typing.SupportsFloat | typing.SupportsIndex, in_sqrt_multiplier: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> None:
        """
                 Converts the :math:`N` to a maximum distance within which to sample.
                 The distance will be calculated as :math:`\\text{std_multiplier} \\times \\sqrt{N \\times \\text{in_sqrt_multiplier}}`.
                 Useful only for performance improvements in large systems.
        
                 :param std_multiplier: The Z-Score, the multiplier of the standard deviation of the end-to-end distribution.
                       E.g., 3. for 99.9994% of all conformations.
                 :param in_sqrt_multiplier: The multiplier with the :math:`N` in the square root. Probably :math:`<b^2>`.
        """
    def validate(self) -> None:
        """
                 Check whether the internal state of the generator is valid.
                 Throws errors if not. 
                 Should in principle always be valid when called from Python; if not, there is a bug in the code.
        """
class MEHPForceEvaluator:
    """
    
          Stateless constitutive model interface for evaluating spring energies, forces,
          stiffnesses, and virial stress tensors in MEHP networks.
    
          To create a custom force evaluator in Python, subclass this class and implement:
    
          - :meth:`evaluate_spring_energy`: Compute strain energy for one spring
          - :meth:`evaluate_spring_force`: Compute force vector for one spring
        
    """
    def __init__(self) -> None:
        ...
    def evaluate_network_energy_and_gradient(self, net: MEHPNetwork, spring_vectors: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], is_2d: bool = False, compute_gradient: bool = True) -> tuple:
        """
                Evaluate total network elastic energy and optional link gradient.
        
                :param net: Network topology and contour lengths.
                :param spring_vectors: Flattened spring vectors (size 3 * nr_of_springs).
                :param is_2d: If True, treat the system as planar (ignore z components).
                :param compute_gradient: If True, compute and return link gradient ndarray.
                :return: Tuple (energy, gradient_or_none).
        """
    def evaluate_spring_energy(self, spring_vec: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], contour_length: typing.SupportsFloat | typing.SupportsIndex) -> float:
        """
        Evaluate strain energy of a single spring.
        """
    def evaluate_spring_force(self, spring_vec: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], contour_length: typing.SupportsFloat | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
        Evaluate spring force vector dE/d(spring_vec) for a single spring.
        """
    def evaluate_spring_stiffness(self, spring_vec: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], contour_length: typing.SupportsFloat | typing.SupportsIndex) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 3]"]:
        """
        Evaluate 3x3 tangent stiffness matrix dF/d(spring_vec) for a single spring.
        """
    def evaluate_stress_contribution(self, spring_vec: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], i: typing.SupportsInt | typing.SupportsIndex, j: typing.SupportsInt | typing.SupportsIndex, contour_length: typing.SupportsFloat | typing.SupportsIndex) -> float:
        """
        Evaluate the (i, j) virial stress contribution of a single spring.
        """
    def evaluate_stress_tensor(self, net: MEHPNetwork, spring_vectors: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], volume: typing.SupportsFloat | typing.SupportsIndex, is_2d: bool = False) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 3]"]:
        """
        Evaluate the 3x3 virial stress tensor for the network.
        """
    def get_kappa(self) -> float:
        """
        Return the spring stiffness parameter kappa.
        """
    def is_linear(self) -> bool:
        """
        Indicate whether the constitutive force law is linear in displacement.
        """
class MEHPLinearSolverConfig:
    """
    Configuration for the sparse linear-system MEHP equilibrium solver.
    """
    simplification_mode: StructureSimplificationMode
    slip_mode: SlipMode
    solver: SLESolver
    def __init__(self, solver: SLESolver = ..., simplification_mode: StructureSimplificationMode = ..., inactive_removal_cutoff: typing.SupportsFloat | typing.SupportsIndex = 1e-06, residual_reduction: typing.SupportsFloat | typing.SupportsIndex = 1e-06, max_iterations: typing.SupportsInt | typing.SupportsIndex = 50000, spring_breaking_length: typing.SupportsFloat | typing.SupportsIndex = -1.0, in_place_rebuild_threshold: typing.SupportsFloat | typing.SupportsIndex = 0.15, cg_fallback_max_iterations: typing.SupportsInt | typing.SupportsIndex = 20000, gradient_descent_learning_rate: typing.SupportsFloat | typing.SupportsIndex = 0.01, slip_mode: SlipMode = ...) -> None:
        ...
    @property
    def cg_fallback_max_iterations(self) -> int:
        ...
    @cg_fallback_max_iterations.setter
    def cg_fallback_max_iterations(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def gradient_descent_learning_rate(self) -> float:
        ...
    @gradient_descent_learning_rate.setter
    def gradient_descent_learning_rate(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def in_place_rebuild_threshold(self) -> float:
        ...
    @in_place_rebuild_threshold.setter
    def in_place_rebuild_threshold(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def inactive_removal_cutoff(self) -> float:
        ...
    @inactive_removal_cutoff.setter
    def inactive_removal_cutoff(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def max_iterations(self) -> int:
        ...
    @max_iterations.setter
    def max_iterations(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def residual_reduction(self) -> float:
        ...
    @residual_reduction.setter
    def residual_reduction(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def spring_breaking_length(self) -> float:
        ...
    @spring_breaking_length.setter
    def spring_breaking_length(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class MEHPNetwork:
    """
    
          Canonical polymer network representation shared across MEHP linear-system,
          slip-link relaxation, and nonlinear optimization solvers.
        
    """
    assume_box_large_enough: bool
    assume_complete: bool
    remove_primary_loops: bool
    @staticmethod
    def from_crosslinker_universe(box: Box, crosslinker_universe: CrosslinkerUniverse, entanglements: StrandPairEntanglements = ..., crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, entanglement_mode: EntanglementMode = ..., remove_dangling_chains: bool = False, include_free_chains: bool = True, remove_primary_loops: bool = False, assume_box_large_enough: bool = False) -> MEHPNetwork:
        """
                Build a canonical MEHPNetwork from a coarse-grained CrosslinkerUniverse and optional StrandPairEntanglements.
        
                :param box: Simulation box.
                :param crosslinker_universe: Crosslinker & strand representation.
                :param entanglements: Optional pre-sampled strand-pair entanglements.
                :param crosslinker_type: Atom type identifying crosslinkers.
                :param entanglement_mode: Entanglement representation mode.
                :param remove_dangling_chains: Omit dangling chains when True.
                :param include_free_chains: Include free/soluble chains when True.
                :param remove_primary_loops: Remove polymer primary loops upon construction.
                :param assume_box_large_enough: Assume box large enough for minimum-image PBC.
                :return: Constructed MEHPNetwork.
        """
    @staticmethod
    def from_universe(universe: Universe, entanglements: AtomPairEntanglements = ..., crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, entanglement_mode: EntanglementMode = ..., keep_2functional_crosslinkers: bool = False, remove_primary_loops: bool = False, assume_box_large_enough: bool = False) -> MEHPNetwork:
        """
                Build a canonical MEHPNetwork from an atomistic Universe using DFS strand contraction.
        
                :param universe: Atomistic polymer Universe.
                :param entanglements: Optional pre-sampled atom-pair entanglements.
                :param crosslinker_type: Atom type identifying crosslinker beads.
                :param entanglement_mode: Entanglement representation mode.
                :param keep_2functional_crosslinkers: Keep degree-2 crosslinkers as explicit network nodes.
                :param remove_primary_loops: Remove polymer primary loops upon construction.
                :param assume_box_large_enough: Whether springs can be evaluated with minimum-image PBC without baked offsets.
                :return: Constructed MEHPNetwork.
        """
    @staticmethod
    def from_universe_with_random_entanglements(universe: Universe, nr_of_entanglements_to_sample: typing.SupportsInt | typing.SupportsIndex, upper_sampling_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_sampling_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, min_nr_of_entanglements: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, seed: str = '', crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, entanglement_mode: EntanglementMode = ..., with_crosslinkers: bool = True, filter_dangling_and_soluble: bool = False, keep_2functional_crosslinkers: bool = False, remove_primary_loops: bool = False, assume_box_large_enough: bool = False) -> MEHPNetwork:
        """
                Sample random entanglements from an atomistic Universe and build a canonical MEHPNetwork.
        
                :param universe: Atomistic polymer Universe.
                :param nr_of_entanglements_to_sample: Target number of entanglements to sample.
                :param upper_sampling_cutoff: Upper distance cutoff for pair sampling.
                :param lower_sampling_cutoff: Lower distance cutoff for pair sampling.
                :param min_nr_of_entanglements: Minimum entanglements required if possible.
                :param same_strand_cutoff: Minimum chemical distance for same-strand pairs.
                :param seed: PRNG seed string (empty for random seed).
                :param crosslinker_type: Atom type identifying crosslinker beads.
                :param entanglement_mode: Entanglement representation mode.
                :param with_crosslinkers: Exclude crosslinker beads as entanglement targets when True.
                :param filter_dangling_and_soluble: Omit entanglements touching dangling/soluble strands during sampling.
                :param keep_2functional_crosslinkers: Keep degree-2 crosslinkers as explicit nodes.
                :param remove_primary_loops: Remove polymer primary loops upon construction.
                :param assume_box_large_enough: Assume box large enough for minimum-image PBC.
                :return: Constructed MEHPNetwork.
        """
    def __getstate__(self) -> tuple[bytes]:
        ...
    def __init__(self) -> None:
        """
        Construct an empty MEHPNetwork.
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def compute_active_weight_fraction(self, u: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001, is_2d: bool = False) -> float:
        """
        Compute the weight fraction of active strands.
        """
    def compute_dangling_weight_fraction(self, u: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001, is_2d: bool = False) -> float:
        """
        Compute the weight fraction of dangling strands.
        """
    def compute_soluble_weight_fraction(self, u: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001, is_2d: bool = False) -> float:
        """
        Compute the weight fraction of soluble strands.
        """
    @typing.overload
    def deform_to(self, target_box: Box) -> None:
        """
        Affinely deform network coordinates and baked box offsets to target_box.
        """
    @typing.overload
    def deform_to(self, target_box: Box, u: typing.Any) -> None:
        """
        Affinely deform network coordinates, displacements u, and baked box offsets to target_box.
        """
    def evaluate_spring_lengths(self, u: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], is_2d: bool = False) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Evaluate Euclidean lengths of all springs.
        """
    def evaluate_spring_vectors(self, u: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], is_2d: bool = False) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Evaluate all flattened spring end-to-end vectors (size 3 * nr_of_springs).
        """
    def find_active_springs(self, u: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001, is_2d: bool = False) -> typing.Annotated[numpy.typing.NDArray[numpy.bool], "[m, 1]"]:
        """
        Return boolean mask of active springs.
        """
    def find_active_strands(self, u: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001, is_2d: bool = False) -> typing.Annotated[numpy.typing.NDArray[numpy.bool], "[m, 1]"]:
        """
        Return boolean mask of active strands.
        """
    def get_crosslinker_verse(self, u: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"] = ...) -> Universe:
        """
        Reconstruct a coarse-grained Universe of crosslinker nodes bonded by polymer strands.
        """
    def get_gamma_factor(self, u: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], b02: typing.SupportsFloat | typing.SupportsIndex = -1.0, nr_of_chains: typing.SupportsInt | typing.SupportsIndex = -1, is_2d: bool = False) -> float:
        """
        Compute the average topological Gamma factor ``sum(gamma_factors) / denom``.
        
        :param u: Displacement vector.
        :param b02: Mean-squared bond length :math:`b_0^2` (defaults to network bond length squared).
        :param nr_of_chains: Denominator for averaging. When ``< 1`` (default ``-1``), divides by
                             ``nr_of_strands`` — note that in ``EntanglementMode.SPRING``,
                             ``nr_of_strands`` includes entanglement coupling springs. Pass ``1`` (or
                             sum ``get_gamma_factors``) for the unnormalized sum ``sum(gamma_factors)``.
        :param is_2d: Whether to evaluate in 2D.
        """
    def simplify(self, u: typing.Any, mode: StructureSimplificationMode, tolerance: typing.SupportsFloat | typing.SupportsIndex = 1e-09, is_2d: bool = False, breaking_length: typing.SupportsFloat | typing.SupportsIndex = -1.0) -> int:
        """
        Run iterative topology simplification until convergence.
        """
    def validate(self, u: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"] = ...) -> None:
        """
        Validate network dimensions, index mappings, and contour invariants.
        """
    @property
    def L(self) -> typing.Annotated[list[float], "FixedSize(3)"]:
        ...
    @L.setter
    def L(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"]) -> None:
        ...
    @property
    def box(self) -> Box:
        ...
    @property
    def coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    @coordinates.setter
    def coordinates(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @property
    def default_bond_length(self) -> float:
        ...
    @default_bond_length.setter
    def default_bond_length(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def default_r02(self) -> float:
        ...
    @default_r02.setter
    def default_r02(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def link_indices_of_strand(self) -> list[list[int]]:
        ...
    @link_indices_of_strand.setter
    def link_indices_of_strand(self, arg0: collections.abc.Sequence[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]]) -> None:
        ...
    @property
    def link_is_entanglement(self) -> typing.Annotated[numpy.typing.NDArray[numpy.bool], "[m, 1]"]:
        ...
    @link_is_entanglement.setter
    def link_is_entanglement(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.bool, "[m, 1]"]) -> None:
        ...
    @property
    def mean_spring_contour_length(self) -> float:
        ...
    @mean_spring_contour_length.setter
    def mean_spring_contour_length(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def mean_strand_contour_length(self) -> float:
        ...
    @mean_strand_contour_length.setter
    def mean_strand_contour_length(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def nr_of_crosslink_swaps_endured(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @nr_of_crosslink_swaps_endured.setter
    def nr_of_crosslink_swaps_endured(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.int32, "[m, 1]"]) -> None:
        ...
    @property
    def nr_of_links(self) -> int:
        ...
    @nr_of_links.setter
    def nr_of_links(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def nr_of_nodes(self) -> int:
        ...
    @nr_of_nodes.setter
    def nr_of_nodes(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def nr_of_springs(self) -> int:
        ...
    @nr_of_springs.setter
    def nr_of_springs(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def nr_of_strands(self) -> int:
        ...
    @nr_of_strands.setter
    def nr_of_strands(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def old_atom_ids(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @old_atom_ids.setter
    def old_atom_ids(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.int32, "[m, 1]"]) -> None:
        ...
    @property
    def old_atom_types(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @old_atom_types.setter
    def old_atom_types(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.int32, "[m, 1]"]) -> None:
        ...
    @property
    def spring_box_offset(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    @spring_box_offset.setter
    def spring_box_offset(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @property
    def spring_contour_length(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    @spring_contour_length.setter
    def spring_contour_length(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @property
    def spring_index_a(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @spring_index_a.setter
    def spring_index_a(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.int32, "[m, 1]"]) -> None:
        ...
    @property
    def spring_index_b(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @spring_index_b.setter
    def spring_index_b(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.int32, "[m, 1]"]) -> None:
        ...
    @property
    def spring_indices_of_strand(self) -> list[list[int]]:
        ...
    @spring_indices_of_strand.setter
    def spring_indices_of_strand(self, arg0: collections.abc.Sequence[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]]) -> None:
        ...
    @property
    def spring_is_entanglement(self) -> typing.Annotated[numpy.typing.NDArray[numpy.bool], "[m, 1]"]:
        ...
    @spring_is_entanglement.setter
    def spring_is_entanglement(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.bool, "[m, 1]"]) -> None:
        ...
    @property
    def spring_partitions(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    @spring_partitions.setter
    def spring_partitions(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @property
    def strand_contour_length(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    @strand_contour_length.setter
    def strand_contour_length(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @property
    def strand_index_of_spring(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @strand_index_of_spring.setter
    def strand_index_of_spring(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.int32, "[m, 1]"]) -> None:
        ...
    @property
    def strand_indices_of_link(self) -> list[list[int]]:
        ...
    @strand_indices_of_link.setter
    def strand_indices_of_link(self, arg0: collections.abc.Sequence[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]]) -> None:
        ...
    @property
    def vol(self) -> float:
        ...
    @vol.setter
    def vol(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class MEHPNonlinearSolverConfig:
    """
    Configuration for general nonlinear optimization MEHP equilibrium solving via NLopt.
    """
    algorithm: NLoptAlgorithm
    simplification_mode: StructureSimplificationMode
    slip_mode: SlipMode
    def __init__(self, algorithm: NLoptAlgorithm = ..., max_nr_of_steps: typing.SupportsInt | typing.SupportsIndex = 250000, xtol: typing.SupportsFloat | typing.SupportsIndex = 1e-12, ftol: typing.SupportsFloat | typing.SupportsIndex = 1e-09, simplification_mode: StructureSimplificationMode = ..., inactive_removal_cutoff: typing.SupportsFloat | typing.SupportsIndex = 1e-06, spring_breaking_length: typing.SupportsFloat | typing.SupportsIndex = -1.0, slip_mode: SlipMode = ...) -> None:
        ...
    @property
    def ftol(self) -> float:
        ...
    @ftol.setter
    def ftol(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def inactive_removal_cutoff(self) -> float:
        ...
    @inactive_removal_cutoff.setter
    def inactive_removal_cutoff(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def max_nr_of_steps(self) -> int:
        ...
    @max_nr_of_steps.setter
    def max_nr_of_steps(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def spring_breaking_length(self) -> float:
        ...
    @spring_breaking_length.setter
    def spring_breaking_length(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def xtol(self) -> float:
        ...
    @xtol.setter
    def xtol(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class MEHPSimulator:
    """
    
          Unified Mechanical Equilibrium of Heterogeneous Polymers (MEHP) simulator
          coordinating a canonical MEHPNetwork, a stateless MEHPForceEvaluator,
          and an interchangeable MEHPSolver strategy.
        
    """
    force_evaluator: MEHPForceEvaluator
    is_2d: bool
    @staticmethod
    def construct_from_string(s: str) -> MEHPSimulator:
        ...
    @staticmethod
    def construct_with_random_entanglements(universe: Universe, nr_of_entanglements_to_sample: typing.SupportsInt | typing.SupportsIndex, upper_sampling_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_sampling_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, min_nr_of_entanglements: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, seed: str = '', crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, entanglement_mode: EntanglementMode = ..., with_crosslinkers: bool = True, filter_dangling_and_soluble: bool = False, force_evaluator: MEHPForceEvaluator = None, is_2d: bool = False, keep_2functional_crosslinkers: bool = False, remove_primary_loops: bool = False, assume_box_large_enough: bool = False) -> MEHPSimulator:
        """
        Construct MEHPSimulator from a Universe with randomly sampled entanglements.
        """
    @staticmethod
    def read_restart_file(file: str) -> MEHPSimulator:
        ...
    def __copy__(self) -> MEHPSimulator:
        ...
    def __getstate__(self) -> tuple[bytes]:
        ...
    @typing.overload
    def __init__(self, network: MEHPNetwork, force_evaluator: MEHPForceEvaluator = None, is_2d: bool = False) -> None:
        """
        Instantiate MEHPSimulator from an existing MEHPNetwork.
        """
    @typing.overload
    def __init__(self, universe: Universe, entanglements: AtomPairEntanglements = ..., crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, entanglement_mode: EntanglementMode = ..., force_evaluator: MEHPForceEvaluator = None, is_2d: bool = False, keep_2functional_crosslinkers: bool = False, remove_primary_loops: bool = False, assume_box_large_enough: bool = False) -> None:
        """
        Instantiate MEHPSimulator from an atomistic Universe.
        """
    @typing.overload
    def __init__(self, box: Box, crosslinker_universe: CrosslinkerUniverse, entanglements: StrandPairEntanglements = ..., crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, entanglement_mode: EntanglementMode = ..., force_evaluator: MEHPForceEvaluator = None, is_2d: bool = False, remove_dangling_chains: bool = False, include_free_chains: bool = True, remove_primary_loops: bool = False, assume_box_large_enough: bool = False) -> None:
        """
        Instantiate MEHPSimulator from a Box and a CrosslinkerUniverse.
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    @typing.overload
    def add_slip_links(self, entanglements: StrandPairEntanglements, clamp_alpha: bool = False) -> int:
        ...
    @typing.overload
    def add_slip_links(self, strand_idx_1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], strand_idx_2: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], y: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], z: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], alpha_1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], alpha_2: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], clamp_alpha: bool = False) -> int:
        ...
    @typing.overload
    def add_slip_links(self, strand_idx_1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], strand_idx_2: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], y: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], z: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], clamp_alpha: bool = False) -> int:
        ...
    @typing.overload
    def add_sliplinks(self, entanglements: StrandPairEntanglements, clamp_alpha: bool = False) -> int:
        ...
    @typing.overload
    def add_sliplinks(self, strand_idx_1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], strand_idx_2: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], y: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], z: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], alpha_1: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], alpha_2: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], clamp_alpha: bool = False) -> int:
        ...
    @typing.overload
    def add_sliplinks(self, strand_idx_1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], strand_idx_2: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], x: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], y: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], z: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], clamp_alpha: bool = False) -> int:
        ...
    def break_too_long_strands(self, breaking_length: typing.SupportsFloat | typing.SupportsIndex) -> int:
        ...
    def config_assume_box_large_enough(self, assume: bool = True) -> None:
        ...
    def config_assume_network_is_complete(self, assume: bool = False) -> None:
        ...
    def config_auto_correlator_output(self, values: collections.abc.Sequence[OutputConfiguration], num_corr_in: typing.SupportsInt | typing.SupportsIndex = 32, p: typing.SupportsInt | typing.SupportsIndex = 16, m: typing.SupportsInt | typing.SupportsIndex = 2) -> None:
        ...
    def config_average_output(self, values: collections.abc.Sequence[OutputConfiguration]) -> None:
        ...
    def config_mean_bond_length(self, mean_bond_length: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def config_mean_squared_bond_length(self, r02: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def config_remove_primary_loops(self, remove: bool = False) -> None:
        ...
    def config_restart_output(self, file: str, output_every: typing.SupportsInt | typing.SupportsIndex = 50000) -> None:
        ...
    def config_step_output(self, output_configuration: collections.abc.Sequence[OutputConfiguration]) -> None:
        ...
    def count_active_clustered_atoms(self, tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001) -> float:
        ...
    def deform_to(self, target_box: Box) -> None:
        ...
    def get_active_weight_fraction(self, tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001) -> float:
        ...
    def get_average_contour_length(self) -> float:
        ...
    def get_average_spring_length(self) -> float:
        ...
    def get_bond_lengths(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_crosslinker_universe(self) -> Universe:
        ...
    def get_crosslinker_verse(self) -> Universe:
        ...
    def get_current_displacements(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_current_time(self, current_step: typing.SupportsFloat | typing.SupportsIndex = 0.0) -> float:
        ...
    def get_dangling_weight_fraction(self, tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001) -> float:
        ...
    def get_default_bond_length(self) -> float:
        ...
    def get_default_r0_square(self) -> float:
        ...
    def get_effective_functionality_of_atoms(self, tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001) -> dict[int, int]:
        ...
    def get_energy(self) -> float:
        ...
    def get_exit_reason(self) -> ExitReason:
        ...
    def get_force_evaluator(self) -> MEHPForceEvaluator:
        ...
    def get_gamma(self) -> float:
        ...
    def get_gamma_factor(self, b02: typing.SupportsFloat | typing.SupportsIndex = -1.0, nr_of_chains: typing.SupportsInt | typing.SupportsIndex = -1) -> float:
        """
        Compute the average topological Gamma factor ``sum(gamma_factors) / denom``.
        
        :param b02: Mean-squared bond length :math:`b_0^2` (defaults to network bond length squared).
        :param nr_of_chains: Denominator for averaging. When ``< 1`` (default ``-1``), divides by
                             ``network.nr_of_strands`` — note that in ``EntanglementMode.SPRING``,
                             ``nr_of_strands`` includes entanglement coupling springs, so comparing
                             ``MERGED_LINKS`` vs. ``SPRING`` with the default ``-1`` is not apples-to-apples.
                             Pass ``1`` or use ``np.sum(get_gamma_factors(b02))`` when computing shear modulus.
        """
    def get_gamma_factors(self, b02: typing.SupportsFloat | typing.SupportsIndex = -1.0) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Compute the per-spring topological Gamma factors :math:`\\Gamma_i = \\mathbf{r}_i^2 / (N_i b_0^2)`.
        
        Summing this vector (``np.sum(sim.get_gamma_factors(b02))``) gives the unnormalized total
        needed to compute the shear modulus :math:`G = (k_B T / V) \\sum_i \\Gamma_i`.
        
        :param b02: Mean-squared bond length :math:`b_0^2` (defaults to network bond length squared).
        """
    def get_gamma_factors_in_dir(self, dir: typing.SupportsInt | typing.SupportsIndex, b02: typing.SupportsFloat | typing.SupportsIndex = -1.0) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Compute the per-spring directional Gamma factors :math:`\\Gamma_{i,\\alpha} = r_{i,\\alpha}^2 / (N_i b_0^2)`.
        
        :param dir: Cartesian direction index (``0``, ``1``, or ``2``).
        :param b02: Mean-squared bond length :math:`b_0^2` (defaults to network bond length squared).
        """
    def get_indices_of_active_nodes(self, tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001) -> list[int]:
        ...
    def get_initial_coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_is_2d(self) -> bool:
        ...
    def get_last_solve_result(self) -> MEHPSolveResult:
        ...
    def get_network(self) -> MEHPNetwork:
        ...
    def get_nr_of_active_nodes(self, tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001) -> int:
        ...
    def get_nr_of_active_springs(self, tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001) -> int:
        ...
    def get_nr_of_active_strands(self, tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001) -> int:
        ...
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
        ...
    def get_nr_of_links(self) -> int:
        ...
    def get_nr_of_nodes(self) -> int:
        ...
    def get_nr_of_particles(self) -> int:
        ...
    def get_nr_of_springs(self) -> int:
        ...
    def get_nr_of_strands(self) -> int:
        ...
    def get_num_bonds_to_form(self) -> int:
        ...
    def get_num_intra_chain_slip_links(self) -> int:
        ...
    def get_num_intra_chain_sliplinks(self) -> int:
        ...
    def get_num_relocations(self) -> int:
        ...
    def get_num_shifts(self) -> int:
        ...
    def get_pressure(self) -> float:
        ...
    def get_residual(self) -> float:
        ...
    def get_residual_norm(self) -> float:
        ...
    def get_residuals(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_simulation_has_run(self) -> bool:
        ...
    def get_soluble_weight_fraction(self, tolerance: typing.SupportsFloat | typing.SupportsIndex = 0.001) -> float:
        ...
    def get_spring_contour_length(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_spring_lengths(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_spring_partitions(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_spring_vectors(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_strand_contour_length(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_stress_tensor(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 3]"]:
        ...
    def get_temperature(self) -> float:
        ...
    def get_timestep(self) -> float:
        ...
    def get_volume(self) -> float:
        ...
    def has_outputs(self) -> bool:
        ...
    def randomly_add_slip_links(self, nr_of_entanglements_to_sample: typing.SupportsInt | typing.SupportsIndex, cutoff: typing.SupportsFloat | typing.SupportsIndex = 2.0, min_nr_of_entanglements: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 2.0, seed: str = '', clamp_alpha: bool = False) -> int:
        ...
    def randomly_add_sliplinks(self, nr_of_entanglements_to_sample: typing.SupportsInt | typing.SupportsIndex, cutoff: typing.SupportsFloat | typing.SupportsIndex = 2.0, min_nr_of_entanglements: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 2.0, seed: str = '', clamp_alpha: bool = False) -> int:
        ...
    def remove_inactive_links(self, tolerance: typing.SupportsFloat | typing.SupportsIndex = 1e-06) -> int:
        ...
    @typing.overload
    def run(self) -> MEHPSolveResult:
        ...
    @typing.overload
    def run(self, solver: MEHPSolver) -> MEHPSolveResult:
        ...
    @typing.overload
    def run(self, config: MEHPLinearSolverConfig) -> MEHPSolveResult:
        ...
    @typing.overload
    def run(self, config: MEHPSlipLinkSolverConfig) -> MEHPSolveResult:
        ...
    @typing.overload
    def run(self, config: MEHPNonlinearSolverConfig) -> MEHPSolveResult:
        ...
    def serialize_to_string(self) -> bytes:
        ...
    def set_current_displacements(self, new_displacements: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    def set_force_evaluator(self, force_evaluator: MEHPForceEvaluator) -> None:
        ...
    def set_is_2d(self, is_2d: bool) -> None:
        ...
    def simplify_network(self, mode: StructureSimplificationMode, tolerance: typing.SupportsFloat | typing.SupportsIndex = 1e-06, breaking_length: typing.SupportsFloat | typing.SupportsIndex = -1.0) -> int:
        ...
    def unlink_bifunctional_links(self, entanglements_only: bool = False) -> int:
        ...
    def validate_network(self) -> None:
        ...
    def write_restart_file(self, file: str) -> None:
        ...
    @property
    def last_solve_result(self) -> MEHPSolveResult:
        ...
    @property
    def network(self) -> MEHPNetwork:
        ...
class MEHPSlipLinkSolverConfig:
    """
    Configuration for alternating coordinate/partition slip-link relaxation.
    """
    do_inner_iterations: bool
    link_swapping_mode: LinkSwappingMode
    simplification_mode: StructureSimplificationMode
    slip_mode: SlipMode
    def __init__(self, max_nr_of_steps: typing.SupportsInt | typing.SupportsIndex = 50000, xtol: typing.SupportsFloat | typing.SupportsIndex = 1e-05, initial_residual_to_use: typing.SupportsFloat | typing.SupportsIndex = -1.0, simplification_mode: StructureSimplificationMode = ..., inactive_removal_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.001, do_inner_iterations: bool = False, link_swapping_mode: LinkSwappingMode = ..., swapping_frequency: typing.SupportsInt | typing.SupportsIndex = 10, one_over_spring_partition_upper_limit: typing.SupportsFloat | typing.SupportsIndex = 1.0, nr_of_crosslink_swaps_allowed_per_sliplink: typing.SupportsInt | typing.SupportsIndex = -1, slip_mode: SlipMode = ..., spring_breaking_length: typing.SupportsFloat | typing.SupportsIndex = -1.0) -> None:
        ...
    @property
    def inactive_removal_cutoff(self) -> float:
        ...
    @inactive_removal_cutoff.setter
    def inactive_removal_cutoff(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def initial_residual_to_use(self) -> float:
        ...
    @initial_residual_to_use.setter
    def initial_residual_to_use(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def max_nr_of_steps(self) -> int:
        ...
    @max_nr_of_steps.setter
    def max_nr_of_steps(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def nr_of_crosslink_swaps_allowed_per_sliplink(self) -> int:
        ...
    @nr_of_crosslink_swaps_allowed_per_sliplink.setter
    def nr_of_crosslink_swaps_allowed_per_sliplink(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def one_over_spring_partition_upper_limit(self) -> float:
        ...
    @one_over_spring_partition_upper_limit.setter
    def one_over_spring_partition_upper_limit(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def spring_breaking_length(self) -> float:
        ...
    @spring_breaking_length.setter
    def spring_breaking_length(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def swapping_frequency(self) -> int:
        ...
    @swapping_frequency.setter
    def swapping_frequency(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def xtol(self) -> float:
        ...
    @xtol.setter
    def xtol(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class MEHPSolveResult:
    """
    Diagnostic summary returned by an MEHPSolver execution.
    """
    def __init__(self) -> None:
        ...
    @property
    def cg_fallback_triggered(self) -> bool:
        ...
    @property
    def effective_sle_solver(self) -> SLESolver:
        ...
    @property
    def energy(self) -> float:
        ...
    @property
    def exit_reason(self) -> ExitReason:
        ...
    @property
    def nr_of_iterations(self) -> int:
        ...
    @property
    def nr_of_relocations(self) -> int:
        ...
    @property
    def nr_of_shifts(self) -> int:
        ...
    @property
    def one_over_spring_partition_upper_limit_used(self) -> float:
        ...
    @property
    def residual(self) -> float:
        ...
class MEHPSolver:
    """
    Abstract strategy interface for solving mechanical equilibrium on an MEHPNetwork.
    """
    def solve(self, net: MEHPNetwork, u: typing.Any, force_evaluator: MEHPForceEvaluator, is_2d: bool = False) -> MEHPSolveResult:
        """
        Solve for equilibrium link displacements on net starting from u.
        """
class MaxDistanceProvider:
    """
    
         A generic implementation of a class, that shall provide a maximum distance for the MC sampling.
         
    """
    def get_max_distance(self, N: typing.SupportsFloat | typing.SupportsIndex) -> float:
        """
                 Get the maximum distance for a given N.
        
                 :param N: Number of segments.
                 :return: Maximum distance for sampling.
        """
class Molecule:
    """
    
           An (ideally) connected series of atoms/beads.
      
    """
    __hash__: typing.ClassVar[None] = None
    def __contains__(self, arg0: Atom) -> bool:
        """
                  Check whether a particular atom is contained in this molecule.
        """
    def __copy__(self) -> Molecule:
        """
                 Create a copy of this molecule.
                 
                 :return: A new Molecule instance that is a copy of this one
        """
    def __eq__(self, arg0: Molecule) -> bool:
        ...
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> Atom:
        """
               Access an atom by its vertex index.
        """
    def __init__(self, other: Molecule) -> None:
        """
                 Construct a copy of a molecule.
        
                 :param other: The molecule to copy
        """
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
                 Compute the length of each bond in the molecule, respecting periodic boundaries.
                 
                 :return: A vector of bond lengths
        """
    def compute_end_to_end_distance(self) -> float:
        """
                    Compute the end-to-end distance (:math:`R_{ee}`) of this molecule. 
        
                    .. warning::
                       Returns 0.0 if the molecule does not have two or more atoms.
                       
                    :return: The end-to-end distance
        """
    def compute_end_to_end_distance_with_derived_image_flags(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2) -> float:
        """
                    Compute the end-to-end distance (:math:`R_{ee}`) of this molecule,
                    but ignoring the image flags attached to the atoms. 
                    This only works for Molecules that can be lined up with 
                    :meth:`~pylimer_tools_cpp.Molecule.get_atoms_lined_up`,
                    as it needs the atoms sorted such that the periodic box can still be respected somewhat.
        
                    .. warning::
                       Returns 0.0 if the molecule does not have two or more atoms.
                       Requires bonds to be shorter than half the box length.
                       
                    :return: The end-to-end distance with derived image flags
        """
    def compute_end_to_end_vector(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                    Compute the end-to-end vector (:math:`\\overrightarrow{R}_{ee}`) of this molecule. 
        
                    .. warning::
                       Returns 0.0 if the molecule does not have two or more atoms.
                       
                    :return: The end-to-end vector
        """
    def compute_end_to_end_vector_with_derived_image_flags(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                    Compute the end-to-end vector (:math:`\\overrightarrow{R}_{ee}`) of this molecule,
                    but ignoring the image flags attached to the atoms. 
                    This only works for Molecules that can be lined up with 
                    :meth:`~pylimer_tools_cpp.Molecule.get_atoms_lined_up`,
                    as it needs the atoms sorted such that the periodic box can still be respected somewhat.
        
                    .. warning::
                       Returns 0.0 if the molecule does not have two or more atoms.
                       Requires bonds to be shorter than half the box length.
                       
                    :return: The end-to-end vector with derived image flags
        """
    def compute_radius_of_gyration(self) -> float:
        """
                    Compute the radius of gyration, :math:`R_g^2` of this molecule.
                    
                    :math:`{R_g}^2 = \\\\frac{1}{M} \\sum_i m_i (r_i - r_{cm})^2`,
                    where :math:`M` is the total mass of the molecule, :math:`r_{cm}`
                    are the coordinates of the center of mass of the molecule and the
                    sum is over all contained atoms.
                    
                    :return: The radius of gyration squared
        """
    def compute_radius_of_gyration_with_derived_image_flags(self) -> float:
        """
                    Compute the radius of gyration, :math:`R_g^2` of this molecule,
                    but ignoring the image flags attached to the atoms.
                    This only works for Molecules that can be lined up with 
                    :meth:`~pylimer_tools_cpp.Molecule.get_atoms_lined_up`,
                    as it needs the atoms sorted such that the periodic box can still be respected somewhat.
                    In other words, this function computes the radius of gyration 
                    assuming the distance between two lined-up beads 
                    is smaller than half the periodic box in each direction.
                    
                    See also: :meth:`~pylimer_tools_cpp.Molecule.compute_radius_of_gyration`.
                    
                    :return: The radius of gyration squared with derived image flags
        """
    def compute_total_length(self) -> float:
        """
                 Compute the sum of the lengths of all bonds.
                 In most cases, this is equal to the contour length.
                 
                 :return: The total contour length of the molecule
        """
    def compute_total_mass(self) -> float:
        """
                    Compute the total mass of this molecule.
                    
                    :return: The total mass of all atoms in this molecule
        """
    def compute_total_vector(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, close_loop: bool = True) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                       Compute the sum of all bond vectors.
                       
                       :param crosslinker_type: The type of crosslinker atoms
                       :param close_loop: Whether to close the loop for calculations
                       :return: The vector sum of all bonds
        """
    def compute_vector_from_to(self, atom_id_from: typing.SupportsInt | typing.SupportsIndex, atom_id_to: typing.SupportsInt | typing.SupportsIndex, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, require_order: bool = True) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                       Compute the sum of all bond vectors between two specified atoms.
                       
                       :param atom_id_from: ID of the starting atom
                       :param atom_id_to: ID of the ending atom
                       :param crosslinker_type: The type of crosslinker atoms
                       :param require_order: Whether to require specific ordering
                       :return: The vector sum from start to end atom
        """
    def get_atom_by_id(self, atom_id: typing.SupportsInt | typing.SupportsIndex) -> Atom:
        """
                    Get an atom by its id.
                    
                    :param atom_id: The atom ID to search for
                    :return: The atom with the specified ID
        """
    def get_atom_by_vertex_idx(self, vertex_idx: typing.SupportsInt | typing.SupportsIndex) -> Atom:
        """
                    Get an atom for a specific vertex.
                    
                    :param vertex_idx: The vertex index to query
                    :return: The atom at the specified vertex
        """
    def get_atom_id_by_vertex_idx(self, vertex_id: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
                 Get the ID of the atom by the vertex ID of the underlying graph.
                 
                 :param vertex_id: The vertex index in the underlying graph
                 :return: The atom ID corresponding to the vertex
        """
    def get_atom_types(self) -> list[int]:
        """
                 Query all types (each one for each atom) ordered by atom vertex id.
                 
                 :return: A vector of atom types in vertex order
        """
    def get_atoms(self) -> list[Atom]:
        """
                    Return all atom objects enclosed in this molecule, ordered by vertex id.
                    
                    :return: List of atoms in vertex order
        """
    def get_atoms_by_degree(self, degree: typing.SupportsInt | typing.SupportsIndex) -> list[Atom]:
        """
                    Get the atoms that have the specified number of bonds.
                    
                    :param degree: The number of bonds (degree/functionality)
                    :return: List of atoms with the specified degree
        """
    def get_atoms_by_type(self, type: typing.SupportsInt | typing.SupportsIndex) -> list[Atom]:
        """
                    Get the atoms with the specified type.
                    
                    :param type: The atom type to search for
                    :return: List of atoms with the specified type
        """
    def get_atoms_connected_to(self, atom: Atom) -> list[Atom]:
        """
                    Get the atoms connected to a specified atom.
        
                    Internally uses :meth:`~pylimer_tools_cpp.Molecule.get_atoms_connected_to`.
                    
                    :param atom: The atom to query connections for
                    :return: List of connected atoms
        """
    def get_atoms_connected_to_vertex(self, vertex_idx: typing.SupportsInt | typing.SupportsIndex) -> list[Atom]:
        """
                    Get the atoms connected to a specified vertex id.
                    
                    :param vertex_idx: The vertex index to query
                    :return: List of connected atoms
        """
    def get_atoms_lined_up(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, assumed_coordinates: bool = False, close_loop: bool = False) -> list[Atom]:
        """
                    Return all atom objects enclosed in this molecule based on the connectivity.
        
                    This method works only for lone chains, atoms and loops, 
                    as it throws an error if the molecule does not allow such a "line-up", 
                    for example because of crosslinkers.
        
                    Use the `crosslinker_type` parameter to force the atoms in a primary loop 
                    to start with the crosslink.
                    
                    :param crosslinker_type: The type of crosslinker atoms
                    :param assumed_coordinates: Whether to assume coordinates are valid
                    :param close_loop: Whether to close loops
                    :return: List of atoms in connected order
        """
    def get_bonds(self) -> dict[str, list[int]]:
        """
                    Get all bonds. Returns a dict with three properties: 'bond_from', 'bond_to' and 'bond_type'.
                    
                    :return: Dictionary with bond information
        """
    def get_edge_ids_from(self, vertex_id: typing.SupportsInt | typing.SupportsIndex) -> list[int]:
        """
                 Get the edge IDs incident to a specific vertex.
                 
                 :param vertex_id: The vertex to query
                 :return: A vector of edge IDs connected to the vertex
        """
    def get_edge_ids_from_to(self, vertex_id_from: typing.SupportsInt | typing.SupportsIndex, vertex_id_to: typing.SupportsInt | typing.SupportsIndex) -> list[int]:
        """
                 Get the edge IDs of the edges between two specific vertices.
                 
                 :param vertex_id_from: Starting vertex ID
                 :param vertex_id_to: Ending vertex ID
                 :return: Vector of edge IDs between the specified vertices
        """
    def get_edges(self) -> dict[str, list[int]]:
        """
                    Get all bonds. Returns a dict with three properties: 'edge_from', 'edge_to' and 'edge_type'.
                    The order is not necessarily related to any structural property.
                    
                    .. note::
                       The integer values returned refer to the vertex ids, not the atom ids.
                       Use :meth:`~pylimer_tools_cpp.Molecule.get_atom_id_by_idx` to translate them to atom ids, or 
                       :meth:`~pylimer_tools_cpp.Molecule.get_bonds` to have that done for you.
                       
                    :return: Dictionary with edge information
        """
    def get_key(self) -> str:
        """
                    Get a unique identifier for this molecule.
                    
                    :return: A unique string identifier for this molecule
        """
    def get_nr_of_atoms(self) -> int:
        """
                 Count and return the number of atoms associated with this molecule.
                 
                 :return: The number of atoms in this molecule
        """
    def get_nr_of_bonds(self) -> int:
        """
                 Count and return the number of bonds associated with this molecule.
                 
                 :return: The number of bonds in this molecule
        """
    def get_nr_of_edges_from_to(self, vertex_id_from: typing.SupportsInt | typing.SupportsIndex, vertex_id_to: typing.SupportsInt | typing.SupportsIndex, max_length: typing.SupportsInt | typing.SupportsIndex = -1) -> int:
        """
                 Get the number of edges in the shortest path between two specific vertices.
        
                 If `max_length` is provided and positive, it will only consider paths up to that length.
                 
                 :param vertex_id_from: Starting vertex ID
                 :param vertex_id_to: Ending vertex ID  
                 :param max_length: Maximum path length to consider (-1 for no limit)
                 :return: Number of edges in the shortest path, or -1 if no path exists
        """
    def get_strand_ends(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, close_loop: bool = False) -> list[Atom]:
        """
                  Get the ends of the given strand (= molecule).
                  In case of a primary loop, the crosslink is returned, if there is one.
                  Use the argument `close_loop` to decide, whether this should be returned once or twice.
        
                  .. note:: 
                       Currently only works for linear strands.
                       
                  :param crosslinker_type: The type of crosslinker atoms
                  :param close_loop: Whether to return the crosslinker twice for loops
                  :return: List of end atoms
        """
    def get_strand_type(self) -> MoleculeType:
        """
                   Get the type of this molecule (see :class:`~pylimer_tools_cpp.MoleculeType` enum).
        
                   .. note:: 
                      This type might be unset; currently, only 
                      :meth:`~pylimer_tools_cpp.Universe.get_chains_with_crosslinker` assigns them automatically.
        """
    def get_vertex_idx_by_atom_id(self, atom_id: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
                 Get the vertex ID of the underlying graph for an atom with a specified ID.
                 
                 :param atom_id: The atom ID to look up
                 :return: The vertex index corresponding to the atom
        """
class MoleculeIterator:
    """
    
           An iterator to iterate through the atoms in :class:`~pylimer_tools_cpp.Molecule`.
      
    """
    def __iter__(self) -> MoleculeIterator:
        ...
    def __next__(self) -> Atom:
        ...
class MoleculeType(enum.Enum):
    """
    An enum representing the type of molecule/chain/strand.
    """
    DANGLING_CHAIN: typing.ClassVar[MoleculeType]  # value = <MoleculeType.DANGLING_CHAIN: 3>
    FREE_CHAIN: typing.ClassVar[MoleculeType]  # value = <MoleculeType.FREE_CHAIN: 4>
    NETWORK_STRAND: typing.ClassVar[MoleculeType]  # value = <MoleculeType.NETWORK_STRAND: 1>
    PRIMARY_LOOP: typing.ClassVar[MoleculeType]  # value = <MoleculeType.PRIMARY_LOOP: 2>
    UNDEFINED: typing.ClassVar[MoleculeType]  # value = <MoleculeType.UNDEFINED: 0>
class NLoptAlgorithm(enum.Enum):
    """
    Supported NLopt optimization algorithms for MEHP force relaxation.
    """
    GD_STOGO: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.GD_STOGO: 18>
    GN_CRS2_LM: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.GN_CRS2_LM: 17>
    GN_DIRECT: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.GN_DIRECT: 16>
    LD_LBFGS: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LD_LBFGS: 1>
    LD_MMA: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LD_MMA: 0>
    LD_SLSQP: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LD_SLSQP: 2>
    LD_TNEWTON: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LD_TNEWTON: 3>
    LD_TNEWTON_PRECOND: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LD_TNEWTON_PRECOND: 5>
    LD_TNEWTON_PRECOND_RESTART: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LD_TNEWTON_PRECOND_RESTART: 6>
    LD_TNEWTON_RESTART: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LD_TNEWTON_RESTART: 4>
    LD_VAR1: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LD_VAR1: 7>
    LD_VAR2: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LD_VAR2: 8>
    LN_BOBYQA: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LN_BOBYQA: 10>
    LN_COBYLA: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LN_COBYLA: 12>
    LN_NELDERMEAD: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LN_NELDERMEAD: 11>
    LN_NEWUOA: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LN_NEWUOA: 13>
    LN_NEWUOA_BOUND: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LN_NEWUOA_BOUND: 14>
    LN_PRAXIS: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LN_PRAXIS: 15>
    LN_SBPLX: typing.ClassVar[NLoptAlgorithm]  # value = <NLoptAlgorithm.LN_SBPLX: 9>
class NoMaxDistanceProvider(MaxDistanceProvider):
    """
    
        For MC generation, to disable the neighbour list usage.
        
    """
    def __init__(self) -> None:
        ...
    def __repr__(self) -> str:
        """
                 Return a string representation of the NoMaxDistanceProvider.
                 
                 :return: String representation
        """
    def get_max_distance(self, N: typing.SupportsFloat | typing.SupportsIndex) -> float:
        """
                 Get the maximum distance for a given N (always returns -1 to disable).
        
                 :param N: Number of segments (ignored).
                 :return: Always returns -1 to disable maximum distance checks.
        """
class NonGaussianSpringForceEvaluator(MEHPForceEvaluator):
    """
    
          Non-Gaussian freely-jointed-chain (inverse Langevin) spring evaluator.
    
          :param kappa: Stiffness parameter (>= 0).
          :param l: Kuhn segment length (> 0).
        
    """
    def __init__(self, kappa: typing.SupportsFloat | typing.SupportsIndex = 1.0, l: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> None:
        ...
    def get_segment_length(self) -> float:
        """
        Return the Kuhn segment length l.
        """
    @property
    def l(self) -> float:
        """
        Kuhn segment length l.
        """
class NonlinearOptimizationSolver(MEHPSolver):
    """
    General nonlinear energy minimization solver for MEHP networks using NLopt.
    """
    config: MEHPNonlinearSolverConfig
    def __init__(self, config: MEHPNonlinearSolverConfig = ...) -> None:
        ...
    def get_config(self) -> MEHPNonlinearSolverConfig:
        ...
    def set_config(self, config: MEHPNonlinearSolverConfig) -> None:
        ...
class NormalModeAnalyzer:
    """
    
        Stateless graph Laplacian builder and normal mode solver.
    
        Please cite :cite:t:`gusev_molecular_2024` if you use this method in your work.
        
    """
    @staticmethod
    def from_matrix(matrix: scipy.sparse.csc_matrix) -> NormalModeAnalyzer:
        """
                        Construct a NormalModeAnalyzer directly from a sparse Laplacian/connectivity matrix.
        
                        :param matrix: Square sparse connectivity matrix
                        :return: NormalModeAnalyzer instance
        """
    @staticmethod
    def from_universe(universe: Universe) -> NormalModeAnalyzer:
        """
                        Construct a NormalModeAnalyzer from a Universe's atom-bond graph.
        
                        :param universe: Universe instance
                        :return: NormalModeAnalyzer instance
        """
    def __getstate__(self) -> tuple[bytes]:
        ...
    def __init__(self, spring_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], spring_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], nr_of_nodes: typing.SupportsInt | typing.SupportsIndex = 0, weights: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex] | None = None) -> None:
        """
                 Initialize the NormalModeAnalyzer from bond/spring edge lists.
        
                 :param spring_from: Starting node indices for springs/bonds
                 :param spring_to: Ending node indices for springs/bonds
                 :param nr_of_nodes: Optional total number of nodes (to include isolated vertices)
                 :param weights: Optional per-spring weights (defaults to unit weight 1.0)
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def chebyshev_moments(self, config: KPMConfig = ...) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                 Compute Jackson-damped Chebyshev moments using KPM.
        
                 :param config: KPMConfig instance
                 :return: Vector of damped Chebyshev moments
        """
    def density_of_states(self, lambda_grid: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], config: KPMConfig = ...) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                 Compute continuous density of states (DOS) on lambda_grid using KPM.
        
                 :param lambda_grid: Eigenvalue points at which to evaluate DOS
                 :param config: KPMConfig instance
                 :return: Density of states values on the grid
        """
    def solve_approx(self, config: pylimer_tools_cpp.SLQConfig | pylimer_tools_cpp.HybridConfig | pylimer_tools_cpp.KPMConfig | None = None) -> NormalModeSpectrum:
        """
              Compute an approximate / stochastic normal mode spectrum.
        
              :param config: SLQConfig, HybridConfig, KPMConfig, or None for AUTO
              :return: NormalModeSpectrum instance
        """
    def solve_exact(self, solver: NormalModeSolver = ..., compute_eigenvectors: bool = False, decompose_components: bool = True, num_modes: typing.SupportsInt | typing.SupportsIndex | None = None) -> NormalModeSpectrum:
        """
                 Compute exact eigenvalues (and optionally eigenvectors).
        
                 :param solver: Solver algorithm (NormalModeSolver enum; default: NormalModeSolver.AUTO). Ignored when ``num_modes < size``.
                 :param compute_eigenvectors: Whether to also compute eigenvectors (default: False)
                 :param decompose_components: Whether to decompose into connected components (default: True). Ignored when ``num_modes < size``.
                 :param num_modes: Optional number of smallest modes to compute via Spectra sparse solver
                 :return: NormalModeSpectrum instance
        """
    @property
    def matrix(self) -> scipy.sparse.csc_matrix:
        """
        Assembled sparse connectivity matrix.
        """
    @property
    def num_soluble_clusters(self) -> int:
        """
        Number of connected components (soluble clusters / zero modes).
        """
    @property
    def size(self) -> int:
        """
        Dimension of the connectivity matrix.
        """
class NormalModeSolver(enum.Enum):
    """
    Exact eigenvalue solver algorithms for normal mode analysis
    """
    AUTO: typing.ClassVar[NormalModeSolver]  # value = <NormalModeSolver.AUTO: 0>
    CUDA: typing.ClassVar[NormalModeSolver]  # value = <NormalModeSolver.CUDA: 5>
    DENSE_EIGEN: typing.ClassVar[NormalModeSolver]  # value = <NormalModeSolver.DENSE_EIGEN: 2>
    DENSE_LAPACK: typing.ClassVar[NormalModeSolver]  # value = <NormalModeSolver.DENSE_LAPACK: 1>
    METAL: typing.ClassVar[NormalModeSolver]  # value = <NormalModeSolver.METAL: 4>
    SPARSE_SLICING: typing.ClassVar[NormalModeSolver]  # value = <NormalModeSolver.SPARSE_SLICING: 3>
class NormalModeSpectrum:
    """
    
        Immutable normal mode spectrum and viscoelastic observable evaluator.
    
        Represents any discrete or quadrature-weighted spectrum :math:`\\{(\\lambda_i, w_i)\\}`.
        
    """
    def __getstate__(self) -> tuple[bytes]:
        ...
    def __init__(self, eigenvalues: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], weights: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"] | None = None, eigenvectors: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"] | None = None, is_complete: bool = True) -> None:
        """
                 Construct a NormalModeSpectrum from eigenvalues and optional weights/eigenvectors.
        
                 :param eigenvalues: 1D array of eigenvalues
                 :param weights: Optional 1D array of quadrature weights (defaults to all 1.0)
                 :param eigenvectors: Optional 2D array of eigenvectors (columns correspond to modes)
                 :param is_complete: Whether this spectrum covers all modes of the system
        """
    def __len__(self) -> int:
        ...
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def loss_modulus(self, omega: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                 Evaluate the loss modulus :math:`G''(\\omega)` (without prefactor).
        
                 :param omega: Angular frequencies
                 :return: Loss modulus values
        """
    def moduli(self, omega: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> tuple[typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"], typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]]:
        """
                 Evaluate both storage and loss moduli :math:`(G'(\\omega), G''(\\omega))`
                 in a single pass.
        
                 :param omega: Angular frequencies
                 :return: Tuple of (storage_modulus, loss_modulus)
        """
    def storage_modulus(self, omega: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                 Evaluate the storage modulus :math:`G'(\\omega)` (without prefactor).
        
                 :param omega: Angular frequencies
                 :return: Storage modulus values
        """
    def stress_autocorrelation(self, t: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                 Evaluate stress autocorrelation :math:`G(t)`.
        
                 :param t: Time points at which to evaluate
                 :return: Stress autocorrelation values
        """
    @property
    def eigenvalues(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Eigenvalues of the spectrum.
        """
    @property
    def eigenvectors(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, n]"] | None:
        """
        Eigenvectors matrix (columns are modes), or None if not computed.
        """
    @property
    def is_complete(self) -> bool:
        """
        Whether this spectrum is complete.
        """
    @property
    def num_modes(self) -> int:
        """
        Number of modes/quadrature nodes in the spectrum.
        """
    @property
    def num_zero_modes(self) -> int:
        """
        Number of zero modes in the spectrum.
        """
    @property
    def weights(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Quadrature weights of the spectrum.
        """
class OutputConfiguration:
    """
    
         A configuration object to configure the output values and frequency
         for simulation classes in this package.
    
         This class specifies which quantities to output and how often to write them
         during simulations.
        
    """
    def __init__(self) -> None:
        """
                  Create a new OutputConfiguration instance.
        
                  :return: A new OutputConfiguration object with default settings
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
    def double_values(self) -> list[ComputedDoubleValues]:
        """
                            List of double-valued quantities to output.
                            
                            Use :class:`~pylimer_tools_cpp.ComputedDoubleValues` enum to specify which floating-point quantities
                            should be computed and written to output.
        """
    @double_values.setter
    def double_values(self, arg0: collections.abc.Sequence[ComputedDoubleValues]) -> None:
        ...
    @property
    def filename(self) -> str:
        """
              The path and name of the file to write to.
              An empty string ("") means standard output (console).
        """
    @filename.setter
    def filename(self, arg0: str) -> None:
        ...
    @property
    def int_values(self) -> list[ComputedIntValues]:
        """
                            List of integer-valued quantities to output.
                            
                            Use :class:`~pylimer_tools_cpp.ComputedIntValues` enum to specify which integer quantities
                            should be computed and written to output.
        """
    @int_values.setter
    def int_values(self, arg0: collections.abc.Sequence[ComputedIntValues]) -> None:
        ...
    @property
    def output_every(self) -> int:
        """
             How often to write the values to the output.
             For averages, this value also says how many values will be averaged.
        """
    @output_every.setter
    def output_every(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def use_every(self) -> int:
        """
             For autocorrelation and averaging, how often to include values.
        
             Use a value of 1 to take average of or autocorrelate, respectively,
             all values encountered during the simulation or optimization procedure.
        """
    @use_every.setter
    def use_every(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
class ResistanceDistanceAnalysis:
    """
    
          Analysis of topological resistance distances, Kirchhoff index,
          and random-walk commute times on graph Laplacians.
        
    """
    @staticmethod
    def create_complete(num_nodes: typing.SupportsInt | typing.SupportsIndex) -> ResistanceDistanceAnalysis:
        """
                Create a ResistanceDistanceAnalysis instance for a complete graph K_N.
        
                :param num_nodes: Number of nodes in the complete graph (at least 2)
                :return: ResistanceDistanceAnalysis instance
        """
    @staticmethod
    def create_linear_chain(num_nodes: typing.SupportsInt | typing.SupportsIndex) -> ResistanceDistanceAnalysis:
        """
                Create a ResistanceDistanceAnalysis instance for a linear chain of num_nodes nodes.
        
                :param num_nodes: Number of nodes in the linear chain
                :return: ResistanceDistanceAnalysis instance
        """
    @staticmethod
    def create_ring(num_nodes: typing.SupportsInt | typing.SupportsIndex) -> ResistanceDistanceAnalysis:
        """
                Create a ResistanceDistanceAnalysis instance for a ring graph of num_nodes nodes.
        
                :param num_nodes: Number of nodes in the ring (at least 3)
                :return: ResistanceDistanceAnalysis instance
        """
    @staticmethod
    def create_star(num_arms: typing.SupportsInt | typing.SupportsIndex, arm_length: typing.SupportsInt | typing.SupportsIndex = 1) -> ResistanceDistanceAnalysis:
        """
                Create a ResistanceDistanceAnalysis instance for a star graph.
        
                :param num_arms: Number of arms (at least 2)
                :param arm_length: Number of edges per arm (defaults to 1)
                :return: ResistanceDistanceAnalysis instance
        """
    @staticmethod
    def from_matrix(L: scipy.sparse.csc_matrix) -> ResistanceDistanceAnalysis:
        """
                Construct a ResistanceDistanceAnalysis instance from a sparse graph Laplacian matrix L.
        
                :param L: Symmetric positive semidefinite graph Laplacian matrix
                :return: ResistanceDistanceAnalysis instance
        """
    @staticmethod
    def from_universe(universe: Universe) -> ResistanceDistanceAnalysis:
        """
                Construct a ResistanceDistanceAnalysis instance from a Universe's atom-bond graph.
        
                :param universe: Universe instance
                :return: ResistanceDistanceAnalysis instance
        """
    @typing.overload
    def __init__(self, spring_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], spring_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], weights: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex] = [], num_nodes: typing.SupportsInt | typing.SupportsIndex = 0) -> None:
        """
                   Construct from edge lists (spring_from, spring_to), optional edge weights,
                   and optional number of nodes.
        
                   :param spring_from: Source node indices of edges
                   :param spring_to: Target node indices of edges
                   :param weights: Optional weights (spring constants) of edges (defaults to 1.0)
                   :param num_nodes: Optional total number of nodes (deduced if 0)
        """
    @typing.overload
    def __init__(self, L: scipy.sparse.csc_matrix) -> None:
        """
                   Construct from a pre-assembled sparse scalar graph Laplacian matrix L.
        
                   :param L: Symmetric positive semidefinite graph Laplacian matrix
        """
    @typing.overload
    def __init__(self, universe: Universe) -> None:
        """
                   Construct from a Universe's atom-bond graph.
        
                   :param universe: Universe instance
        """
    def compute_all_pairs_commute_time(self, chunk_size: typing.SupportsInt | typing.SupportsIndex = 128) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, n]"]:
        """
                   Compute the N x N matrix of all-pairs random-walk commute times.
        
                   :param chunk_size: Number of columns per solve chunk (defaults to 128)
                   :return: N x N matrix of commute times
        """
    def compute_all_pairs_effective_resistance(self, chunk_size: typing.SupportsInt | typing.SupportsIndex = 128) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, n]"]:
        """
                   Compute the N x N symmetric matrix of all-pairs effective resistances.
        
                   :param chunk_size: Number of columns per solve chunk (defaults to 128)
                   :return: N x N matrix of effective resistances
        """
    def compute_commute_time(self, i: typing.SupportsInt | typing.SupportsIndex, j: typing.SupportsInt | typing.SupportsIndex) -> float:
        """
                   Compute the random-walk commute time C_ij = 2 * M * R_ij between nodes i and j.
        
                   :param i: First node index
                   :param j: Second node index
                   :return: Commute time C_ij
        """
    def compute_effective_resistance(self, i: typing.SupportsInt | typing.SupportsIndex, j: typing.SupportsInt | typing.SupportsIndex) -> float:
        """
                   Compute the effective resistance between nodes i and j.
        
                   :param i: First node index
                   :param j: Second node index
                   :return: Effective resistance R_ij (or inf if in different components)
        """
    def compute_effective_resistance_batched(self, pairs: collections.abc.Sequence[tuple[typing.SupportsInt | typing.SupportsIndex, typing.SupportsInt | typing.SupportsIndex]], chunk_size: typing.SupportsInt | typing.SupportsIndex = 128) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                   Compute the effective resistance for a batch of node pairs using chunked multi-RHS solves.
        
                   :param pairs: List of (i, j) node index pairs
                   :param chunk_size: Number of columns per solve chunk (defaults to 128)
                   :return: Vector of effective resistances for each pair
        """
    def compute_kirchhoff_index(self, chunk_size: typing.SupportsInt | typing.SupportsIndex = 128) -> float:
        """
                Compute the Kirchhoff index K(G) = sum_{i < j} R_ij without allocating the full N x N matrix.
        
                :param chunk_size: Number of columns per solve chunk (defaults to 128)
                :return: Kirchhoff index (or inf if disconnected)
        """
    def compute_kirchhoff_index_spectral(self) -> float:
        """
                   Compute the Kirchhoff index from the non-zero eigenvalues of the Laplacian.
        
                   :return: Kirchhoff index (or inf if disconnected)
        """
    def compute_moore_penrose_pseudo_inverse(self, chunk_size: typing.SupportsInt | typing.SupportsIndex = 128) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, n]"]:
        """
                   Compute the Moore-Penrose pseudo-inverse L^+ of the graph Laplacian.
        
                   :param chunk_size: Number of columns per solve chunk (defaults to 128)
                   :return: N x N Moore-Penrose pseudo-inverse matrix
        """
    def compute_resistance_centrality(self, chunk_size: typing.SupportsInt | typing.SupportsIndex = 128) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                   Compute the resistance centrality (N - 1) / sum_{j != i} R_ij for each node.
        
                   :param chunk_size: Number of columns per solve chunk (defaults to 128)
                   :return: Vector of resistance centralities
        """
    def compute_total_resistance_per_node(self, chunk_size: typing.SupportsInt | typing.SupportsIndex = 128) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                   Compute the sum of effective resistances from each node to all other nodes.
        
                   :param chunk_size: Number of columns per solve chunk (defaults to 128)
                   :return: Vector of total resistance per node
        """
    def get_connected_components(self) -> list[list[int]]:
        """
                Get the list of node indices for each connected component.
        
                :return: List of node index lists per connected component
        """
    def get_laplacian(self) -> scipy.sparse.csc_matrix:
        """
                   Get the sparse graph Laplacian matrix L.
        
                   :return: Sparse Laplacian matrix
        """
    def get_nr_of_connected_components(self) -> int:
        """
                   Get the number of connected components in the graph.
        
                   :return: Number of connected components
        """
    def get_nr_of_edges(self) -> int:
        """
                   Get the number of non-self-loop edges in the graph.
        
                   :return: Number of edges
        """
    def get_nr_of_nodes(self) -> int:
        """
                   Get the number of nodes in the graph.
        
                   :return: Number of nodes
        """
    def get_pinned_laplacian(self) -> scipy.sparse.csc_matrix:
        """
                   Get the sparse pinned Laplacian matrix with reference nodes pinned to 1 on the diagonal.
        
                   :return: Sparse pinned Laplacian matrix
        """
    def get_total_edge_weight(self) -> float:
        """
                   Get the total weight of all edges in the graph.
        
                   :return: Total edge weight
        """
    def is_connected(self) -> bool:
        """
                   Check whether the graph consists of a single connected component.
        
                   :return: True if connected, False otherwise
        """
    def solve_identity_chunked(self, chunk_size: typing.SupportsInt | typing.SupportsIndex = 128) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, n]"]:
        """
                   Solve L_pinned * X = I in multi-column chunks.
        
                   :param chunk_size: Number of columns per solve chunk (defaults to 128)
                   :return: N x N solution matrix X
        """
class RouseBrownianSimulator:
    """
    
          Implicit time-stepping Rouse and Brownian dynamics simulator for polymer chains
          and general spring networks.
        
    """
    @staticmethod
    def compute_center_of_mass_msd(trajectory: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"]], time_averaged: bool = False) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Compute center-of-mass mean-squared displacement from a trajectory of (N, 3) position matrices.
        """
    @staticmethod
    def compute_end_to_end_autocorrelation(trajectory: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"]], node_a: typing.SupportsInt | typing.SupportsIndex = 0, node_b: typing.SupportsInt | typing.SupportsIndex | None = None, time_averaged: bool = False) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Compute normalized end-to-end vector autocorrelation from a trajectory of (N, 3) position matrices (node_b=None defaults to N - 1).
        """
    @staticmethod
    def compute_msd(trajectory: collections.abc.Sequence[typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"]], time_averaged: bool = False) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Compute bead mean-squared displacement from a trajectory of (N, 3) position matrices.
        """
    @staticmethod
    def create_complete(num_nodes: typing.SupportsInt | typing.SupportsIndex, gamma: typing.SupportsFloat | typing.SupportsIndex = 1.0, k: typing.SupportsFloat | typing.SupportsIndex = 1.0, T: typing.SupportsFloat | typing.SupportsIndex = 1.0, dt: typing.SupportsFloat | typing.SupportsIndex = 0.01) -> RouseBrownianSimulator:
        """
                Create a RouseBrownianSimulator for a complete graph K_N of num_nodes beads.
        
                :param num_nodes: Number of nodes in the complete graph (at least 2)
                :param gamma: Bead friction coefficient
                :param k: Hookean spring constant
                :param T: Temperature
                :param dt: Integration time step
                :return: RouseBrownianSimulator instance
        """
    @staticmethod
    def create_linear_chain(num_nodes: typing.SupportsInt | typing.SupportsIndex, gamma: typing.SupportsFloat | typing.SupportsIndex = 1.0, k: typing.SupportsFloat | typing.SupportsIndex = 1.0, T: typing.SupportsFloat | typing.SupportsIndex = 1.0, dt: typing.SupportsFloat | typing.SupportsIndex = 0.01) -> RouseBrownianSimulator:
        """
                Create a RouseBrownianSimulator for a linear chain of num_nodes beads.
        """
    @staticmethod
    def create_ring(num_nodes: typing.SupportsInt | typing.SupportsIndex, gamma: typing.SupportsFloat | typing.SupportsIndex = 1.0, k: typing.SupportsFloat | typing.SupportsIndex = 1.0, T: typing.SupportsFloat | typing.SupportsIndex = 1.0, dt: typing.SupportsFloat | typing.SupportsIndex = 0.01) -> RouseBrownianSimulator:
        """
                Create a RouseBrownianSimulator for a ring polymer of num_nodes beads.
        """
    @staticmethod
    def create_star(num_arms: typing.SupportsInt | typing.SupportsIndex, arm_length: typing.SupportsInt | typing.SupportsIndex = 1, gamma: typing.SupportsFloat | typing.SupportsIndex = 1.0, k: typing.SupportsFloat | typing.SupportsIndex = 1.0, T: typing.SupportsFloat | typing.SupportsIndex = 1.0, dt: typing.SupportsFloat | typing.SupportsIndex = 0.01) -> RouseBrownianSimulator:
        """
                Create a RouseBrownianSimulator for a star polymer with num_arms arms of arm_length bonds each.
        """
    @staticmethod
    def from_matrix(L: scipy.sparse.csc_matrix, gamma: typing.SupportsFloat | typing.SupportsIndex = 1.0, k: typing.SupportsFloat | typing.SupportsIndex = 1.0, T: typing.SupportsFloat | typing.SupportsIndex = 1.0, dt: typing.SupportsFloat | typing.SupportsIndex = 0.01) -> RouseBrownianSimulator:
        """
                Create a RouseBrownianSimulator from a sparse graph Laplacian matrix L.
        
                :param L: Symmetric positive semidefinite graph Laplacian matrix
                :param gamma: Bead friction coefficient
                :param k: Hookean spring constant
                :param T: Temperature
                :param dt: Integration time step
                :return: RouseBrownianSimulator instance
        """
    @staticmethod
    def from_universe(universe: Universe, gamma: typing.SupportsFloat | typing.SupportsIndex = 1.0, k: typing.SupportsFloat | typing.SupportsIndex = 1.0, T: typing.SupportsFloat | typing.SupportsIndex = 1.0, dt: typing.SupportsFloat | typing.SupportsIndex = 0.01) -> RouseBrownianSimulator:
        """
                Create a RouseBrownianSimulator from a Universe's atom-bond graph and coordinates.
        
                :param universe: Universe instance
                :param gamma: Bead friction coefficient
                :param k: Hookean spring constant
                :param T: Temperature
                :param dt: Integration time step
                :return: RouseBrownianSimulator instance
        """
    @typing.overload
    def __init__(self, spring_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], spring_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], gamma: typing.SupportsFloat | typing.SupportsIndex = 1.0, k: typing.SupportsFloat | typing.SupportsIndex = 1.0, T: typing.SupportsFloat | typing.SupportsIndex = 1.0, dt: typing.SupportsFloat | typing.SupportsIndex = 0.01, num_nodes: typing.SupportsInt | typing.SupportsIndex = 0, weights: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex] = []) -> None:
        """
                   Construct from edge lists (spring_from, spring_to), friction gamma,
                   spring constant k, temperature T, time step dt, optional num_nodes,
                   and optional edge weights.
        
                   :param spring_from: Source node indices of springs
                   :param spring_to: Target node indices of springs
                   :param gamma: Bead friction coefficient
                   :param k: Hookean spring constant
                   :param T: Temperature
                   :param dt: Integration time step
                   :param num_nodes: Total number of nodes (deduced if 0)
                   :param weights: Optional positive edge weights (defaults to 1.0)
        """
    @typing.overload
    def __init__(self, L: scipy.sparse.csc_matrix, gamma: typing.SupportsFloat | typing.SupportsIndex = 1.0, k: typing.SupportsFloat | typing.SupportsIndex = 1.0, T: typing.SupportsFloat | typing.SupportsIndex = 1.0, dt: typing.SupportsFloat | typing.SupportsIndex = 0.01) -> None:
        """
                   Construct from a pre-assembled sparse scalar graph Laplacian matrix L.
        
                   :param L: Symmetric positive semidefinite graph Laplacian matrix
                   :param gamma: Bead friction coefficient
                   :param k: Hookean spring constant
                   :param T: Temperature
                   :param dt: Integration time step
        """
    @typing.overload
    def __init__(self, universe: Universe, gamma: typing.SupportsFloat | typing.SupportsIndex = 1.0, k: typing.SupportsFloat | typing.SupportsIndex = 1.0, T: typing.SupportsFloat | typing.SupportsIndex = 1.0, dt: typing.SupportsFloat | typing.SupportsIndex = 0.01) -> None:
        """
                   Construct from a Universe's atom-bond graph and initial coordinates.
        
                   :param universe: Universe instance
                   :param gamma: Bead friction coefficient
                   :param k: Hookean spring constant
                   :param T: Temperature
                   :param dt: Integration time step
        """
    def get_center_of_mass(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
        Compute the 3D center of mass of all nodes.
        """
    def get_diffusion_coefficient(self) -> float:
        """
        Get the bead diffusion coefficient D = kB * T / gamma.
        """
    def get_dt(self) -> float:
        """
        Get the integration time step dt.
        """
    def get_end_to_end_distance(self, node_a: typing.SupportsInt | typing.SupportsIndex = 0, node_b: typing.SupportsInt | typing.SupportsIndex | None = None) -> float:
        """
        Compute the Euclidean distance between node_a and node_b (None defaults to N - 1).
        """
    def get_end_to_end_vector(self, node_a: typing.SupportsInt | typing.SupportsIndex = 0, node_b: typing.SupportsInt | typing.SupportsIndex | None = None) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
        Compute the 3D vector from node_a to node_b (None defaults to N - 1).
        """
    def get_gamma(self) -> float:
        """
        Get the friction coefficient gamma.
        """
    def get_k(self) -> float:
        """
        Get the Hookean spring constant k.
        """
    def get_kb(self) -> float:
        """
        Get the Boltzmann constant kB.
        """
    def get_laplacian(self) -> scipy.sparse.csc_matrix:
        """
        Get the sparse graph Laplacian matrix L.
        """
    def get_nr_of_nodes(self) -> int:
        """
        Get the number of nodes (beads) in the simulator.
        """
    def get_positions(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, n]"]:
        """
        Get the (N, 3) bead positions matrix.
        """
    def get_positions_flattened(self, interleaved: bool = True) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Get bead positions as a flattened 3N vector.
        """
    def get_potential_energy(self) -> float:
        """
        Compute the harmonic spring potential energy of the current configuration.
        """
    def get_radius_of_gyration(self) -> float:
        """
        Compute the radius of gyration of the current bead configuration.
        """
    def get_shifted_matrix(self) -> scipy.sparse.csc_matrix:
        """
        Get the shifted implicit stepping matrix A = I + dt * (k / gamma) * L.
        """
    def get_solver_type(self) -> RouseSolverType:
        """
        Get the active linear solver type.
        """
    def get_temperature(self) -> float:
        """
        Get the temperature T.
        """
    def initialize_linear_positions(self, bond_length: typing.SupportsFloat | typing.SupportsIndex = 1.0, axis: typing.SupportsInt | typing.SupportsIndex = 0) -> None:
        """
        Initialize bead positions spaced along a coordinate axis.
        """
    def initialize_random_walk(self, bond_length: typing.SupportsFloat | typing.SupportsIndex = 1.0, seed: typing.SupportsInt | typing.SupportsIndex = 42) -> None:
        """
        Initialize bead positions as a 3D random walk.
        """
    def run(self, total_steps: typing.SupportsInt | typing.SupportsIndex, sample_interval: typing.SupportsInt | typing.SupportsIndex = 1, save_trajectory: bool = False, end_to_end_node_a: typing.SupportsInt | typing.SupportsIndex = 0, end_to_end_node_b: typing.SupportsInt | typing.SupportsIndex | None = None) -> RouseSimulationResult:
        """
        Run the simulation for total_steps and record dynamical observables every sample_interval steps (end_to_end_node_b=None defaults to N - 1).
        """
    def set_dt(self, dt: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Set the integration time step dt.
        """
    def set_gamma(self, gamma: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Set the friction coefficient gamma.
        """
    def set_k(self, k: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Set the Hookean spring constant k.
        """
    def set_kb(self, kb: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Set the Boltzmann constant kB.
        """
    def set_parameters(self, gamma: typing.SupportsFloat | typing.SupportsIndex, k: typing.SupportsFloat | typing.SupportsIndex, T: typing.SupportsFloat | typing.SupportsIndex, dt: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Set gamma, k, T, and dt simultaneously.
        """
    def set_positions(self, positions: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"]) -> None:
        """
        Set the (N, 3) bead positions matrix.
        """
    def set_positions_flattened(self, positions: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], interleaved: bool = True) -> None:
        """
        Set bead positions from a flattened 3N vector.
        """
    def set_seed(self, seed: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Set the random number generator seed for thermal noise.
        """
    def set_solver_type(self, solver_type: RouseSolverType) -> None:
        """
        Set the linear solver algorithm used for implicit time-stepping.
        """
    def set_temperature(self, T: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Set the temperature T.
        """
    def step(self, num_steps: typing.SupportsInt | typing.SupportsIndex = 1) -> None:
        """
        Advance the simulation by num_steps implicit Brownian steps.
        """
    def step_deterministic(self, num_steps: typing.SupportsInt | typing.SupportsIndex = 1) -> None:
        """
        Advance the simulation by num_steps deterministic relaxation steps (without thermal noise).
        """
    def step_deterministic_with_external_forces(self, forces: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"], num_steps: typing.SupportsInt | typing.SupportsIndex = 1) -> None:
        """
        Advance the simulation by num_steps deterministic steps under external forces.
        """
    def step_with_external_forces(self, forces: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"], num_steps: typing.SupportsInt | typing.SupportsIndex = 1) -> None:
        """
        Advance the simulation by num_steps under external forces and thermal noise.
        """
class RouseSimulationResult:
    """
    
          Output container for Rouse dynamics simulation runs.
        
    """
    @property
    def end_to_end_correlation(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Normalized end-to-end vector autocorrelation at sampled times.
        """
    @property
    def msd(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Bead mean-squared displacement at sampled times.
        """
    @property
    def msd_center_of_mass(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Center-of-mass mean-squared displacement at sampled times.
        """
    @property
    def times(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
        Sampled simulation times.
        """
    @property
    def trajectory(self) -> list[typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, n]"]]:
        """
        List of (N, 3) position matrices if save_trajectory was enabled.
        """
class RouseSolverType(enum.IntEnum):
    """
    Linear solver algorithm for implicit Rouse/Brownian time-stepping.
    """
    CONJUGATE_GRADIENT: typing.ClassVar[RouseSolverType]  # value = <RouseSolverType.CONJUGATE_GRADIENT: 2>
    DIRECT_LDLT: typing.ClassVar[RouseSolverType]  # value = <RouseSolverType.DIRECT_LDLT: 0>
    DIRECT_LLT: typing.ClassVar[RouseSolverType]  # value = <RouseSolverType.DIRECT_LLT: 1>
    @classmethod
    def __new__(cls, value):
        ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """
class SLESolver(enum.Enum):
    """
    Solver for sparse linear equation systems
    """
    ACCELERATE_LDLT: typing.ClassVar[SLESolver]  # value = <SLESolver.ACCELERATE_LDLT: 7>
    ACCELERATE_LLT: typing.ClassVar[SLESolver]  # value = <SLESolver.ACCELERATE_LLT: 6>
    AUTO: typing.ClassVar[SLESolver]  # value = <SLESolver.AUTO: 1>
    BICGSTAB: typing.ClassVar[SLESolver]  # value = <SLESolver.BICGSTAB: 15>
    BICGSTAB_DIAGONALIZED: typing.ClassVar[SLESolver]  # value = <SLESolver.BICGSTAB_DIAGONALIZED: 16>
    BICGSTAB_IDENTITY: typing.ClassVar[SLESolver]  # value = <SLESolver.BICGSTAB_IDENTITY: 17>
    BICGSTAB_INCOMPLETE_LU: typing.ClassVar[SLESolver]  # value = <SLESolver.BICGSTAB_INCOMPLETE_LU: 18>
    CONJUGATE_GRADIENT: typing.ClassVar[SLESolver]  # value = <SLESolver.CONJUGATE_GRADIENT: 8>
    CONJUGATE_GRADIENT_DIAGONALIZED: typing.ClassVar[SLESolver]  # value = <SLESolver.CONJUGATE_GRADIENT_DIAGONALIZED: 9>
    CONJUGATE_GRADIENT_IDENTITY: typing.ClassVar[SLESolver]  # value = <SLESolver.CONJUGATE_GRADIENT_IDENTITY: 10>
    CONJUGATE_GRADIENT_INCOMPLETE_CHOLESKY: typing.ClassVar[SLESolver]  # value = <SLESolver.CONJUGATE_GRADIENT_INCOMPLETE_CHOLESKY: 11>
    DEFAULT: typing.ClassVar[SLESolver]  # value = <SLESolver.DEFAULT: 0>
    GRADIENT_DESCENT: typing.ClassVar[SLESolver]  # value = <SLESolver.GRADIENT_DESCENT: 19>
    GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG: typing.ClassVar[SLESolver]  # value = <SLESolver.GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG: 21>
    GRADIENT_DESCENT_BARZILAI_BORWEIN_MOMENTUM: typing.ClassVar[SLESolver]  # value = <SLESolver.GRADIENT_DESCENT_BARZILAI_BORWEIN_MOMENTUM: 22>
    GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT: typing.ClassVar[SLESolver]  # value = <SLESolver.GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT: 20>
    LEAST_SQUARES_CONJUGATE_GRADIENT: typing.ClassVar[SLESolver]  # value = <SLESolver.LEAST_SQUARES_CONJUGATE_GRADIENT: 12>
    LEAST_SQUARES_CONJUGATE_GRADIENT_DIAGONALIZED: typing.ClassVar[SLESolver]  # value = <SLESolver.LEAST_SQUARES_CONJUGATE_GRADIENT_DIAGONALIZED: 13>
    LEAST_SQUARES_CONJUGATE_GRADIENT_IDENTITY: typing.ClassVar[SLESolver]  # value = <SLESolver.LEAST_SQUARES_CONJUGATE_GRADIENT_IDENTITY: 14>
    METAL_BARZILAI_BORWEIN: typing.ClassVar[SLESolver]  # value = <SLESolver.METAL_BARZILAI_BORWEIN: 24>
    METAL_CG: typing.ClassVar[SLESolver]  # value = <SLESolver.METAL_CG: 23>
    SIMPLICIAL_LDLT: typing.ClassVar[SLESolver]  # value = <SLESolver.SIMPLICIAL_LDLT: 3>
    SIMPLICIAL_LLT: typing.ClassVar[SLESolver]  # value = <SLESolver.SIMPLICIAL_LLT: 2>
    SPARSE_LU: typing.ClassVar[SLESolver]  # value = <SLESolver.SPARSE_LU: 4>
    SPARSE_QR: typing.ClassVar[SLESolver]  # value = <SLESolver.SPARSE_QR: 5>
class SLQConfig:
    """
    Configuration for Stochastic Lanczos Quadrature (SLQ).
    """
    backend: ComputeBackend
    def __init__(self, num_steps: typing.SupportsInt | typing.SupportsIndex = 200, num_random_vectors: typing.SupportsInt | typing.SupportsIndex = 20, backend: ComputeBackend = ...) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def num_random_vectors(self) -> int:
        ...
    @num_random_vectors.setter
    def num_random_vectors(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_steps(self) -> int:
        ...
    @num_steps.setter
    def num_steps(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
class SlipLinkRelaxationSolver(MEHPSolver):
    """
    Iterative alternating coordinate and contour-partition relaxation solver for slip-link networks.
    """
    config: MEHPSlipLinkSolverConfig
    def __init__(self, config: MEHPSlipLinkSolverConfig = ...) -> None:
        ...
    def get_config(self) -> MEHPSlipLinkSolverConfig:
        ...
    def set_config(self, config: MEHPSlipLinkSolverConfig) -> None:
        ...
class SlipMode(enum.Enum):
    """
    Mobility of entanglement links along polymer strands.
    """
    FIXED: typing.ClassVar[SlipMode]  # value = <SlipMode.FIXED: 0>
    SLIDING: typing.ClassVar[SlipMode]  # value = <SlipMode.SLIDING: 1>
class StrandCategory(enum.IntEnum):
    """
    Topological category of a strand in a CrosslinkerUniverse.
    """
    ALL: typing.ClassVar[StrandCategory]  # value = <StrandCategory.ALL: 0>
    DANGLING_CHAIN: typing.ClassVar[StrandCategory]  # value = <StrandCategory.DANGLING_CHAIN: 3>
    FREE_CHAIN: typing.ClassVar[StrandCategory]  # value = <StrandCategory.FREE_CHAIN: 4>
    NETWORK_STRAND: typing.ClassVar[StrandCategory]  # value = <StrandCategory.NETWORK_STRAND: 1>
    PRIMARY_LOOP: typing.ClassVar[StrandCategory]  # value = <StrandCategory.PRIMARY_LOOP: 5>
    TWO_CROSSLINKERS_INCLUDING_LOOPS: typing.ClassVar[StrandCategory]  # value = <StrandCategory.TWO_CROSSLINKERS_INCLUDING_LOOPS: 2>
    @classmethod
    def __new__(cls, value):
        ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """
class StrandPairEntanglements:
    """
    
          A struct to store sampled entanglements on strands (without requiring an explicit Universe of all beads).
        
    """
    def __init__(self) -> None:
        """
        Get an instance of this struct
        """
    def __len__(self) -> int:
        ...
    def size(self) -> int:
        ...
    def to_atom_pair_entanglements(self, crosslinker_universe: CrosslinkerUniverse) -> AtomPairEntanglements:
        """
              Convert the strand pair entanglements to AtomPairEntanglements using the deterministic atom indexing of CrosslinkerUniverse.
        """
    @property
    def alpha_1(self) -> list[float]:
        """
              Contour fraction along the first strand in (0, 1).
        """
    @alpha_1.setter
    def alpha_1(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def alpha_2(self) -> list[float]:
        """
              Contour fraction along the second strand in (0, 1).
        """
    @alpha_2.setter
    def alpha_2(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def bead_idx_1(self) -> list[int]:
        """
              Bead index along the first strand in [0, beads_in_strand - 1].
        """
    @bead_idx_1.setter
    def bead_idx_1(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def bead_idx_2(self) -> list[int]:
        """
              Bead index along the second strand in [0, beads_in_strand - 1].
        """
    @bead_idx_2.setter
    def bead_idx_2(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def position_1(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
              Sampled 3D positions of the first entangled beads.
        """
    @position_1.setter
    def position_1(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @property
    def position_2(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
              Sampled 3D positions of the second entangled beads.
        """
    @position_2.setter
    def position_2(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @property
    def strand_idx_1(self) -> list[int]:
        """
              First strand index of each entanglement pair.
        """
    @strand_idx_1.setter
    def strand_idx_1(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def strand_idx_2(self) -> list[int]:
        """
              Second strand index of each entanglement pair.
        """
    @strand_idx_2.setter
    def strand_idx_2(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def x(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
              X-coordinates of the sampled entanglement midpoints.
        """
    @x.setter
    def x(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @property
    def y(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
              Y-coordinates of the sampled entanglement midpoints.
        """
    @y.setter
    def y(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
    @property
    def z(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
              Z-coordinates of the sampled entanglement midpoints.
        """
    @z.setter
    def z(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        ...
class StructureSimplificationMode(enum.Enum):
    """
    How the structure shall be simplified during the optimization in order to remove non-trapped entanglement links.
    """
    INACTIVE_ONLY: typing.ClassVar[StructureSimplificationMode]  # value = <StructureSimplificationMode.INACTIVE_ONLY: 2>
    INACTIVE_THEN_X2F: typing.ClassVar[StructureSimplificationMode]  # value = <StructureSimplificationMode.INACTIVE_THEN_X2F: 3>
    NO_SIMPLIFICATION: typing.ClassVar[StructureSimplificationMode]  # value = <StructureSimplificationMode.NO_SIMPLIFICATION: 0>
    X1F_X2F_THEN_INACTIVE: typing.ClassVar[StructureSimplificationMode]  # value = <StructureSimplificationMode.X1F_X2F_THEN_INACTIVE: 4>
    X2F_ONLY: typing.ClassVar[StructureSimplificationMode]  # value = <StructureSimplificationMode.X2F_ONLY: 1>
class Universe:
    """
    
        Represents a full Polymer Network structure, a collection of molecules.
        
        This is the main class for representing molecular systems, containing
        atoms, bonds, angles, and the simulation box.
        
    """
    def __contains__(self, arg0: Atom) -> bool:
        """
                  Check whether a particular atom is contained in this universe.
        """
    def __copy__(self) -> Universe:
        """
                 Create a copy of this universe.
                 
                 :return: A new Universe instance that is a copy of this one
        """
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> Atom:
        """
               Access an atom by its vertex index.
        """
    def __getstate__(self) -> tuple[str]:
        ...
    def __init__(self, Lx: typing.SupportsFloat | typing.SupportsIndex, Ly: typing.SupportsFloat | typing.SupportsIndex, Lz: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
                 Instantiate this Universe (Collection of Molecules) providing the box lengths.
                 
                 :param Lx: Box length in x direction
                 :param Ly: Box length in y direction
                 :param Lz: Box length in z direction
        """
    def __len__(self) -> int:
        """
               Get the number of atoms in this universe.
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def add_angles(self, angles_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], angles_via: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], angles_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], angle_types: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        """
                 Add angles to the Universe. No relation to the underlying graph, 
                 just a method to preserve read & write capabilities.
                 
                 :param angles_from: Vector of atom IDs for angle start points
                 :param angles_via: Vector of atom IDs for angle middle points
                 :param angles_to: Vector of atom IDs for angle end points
                 :param angle_types: Vector of angle types
        """
    def add_atoms(self, ids: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], types: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], x: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], y: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], z: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], nx: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], ny: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], nz: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        """
                 Add atoms to the Universe, vertices to the underlying graph.
                 
                 :param ids: Vector of atom IDs
                 :param types: Vector of atom types
                 :param x: Vector of x coordinates
                 :param y: Vector of y coordinates
                 :param z: Vector of z coordinates
                 :param nx: Vector of periodic image flags in x direction
                 :param ny: Vector of periodic image flags in y direction
                 :param nz: Vector of periodic image flags in z direction
        """
    @typing.overload
    def add_bonds(self, bonds_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], bonds_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        """
                 Add bonds to the underlying atoms, edges to the underlying graph. 
                 If the connected atoms are not found, the bonds are silently skipped.
                 
                 :param bonds_from: Vector of atom IDs for bond start points
                 :param bonds_to: Vector of atom IDs for bond end points
        """
    @typing.overload
    def add_bonds(self, nr_of_bonds: typing.SupportsInt | typing.SupportsIndex, bonds_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], bonds_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], bond_types: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], ignore_non_existent_atoms: bool = False, simplify_universe: bool = True) -> None:
        """
                 Add bonds to the underlying atoms, edges to the underlying graph.
                 
                 :param nr_of_bonds: Number of bonds to add
                 :param bonds_from: Vector of atom IDs for bond start points
                 :param bonds_to: Vector of atom IDs for bond end points
                 :param bond_types: Vector of bond types
                 :param ignore_non_existent_atoms: Whether to skip bonds to non-existent atoms
                 :param simplify_universe: Whether to simplify the universe after adding bonds
        """
    def add_bonds_with_types(self, bonds_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], bonds_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], bond_types: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        """
                 Add bonds to the underlying atoms, edges to the underlying graph. 
                 If the connected atoms are not found, the bonds are silently skipped.
                 
                 :param bonds_from: Vector of atom IDs for bond start points
                 :param bonds_to: Vector of atom IDs for bond end points
                 :param bond_types: Vector of bond types
        """
    def add_dihedral_angles(self, angles_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], angles_via1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], angles_via2: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], angles_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], angle_types: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        """
                 Add dihedral angles to the Universe. No relation to the underlying graph, 
                 just a method to preserve read & write capabilities.
                 
                 :param angles_from: Vector of atom IDs for dihedral start points
                 :param angles_via1: Vector of atom IDs for first middle points
                 :param angles_via2: Vector of atom IDs for second middle points
                 :param angles_to: Vector of atom IDs for dihedral end points
                 :param angle_types: Vector of dihedral angle types
        """
    def clear_caches(self) -> None:
        """
                 Clear internal caches of decomposed molecules and chains with crosslinkers to reclaim memory.
        """
    def compute_angles(self) -> list[float]:
        """
                 Compute the angle of each angle in the molecule, respecting periodic boundaries.
                 
                 :return: A list of angle values in radians
        """
    def compute_bond_length_stats(self) -> BondLengthStats:
        """
                  Compute summary bond length statistics (count, mean, mean squared, min, max)
                  from the Universe's actual bond coordinates.
        
                  :return: BondLengthStats object.
        """
    def compute_bond_lengths(self) -> list[float]:
        """
                 Compute the length of each bond in the molecule, respecting periodic boundaries.
                 
                 :return: A list of bond lengths
        """
    def compute_bond_vectors(self) -> list[typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]]:
        """
                 Compute the vectors of each bond in the molecule, respecting periodic boundaries.
                 
                 :return: A list of bond vectors
        """
    def compute_dxs(self, atom_ids_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], atom_ids_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> list[float]:
        """
                 Compute the dx distance for certain bonds (length in x direction).
                 
                 :param atom_ids_to: Vector of destination atom IDs
                 :param atom_ids_from: Vector of source atom IDs
                 :return: Vector of x-direction distances
        """
    def compute_dys(self, atom_ids_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], atom_ids_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> list[float]:
        """
                 Compute the dy distance for certain bonds (length in y direction).
                 
                 :param atom_ids_to: Vector of destination atom IDs
                 :param atom_ids_from: Vector of source atom IDs
                 :return: Vector of y-direction distances
        """
    def compute_dzs(self, atom_ids_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], atom_ids_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> list[float]:
        """
                 Compute the dz distance for certain bonds (length in z direction).
                 
                 :param atom_ids_to: Vector of destination atom IDs
                 :param atom_ids_from: Vector of source atom IDs
                 :return: Vector of z-direction distances
        """
    def compute_end_to_end_distances(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex, derive_image_flags: bool = False) -> list[float]:
        """
                  Compute the end-to-end distance of each strand in the network.
        
                  .. note::
                       Internally, this uses either :meth:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance` 
                       or :meth:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags`, 
                       depending on `derive_image_flags`.
                       Invalid strands (where said function returns 0.0 or -1.0) are ignored.
                       
                  :param crosslinker_type: The type of crosslinker atoms
                  :param derive_image_flags: Whether to derive image flags from connectivity
                  :return: List of end-to-end distances
        """
    def compute_end_to_end_stats(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, category: StrandCategory = ..., derive_image_flags: bool = True) -> EndToEndStats:
        """
                  Compute summary end-to-end distance statistics (count, excluded zero-length count,
                  mean, mean squared) for strands matching ``category`` from the Universe's actual coordinates.
        
                  :param crosslinker_type: The type of crosslinker atoms (default: 2)
                  :param category: StrandCategory filter (default: NETWORK_STRAND)
                  :param derive_image_flags: Whether to derive image flags from connectivity (default: True)
                  :return: EndToEndStats object.
        """
    def compute_mean_end_to_end_distance(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex, derive_image_flags: bool = False) -> float:
        """
                  Compute the mean of the end-to-end distances of each strand in the network.
        
                  .. note::
                       Internally, this uses either :meth:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance` 
                       or :meth:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags`, 
                       depending on `derive_image_flags`.
                       Invalid strands (where said function returns 0.0 or -1.0) are ignored.
                       
                  :param crosslinker_type: The type of crosslinker atoms
                  :param derive_image_flags: Whether to derive image flags from connectivity
                  :return: Mean end-to-end distance
        """
    def compute_mean_squared_end_to_end_distance(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex, only_those_with_two_crosslinkers: bool = False, derive_image_flags: bool = False) -> float:
        """
                  Compute the mean square of the end-to-end distances of each strand (incl. crosslinks) in the network.
        
                  .. note::
                       Internally, this uses either :meth:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance` 
                       or :meth:`~pylimer_tools_cpp.Molecule.compute_end_to_end_distance_with_derived_image_flags`, 
                       depending on `derive_image_flags`.
                       Invalid strands (where said function returns 0.0 or -1.0) are ignored.
                       
                  :param crosslinker_type: The type of crosslinker atoms
                  :param only_those_with_two_crosslinkers: Whether to only consider strands with two crosslinkers
                  :param derive_image_flags: Whether to derive image flags from connectivity
                  :return: Mean squared end-to-end distance
        """
    def compute_mean_strand_length(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex) -> float:
        """
                      Compute the mean number of beads per strand.
                      
                      :param crosslinker_type: The type of crosslinker atoms
                      :return: Mean strand length
        """
    def compute_number_average_molecular_weight(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex) -> float:
        """
                      Compute the number average molecular weight.
        
                      .. note:: 
                            Crosslinkers are ignored completely.
                            
                      :param crosslinker_type: The type of crosslinker atoms
                      :return: Number average molecular weight
        """
    def compute_polydispersity_index(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex) -> float:
        """
                      Compute the polydispersity index: 
                      the weight average molecular weight over the number average molecular weight.
        """
    def compute_temperature(self, dimensions: typing.SupportsInt | typing.SupportsIndex = 3, k_b: typing.SupportsFloat | typing.SupportsIndex = 1.0) -> float:
        """
              Use the velocities per atom to compute the temperature from the kinetic energy of the system.
              
              :param dimensions: Number of dimensions (typically 3)
              :param k_b: Boltzmann constant value
              :return: Computed temperature
        """
    def compute_total_mass(self) -> float:
        """
                  Compute the total mass of this network/universe in whatever mass unit was used when 
                  :meth:`~pylimer_tools_cpp.Universe.set_masses` was called.
                  
                  :return: Total mass of the universe
        """
    def compute_weight_average_molecular_weight(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex) -> float:
        """
                      Compute the weight average molecular weight.
        
                      .. note:: 
                            Crosslinkers are ignored completely.
                            
                      :param crosslinker_type: The type of crosslinker atoms
                      :return: Weight average molecular weight
        """
    def compute_weight_fractions(self) -> dict[int, float]:
        """
                    Compute the weight fractions of each atom type in the network.
        
                    If no masses are stored, assumes a mass of 1 for each atom.
        
                    If the total mass is 0., returns the total mass per atom type.
        """
    def contract_vertices_along_bond_type(self, bond_type: typing.SupportsInt | typing.SupportsIndex) -> Universe:
        """
                  Merge vertices along a specific bond type.
        
                  May result in new self-loops; use :meth:`~pylimer_tools_cpp.Universe.simplify` to remove them.
                  
                  :param bond_type: The bond type to contract along
        """
    def count_atom_types(self) -> dict[int, int]:
        """
                  Count how often each atom type is present.
                  
                  :return: Dictionary mapping atom types to their counts
        """
    def count_atoms_in_skin_distance(self, distances: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], unwrapped: bool = False) -> list[int]:
        """
                  This is a function that may help you to compute the radial distribution function.
                  It loops through all atoms and counts neighbors within specified distance ranges.
        
                  :param distances: The edges of the bins for distance counting
                  :param unwrapped: Whether to measure the distance in unwrapped coordinates or as PBC-corrected distance
                  :return: Array of counts for each distance bin
        """
    def count_loop_lengths(self, max_length: typing.SupportsInt | typing.SupportsIndex = -1) -> dict[int, int]:
        """
                  Find all loops (below a specific length) and count the number of atoms involved in them.
                  Returns the count, how many loops per length are found.
        """
    def detect_angles(self) -> dict[str, list[int]]:
        """
                  Detect angles in the network based on the current bonds.
                  Return the result in the same format as :meth:`~pylimer_tools_cpp.Universe.get_angles`, 
                  but all angles that are detected in the network, rather than the ones already set.
                  Note that the angle types are determined by 
                  :meth:`~pylimer_tools_cpp.Universe.hash_angle_type`,
                  which serves angle types that should be mapped by you back to smaller numbers, 
                  before serving them again to :meth:`~pylimer_tools_cpp.Universe.add_angles`,
                  if you want to have them written e.g. for LAMMPS.
                  
                  :return: Dictionary with detected angle information
        """
    def detect_dihedral_angles(self) -> dict[str, list[int]]:
        """
                  Detect dihedral angles in the network based on the current bonds.
                  Return the result in the same format as :meth:`~pylimer_tools_cpp.Universe.get_dihedral_angles`, 
                  but all dihedral angles that are detected in the network, rather than the ones already set.
                  Note that the angle types are determined by 
                  :meth:`~pylimer_tools_cpp.Universe.hash_dihedral_angle_type`,
                  which serves angle types that should be mapped by you back to smaller numbers, 
                  before serving them to :meth:`~pylimer_tools_cpp.Universe.add_dihedral_angles`,
                  if you want to have them written e.g. for LAMMPS.
        
                  :return: Dictionary with detected dihedral angle information
        """
    def determine_effective_functionality_per_type(self) -> dict[int, float]:
        """
                    Find the average functionality of each atom type in the network.
                    
                    :return: Dictionary mapping atom types to average functionality
        """
    def determine_functionality_per_type(self) -> dict[int, int]:
        """
                    Find the maximum functionality of each atom type in the network.
                    
                    :return: Dictionary mapping atom types to maximum functionality
        """
    def find_loops(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex, max_length: typing.SupportsInt | typing.SupportsIndex = -1, skip_self_loops: bool = False) -> dict[int, list[list[Atom]]]:
        """
                    Decompose the Universe into loops.
                    The primary index specifies the degree of the loop.
        
                    CAUTION:
                       There are exponentially many paths between two crosslinkers of a network,
                       and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
                       You can use the `max_length` parameter to restrict the algorithm to only search for loops up to a certain length.
                       Use a negative value to find all loops and paths.
        """
    def find_minimal_order_loop_from(self, loop_start: typing.SupportsInt | typing.SupportsIndex, loop_step1: typing.SupportsInt | typing.SupportsIndex, max_length: typing.SupportsInt | typing.SupportsIndex = -1, skip_self_loops: bool = False) -> list[Atom]:
        """
                    Find the loops in the network starting with one connection.
        
                    .. warning::
                       There are exponentially many paths between two crosslinkers of a network,
                       and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
                       You can use the `max_length` parameter to restrict the algorithm to only search for loops up to a certain length.
                       Use a negative value to find all loops and paths.
        
                    :param loop_start: The atom ID to start the search from
                    :param loop_step1: The first step to take, i.e., the first bond to follow
                    :param max_length: The maximum length of the loop to find, or -1 for no limit
                    :param skip_self_loops: Whether to skip self-loops (i.e., loops that start and end at the same atom; only relevant if `loop_start` is equal to `loop_step1`)
                    :return: List of loops found
        """
    def get_angles(self) -> dict[str, list[int]]:
        """
                   Get all angles added to this network.
        
                   Returns a dict with three properties: 'angle_from', 'angle_via' and 'angle_to'.
        
                   .. note::
                       The integer values returned refer to the atom IDs, not the vertex IDs.
                       Use :meth:`~pylimer_tools_cpp.Universe.get_idx_by_atom_id` to translate them to vertex IDs.
                       
                   :return: Dictionary with angle information
        """
    def get_atom(self, atom_id: typing.SupportsInt | typing.SupportsIndex) -> Atom:
        """
                 Find an atom by its ID.
                 
                 :param atom_id: The ID of the atom to find
                 :return: The atom with the specified ID
        """
    def get_atom_by_vertex_id(self, vertex_id: typing.SupportsInt | typing.SupportsIndex) -> Atom:
        """
              Find an atom by the ID of the vertex of the underlying graph.
              
              :param vertex_id: The vertex ID to query
              :return: The atom at the specified vertex
        """
    def get_atom_id_by_vertex_idx(self, vertex_id: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
                 Get the ID of the atom by the vertex ID of the underlying graph.
                 
                 :param vertex_id: The vertex index in the underlying graph
                 :return: The atom ID corresponding to the vertex
        """
    def get_atom_types(self) -> list[int]:
        """
                  Get all types (each one for each atom) ordered by atom vertex ID.
                  
                  :return: Vector of atom types in vertex order
        """
    def get_atoms(self) -> list[Atom]:
        """
                    Get all atoms in the universe.
                    
                    :return: List of all atoms
        """
    def get_atoms_by_degree(self, functionality: typing.SupportsInt | typing.SupportsIndex) -> list[Atom]:
        """
                    Get the atoms that have the specified number of bonds.
        """
    def get_atoms_by_type(self, atom_type: typing.SupportsInt | typing.SupportsIndex) -> list[Atom]:
        """
                    Query all atoms by their type.
        """
    def get_atoms_connected_to(self, atom: Atom) -> list[Atom]:
        """
                    Get the atoms connected to a specified atom.
        
                    Internally uses :meth:`~pylimer_tools_cpp.Universe.get_atoms_connected_to`.
                    
                    :param atom: The atom to query connections for
                    :return: List of connected atoms
        """
    def get_atoms_connected_to_vertex(self, vertex_idx: typing.SupportsInt | typing.SupportsIndex) -> list[Atom]:
        """
                    Get the atoms connected to a specified vertex ID.
                    
                    :param vertex_idx: The vertex index to query
                    :return: List of connected atoms
        """
    def get_bonds(self) -> dict[str, list[int]]:
        """
                    Get all bonds. Returns a dict with three properties: 'bond_from', 'bond_to' and 'bond_type'.
                    The order is not necessarily related to any structural characteristic.
                    
                    :return: Dictionary with bond information
        """
    def get_box(self) -> Box:
        """
                    Get the underlying bounding box object.
                    
                    :return: The simulation box
        """
    def get_chains_with_crosslinker(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex) -> list[Molecule]:
        """
                    Decompose the Universe into strands (molecules, which could be either chains, or even lonely atoms), without omitting the crosslinkers
                    (as in :meth:`~pylimer_tools_cpp.Universe.get_molecules`).
                    In turn, e.g. for a tetrafunctional crosslinker, it will be 4 times in the resulting molecules.
                    
                    .. note::
                       Crosslinkers without bonds to non-crosslinkers are not returned 
                       (i.e., single crosslinkers, are not counted as strands).
                       
                    :param crosslinker_type: The type of crosslinker atoms
                    :return: List of chains including crosslinkers
        """
    def get_clusters(self) -> list[Universe]:
        """
                  Get the components of the universe that are not connected to each other.
                  Returns a list of :class:`~pylimer_tools_cpp.Universe` objects.
                  Unconnected, free atoms/beads become their own :class:`~pylimer_tools_cpp.Universe`.
                  
                  :return: List of disconnected Universe components
        """
    def get_edge_ids_from(self, vertex_id: typing.SupportsInt | typing.SupportsIndex) -> list[int]:
        """
                  Get the edge IDs incident to a specific vertex.
        
                  :param vertex_id: The ID of the vertex to query
                  :return: List of edge IDs connected to the vertex
        """
    def get_edge_ids_from_to(self, vertex_id_from: typing.SupportsInt | typing.SupportsIndex, vertex_id_to: typing.SupportsInt | typing.SupportsIndex) -> list[int]:
        """
                  Get the edge IDs of all edges between two specific vertices.
        
                  :param vertex_id_from: The starting vertex ID
                  :param vertex_id_to: The ending vertex ID
                  :return: List of edge IDs connecting the two vertices
        """
    def get_edges(self) -> dict[str, list[int]]:
        """
                    Get all edges. Returns a dict with three properties: 'edge_from', 'edge_to' and 'edge_type'.
                    The order is not necessarily related to any structural characteristic.
                    
                    .. note::
                       The integer values returned refer to the vertex IDs, not the atom IDs.
                       Use :meth:`~pylimer_tools_cpp.Universe.get_atom_id_by_idx` to translate them to atom IDs, or
                       :meth:`~pylimer_tools_cpp.Universe.get_bonds` to have that done for you.
                       
                    :return: Dictionary with edge information
        """
    def get_masses(self) -> dict[int, float]:
        """
                    Get the mass of one atom per type.
                    
                    :return: Dictionary mapping atom types to masses
        """
    def get_molecules(self, atom_type_to_omit: typing.SupportsInt | typing.SupportsIndex) -> list[Molecule]:
        """
                  Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms.
                  
                  Reduces the Universe to a list of molecules. 
                  Specify the crosslinker_type to an existing type ID, 
                  then those atoms will be omitted, and this function returns chains instead.
        
                  :param atom_type_to_omit: The type of atom to omit from the universe to end up with the desired molecules (e.g., the type of the crosslinkers).
                  :return: List of molecules
        """
    def get_network_of_crosslinker(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex) -> Universe:
        """
                    Reduce the network to contain only crosslinkers, replacing all the strands with a single bond.
                    Useful e.g. to reduce the memory usage and runtime of 
                    :meth:`~pylimer_tools_cpp.Universe.find_loops` or 
                    :meth:`~pylimer_tools_cpp.Universe.has_infinite_strand`.
                    
                    Further use :meth:`~pylimer_tools_cpp.Universe.simplify` to remove primary loops.
                    
                    :param crosslinker_type: The type of crosslinker atoms
                    :return: Reduced network containing only crosslinkers
        """
    def get_nr_of_angles(self) -> int:
        """
                    Query the number of angles that have been added to this universe.
        """
    def get_nr_of_atoms(self) -> int:
        """
                    Query the number of atoms in this universe.
                    
                    :return: Number of atoms
        """
    def get_nr_of_bonds(self) -> int:
        """
                    Query the number of bonds associated with this universe.
        """
    def get_nr_of_bonds_of_atom(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
               Count the number of immediate neighbors of an atom, specified by its ID.
        
               :param atom_id: The ID of the atom to query
               :return: Number of bonds connected to the atom
        """
    def get_nr_of_bonds_of_vertex(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
               Count the number of immediate neighbors of an atom, specified by its vertex ID.
        
               :param vertex_id: The vertex ID of the atom to query
               :return: Number of bonds connected to the vertex
        """
    def get_nr_of_dihedral_angles(self) -> int:
        """
                    Query the number of dihedral angles that have been added to this universe.
        """
    def get_nr_of_edges_from_to(self, vertex_id_from: typing.SupportsInt | typing.SupportsIndex, vertex_id_to: typing.SupportsInt | typing.SupportsIndex, max_length: typing.SupportsInt | typing.SupportsIndex = -1) -> int:
        """
                  Get the number of edges in the shortest path between two specific vertices.
        
                  If `max_length` is provided and positive, it will only consider paths up to that length.
        
                  :param vertex_id_from: The starting vertex ID
                  :param vertex_id_to: The ending vertex ID  
                  :param max_length: Maximum path length to consider (default: -1 for no limit)
                  :return: Number of edges in shortest path, or -1 if no path exists
        """
    def get_timestep(self) -> int:
        """
                    Query the timestep when this universe was captured.
        """
    def get_vertex_degrees(self) -> list[int]:
        """
                  Get the degree (functionality) of each vertex.
        """
    def get_vertex_idx_by_atom_id(self, atom_id: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
                 Get the vertex ID of the underlying graph for an atom with a specified ID.
                 
                 :param atom_id: The atom ID to look up
                 :return: The vertex index corresponding to the atom
        """
    def get_volume(self) -> float:
        """
                    Query the volume of the underlying bounding box.
                    
                    :return: The volume of the simulation box
        """
    def has_atom_with_id(self, atom_id: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
                 Check whether this universe contains an atom with the specified ID.
        
                 :param atom_id: The atom ID to check
                 :return: True if the atom exists in this universe, False otherwise
        """
    def has_infinite_strand(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
                   Check whether there is a strand (with crosslinker) in the universe that loops through periodic images without coming back.
                   
                    .. warning::
                       There are exponentially many paths between two crosslinkers of a network,
                       and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
                       
                   :return: True if infinite strands are detected, False otherwise
        """
    def hash_angle_type(self, angle_from: typing.SupportsInt | typing.SupportsIndex, angle_via: typing.SupportsInt | typing.SupportsIndex, angle_to: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
                  Convert the three integers to one long number/hash.
                  Used internally for duplicate detection.
                  
                  :param angle_from: First atom ID in the angle
                  :param angle_via: Middle atom ID in the angle
                  :param angle_to: Last atom ID in the angle
                  :return: Hash value for the angle type
        """
    def hash_dihedral_angle_type(self, angle_from: typing.SupportsInt | typing.SupportsIndex, angle_via1: typing.SupportsInt | typing.SupportsIndex, angle_via2: typing.SupportsInt | typing.SupportsIndex, angle_to: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
                  Convert the four integers to one long number/hash.
                  Used internally for duplicate detection.
                  
                  :param angle_from: First atom ID in the dihedral
                  :param angle_via1: Second atom ID in the dihedral
                  :param angle_via2: Third atom ID in the dihedral
                  :param angle_to: Fourth atom ID in the dihedral
                  :return: Hash value for the dihedral angle type
        """
    def identify_active_atoms(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2) -> list[bool]:
        """
                    Identify which atoms belong to the topologically active 2-core of the network.
        
                    Atoms on pruned leaves, dangling branches, and soluble components will have false.
        
                    .. note::
                       Uses :func:`peel_topological_leaves` and is therefore purely topological
                       (ignores coordinates and PBC). Local (non-PBC) primary loops and isolated
                       multi-crosslinker cycles are still marked ``True``, while crosslinker-less
                       PBC-spanning rings are marked ``False``.
        
                    :param crosslinker_type: Type ID of crosslinker atoms (default: 2)
                    :return: Boolean vector of length get_nr_of_atoms() (indexed by vertex index)
        """
    def interpolate_edges(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex, interpolation_factor: typing.SupportsFloat | typing.SupportsIndex) -> list[tuple[int, int]]:
        """
                  Get more or less edges than currently present, interpolating between junctions.
                  
                  :param crosslinker_type: The type of crosslinker atoms
                  :param interpolation_factor: Factor for edge interpolation
                  :return: Interpolated edge structure
        """
    def peel_topological_leaves(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2) -> list[bool]:
        """
                    Peel 1-core dangling trees and subgraphs from the network topology.
        
                    Returns a boolean vector indicating for each strand whether it is ignored (pruned).
        
                    .. note::
                       Purely topological (2-core peeling on crosslinker connectivity); it does not
                       inspect coordinates or periodic boundary conditions (PBC). Crosslinker-less
                       solvent rings are always pruned (even if PBC-spanning), whereas primary loops
                       and multi-crosslinker cycles are always kept (even if local / non-PBC, since
                       a primary loop can still become elastically active when an entanglement or
                       slip-link is associated with it).
                       Use :class:`~pylimer_tools_cpp.MEHPSimulator` with
                       ``assume_box_large_enough=False`` to separate PBC-spanning (or entangled)
                       active loops from local zero-tension loops.
        
                    :param crosslinker_type: Type ID of crosslinker atoms (default: 2)
                    :return: Vector of booleans where true means the strand is pruned/ignored
        """
    def remove_all_angles(self) -> None:
        """
                  Remove all angles from the Universe. 
                  This will not remove the atoms or bonds, just the angles.
                  
                  .. warning::
                    This will not remove dihedral angles, use :meth:`~pylimer_tools_cpp.Universe.remove_all_dihedral_angles` for that.
        """
    def remove_all_dihedral_angles(self) -> None:
        """
                  Remove all dihedral angles from the Universe. 
                  This will not remove the atoms or bonds, just the stored dihedral angles.
        """
    def remove_atoms(self, atom_ids: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        """
                  Remove atoms and all associated bonds by their atom IDs. 
                  
                  :param atom_ids: Vector of atom IDs to remove
        """
    def remove_bonds(self, bonds_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], bonds_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        """
                  Remove bonds by their connected atom IDs. 
                  
                  :param bonds_from: Vector of starting atom IDs
                  :param bonds_to: Vector of ending atom IDs
        """
    def remove_bonds_by_type(self, bond_type: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
                  Remove bonds with a specific type. 
                  
                  :param bond_type: The bond type to remove
        """
    def replace_atom(self, atom_id: typing.SupportsInt | typing.SupportsIndex, replacement_atom: Atom) -> None:
        """
                  Replace the properties of an atom with the properties of another given atom.
                  
                  :param atom_id: ID of the atom to replace
                  :param replacement_atom: The atom with new properties
        """
    def replace_atom_type(self, atom_id: typing.SupportsInt | typing.SupportsIndex, new_type: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
                  Replace the type of an atom with another type.
                  
                  :param atom_id: ID of the atom to modify
                  :param new_type: The new atom type
        """
    def resample_velocities(self, mean: typing.SupportsFloat | typing.SupportsIndex, variance: typing.SupportsFloat | typing.SupportsIndex, seed: str = '', is_2d: bool = False) -> None:
        """
                 Resample velocities for atoms in the universe.
                 
                 :param mean: Mean velocity
                 :param variance: Velocity variance
                 :param seed: Random seed for velocity generation
                 :param is_2d: Whether to omit sampling the z direction
        """
    def set_box(self, box: Box, rescale_atoms: bool = False) -> None:
        """
                  Override the currently assigned box with the one specified.
                  
                  :param box: The new box to assign
                  :param rescale_atoms: Whether to rescale atom positions
        """
    def set_box_lengths(self, lx: typing.SupportsFloat | typing.SupportsIndex, ly: typing.SupportsFloat | typing.SupportsIndex, lz: typing.SupportsFloat | typing.SupportsIndex, rescale_atoms: bool = False) -> None:
        """
                  Override the currently assigned box with one with the side lengths specified.
                  
                  :param lx: Length in x direction
                  :param ly: Length in y direction
                  :param lz: Length in z direction
                  :param rescale_atoms: Whether to rescale atom positions
        """
    def set_mass(self, atom_type: typing.SupportsInt | typing.SupportsIndex, mass: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
                 Set the mass for a specific atom type.
                 
                 :param atom_type: The atom type to set mass for
                 :param mass: The mass value to assign
        """
    def set_masses(self, mass_per_type: collections.abc.Mapping[typing.SupportsInt | typing.SupportsIndex, typing.SupportsFloat | typing.SupportsIndex]) -> None:
        """
                 Set the mass per type of atom.
                 
                 :param mass_per_type: Map of atom types to their masses
        """
    def set_timestep(self, timestep: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
                 Set the timestep when this Universe was captured.
                 
                 :param timestep: The timestep value
        """
    def set_vertex_property(self, vertex_id: typing.SupportsInt | typing.SupportsIndex, property_name: str, value: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
                  Set a specific property for a specific vertex.
                  
                  :param vertex_id: The vertex ID to modify
                  :param property_name: Name of the property to set
                  :param value: Value to assign to the property
        """
    def simplify(self) -> None:
        """
                 Remove self links and double bonds. This function is called 
                 automatically after adding bonds.
                 
                 This operation cleans up the graph structure by removing
                 redundant connections.
        """
    def to_crosslinker_universe(self, crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2) -> CrosslinkerUniverse:
        """
                    Convert this Universe into a CrosslinkerUniverse representation.
        
                    Extracts all crosslinkers and strands, lining up beads deterministically.
        
                    :param crosslinker_type: Type ID of crosslinker atoms (default: 2)
                    :return: CrosslinkerUniverse struct
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
         Alternatively, use :meth:`~pylimer_tools_cpp.UniverseSequence.forget_at_index`
         to have the UniverseSequence forget about already read universes.
         
    """
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> Universe:
        """
              Get a universe by its index.
              
              :param index: The index of the universe to retrieve
              :return: The Universe at the specified index
        """
    def __init__(self) -> None:
        """
                 Construct an empty UniverseSequence.
                 
                 Use initialization methods to populate it with data.
        """
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
                 Get the number of universes in this sequence.
                 
                 :return: The total number of universes available
        """
    def at_index(self, index: typing.SupportsInt | typing.SupportsIndex) -> Universe:
        """
                 Get the Universe at the given index (as of in the sequence given 
                 by the dump file).
                 
                 :param index: The index of the universe to retrieve
                 :return: The Universe at the specified index
        """
    def compute_distance_autocorrelation_from_to(self, atom_ids_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], atom_ids_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], nr_of_origins: typing.SupportsInt | typing.SupportsIndex = 25, reduce_memory: bool = False) -> dict[int, float]:
        """
                  Compute the autocorrelation of the dot product of the distance vector from certain to other atoms.
        
                  For example, this can be used to compute Eq. 4.51 from Masao Doi, Introduction to Polymer Physics, p. 74.
        """
    def compute_distance_from_to_atoms(self, atom_ids_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], atom_ids_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], reduce_memory: bool = False) -> list[float]:
        """
                  Compute the root square norm of all the (unwrapped!) distances for the given pair of atoms.
                  
                  Can be used to somewhat faster compute e.g. all the end-to-end or bond distances.
                  Pay attention that the image flags are correct, otherwise, this data may not be useable.
        """
    def compute_msd_for_atom_properties(self, atom_ids: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], x_property: str, y_property: str, z_property: str, nr_of_origins: typing.SupportsInt | typing.SupportsIndex = 25, reduce_memory: bool = False, max_tau: typing.SupportsInt | typing.SupportsIndex = -1) -> dict[int, float]:
        """
                  Compute the mean square displacement for atoms using specified property names.
                  
                  :param atom_ids: List of atom IDs for which to compute the MSD
                  :param x_property: Name of the x-coordinate property in the dump file (e.g., "x", "xu", "xs")
                  :param y_property: Name of the y-coordinate property in the dump file (e.g., "y", "yu", "ys")
                  :param z_property: Name of the z-coordinate property in the dump file (e.g., "z", "zu", "zs")
                  :param nr_of_origins: Number of time origins to use for averaging. Higher values provide better statistics but increase computation time (default: 25)
                  :param reduce_memory: If True, reduces memory usage by forgetting universes after processing them (default: False)
                  :param max_tau: Maximum time lag (tau) to compute. If -1, computes for all possible tau values. For better statistics, consider setting this to approximately half the sequence length (default: -1)
                  :return: Dictionary mapping time lag (tau) to mean square displacement values
        """
    def compute_msd_for_atoms(self, atom_ids: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], nr_of_origins: typing.SupportsInt | typing.SupportsIndex = 25, reduce_memory: bool = False, max_tau: typing.SupportsInt | typing.SupportsIndex = -1) -> dict[int, float]:
        """
                  Compute the mean square displacement for atoms with the specified IDs.
                  
                  :param atom_ids: List of atom IDs for which to compute the MSD
                  :param nr_of_origins: Number of time origins to use for averaging. Higher values provide better statistics but increase computation time (default: 25)
                  :param reduce_memory: If True, reduces memory usage by forgetting universes after processing them (default: False)
                  :param max_tau: Maximum time lag (tau) to compute. If -1, computes for all possible tau values. For better statistics, consider setting this to approximately half the sequence length (default: -1)
                  :return: Dictionary mapping time lag (tau) to mean square displacement values
        """
    def compute_vector_from_to_atoms(self, atom_ids_from: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], atom_ids_to: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], reduce_memory: bool = False) -> list[typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]]:
        """
                  Compute the (unwrapped!) distances for the given pair of atoms.
                  
                  Can be used to somewhat faster compute e.g. all the end-to-end or bond vectors.
                  Pay attention that the image flags are correct, otherwise, this data may not be useable.
        """
    def forget_at_index(self, index: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
              Clear the memory of the Universe at the given index (as of in the
              sequence given by the dump file).
        """
    def get_all(self) -> list[Universe]:
        """
                    Get all universes initialized back in a list.
                    For big dump files or lots of data files, this might lead to memory issues.
                    Use :meth:`~pylimer_tools_cpp.UniverseSequence.__iter__`
                    or :meth:`~pylimer_tools_cpp.UniverseSequence.at_index`
                    and :meth:`~pylimer_tools_cpp.UniverseSequence.forget_at_index`
                    to craft a more memory-efficient retrieval mechanism.
                    
                    Returns:
                        A list of all Universe objects in the sequence
        """
    def get_length(self) -> int:
        """
                    Get the number of universes in this sequence.
        """
    def initialize_from_data_sequence(self, data_files: collections.abc.Sequence[str]) -> None:
        """
        Reset and initialize the Universes from an ordered list of Lammps data (:code:`write_data`) files.
        """
    def initialize_from_dump_file(self, initial_data_file: str, dump_file: str) -> None:
        """
                  Reset and initialize the Universes from a Lammps :code:`dump` output. 
                
                  NOTE:
                       If you have not output the id of the atoms in the dump file, they will be assigned sequentially. 
                       If you have not output the type of the atoms in the dump file, they will be set to -1 if they cannot be infered from the data file.
        """
    def next(self) -> Universe:
        """
                 Get the Universe that's next in the sequence.
                 
                 :return: The next Universe in the sequence
        """
    def reset_iterator(self) -> None:
        """
                  Reset the internal iterator, such that a subsequent call to 
                  :meth:`~pylimer_tools_cpp.UniverseSequence.next` returns the first one again.
        """
    def set_data_file_atom_style(self, atom_styles: collections.abc.Sequence[AtomStyle]) -> None:
        """
                  Set the format of the data files to be read. See :obj:`~pylimer_tools_cpp.AtomStyle`.
        """
class ZScoreMaxDistanceProvider(MaxDistanceProvider):
    """
    
         For MC generation, converts the :math:`N` to a maximum distance within which to sample.
         The distance will be calculated as :math:`\\text{std_multiplier} \\times \\sqrt{N \\times \\text{in_sqrt_multiplier}}`.
         Useful only for performance improvements in large systems.
    
         :param std_multiplier: The Z-Score, the multiplier of the standard deviation of the end-to-end distribution.
             E.g., 3.29 for 99.9% of all conformations.
         :param in_sqrt_multiplier: The multiplier with the :math:`N` in the square root. Probably :math:`<b^2>`.
        
    """
    def __init__(self, std_multiplier: typing.SupportsFloat | typing.SupportsIndex, in_sqrt_multiplier: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    def __repr__(self) -> str:
        """
                 Return a string representation of the ZScoreMaxDistanceProvider.
                 
                 :return: String representation
        """
    def get_max_distance(self, N: typing.SupportsFloat | typing.SupportsIndex) -> float:
        """
                 Get the maximum distance for a given N.
        
                 :param N: Number of segments.
        """
def do_linear_walk_chain_from_to(box: Box, from_coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], to_coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], chain_len: typing.SupportsInt | typing.SupportsIndex, include_ends: bool = False) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
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
def do_random_walk(chain_len: typing.SupportsInt | typing.SupportsIndex, bead_distance: typing.SupportsFloat | typing.SupportsIndex = 1.0, mean_squared_bead_distance: typing.SupportsFloat | typing.SupportsIndex = 1.0, seed: str = '') -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
    """
                Do a random walk, return the coordinates of each point visited.
    
                :param chain_len: Length of the chain to generate.
                :param bead_distance: Mean distance between consecutive beads (default: 1.0).
                :param mean_squared_bead_distance: Mean squared distance between consecutive beads (default: 1.0).
                :param seed: Random seed for reproducibility (default: empty string for random seed).
                :return: Coordinates of each point as a flat array (i.e., [x1, y1, z1, x2, y2, z2, ...]).
    """
def do_random_walk_chain_from_to(box: Box, from_coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], to_coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], chain_len: typing.SupportsInt | typing.SupportsIndex, bead_distance: typing.SupportsFloat | typing.SupportsIndex = 1.0, mean_squared_bead_distance: typing.SupportsFloat | typing.SupportsIndex = 1.0, seed: str = '') -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
    """
                Do a random walk from one point to another.
    
                :param box: Simulation box for periodic boundary conditions.
                :param from_coordinates: Starting coordinates as 3D vector.
                :param to_coordinates: Target coordinates as 3D vector.
                :param chain_len: Number of beads to place between start and end points.
                :param bead_distance: Mean distance between consecutive beads (default: 1.0).
                :param mean_squared_bead_distance: Mean squared distance between consecutive beads (default: 1.0).
                :param seed: Random seed for reproducibility (default: empty string for no seed).
                :return: Coordinates of the chain as a flat array (i.e., [x1, y1, z1, x2, y2, z2, ...]).
    """
def do_random_walk_chain_from_to_mc(box: Box, from_coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], to_coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], chain_len: typing.SupportsInt | typing.SupportsIndex, bead_distance: typing.SupportsFloat | typing.SupportsIndex = 1.0, mean_squared_bead_distance: typing.SupportsFloat | typing.SupportsIndex = 1.0, seed: str = '', n_iterations: typing.SupportsInt | typing.SupportsIndex = 10000) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
    """
                Do a random walk from one point to another.
                Then, relax the points in between using a Metropolis-Monte Carlo simulation.
    
                :param box: Simulation box for periodic boundary conditions.
                :param from_coordinates: Starting coordinates as 3D vector.
                :param to_coordinates: Target coordinates as 3D vector.
                :param chain_len: Number of beads to place between start and end points.
                :param bead_distance: Mean distance between consecutive beads (default: 1.0).
                :param mean_squared_bead_distance: Mean squared distance between consecutive beads (default: 1.0).
                :param seed: Random seed for reproducibility (default: empty string for no seed).
                :param n_iterations: Number of Monte Carlo iterations for relaxation (default: 10000).
                :return: Coordinates of the relaxed chain as a flat array (i.e., [x1, y1, z1, x2, y2, z2, ...]).
    """
def inverse_langevin(x: typing.SupportsFloat | typing.SupportsIndex) -> float:
    """
         A somewhat accurate (for :math:`x \\in (-1, 1)`) implementation of the inverse Langevin function.
    
         Source: https://scicomp.stackexchange.com/a/30251
    
         :param x: Input value in the range (-1, 1)
         :return: Inverse Langevin function value
    """
def peel_topological_leaves(crosslinker_universe: CrosslinkerUniverse) -> list[bool]:
    """
        Peel 1-core dangling trees and subgraphs iteratively from the CrosslinkerUniverse topology.
        Returns a boolean vector indicating for each strand whether it is ignored (pruned).
    
        .. note::
           This is a purely topological 2-core peeling on ``strand_from``/``strand_to`` and
           does not account for coordinates or periodic boundary conditions (PBC):
    
           - Free/dangling strands (including crosslinker-less solvent rings) are always pruned,
             even if they wrap across the PBC.
           - Primary loops (``strand_from == strand_to >= 0``) and multi-crosslinker cycles are
             always kept in the 2-core, even if they do not span the PBC (note that a
             non-PBC primary loop can still become elastically active if an entanglement or
             slip-link is associated with it, whereas in a phantom network without
             entanglements it relaxes to zero tension).
    
           To distinguish PBC-spanning (or entangled) active loops from non-PBC
           zero-tension loops, use :class:`~pylimer_tools_cpp.MEHPSimulator` with
           ``assume_box_large_enough=False``.
    
        :param crosslinker_universe: The CrosslinkerUniverse representation.
        :return: Boolean vector of ignored/pruned strands.
    """
@typing.overload
def randomly_sample_entanglements(universe: Universe, nr_of_samples: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_samples: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = '', crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, ignore_crosslinks: bool = True, filtering_mode: EntanglementFilteringMode = ...) -> AtomPairEntanglements:
    """
        Randomly find pairs of atoms that are close together and could be entanglements.
    
        :param universe: The universe of atoms from which to sample entanglements from.
        :param nr_of_samples: The number of pairs of atoms to randomly sample.
        :param upper_cutoff: The maximum distance between atoms for a pair to be considered a potential entanglement.
        :param lower_cutoff: The minimum distance between atoms for a pair to be considered a potential entanglement.
        :param minimum_nr_of_samples: The minimum number of entanglements to be found.
        :param same_strand_cutoff: The minimum shortest path distance on the graph between atoms for an entanglement.
        :param seed: A seed for the random number generator.
        :param crosslinker_type: Atom type ID representing crosslinkers.
        :param ignore_crosslinks: Whether to exclude crosslink atoms from forming entanglements.
        :param filtering_mode: Filtering mode for entanglements (default: EntanglementFilteringMode.NONE).
    """
@typing.overload
def randomly_sample_entanglements(box: Box, crosslinker_universe: CrosslinkerUniverse, nr_of_samples: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_samples: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = '', filtering_mode: EntanglementFilteringMode = ..., n_mc_steps: typing.SupportsInt | typing.SupportsIndex = 0) -> StrandPairEntanglements:
    ...
@typing.overload
def randomly_sample_entanglements(crosslinker_universe: CrosslinkerUniverse, box: Box, nr_of_samples: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_samples: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = '', filtering_mode: EntanglementFilteringMode = ..., n_mc_steps: typing.SupportsInt | typing.SupportsIndex = 0) -> StrandPairEntanglements:
    ...
@typing.overload
def randomly_sample_entanglements_on_demand(box: Box, crosslinker_universe: CrosslinkerUniverse, nr_of_samples: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_samples: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = '', filtering_mode: EntanglementFilteringMode = ..., n_mc_steps: typing.SupportsInt | typing.SupportsIndex = 0) -> StrandPairEntanglements:
    """
        Sample pairs of strand beads that could form entanglements on demand from a CrosslinkerUniverse
        using per-strand bounding boxes and on-demand Brownian bridge bead resampling.
        Only strand seeds and bounding boxes are stored, achieving O(N_strands) memory scaling
        while sampling from the exact discrete multi-bead chain conformations.
    
        :param box: The simulation box for periodic boundary conditions.
        :param crosslinker_universe: The CrosslinkerUniverse of crosslinkers and strands.
        :param nr_of_samples: The number of entanglement pairs to randomly sample.
        :param upper_cutoff: The maximum distance between two beads to be considered an entanglement.
        :param lower_cutoff: The minimum distance between two beads to be considered an entanglement.
        :param minimum_nr_of_samples: The minimum number of entanglements to be found.
        :param same_strand_cutoff: The minimum contour separation (in beads) between two beads on the same strand.
        :param seed: A seed for the random number generator.
        :param filtering_mode: Filtering mode for entanglements (default: EntanglementFilteringMode.OMIT_INACTIVE).
        :param n_mc_steps: Number of MC equilibration steps per chain (default: 0).
    """
@typing.overload
def randomly_sample_entanglements_on_demand(crosslinker_universe: CrosslinkerUniverse, box: Box, nr_of_samples: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_samples: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = '', filtering_mode: EntanglementFilteringMode = ..., n_mc_steps: typing.SupportsInt | typing.SupportsIndex = 0) -> StrandPairEntanglements:
    ...
def randomly_sample_entanglements_v2(universe: Universe, nr_of_samples: typing.SupportsInt | typing.SupportsIndex, upper_cutoff: typing.SupportsFloat | typing.SupportsIndex, lower_cutoff: typing.SupportsFloat | typing.SupportsIndex = 0.0, minimum_nr_of_samples: typing.SupportsInt | typing.SupportsIndex = 0, same_strand_cutoff: typing.SupportsFloat | typing.SupportsIndex = 3.0, seed: str = '', crosslinker_type: typing.SupportsInt | typing.SupportsIndex = 2, ignore_crosslinks: bool = True, filtering_mode: EntanglementFilteringMode = ...) -> AtomPairEntanglements:
    """
    Alias for randomly_sample_entanglements (v2 has been unified into standard method).
    """
def sample_strand_beads_on_demand(box: Box, crosslinker_universe: CrosslinkerUniverse, strand_idx: typing.SupportsInt | typing.SupportsIndex, strand_seed: typing.SupportsInt | typing.SupportsIndex, n_mc_steps: typing.SupportsInt | typing.SupportsIndex = 0) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
    """
        Deterministically realize the internal bead coordinates of a single strand
        in a CrosslinkerUniverse using a per-strand PRNG seed.
    
        :param box: The simulation box.
        :param crosslinker_universe: The CrosslinkerUniverse of crosslinkers and strands.
        :param strand_idx: Index of the strand in `crosslinker_universe`.
        :param strand_seed: Deterministic per-strand PRNG seed.
        :param n_mc_steps: Number of Monte Carlo equilibration steps (default: 0).
        :return: Vector of length `3 * beads_in_strand[strand_idx]` with bead coordinates.
    """
def split_csv(file_path: str, delimiter: str = ',') -> list[str]:
    """
            Read a file containing a number of CSVs. Returns them split up.
            
            :param file_path: Path to the file containing CSV data
            :param delimiter: Delimiter used in the CSV file (default: ',')
            :return: Split CSV data structures
    """
def version_information() -> str:
    """
        Returns  a string of the the current version, incl. git hash and date of compilation.
    """
STOP: BackTrackStatus  # value = <BackTrackStatus.STOP: 0>
TRACK_BACKWARD: BackTrackStatus  # value = <BackTrackStatus.TRACK_BACKWARD: 2>
TRACK_FORWARD: BackTrackStatus  # value = <BackTrackStatus.TRACK_FORWARD: 1>
__version__: str = '0.7.0'
