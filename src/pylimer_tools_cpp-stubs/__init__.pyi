"""

    pylimer_tools_cpp
    -----------------

    A collection of utility python functions for handling LAMMPS output and polymers in Python.

    .. autosummary::
        :toctree: _generate

    
"""
from __future__ import annotations
import collections.abc
import numpy
import numpy.typing
import scipy.sparse
import typing
__all__: list[str] = ['Atom', 'AtomPairEntanglements', 'AtomStyle', 'AveFileReader', 'BackTrackStatus', 'Box', 'ComputedDoubleValues', 'ComputedIntValues', 'DPDSimulator', 'DataFileReader', 'DataFileWriter', 'DumpFileReader', 'ExitReason', 'LazyUniverseSequenceIterator', 'LinearMaxDistanceProvider', 'LinkSwappingMode', 'MCUniverseGenerator', 'MEHPForceBalance', 'MEHPForceBalance2', 'MEHPForceEvaluator', 'MEHPForceRelaxation', 'MaxDistanceProvider', 'Molecule', 'MoleculeIterator', 'MoleculeType', 'NeighbourList', 'NoMaxDistanceProvider', 'NonGaussianSpringForceEvaluator', 'NormalModeAnalyzer', 'OutputConfiguration', 'SLESolver', 'STOP', 'SimpleSpringMEHPForceEvaluator', 'SimplifiedBalance2Network', 'SimplifiedBalanceNetwork', 'SimplifiedNetwork', 'StructureSimplificationMode', 'TRACK_BACKWARD', 'TRACK_FORWARD', 'Universe', 'UniverseSequence', 'ZScoreMaxDistanceProvider', 'do_linear_walk_chain_from_to', 'do_random_walk', 'do_random_walk_chain_from_to', 'do_random_walk_chain_from_to_mc', 'inverse_langevin', 'randomly_sample_entanglements', 'split_csv', 'version_information']
class Atom:
    """
    
           A single bead or atom.
      
    """
    __hash__: typing.ClassVar[None] = None
    def __eq__(self, arg0: Atom) -> bool:
        ...
    def __getstate__(self) -> tuple:
        ...
    @typing.overload
    def __init__(self, id: typing.SupportsInt, type: typing.SupportsInt, x: typing.SupportsFloat, y: typing.SupportsFloat, z: typing.SupportsFloat, nx: typing.SupportsInt, ny: typing.SupportsInt, nz: typing.SupportsInt) -> None:
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
    def __init__(self, properties: collections.abc.Mapping[str, typing.SupportsFloat]) -> None:
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
    def pair_of_atom(self, arg0: collections.abc.Sequence[typing.SupportsInt]) -> None:
        ...
    @property
    def pairs_of_atoms(self) -> list[tuple[int, int]]:
        """
              A list of pairs of atom ids that are close together and could be entanglements
        """
    @pairs_of_atoms.setter
    def pairs_of_atoms(self, arg0: collections.abc.Sequence[tuple[typing.SupportsInt, typing.SupportsInt]]) -> None:
        ...
class AtomStyle:
    """
    An enumeration of the LAMMPS atom styles.
    
    Members:
    
      NONE : LAMMPS atom style 'none'
    
      ANGLE : LAMMPS atom style 'angle'
    
      ATOMIC : LAMMPS atom style 'atomic'
    
      BODY : LAMMPS atom style 'body'
    
      BOND : LAMMPS atom style 'bond'
    
      BPM_SPHERE : LAMMPS atom style 'bpm/sphere'
    
      CHARGE : LAMMPS atom style 'charge'
    
      DIELECTRIC : LAMMPS atom style 'dielectric'
    
      DIPOLE : LAMMPS atom style 'dipole'
    
      DPD : LAMMPS atom style 'dpd'
    
      EDPD : LAMMPS atom style 'edpd'
    
      ELECTRON : LAMMPS atom style 'electron'
    
      ELLIPSOID : LAMMPS atom style 'ellipsoid'
    
      FULL : LAMMPS atom style 'full'
    
      LINE : LAMMPS atom style 'line'
    
      MDPD : LAMMPS atom style 'mdpd'
    
      MOLECULAR : LAMMPS atom style 'molecular'
    
      PERI : LAMMPS atom style 'peri'
    
      RHEO : LAMMPS atom style 'rheo'
    
      RHEO_THERMAL : LAMMPS atom style 'rheo/thermal'
    
      SMD : LAMMPS atom style 'smd'
    
      SPH : LAMMPS atom style 'sph'
    
      SPHERE : LAMMPS atom style 'sphere'
    
      SPIN : LAMMPS atom style 'spin'
    
      TDPD : LAMMPS atom style 'tdpd'
    
      TEMPLATE : LAMMPS atom style 'template'
    
      TRI : LAMMPS atom style 'tri'
    
      WAVEPACKET : LAMMPS atom style 'wavepacket'
    
      HYBRID : LAMMPS atom style 'hybrid'
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
    __members__: typing.ClassVar[dict[str, AtomStyle]]  # value = {'NONE': <AtomStyle.NONE: 0>, 'ANGLE': <AtomStyle.ANGLE: 1>, 'ATOMIC': <AtomStyle.ATOMIC: 2>, 'BODY': <AtomStyle.BODY: 3>, 'BOND': <AtomStyle.BOND: 4>, 'BPM_SPHERE': <AtomStyle.BPM_SPHERE: 5>, 'CHARGE': <AtomStyle.CHARGE: 6>, 'DIELECTRIC': <AtomStyle.DIELECTRIC: 7>, 'DIPOLE': <AtomStyle.DIPOLE: 8>, 'DPD': <AtomStyle.DPD: 9>, 'EDPD': <AtomStyle.EDPD: 10>, 'ELECTRON': <AtomStyle.ELECTRON: 11>, 'ELLIPSOID': <AtomStyle.ELLIPSOID: 12>, 'FULL': <AtomStyle.FULL: 13>, 'LINE': <AtomStyle.LINE: 14>, 'MDPD': <AtomStyle.MDPD: 15>, 'MOLECULAR': <AtomStyle.MOLECULAR: 16>, 'PERI': <AtomStyle.PERI: 17>, 'RHEO': <AtomStyle.RHEO: 18>, 'RHEO_THERMAL': <AtomStyle.RHEO_THERMAL: 19>, 'SMD': <AtomStyle.SMD: 20>, 'SPH': <AtomStyle.SPH: 21>, 'SPHERE': <AtomStyle.SPHERE: 22>, 'SPIN': <AtomStyle.SPIN: 23>, 'TDPD': <AtomStyle.TDPD: 24>, 'TEMPLATE': <AtomStyle.TEMPLATE: 25>, 'TRI': <AtomStyle.TRI: 26>, 'WAVEPACKET': <AtomStyle.WAVEPACKET: 27>, 'HYBRID': <AtomStyle.HYBRID: 28>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
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
    def __init__(self, file_path: str) -> None:
        """
                 Initialize the AveFileReader with a file path.
                 
                 :param file_path: Path to the averages file to read
        """
    def autocorrelate_column(self, column_index: typing.SupportsInt, delta_indices: collections.abc.Sequence[typing.SupportsInt]) -> list[float]:
        """
                  Do autocorrelation on one particular column for a specified set of delta indices.
        
                  Assumes the data is equally spaced.
                  
                  :param column_index: Index of the column to autocorrelate
                  :param delta_indices: List of delta indices for the autocorrelation
                  :return: Autocorrelation values for the specified deltas
        """
    def autocorrelate_column_difference(self, column_index1: typing.SupportsInt, column_index2: typing.SupportsInt, delta_indices: collections.abc.Sequence[typing.SupportsInt]) -> list[float]:
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
class BackTrackStatus:
    """
    
         Enum for controlling the strand linking process in linkStrandsCallback.
         
    
    Members:
    
      STOP : Stop the linking process
    
      TRACK_FORWARD : Continue linking forward
    
      TRACK_BACKWARD : Track backward in the linking process
    """
    STOP: typing.ClassVar[BackTrackStatus]  # value = <BackTrackStatus.STOP: 0>
    TRACK_BACKWARD: typing.ClassVar[BackTrackStatus]  # value = <BackTrackStatus.TRACK_BACKWARD: 2>
    TRACK_FORWARD: typing.ClassVar[BackTrackStatus]  # value = <BackTrackStatus.TRACK_FORWARD: 1>
    __members__: typing.ClassVar[dict[str, BackTrackStatus]]  # value = {'STOP': <BackTrackStatus.STOP: 0>, 'TRACK_FORWARD': <BackTrackStatus.TRACK_FORWARD: 1>, 'TRACK_BACKWARD': <BackTrackStatus.TRACK_BACKWARD: 2>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class Box:
    """
    
            The box that the simulation is run in.
    
            .. note:: 
              Currently, only rectangular boxes are supported.
            
    """
    def __getstate__(self) -> tuple:
        ...
    @typing.overload
    def __init__(self, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat, arg3: typing.SupportsFloat, arg4: typing.SupportsFloat, arg5: typing.SupportsFloat) -> None:
        ...
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
    def apply_simple_shear(self, shear_magnitude: typing.SupportsFloat, shear_direction: typing.SupportsInt = 0) -> None:
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
    def is_valid_offset(self, potential_offset: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], abs_precision: typing.SupportsFloat = 1e-05) -> bool:
        """
                  Check whether the passed offset is a valid one in this box.
                  
                  :param potential_offset: The offset to validate
                  :param abs_precision: Absolute precision for the validation
                  :return: True if the offset is valid, False otherwise
        """
class ComputedDoubleValues:
    """
    Floating point output quantities
    
    Members:
    
      TIMESTEP : Results in the output column "TimeStep".
    
      TIME : Results in the output column "Time".
    
      VOLUME : Results in the output column "Volume".
    
      PRESSURE : Results in the output column "Pressure".
    
      TEMPERATURE : Results in the output column "Temperature".
    
      STRESS_XX : Results in the output column "Stress[0,0]".
    
      STRESS_YY : Results in the output column "Stress[1,1]".
    
      STRESS_ZZ : Results in the output column "Stress[2,2]".
    
      STRESS_XY : Results in the output column "Stress[0,1]".
    
      STRESS_YZ : Results in the output column "Stress[1,2]".
    
      STRESS_XZ : Results in the output column "Stress[0,2]".
    
      STRESS_NXY : Results in the output column "Stress[0,0]-Stress[1,1]".
    
      STRESS_NYZ : Results in the output column "Stress[1,1]-Stress[2,2]".
    
      STRESS_NXZ : Results in the output column "Stress[0,0]-Stress[2,2]".
    
      GAMMA : Results in the output column "Gamma".
    
      RESIDUAL : Results in the output column "Residual".
    
      MEAN_B : Results in the output column "<b>".
    
      MAX_B : Results in the output column "max(b)".
    
      MSD : Results in the output column "MSD".
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
    __members__: typing.ClassVar[dict[str, ComputedDoubleValues]]  # value = {'TIMESTEP': <ComputedDoubleValues.TIMESTEP: 0>, 'TIME': <ComputedDoubleValues.TIME: 1>, 'VOLUME': <ComputedDoubleValues.VOLUME: 2>, 'PRESSURE': <ComputedDoubleValues.PRESSURE: 3>, 'TEMPERATURE': <ComputedDoubleValues.TEMPERATURE: 4>, 'STRESS_XX': <ComputedDoubleValues.STRESS_XX: 5>, 'STRESS_YY': <ComputedDoubleValues.STRESS_YY: 6>, 'STRESS_ZZ': <ComputedDoubleValues.STRESS_ZZ: 7>, 'STRESS_XY': <ComputedDoubleValues.STRESS_XY: 8>, 'STRESS_YZ': <ComputedDoubleValues.STRESS_YZ: 9>, 'STRESS_XZ': <ComputedDoubleValues.STRESS_XZ: 10>, 'STRESS_NXY': <ComputedDoubleValues.STRESS_NXY: 11>, 'STRESS_NYZ': <ComputedDoubleValues.STRESS_NYZ: 12>, 'STRESS_NXZ': <ComputedDoubleValues.STRESS_NXZ: 13>, 'GAMMA': <ComputedDoubleValues.GAMMA: 14>, 'RESIDUAL': <ComputedDoubleValues.RESIDUAL: 15>, 'MEAN_B': <ComputedDoubleValues.MEAN_B: 16>, 'MAX_B': <ComputedDoubleValues.MAX_B: 17>, 'MSD': <ComputedDoubleValues.MSD: 18>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
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
    Integer output quantities
    
    Members:
    
      STEP : Results in the output column "Step".
    
      NUM_SHIFT : Results in the output column "numShift".
    
      NUM_RELOC : Results in the output column "numReloc".
    
      NUM_ATOMS : Results in the output column "numAtoms".
    
      NUM_EXTRA_ATOMS : Results in the output column "numExtraAtoms".
    
      NUM_BONDS : Results in the output column "numBonds".
    
      NUM_EXTRA_BONDS : Results in the output column "numExtraBonds".
    
      NUM_BONDS_TO_FORM : Results in the output column "numBondsToForm".
    """
    NUM_ATOMS: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_ATOMS: 3>
    NUM_BONDS: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_BONDS: 5>
    NUM_BONDS_TO_FORM: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_BONDS_TO_FORM: 7>
    NUM_EXTRA_ATOMS: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_EXTRA_ATOMS: 4>
    NUM_EXTRA_BONDS: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_EXTRA_BONDS: 6>
    NUM_RELOC: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_RELOC: 2>
    NUM_SHIFT: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.NUM_SHIFT: 1>
    STEP: typing.ClassVar[ComputedIntValues]  # value = <ComputedIntValues.STEP: 0>
    __members__: typing.ClassVar[dict[str, ComputedIntValues]]  # value = {'STEP': <ComputedIntValues.STEP: 0>, 'NUM_SHIFT': <ComputedIntValues.NUM_SHIFT: 1>, 'NUM_RELOC': <ComputedIntValues.NUM_RELOC: 2>, 'NUM_ATOMS': <ComputedIntValues.NUM_ATOMS: 3>, 'NUM_EXTRA_ATOMS': <ComputedIntValues.NUM_EXTRA_ATOMS: 4>, 'NUM_BONDS': <ComputedIntValues.NUM_BONDS: 5>, 'NUM_EXTRA_BONDS': <ComputedIntValues.NUM_EXTRA_BONDS: 6>, 'NUM_BONDS_TO_FORM': <ComputedIntValues.NUM_BONDS_TO_FORM: 7>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
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
    def __init__(self, universe: Universe, crosslinker_type: typing.SupportsInt = 2, slipspring_bond_type: typing.SupportsInt = 9, is_2d: bool = False, seed: str = '') -> None:
        """
        Get an instance of this class
        """
    def assume_box_large_enough(self) -> None:
        """
                  Configure whether to run PBC on the bonds or not.
        
                  If your bonds could get larger than half the box length, this must be kept false (default).
                  Otherwise, you can set it to true and therewith get some securities.
        """
    def config_a(self, A: typing.SupportsFloat = 25.0) -> None:
        """
                  Configure the force-field (pair-style) parameter `A`.
        """
    def config_allow_relocation_in_network(self, allow_relocation_in_network: bool = False) -> None:
        """
                  Configure whether a relocation step may happen when a slip-spring has ended at a crosslink.
        
                  Side-effect: if true, the relocations may also happen *to* a slip-spring next to a crosslink.
        
                  :param allow_relocation_in_network (bool): Whether to allow relocation in the network or not.
        """
    def config_auto_correlator_output(self, values: collections.abc.Sequence[OutputConfiguration], num_corr_in: typing.SupportsInt = 32, p: typing.SupportsInt = 16, m: typing.SupportsInt = 2) -> None:
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
    def config_bond_formation(self, num_bonds_to_form: typing.SupportsInt, max_bonds_per_atom_type: collections.abc.Mapping[typing.SupportsInt, typing.SupportsInt], bond_formation_dist: typing.SupportsFloat = 1.0, attempt_bond_formation_every: typing.SupportsInt = 50, atom_type_form_from: typing.SupportsInt = 2, atom_type_form_to: typing.SupportsInt = 1) -> None:
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
    def config_lambda(self, l: typing.SupportsFloat = 0.65) -> None:
        """
                  Configure the modified velocity verlet integration parameter `\\lambda`.
        """
    def config_num_steps_dpd(self, num_steps: typing.SupportsInt = 500) -> None:
        """
                  Configure the number of steps to do in one DPD sequence.
        
                  :param num_steps: Number of DPD steps per sequence
        """
    def config_num_steps_mc(self, num_steps: typing.SupportsInt = 500) -> None:
        """
                  Configure the number of steps to do in one MC sequence.
        
                  :param num_steps: Number of Monte Carlo steps per sequence
        """
    def config_restart_output(self, file: str, output_every: typing.SupportsInt = 50000) -> None:
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
    def config_sigma(self, sigma: typing.SupportsFloat = 3.0) -> None:
        """
                  Configure the force-field (pair-style) parameter `\\sigma`.
        """
    def config_slipspring_high_cutoff(self, cutoff: typing.SupportsFloat = 2.0) -> None:
        """
                  Configure the lower cut-off of how far a pair may be distanced for a slip-spring to be created.
        """
    def config_slipspring_low_cutoff(self, cutoff: typing.SupportsFloat = 0.5) -> None:
        """
                  Configure the higher cut-off of how far a pair may be distanced for a slip-spring to be created.
        """
    def config_spring_constant(self, k: typing.SupportsFloat = 2.0) -> None:
        """
                  Configure the force-field (bond-style) parameter `k`, the spring constant.
        """
    def config_step_output(self, values: collections.abc.Sequence[OutputConfiguration]) -> None:
        """
                  Set which values to log.
        
                  :param values: a list of OutputConfiguration structs
        """
    def create_slip_springs(self, num: typing.SupportsInt, bond_type: typing.SupportsInt = 9) -> int:
        """
                  Randomly add the specified number of slip-springs to neighbours within the specified cut-offs.
        """
    def get_bond_lengths(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_current_timestep(self) -> int:
        """
                  Get the current timestep number.
        
                  :return: The current timestep index
        """
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
    def get_stress_tensor(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 3]"]:
        ...
    def get_temperature(self) -> float:
        ...
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
    def get_volume(self) -> float:
        ...
    def refresh_current_state(self) -> None:
        """
                  After re-configuring the force-field parameters,
                  this method should be called to update the current stress tensor etc.
        """
    def run_simulation(self, n_steps: typing.SupportsInt, dt: typing.SupportsFloat = 0.06, with_MC: bool = False) -> None:
        ...
    def start_measuring_msd_for_atoms(self, atom_ids: collections.abc.Sequence[typing.SupportsInt]) -> None:
        """
                  Set a new origin for measuing the mean square displacement for a specified set of atoms
        """
    def validate_neighbour_list(self, arg0: typing.SupportsFloat) -> None:
        """
                  Validate the neighbor list consistency.
        
                  :return: True if the neighbor list is valid
        """
    def validate_state(self) -> None:
        """
                  Validate the current simulation state.
        
                  :return: True if the simulation state is valid
        """
    @typing.overload
    def write_restart_file(self, file: str) -> None:
        """
                  Explicitly force the writing of a restart file, now!
        
                  :param file: The file path and name of the restart file to be written.
                               Can end in .xml, .json or anything else (-> binary)
        """
    @typing.overload
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
    def read(self, path_of_file_to_read: str, atom_style: AtomStyle = ..., atom_style2: AtomStyle = ..., atom_style_3: AtomStyle = ...) -> None:
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
    def config_crosslinker_type(self, crosslinker_type: typing.SupportsInt = 2) -> None:
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
    def get_length(self) -> int:
        """
                 Get the number of sections (time-steps) in the file.
                 
                 :return: Number of timesteps/sections in the dump file
        """
    def get_numeric_values_for_at(self, rowIndex: typing.SupportsInt, headerKey: str, columnIndex: str) -> list[float]:
        """
                 Get numeric values for the section at index, the main header headerKey and the column columnIndex.
                 
                 :param rowIndex: Index of the section/timestep to query
                 :param headerKey: Name of the header section to query  
                 :param columnIndex: Index of the column within the header
                 :return: Vector of numeric values for the specified location
        """
    def get_string_values_for_at(self, rowIndex: typing.SupportsInt, headerKey: str, columnIndex: str) -> list[str]:
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
class ExitReason:
    """
    
    An enum representing the reason for exiting
    the simulation or optimization procedure.
    
    Members:
    
      UNSET : Exit reason: "Unset".
    
      F_TOLERANCE : Exit reason: "F (force) tolerance reached".
    
      X_TOLERANCE : Exit reason: "X (displacement) tolerance reached".
    
      MAX_STEPS : Exit reason: "Maximum number of steps reached".
    
      NO_STEPS_POSSIBLE : Exit reason: "No (more) steps possible".
    
      FAILURE : Exit reason: "Failure".
    
      INTERRUPT : Exit reason: "Interrupt".
    
      OTHER : Exit reason: "Other".
    """
    FAILURE: typing.ClassVar[ExitReason]  # value = <ExitReason.FAILURE: 5>
    F_TOLERANCE: typing.ClassVar[ExitReason]  # value = <ExitReason.F_TOLERANCE: 1>
    INTERRUPT: typing.ClassVar[ExitReason]  # value = <ExitReason.INTERRUPT: 6>
    MAX_STEPS: typing.ClassVar[ExitReason]  # value = <ExitReason.MAX_STEPS: 3>
    NO_STEPS_POSSIBLE: typing.ClassVar[ExitReason]  # value = <ExitReason.NO_STEPS_POSSIBLE: 4>
    OTHER: typing.ClassVar[ExitReason]  # value = <ExitReason.OTHER: 7>
    UNSET: typing.ClassVar[ExitReason]  # value = <ExitReason.UNSET: 0>
    X_TOLERANCE: typing.ClassVar[ExitReason]  # value = <ExitReason.X_TOLERANCE: 2>
    __members__: typing.ClassVar[dict[str, ExitReason]]  # value = {'UNSET': <ExitReason.UNSET: 0>, 'F_TOLERANCE': <ExitReason.F_TOLERANCE: 1>, 'X_TOLERANCE': <ExitReason.X_TOLERANCE: 2>, 'MAX_STEPS': <ExitReason.MAX_STEPS: 3>, 'NO_STEPS_POSSIBLE': <ExitReason.NO_STEPS_POSSIBLE: 4>, 'FAILURE': <ExitReason.FAILURE: 5>, 'INTERRUPT': <ExitReason.INTERRUPT: 6>, 'OTHER': <ExitReason.OTHER: 7>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
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
    def __init__(self, max_distance_multiplier: typing.SupportsFloat) -> None:
        ...
    def get_max_distance(self, N: typing.SupportsFloat) -> float:
        """
                 Get the maximum distance for a given N.
        
                 :param N: Number of segments.
        """
class LinkSwappingMode:
    """
    How slip-links may act when they reach each-other or even a crosslink.
    
    Members:
    
      NO_SWAPPING : No Swapping
    
      SLIPLINKS_ONLY : Sliplinks only
    
      ALL : All
    
      ALL_CYCLE : All, restrict to cycles
    
      ALL_MC : All MC
    
      ALL_MC_CYCLE : All MC, restrict to cycles
    
      ALL_MC_TRY : All MC, attempt the move
    
      ALL_MC_TRY_CYCLE : All MC, restrict to cycles, attempt the move
    """
    ALL: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL: 2>
    ALL_CYCLE: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL_CYCLE: 3>
    ALL_MC: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL_MC: 4>
    ALL_MC_CYCLE: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL_MC_CYCLE: 5>
    ALL_MC_TRY: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL_MC_TRY: 6>
    ALL_MC_TRY_CYCLE: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.ALL_MC_TRY_CYCLE: 7>
    NO_SWAPPING: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.NO_SWAPPING: 0>
    SLIPLINKS_ONLY: typing.ClassVar[LinkSwappingMode]  # value = <LinkSwappingMode.SLIPLINKS_ONLY: 1>
    __members__: typing.ClassVar[dict[str, LinkSwappingMode]]  # value = {'NO_SWAPPING': <LinkSwappingMode.NO_SWAPPING: 0>, 'SLIPLINKS_ONLY': <LinkSwappingMode.SLIPLINKS_ONLY: 1>, 'ALL': <LinkSwappingMode.ALL: 2>, 'ALL_CYCLE': <LinkSwappingMode.ALL_CYCLE: 3>, 'ALL_MC': <LinkSwappingMode.ALL_MC: 4>, 'ALL_MC_CYCLE': <LinkSwappingMode.ALL_MC_CYCLE: 5>, 'ALL_MC_TRY': <LinkSwappingMode.ALL_MC_TRY: 6>, 'ALL_MC_TRY_CYCLE': <LinkSwappingMode.ALL_MC_TRY_CYCLE: 7>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
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
    
           Please cite :cite:t:`gusev_molecular_2024` and/or :cite:t:`bernhard_phantom_2025` if you use this method in your work.
      
    """
    def __copy__(self) -> MCUniverseGenerator:
        ...
    def __getstate__(self) -> tuple:
        ...
    def __init__(self, lx: typing.SupportsFloat = 10.0, ly: typing.SupportsFloat = 10.0, lz: typing.SupportsFloat = 10.0) -> None:
        """
                 Initialize a new MCUniverseGenerator with specified box dimensions.
        
                 :param lx: Box length in x-direction (default: 10.0).
                 :param ly: Box length in y-direction (default: 10.0).
                 :param lz: Box length in z-direction (default: 10.0).
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def add_crosslinkers(self, nr_of_crosslinkers: typing.SupportsInt, crosslinker_functionality: typing.SupportsInt = 4, crosslinker_type: typing.SupportsInt = 2, white_noise: bool = True) -> None:
        """
                    Add crosslinkers at random positions.
        
                    :param nr_of_crosslinkers: Number of crosslinkers to add.
                    :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
                    :param crosslinker_type: Atom type of the crosslinkers (default: 2).
                    :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).
        """
    def add_crosslinkers_at(self, coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], crosslinker_functionality: typing.SupportsInt = 4, crosslinker_type: typing.SupportsInt = 2) -> None:
        """
                  Add crosslinkers at specific coordinates.
        
                  :param coordinates: Coordinates of the crosslinkers as a flat array [x1, y1, z1, x2, y2, z2, ...].
                  :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
                  :param crosslinker_type: Atom type of the crosslinkers (default: 2).
        """
    def add_end_functionalized_strands(self, nr_of_strands: typing.SupportsInt, strand_length: collections.abc.Sequence[typing.SupportsInt], crosslinker_functionality: typing.SupportsInt = 4, crosslinker_type: typing.SupportsInt = 2, strand_atom_type: typing.SupportsInt = 1, white_noise: bool = True) -> None:
        """
                    Add strands with crosslinkers at the ends.
        
                    :param nr_of_strands: Number of strands to add.
                    :param strand_length: Length of each strand.
                    :param crosslinker_functionality: Functionality of the crosslinkers (default: 4).
                    :param crosslinker_type: Atom type of the crosslinkers (default: 2).
                    :param strand_atom_type: Atom type of the beads that are not at the ends (default: 1).
                    :param white_noise: Whether to use white noise (true) or blue noise (false) for the positions of the crosslinkers (default: true).
        """
    def add_functionalized_strands(self, nr_of_strands: typing.SupportsInt, strand_length: collections.abc.Sequence[typing.SupportsInt], crosslink_selector: collections.abc.Callable[[typing.SupportsInt, typing.SupportsInt, typing.SupportsInt], tuple[bool, int]], default_crosslinker_functionality: typing.SupportsInt, crosslinker_type: typing.SupportsInt = 2, strand_atom_type: typing.SupportsInt = 1, white_noise: bool = True) -> None:
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
    def add_monofunctional_strands(self, nr_of_monofunctional_strands: typing.SupportsInt, monofunctional_strand_length: collections.abc.Sequence[typing.SupportsInt], monofunctional_strand_atom_type: typing.SupportsInt = 4) -> None:
        """
                 Add multiple monofunctional strands with specified bead types.
        
                 :param nr_of_monofunctional_strands: Number of monofunctional strands to add.
                 :param monofunctional_strand_length: Vector specifying the length of each strand in beads.
                 :param monofunctional_strand_atom_type: Atom type for the strand beads (default: 4).
        """
    def add_randomly_functionalized_strands(self, nr_of_strands: typing.SupportsInt, strand_length: collections.abc.Sequence[typing.SupportsInt], functionalization_probability: typing.SupportsFloat, crosslinker_functionality: typing.SupportsInt = 4, crosslinker_type: typing.SupportsInt = 2, strand_atom_type: typing.SupportsInt = 1, white_noise: bool = True) -> None:
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
    def add_regularly_spaced_functionalized_strands(self, nr_of_strands: typing.SupportsInt, strand_length: collections.abc.Sequence[typing.SupportsInt], spacing_between_crosslinks: typing.SupportsInt, offset_to_first_crosslink: typing.SupportsInt = 0, crosslinker_functionality: typing.SupportsInt = 4, crosslinker_type: typing.SupportsInt = 2, strand_atom_type: typing.SupportsInt = 1, white_noise: bool = True) -> None:
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
    def add_solvent_chains(self, nr_of_solvent_chains: typing.SupportsInt, solvent_chain_length: typing.SupportsInt, solvent_atom_type: typing.SupportsInt = 3, white_noise: bool = True) -> None:
        """
                    Randomly distribute additional, free chains.
        
                    :param nr_of_solvent_chains: Number of solvent chains to add.
                    :param solvent_chain_length: Length of each solvent chain in beads.
                    :param solvent_atom_type: Atom type for solvent chain beads (default: 3).
                    :param white_noise: Whether to use white noise (true) or blue noise (false) for positioning (default: true).
        """
    def add_star_crosslinkers(self, nr_of_stars: typing.SupportsInt, functionality: typing.SupportsInt, beads_per_strand: typing.SupportsInt, crosslinker_atom_type: typing.SupportsInt = 2, strand_atom_type: typing.SupportsInt = 1, white_noise: bool = True) -> None:
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
    def add_strands(self, nr_of_strands: typing.SupportsInt, strand_lengths: collections.abc.Sequence[typing.SupportsInt], strand_atom_type: typing.SupportsInt = 1) -> None:
        """
                    Add strands.
                    Adds them unconnected at first.
                    To link them to crosslinkers, use some of the :code:`link_strand_` methods.
        
                    :param nr_of_strands: Number of strands to add.
                    :param strand_lengths: A list of integers representing the number of beads of each of the strands.
                    :param strand_atom_type: Type of atoms for the strands (default: 1).
        """
    def config_nr_of_mc_steps(self, n_steps: typing.SupportsInt = 2000) -> None:
        """
               Set the number of Monte-Carlo steps during bond length equilibration.
        
               :param n_steps: Number of MC steps to perform (default: 2000)
        """
    def config_primary_loop_probability(self, probability: typing.SupportsFloat = 1.0) -> None:
        """
                 Configure an additional weight reduction for the primary loop formation.
        
                 Defaults to 1., which means the general :math:`P(\\vec{R})` is used without any bias.
                 This results in more primary loops for shorter chains than for longer ones.
        
                 Set to 0. to disable the formation of primary loops.
        """
    def config_secondary_loop_probability(self, probability: typing.SupportsFloat = 1.0) -> None:
        """
                 Configure an additional weight reduction for the secondary loop formation.
        
                 Defaults to 1., which means the general :math:`P(\\vec{R})` is used without any bias.
                 This results in more secondary loops for shorter chains than for longer ones.
        
                 Set to 0. to disable the formation of secondary loops.
        """
    def disable_max_distance(self) -> None:
        """
                 Disable the maximum distance provider to allow unlimited distance sampling.
                 This may slow down performance for large systems but ensures complete sampling.
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
    def get_force_balance(self) -> ...:
        """
                 Get an instance of the force balance procedure.
                 This is a useful shorthand e.g. to skip the sampling of beads within in the strands.
        
                 :return: Configured MEHPForceBalance instance.
        """
    def get_force_balance2(self) -> ...:
        """
                      Get an instance of the force balance procedure.
                      This is a useful shorthand e.g. to skip the sampling of beads within in the strands.
        
                      :return: Configured MEHPForceBalance2 instance.
        """
    def get_force_relaxation(self) -> ...:
        """
                 Get an instance of the force relaxation procedure.
                 This is a useful shorthand e.g. to skip the sampling of beads within in the strands.
        
                 :return: Configured MEHPForceRelaxation instance.
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
    def link_strand(self, strand_idx: typing.SupportsInt, c_infinity: typing.SupportsFloat = 1.0) -> None:
        """
                  Link one particular strand to a crosslinker.
                  This strand will have one bond made, to an appropriate crosslinker, 
                  as chosen by the parameters associated with the strand.
        
                  :param strand_idx: Index of the strand to be linked.
                  :param c_infinity: As needed for the end-to-end distribution, given by :math:`\\langle R^2\\rangle_0 = C_{\\infty} N b^2` (default: 1.0).
        """
    def link_strand_to(self, strand_idx: typing.SupportsInt, link_idx: typing.SupportsInt) -> None:
        """
                  Link a strand to a specific crosslinker.
                  This assumes that you keep track of the order in which you added the crosslinkers and strands,
                  as those will determine the indices.
        
                  .. caution::
        
                       Be aware that some few methods may change the cross-link or strand indices.
        
                  :param strand_idx: Index of the strand to be linked.
                  :param link_idx: Index of the crosslinker to be linked.
        """
    def link_strand_to_strand(self, c_infinity: typing.SupportsFloat = 1.0) -> bool:
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
    def link_strands_callback(self, linking_controller: collections.abc.Callable[[MCUniverseGenerator, typing.SupportsInt], BackTrackStatus], c_infinity: typing.SupportsFloat = 1.0) -> None:
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
    def link_strands_to_conversion(self, crosslinker_conversion: typing.SupportsFloat, c_infinity: typing.SupportsFloat = 1.0) -> None:
        """
                    Actually link the previously added strands to the previously added crosslinkers,
                    until a certain crosslink conversion is reached.
        
                    :param crosslinker_conversion: Target conversion of crosslinkers (0: no connections to crosslinks; 1: all crosslinkers fully connected).
                    :param c_infinity: As needed for the end-to-end distribution, given by :math:`\\langle R^2\\rangle_0 = C_{\\infty} N b^2`.
        """
    def link_strands_to_soluble_fraction(self, soluble_fraction: typing.SupportsFloat, c_infinity: typing.SupportsFloat = 1.0) -> None:
        """
                    Actually link the previously added strands to the previously added crosslinkers,
                    until a certain soluble fraction is reached.
        
                    :param soluble_fraction: Target soluble fraction (0: all material is soluble; 1: no material is soluble).
                    :param c_infinity: As needed for the end-to-end distribution, given by :math:`\\langle R^2\\rangle_0 = C_{\\infty} N b^2` (default: 1.0).
        """
    def link_strands_to_strands_to_conversion(self, target_strand_conversion: typing.SupportsFloat, c_infinity: typing.SupportsFloat = 1.0) -> None:
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
    def set_bead_distance(self, distance: typing.SupportsFloat, update_mean_squared: bool = True) -> None:
        """
                 Set the mean distance between beads when doing MC stepping.
                 Also used for the target crosslinker partner sampling.
        
                 NOTE: Mainly the mean squared bead distance is effectively used in the Monte-Carlo simulation.
        
                 :param distance: Mean distance between beads.
                 :param update_mean_squared: Whether to update the mean squared distance as well, deduced from the assumed gaussian distribution in 3D (default: true).
        """
    def set_mean_squared_bead_distance(self, mean_squared_distance: typing.SupportsFloat, update_mean: bool = True) -> None:
        """
                  Set the mean squared distance between beads.
        
                  :param mean_squared_distance: Mean squared distance between beads
                  :param update_mean: Whether to update the mean bead distance as well, deduced from the assumed gaussian distribution in 3D (default: true)
        """
    def set_seed(self, seed: typing.SupportsInt) -> None:
        """
                 Set the seed for the random generator.
        
                 :param seed: Random seed value for reproducible results.
        """
    def use_linear_max_distance(self, multiplier: typing.SupportsFloat) -> None:
        """
                  Converts the :math:`N` to a maximum distance within which to sample.
                  The distance will be calculated as :math:`N \\times \\text{multiplier}`.
                  For example, commonly, a maximum distance is :math:`N \\times <b>`, 
                  in which case the multiplier is :math:`<b>`.
        
                  Useful only for performance improvements in large systems.
        """
    def use_zscore_max_distance(self, std_multiplier: typing.SupportsFloat, in_sqrt_multiplier: typing.SupportsFloat = 1.0) -> None:
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
class MEHPForceBalance:
    """
    
        A small simulation tool for quickly minimizing the force between the crosslinker beads.
    
        This is the second implementation in the group of MEHP provided by this package.
        The distinct feature here is the slip-links: a form of entanglement,
        represented as an entanglement link, just like a four-functional crosslink,
        but with the ability to slip along the two associated strands, therewith
        adjusting the fraction of the contour length on both sides of the link.
    
        Please cite :cite:t:`bernhard_phantom_2025` if you use this method.
        
    """
    @staticmethod
    def construct_with_random_sliplinks(universe: Universe, nr_of_sliplinks_to_sample: typing.SupportsInt, upper_sampling_cutoff: typing.SupportsFloat = 1.2, lower_sampling_cutoff: typing.SupportsFloat = 0.0, minimum_nr_of_sliplinks: typing.SupportsInt = 0, same_strand_cutoff: typing.SupportsFloat = 3, seed: str = '', crosslinker_type: typing.SupportsInt = 2, is_2d: bool = False, skip_dangling_soluble_entanglements: bool = True) -> MEHPForceBalance:
        """
                  Instantiate this simulator with randomly chosen slip-links.
        """
    @staticmethod
    def get_neighbour_link_indices(network: SimplifiedBalanceNetwork, link_idx: typing.SupportsInt) -> list[int]:
        """
                         Get the indices of neighboring links for a given link.
        
                         :param network: The force balance network
                         :param link_idx: Index of the link to find neighbors for
                         :return: Vector of connected link indices
        """
    def __copy__(self) -> MEHPForceBalance:
        ...
    def __getstate__(self) -> tuple:
        ...
    def __init__(self, universe: Universe, crosslinker_type: typing.SupportsInt = 2, is_2d: bool = False, remove_2functional_crosslinkers: bool = False, remove_dangling_chains: bool = False) -> None:
        """
                  Instantiate the simulator for a certain universe.
        
                  :param universe: The universe to simulate with
                  :param crosslinker_type: The atom type of the crosslinkers. Needed to reduce the network.
                  :param is2D: Whether to ignore the z direction.
                  :param kappa: The spring constant
                  :param remove_2functionalCrosslinkers: whether to keep or remove the 2-functional crosslinkers when setting up the network
                  :param remove_dangling_chains: whether to keep or remove obviously dangling chains when setting up the network
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def add_sliplinks(self, strand_idx_1: collections.abc.Sequence[typing.SupportsInt], strand_idx_2: collections.abc.Sequence[typing.SupportsInt], x: collections.abc.Sequence[typing.SupportsFloat], y: collections.abc.Sequence[typing.SupportsFloat], z: collections.abc.Sequence[typing.SupportsFloat], alpha_1: collections.abc.Sequence[typing.SupportsFloat], alpha_2: collections.abc.Sequence[typing.SupportsFloat], clamp_alpha: bool = False) -> None:
        """
                  Add new slip-links
        """
    def add_sliplinks_based_on_cycles(self, maxLoopLength: typing.SupportsInt = -1) -> int:
        """
                  Detect and add slip-links based on detected entanglements.
        
                  WARNING:
                       Does not work yet.
        """
    def config_assume_box_large_enough(self, box_large_enough: bool = False) -> None:
        """
                  Configure whether to run PBC on the bonds or not.
        
                  :param box_large_enough: If True, assume the box is large enough and don't apply PBC on bonds.
                                          If your bonds could get larger than half the box length, 
                                          this must be kept False (default).
        """
    def config_entanglement_type(self, type: typing.SupportsInt = -1) -> None:
        """
                 To have certain crosslinks behave as entanglements in the removal process,
                 you can specify here a type, that you have used in the universe to specify:
                 - the type of entanglement atoms (expected with functionality f = 3),
                 - and the entanglement-bonds between the entanglement atoms.
        
                 I.e., say you want to model some entanglements as non-slipping,
                 bonds between two strand beads resulting in f = 3 beads, for example,
                 you can call this method to have the "StructureSimplificationMode" also remove these atoms,
                 if they have a functionality of 2 or less while still being connected to its partner bead.
        
                 :param type: The atom type to treat as entanglement.
        """
    def config_mean_bond_length(self, b: typing.SupportsFloat = 1.0) -> None:
        """
                  Configure the :math:`b` used e.g. for the topological Gamma-factor.
        
                  :param b: The mean bond length to use.
        """
    def config_simplification_frequency(self, frequency: typing.SupportsInt = 10) -> None:
        """
                  Configure every how many steps to simplify the structure.
        
                  :param frequency: The number of steps between simplification. Default: 10.
        """
    def config_spring_breaking_distance(self, distance_over_contour_length: typing.SupportsFloat = -1) -> None:
        """
                  Configure the "force" (distance over contour length) at which the bonds break.
                  Can be used to model the effect of fracture, to reduce the stiffening happening upon deformation.
                  Springs breaking will happen before the simplification procedure is run.
                  Negative values will disable spring breaking.
                  Default: -1.
        
                  :param distance_over_contour_length: The threshold for breaking springs.
        """
    def config_spring_constant(self, kappa: typing.SupportsFloat = 1.0) -> None:
        """
                  Configure the spring constant used in the simulation.
        
                  :param kappa: The spring constant.
        """
    def config_step_output(self, output_configuration: collections.abc.Sequence[OutputConfiguration]) -> None:
        """
                  Set which values to log during the simulation.
        
                  :param output_configuration: A list of OutputConfiguration structs
                                              specifying what values to log and how often
        """
    def deform_to(self, new_box: Box) -> None:
        """
                  Perform a deformation of the system box to a different box.
                  All coordinates etc. will be scaled as needed.
        """
    def evaluate_spring_distance(self, network: SimplifiedBalanceNetwork, displacements: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], spring_idx: typing.SupportsInt) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
               Evaluate the distance vector for a specific spring.
        
               :param network: The force balance network
               :param displacements: Current displacement vector
               :param spring_idx: Index of the spring to evaluate
               :return: 3D distance vector for the spring
        """
    def evaluate_spring_distance_from(self, network: SimplifiedBalanceNetwork, displacements: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], spring_idx: typing.SupportsInt, link_idx: typing.SupportsInt) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
               Evaluate the spring distance from a specific link.
        
               :param network: The force balance network
               :param displacements: Current displacement vector
               :param spring_idx: Index of the spring to evaluate
               :param link_idx: Index of the starting link
               :return: 3D distance vector from the specified link
        """
    def evaluate_spring_distance_to(self, network: SimplifiedBalanceNetwork, displacements: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], spring_idx: typing.SupportsInt, link_idx: typing.SupportsInt) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
               Evaluate the spring distance to a specific link.
        
               :param network: The force balance network
               :param displacements: Current displacement vector
               :param spring_idx: Index of the spring to evaluate
               :param link_idx: Index of the target link
               :return: 3D distance vector to the specified link
        """
    def get_average_strand_length(self) -> float:
        """
                   Get the average length of the strands. Note that in contrast to :meth:`~pylimer_tools_cpp.MEHPForceBalance.get_gamma_factor`,
                   this value is normalized by the number of strands rather than the number of chains.
        """
    def get_coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                     Get the current coordinates of the crosslinkers and entanglement links.
        """
    def get_crosslinker_universe(self) -> Universe:
        """
                   Returns the universe [of crosslinkers] with the positions of the current state of the simulation.
        """
    def get_current_spring_lengths(self) -> list[float]:
        """
                  Get the spring distances.
        """
    def get_current_spring_vectors(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Get the spring vectors.
        """
    def get_current_strand_vectors(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    def get_dangling_weight_fraction(self, tolerance: typing.SupportsFloat = 0.001) -> float:
        """
                  Compute the weight fraction of non-active strands
        
                  Caution: ignores atom masses.
        """
    def get_default_mean_bond_length(self) -> float:
        """
                   Returns the value effectively used in :meth:`~pylimer_tools_cpp.MEHPForceBalance.get_gamma_factor` for
                   :math:`b` in :math:`\\langle R_{0,\\eta}^2 = N_{\\eta} b^2\\rangle`.
        """
    def get_displacement_residual_norm(self, one_over_strand_partition_upper_limit: typing.SupportsFloat = 1.0) -> float:
        """
                   Get the current link displacement residual norm.
        """
    def get_displacements(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                   Get the current link displacements.
        """
    def get_effective_functionality_of_atoms(self, tolerance: typing.SupportsFloat = 0.001) -> dict[int, int]:
        """
                  Returns the number of active strands connected to each atom, atomId used as index
        
                  :param tolerance: strands under this length are considered inactive
        """
    def get_exit_reason(self) -> ExitReason:
        """
                   Returns the reason for termination of the simulation
        """
    def get_force_magnitude_vector(self, arg0: typing.SupportsFloat) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Evaluate the norm of the force on each (slip- or cross-) link.
        """
    def get_force_on(self, link_idx: typing.SupportsInt, one_over_strand_partition_upper_limit: typing.SupportsFloat = 1.0) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                  Evaluate the force on a particular (slip- or cross-) link.
        """
    def get_gamma_factor(self, b02: typing.SupportsFloat = -1.0, nr_of_chains: typing.SupportsInt = -1, one_over_strand_partition_upper_limit: typing.SupportsFloat = 1.0) -> float:
        """
                  Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:
        
                  :math:`\\Gamma = \\langle\\gamma_{\\eta}\\rangle`, with :math:`\\gamma_{\\eta} = \\\\frac{\\bar{r_{\\eta}}^2}{R_{0,\\eta}^2}`,
                  which you can use as :math:`G_{\\mathrm{ANT}} = \\Gamma \\nu k_B T`,
                  where :math:`\\eta` is the index of a particular strand,
                  :math:`R_{0}^2` is the melt mean square end to end distance, in phantom systems :math:`$= N_{\\eta}*b^2$`
                  :math:`N_{\\eta}` is the number of atoms in this strand :math:`\\eta`,
                  :math:`b` its mean square bond length,
                  :math:`\\nu` the number density of network strands,
                  :math:`T` the temperature and
                  :math:`k_B` Boltzmann's constant.
        
                  :param b02: The melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\\infinity b^2`.
                  :param nr_of_chains: The value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of springs.
        """
    def get_gamma_factors(self, b02: typing.SupportsFloat, one_over_strand_partition_upper_limit: typing.SupportsFloat = 1.0) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Evaluates the gamma factor for each strand (i.e., the squared distance divided by the contour length multiplied by b02)
        """
    def get_gamma_factors_in_dir(self, b02: typing.SupportsFloat, direction: typing.SupportsInt, one_over_strand_partition_upper_limit: typing.SupportsFloat = 1.0) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                       Evaluates the gamma factor for each strand in the specified direction (i.e., the squared distance divided by the contour length multiplied by b02)
        
                       :param b02: The melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\\infinity b^2`.
                       :param direction: The direction in which to compute the gamma factors (0: x, 1: y, 2: z)
        """
    def get_ids_of_active_nodes(self, tolerance: typing.SupportsFloat = 0.001) -> list[int]:
        """
                  Get the atom ids of the nodes that are considered active.
                  Only crosslink ids are returned (not e.g. entanglement links).
        
                  :param tolerance: strands under this length are considered inactive. A node is active if it has > 1 active strands.
        """
    def get_initial_coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                     Get the initial coordinates of the (remaining) crosslinkers and entanglement links.
        """
    def get_nr_of_active_nodes(self, tolerance: typing.SupportsFloat = 0.001) -> int:
        """
                  Get the number of active nodes (incl. entanglement nodes [atoms with type = entanglementType, present in the universe when creating this simulator],
                  excl. entanglement links [the slip-links created internally when e.g. constructing the simulator with random slip-links]).
        
                  :param tolerance: strands under this length are considered inactive. A node is active if it has > 1 active strands.
        """
    def get_nr_of_active_springs(self, tolerance: typing.SupportsFloat = 0.001) -> int:
        """
                   Get the number of active springs remaining after running the simulation.
        
                  :param tolerance: springs under this length are considered inactive
        """
    def get_nr_of_active_strands(self, tolerance: typing.SupportsFloat = 0.001) -> int:
        """
                   Get the number of active strands remaining after running the simulation.
        
                  :param tolerance: strands under this length are considered inactive
        """
    def get_nr_of_active_strands_in_dir(self, direction: typing.SupportsInt, tolerance: typing.SupportsFloat = 0.001) -> int:
        """
                        Get the number of active strands remaining after running the simulation.
        
                       :param direction: The direction in which to compute the active strands (0: x, 1: y, 2: z)
                       :param tolerance: strands under this length are considered inactive
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
    def get_nr_of_strands(self) -> int:
        """
                  Get the number of strands considered in this simulation.
        """
    def get_overall_strand_lengths(self) -> list[float]:
        """
                  Get the sum of the lengths of the springs of each strand.
        """
    def get_pressure(self) -> float:
        """
                  Returns the pressure at the current state of the simulation.
        """
    def get_soluble_weight_fraction(self, tolerance: typing.SupportsFloat = 0.001) -> float:
        """
                  Compute the weight fraction of strands connected to active
                  strands (any depth).
        
                  Caution: ignores atom masses.
        """
    def get_strand_partition_indices_of_sliplink(self, network: SimplifiedBalanceNetwork, link_idx: typing.SupportsInt) -> list[int]:
        """
                  Get the indices of strand partitions associated with a slip-link.
        
                  :param network: The force balance network
                  :param link_idx: Index of the slip-link
                  :return: Vector of strand partition indices
        """
    def get_strand_partitions(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Get the current strand partitions (the fraction of the contour length associated with each spring).
        """
    def get_stress_on(self, link_idx: typing.SupportsInt, one_over_strand_partition_upper_limit: typing.SupportsFloat = 1.0) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 3]"]:
        """
                  Evaluate the stress on a particular (slip- or cross-) link.
        """
    def get_stress_tensor(self, one_over_strand_partition_upper_limit: typing.SupportsFloat = 1.0) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 3]"]:
        """
                  Returns the stress tensor at the current state of the simulation.
        
                  The units are in :math:`[\\text{units of }\\kappa]/[\\text{distance units}]`,
                  where the units of :math:`\\kappa` should be :math:`[\\text{force}]/[\\text{distance units}]^2`.
                  Make sure to multiply by :math:`\\kappa` or configure it appropriately.
        """
    def get_stress_tensor_link_based(self, one_over_strand_partition_upper_limit: typing.SupportsFloat = 1.0, xlinks_only: bool = False) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 3]"]:
        """
                  Returns the stress tensor at the current state of the simulation.
        """
    def get_weighted_spring_lengths(self, one_over_strand_partition_upper_limit: typing.SupportsFloat = 1.0) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Get the current spring lengths (norm of vector) divided by the strand partition times the contour length.
        """
    def inspect_displacement_to_mean_position_update(self, link_idx: typing.SupportsInt, one_over_strand_partition_upper_limit: typing.SupportsFloat = 1.0) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Helper method to debug and/or understand what happens to certain links when being displaced.
        """
    def inspect_parametrisation_optimsation_for_link(self, link_idx: typing.SupportsInt, displacements: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], strand_partitions: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], max_nr_of_steps: typing.SupportsInt = 100, alpha_tol: typing.SupportsFloat = 1e-09, min_nr_of_steps: typing.SupportsInt = 1, one_over_strand_partition_upper_limit: typing.SupportsFloat = 1.0) -> tuple[typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"], typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"], int, float, float, float, float]:
        """
                  Helper method to debug and/or understand what happens to certain links
                  when being displaced and their partition updated.
        """
    def inspect_strand_partition_update(self, link_idx: typing.SupportsInt) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Helper method to debug and/or understand what happens to certain links
                  when the strand partition is being updated.
        """
    def move_sliplinks_to_their_best_branch(self, arg0: SimplifiedBalanceNetwork, arg1: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], arg2: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], arg3: typing.SupportsFloat, arg4: typing.SupportsInt, arg5: bool, arg6: bool) -> None:
        ...
    def randomly_add_sliplinks(self, nr_of_sliplinks_to_sample: typing.SupportsInt, cutoff: typing.SupportsFloat = 2.0, minimum_nr_of_sliplinks: typing.SupportsInt = 0, same_strand_cutoff: typing.SupportsFloat = 2, exclude_crosslinks: bool = False, seed: typing.SupportsInt = -1) -> int:
        """
                  Randomly sample and add slip-links based on certain criteria.
        """
    def run_force_relaxation(self, max_nr_of_steps: typing.SupportsInt = 250000, x_tolerance: typing.SupportsFloat = 1e-12, initial_residual_norm: typing.SupportsFloat = -1.0, simplification_mode: StructureSimplificationMode = ..., inactive_removal_cutoff: typing.SupportsFloat = 0.001, do_inner_iterations: bool = False, allow_sliplinks_to_pass_each_other: LinkSwappingMode = ..., swapping_frequency: typing.SupportsInt = 10, one_over_strand_partition_upper_limit: typing.SupportsFloat = 1.0, nr_of_crosslink_swaps_allowed_per_sliplink: typing.SupportsInt = -1, disable_slipping: bool = False) -> None:
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
                  :param one_over_strand_partition_upper_limit: Super-secret parameter. Use 1, gradually increase (and then -1) if you want to publish.
                  :param nr_of_crosslink_swaps_allowed_per_sliplink: Use to steer whether slip-links can cross crosslinks when swapping is enabled.
                  :param disable_slipping: Whether slip-links should be prohibited from slipping.
        """
    def set_displacements(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        """
                   Set the current link displacements.
        """
    def set_strand_contour_lengths(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        """
                   Set/overwrite the contour lengths.
        """
    def set_strand_partitions(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        """
                  Set the current strand partitions.
        """
    def swap_sliplinks_incl_xlinks(self, arg0: SimplifiedBalanceNetwork, arg1: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], arg2: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], arg3: typing.SupportsFloat, arg4: bool) -> None:
        ...
    def validate_network(self) -> bool:
        """
                  Validates the internal structures.
        
                  Throws an error if something is not ok.
                  Otherwise, it returns true.
        
                  Can be used e.g. as :code:`assert fb.validate_network()`.
        
                  :return: True if validation passes, raises exception otherwise
        """
    @property
    def network(self) -> SimplifiedBalanceNetwork:
        ...
class MEHPForceBalance2:
    """
    
         A small simulation tool for quickly minimizing the force between the crosslinker beads.
    
         This is the third implementation of the MEHP. 
         It's the fastest implementation by using a simple spring model only, disabling the non-linear
         periodic boundary conditions, and instead builds a sparse linear system of equations that's readily solved.
         However, it allows entanglements to be represented as additional links or/and springs,
         although without slipping along the chain.
    
         Please cite :cite:t:`bernhard_phantom_2025` if you use this method.
          
    """
    @staticmethod
    def get_neighbour_link_indices(network: SimplifiedBalance2Network, link_idx: typing.SupportsInt) -> list[int]:
        """
                         Get the indices of neighboring links for a given link.
        
                         :param network: The force balance 2 network
                         :param link_idx: Index of the link to find neighbors for
                         :return: Vector of connected link indices
        """
    def __copy__(self) -> MEHPForceBalance2:
        ...
    def __getstate__(self) -> tuple:
        ...
    @typing.overload
    def __init__(self, universe: Universe, crosslinker_type: typing.SupportsInt = 2, is_2d: bool = False) -> None:
        """
                   Instantiate the simulator for a certain universe.
        
                   :param universe: The universe giving the basic connectivity to compute with
                   :param crosslinker_type: The atom type of the crosslinkers. Needed to reduce the network.
                   :param is_2d: Whether to ignore the z direction.
                   :param kappa: The spring constant
                   :param remove_2functionalCrosslinkers: whether to keep or remove the 2-functional crosslinkers when setting up the network
                   :param remove_dangling_chains: whether to keep or remove obviously dangling chains when setting up the network
        """
    @typing.overload
    def __init__(self, universe: Universe, entanglements: AtomPairEntanglements, crosslinker_type: typing.SupportsInt = 2, is_2d: bool = False, entanglements_as_springs: bool = False) -> None:
        """
                   Instantiate the simulator for a certain universe with the given entanglements.
        
                   :param universe: The universe giving the basic connectivity to compute with
                   :param entanglements: The entanglements to use in the computation
                   :param crosslinker_type: The atom type of the crosslinkers. Needed to reduce the network.
                   :param is_2d: Whether to ignore the z direction.
                   :param entanglements_as_springs: whether to use the entanglements as springs instead of links
        """
    @typing.overload
    def __init__(self, universe: Universe, nr_of_entanglements_to_sample: typing.SupportsInt, upper_sampling_cutoff: typing.SupportsFloat = 1.2, lower_sampling_cutoff: typing.SupportsFloat = 0.0, minimum_nr_of_sliplinks: typing.SupportsInt = 0, same_strand_cutoff: typing.SupportsFloat = 3, seed: str = '', crosslinker_type: typing.SupportsInt = 2, is_2d: bool = False, skip_dangling_soluble_entanglements: bool = True, entanglements_as_springs: bool = True) -> None:
        """
             Instantiate this simulator with randomly chosen slip-links.
        
             :param universe: The universe containing the basic atoms and connectivity
             :param nr_of_entanglements_to_sample: The number of entanglements to sample
             :param upper_cutoff: maximum distance from one sampled bead to its partner
             :param lower_cutoff: minimum distance from one sampled bead to its partner
             :param minimum_nr_of_entanglements: The minimum number of entanglements that should be sampled
             :param same_strand_cutoff: distance from one sampled bead to its pair within the same strand
             :param seed: The seed for the random number generator
             :param cross_linker_type:
             :param is_2d:
             :param filter_entanglements:
             :param entanglements_as_springs: whether to model the entanglements as merged beads or beads with 1 spring in between
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def config_assume_network_is_complete(self, assume_network_is_complete: bool = False) -> None:
        """
                   Configure whether the network is assumed to be complete.
        
                   This assumption means, that the `universe` instance is not queried for clusters when
                   computing the fraction of e.g. the soluble or dangling strands.
                   Do not set this to true if inactive
                   (dangling or free) strands have been removed from the network.
        
                   This is useful to reduce memory between MC generator and force balance,
                   since the universe representation with all the in-between beads does not need to be stored.
                   However, currently, the removal procedures are not loss-free.
                   Therefore, use this with care.
        """
    def config_mean_bond_length(self, b: typing.SupportsFloat = 1.0) -> None:
        """
              Configure the :math:`b` used e.g. for the topological Gamma-factor.
        
              :param b: The mean bond length to use.
        """
    def config_spring_breaking_distance(self, distance_over_contour_length: typing.SupportsFloat = -1) -> None:
        """
                   Configure the "force" (distance over contour length) at which the bonds break.
                   Can be used to model the effect of fracture, to reduce the stiffening happening upon deformation.
                   Springs breaking will happen before the simplification procedure is run.
                   Negative values will disable spring breaking.
                   Default: -1.
        
                   :param distance_over_contour_length: The threshold for breaking springs.
        """
    def config_spring_constant(self, kappa: typing.SupportsFloat = 1.0) -> None:
        """
                   Configure the spring constant used in the simulation.
        
                   :param kappa: The spring constant.
        """
    def config_step_output(self, values: collections.abc.Sequence[OutputConfiguration]) -> None:
        """
                   Set which values to log.
        
                   :param values: a list of OutputConfiguration structs
        """
    def deform_to(self, new_box: Box) -> None:
        """
                   Perform a deformation of the system box to a different box.
                   All coordinates etc. will be scaled as needed.
        """
    def evaluate_spring_vector(self, network: SimplifiedBalance2Network, displacements: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], spring_idx: typing.SupportsInt) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
               Evaluate the distance vector for a specific spring.
        
               :param network: The force balance 2 network
               :param displacements: Current displacement vector
               :param spring_idx: Index of the spring to evaluate
               :return: 3D distance vector for the spring
        """
    def evaluate_spring_vector_from(self, network: SimplifiedBalance2Network, displacements: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], spring_idx: typing.SupportsInt, link_idx: typing.SupportsInt) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
               Evaluate the spring vector in the direction from a specific link.
        
               :param network: The force balance 2 network
               :param displacements: Current displacement vector
               :param spring_idx: Index of the spring to evaluate
               :param link_idx: Index of the starting link
               :return: 3D distance vector from the specified link
        """
    def evaluate_spring_vector_to(self, network: SimplifiedBalance2Network, displacements: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"], spring_idx: typing.SupportsInt, link_idx: typing.SupportsInt) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
               Evaluate the spring vector in the direction to a specific link.
        
               :param network: The force balance 2 network
               :param displacements: Current displacement vector
               :param spring_idx: Index of the spring to evaluate
               :param link_idx: Index of the target link
               :return: 3D distance vector to the specified link
        """
    def get_average_spring_length(self) -> float:
        """
                    Get the average length of the springs. Note that in contrast to :func:`~pylimer_tools_cpp.MEHPForceBalance2.getGammaFactor()`,
                    this value is normalized by the number of springs rather than the number of chains.
        """
    def get_coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                     Get the current coordinates of the crosslinkers and entanglement links.
        """
    def get_crosslinker_universe(self) -> Universe:
        """
                   Returns the universe [of crosslinkers] with the positions of the current state of the simulation.
        """
    def get_current_spring_lengths(self) -> list[float]:
        """
                   Get the spring lengths (Euclidean distances of the spring vectors).
        """
    def get_current_spring_vectors(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                   Get the spring vectors.
        """
    def get_dangling_weight_fraction(self, tolerance: typing.SupportsFloat = 0.001) -> float:
        """
                   Compute the weight fraction of non-active springs
        
                   Caution: ignores atom masses.
        """
    def get_default_mean_bond_length(self) -> float:
        """
                    Returns the value effectively used in :func:`~pylimer_tools_cpp.MEHPForceBalance2.getGammaFactor()` for
                    :math:`b` in :math:`\\langle R_{0,\\eta}^2 = N_{\\eta} b^2\\rangle`.
        """
    def get_displacement_residual_norm(self) -> float:
        """
                   Get the current link displacement residual norm.
        """
    def get_displacements(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                   Get the current link displacements.
        """
    def get_effective_functionality_of_atoms(self, tolerance: typing.SupportsFloat = 0.001) -> dict[int, int]:
        """
                   Returns the number of active springs connected to each atom, atomId used as index
        
                   :param tolerance: springs under this length are considered inactive
        """
    def get_exit_reason(self) -> ExitReason:
        """
                    Returns the reason for termination of the simulation
        """
    def get_force_magnitude_vector(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                   Evaluate the norm of the force on each (slip- or cross-) link.
        """
    def get_force_on(self, link_idx: typing.SupportsInt) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                   Evaluate the force on a particular (slip- or cross-) link.
        """
    def get_gamma_factor(self, b02: typing.SupportsFloat = -1.0, nr_of_chains: typing.SupportsInt = -1) -> float:
        """
                   Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:
        
                   :math:`\\Gamma = \\langle\\gamma_{\\eta}\\rangle`, with :math:`\\gamma_{\\eta} = \\\\frac{\\bar{r_{\\eta}}^2}{R_{0,\\eta}^2}`,
                   which you can use as :math:`G_{\\mathrm{ANT}} = \\Gamma \\nu k_B T`,
                   where :math:`\\eta` is the index of a particular strand,
                   :math:`R_{0}^2` is the melt mean square end to end distance, in phantom systems :math:`$= N_{\\eta}*b^2$`
                   :math:`N_{\\eta}` is the number of atoms in this strand :math:`\\eta`,
                   :math:`b` its mean square bond length,
                   :math:`T` the temperature and
                   :math:`k_B` Boltzmann's constant.
        
                   :param b02: The melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\\infinity b^2`.
                   :param nr_of_chains: The value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of springs.
        """
    def get_gamma_factors(self, b02: typing.SupportsFloat) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                   Evaluates the gamma factor for each strand (i.e., the squared distance divided by the contour length multiplied by b02)
        """
    def get_gamma_factors_in_dir(self, b02: typing.SupportsFloat, direction: typing.SupportsInt) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                        Evaluates the gamma factor for each strand in the specified direction (i.e., the squared distance divided by the contour length multiplied by b02)
        
                        :param b02: The melt :math:`<b>_0^2`: mean bond length squared; vgl. the required <R_0^2>, computed as phantom = N<b>^2; otherwise, it's the slope in a <R_0^2> vs. N plot, also sometimes labelled :math:`C_\\infinity b^2`.
                        :param direction: The direction in which to compute the gamma factors (0: x, 1: y, 2: z)
        """
    def get_ids_of_active_nodes(self, tolerance: typing.SupportsFloat = 0.001) -> list[int]:
        """
                   Get the atom ids of the nodes that are considered active.
                   Only crosslink ids are returned (not e.g. entanglement links).
        
                   :param tolerance: springs under this length are considered inactive. A node is active if it has > 1 active springs.
        """
    def get_initial_coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                     Get the initial coordinates of the crosslinkers and entanglement links.
        """
    def get_nr_of_active_nodes(self, tolerance: typing.SupportsFloat = 0.001) -> int:
        """
                   Get the number of active nodes (incl. entanglement nodes [atoms with type = entanglementType, present in the universe when creating this simulator],
                   excl. entanglement links [the slip-links created internally when e.g. constructing the simulator with random slip-links]).
        
                   :param tolerance: springs under this length are considered inactive. A node is active if it has > 1 active springs.
        """
    def get_nr_of_active_spring(self, tolerance: typing.SupportsFloat = 0.001) -> int:
        """
                   Get the number of active springs remaining after running the simulation.
        
                  :param tolerance: springs under this length are considered inactive
        """
    def get_nr_of_active_springs(self, tolerance: typing.SupportsFloat = 0.001) -> int:
        """
                   Get the number of active springs remaining after running the simulation.
        
                  :param tolerance: springs under this length are considered inactive
        """
    def get_nr_of_active_springs_in_dir(self, direction: typing.SupportsInt, tolerance: typing.SupportsFloat = 0.001) -> int:
        """
                         Get the number of active springs remaining after running the simulation.
        
                        :param direction: The direction in which to compute the active springs (0: x, 1: y, 2: z)
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
    def get_overall_strand_lengths(self) -> list[float]:
        """
                   Get the sum of the lengths of the springs of each strand.
        """
    def get_pressure(self) -> float:
        """
                   Returns the pressure at the current state of the simulation.
        """
    def get_soluble_weight_fraction(self, tolerance: typing.SupportsFloat = 0.001) -> float:
        """
                   Compute the weight fraction of springs connected to active
                   springs (any depth).
        
                   Caution: ignores atom masses.
        """
    def get_stress_on(self, link_idx: typing.SupportsInt) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 3]"]:
        """
                   Evaluate the stress on a particular (slip- or cross-) link.
        """
    def get_stress_tensor(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 3]"]:
        """
                   Returns the stress tensor at the current state of the simulation.
        
                   The units are in :math:`[\\text{units of }\\kappa]/[\\text{distance units}]`,
                   where the units of :math:`\\kappa` should be :math:`[\\text{force}]/[\\text{distance units}]^2`.
                   Make sure to multiply by :math:`\\kappa` or configure it appropriately.
        """
    def get_stress_tensor_link_based(self, xlinks_only: bool = False) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 3]"]:
        """
                   Returns the stress tensor at the current state of the simulation.
        """
    def get_weighted_spring_lengths(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                   Get the current spring lengths (norm of vector) divided by the spring partition times the contour length.
        """
    def inspect_displacement_to_mean_position_update(self, link_idx: typing.SupportsInt) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                   Helper method to debug and/or understand what happens to certain links when being displaced.
        """
    def run_force_relaxation(self, simplification_mode: StructureSimplificationMode = ..., inactive_removal_cutoff: typing.SupportsFloat = 1e-06, sle_solver: SLESolver = ..., tolerance: typing.SupportsFloat = 1e-15, max_iterations: typing.SupportsInt = 10000) -> None:
        """
                   Run the simulation.
        
                   :param simplification_mode: How to simplify the structure during the minimization.
                   :param inactive_removal_cutoff: The tolerance in distance units of the partial spring length to count as active, relevant if simplification mode is specified to be something other than NO_SIMPLIFICATION.
                   :param sle_solver: The solver to use for the system of linear equations (SLE).
                   :param tolerance: The stopping condition/tolerance for the SLE solver if it's an iterative one.
                   :param max_iterations: The maximum number of iterations to perform if the SLE solver is an iterative one.
        """
    def set_displacements(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        """
                   Set the current link displacements.
        """
    def set_spring_contour_lengths(self, arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        """
                   Set/overwrite the contour lengths.
        """
    def validate_network(self) -> bool:
        """
                   Validates the internal structures.
        
                   Throws an error if something is not ok.
                   Otherwise, it returns true.
        
                   Can be used e.g. as :code:`assert fb.validate_network()`.
        
                   :return: True if validation passes, raises exception otherwise
        """
    @property
    def network(self) -> SimplifiedBalance2Network:
        ...
class MEHPForceEvaluator:
    """
    
         The base interface to change the way the force is evaluated during a MEHP run.
        
    """
    is_2d: bool
    def __init__(self) -> None:
        ...
    def evaluate_stress_contribution(self, spring_distances: typing.SupportsFloat, i: typing.SupportsInt, j: typing.SupportsInt, spring_index: typing.SupportsInt) -> float:
        """
                  An evaluation of the stress-contribution.
        
                  :param springDistances: The three coordinate differences for one spring.
                  :param i: The row index of the stress tensor
                  :param j: The column index of the stress tensor
        """
    @property
    def network(self) -> SimplifiedNetwork:
        ...
class MEHPForceRelaxation:
    """
    
        A small simulation tool for quickly minimizing the force between the crosslinker beads.
    
        This is the first of three force relaxation methods available in this library.
        The relevant feature of this implementation is the configurable spring potential.
        Consequently, it offers a variety of configurable non-linear solvers using NLoptLib.
    
        Please cite :cite:t:`gusev_numerical_2019` if you use this method in your work.
        
    """
    def __getstate__(self) -> tuple:
        ...
    def __init__(self, universe: Universe, crosslinker_type: typing.SupportsInt = 2, is_2d: bool = False, force_evaluator: MEHPForceEvaluator = None, kappa: typing.SupportsFloat = 1.0, remove_2functional_crosslinkers: bool = False, remove_dangling_chains: bool = False) -> None:
        """
                  Instantiate the simulator for a certain universe.
        
                  :param universe: The universe to simulate with
                  :param crosslinker_type: The atom type of the crosslinkers. Needed to reduce the network.
                  :param is2d: Whether to ignore the z direction.
                  :param force_evaluator: The force evaluator to use
                  :param kappa: The spring constant
                  :param remove_2functional_crosslinkers: Whether to replace two-functional crosslinkers with a "normal" chain bead
                  :param remove_dangling_chains: Whether to remove dangling chains before running the simulation.
                       **Caution**: Removing the dangling chains will result in incorrect results fo the computation of
                       :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.get_soluble_weight_fraction` and
                       :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.get_dangling_weight_fraction`
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def assume_box_large_enough(self, box_large_enough: bool = False) -> None:
        """
                  Configure whether to run PBC on the bonds or not.
        
                  :param box_large_enough: If True, assume the box is large enough and don't apply PBC on bonds.
                                          If your bonds could get larger than half the box length, 
                                          this must be kept False (default).
        """
    def config_rerun_epsilon(self, epsilon: typing.SupportsFloat = 0.001) -> None:
        """
                  Configure the offset from the lower and upper bounds for the simulation to suggest another run.
                  
                  :param epsilon: The epsilon value to use for the rerun check
                       (See: :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.requires_another_run`)
        """
    @typing.overload
    def config_step_output(self, output_configuration: collections.abc.Sequence[OutputConfiguration]) -> None:
        """
                  Set which values to log during the simulation.
        
                  :param output_configuration: An OutputConfiguration struct or list of OutputConfiguration structs
                                              specifying what values to log and how often
        """
    @typing.overload
    def config_step_output(self, output_configuration: collections.abc.Sequence[OutputConfiguration]) -> None:
        """
                  Set which values to log during the simulation.
        
                  :param output_configuration: An OutputConfiguration struct or list of OutputConfiguration structs
                                              specifying what values to log and how often
        """
    def count_active_clustered_atoms(self, tolerance: typing.SupportsFloat = 0.001) -> float:
        """
                  Counts the active clustered atoms in the system.
        
                  :param tolerance: Springs under this length are considered inactive
                  :return: The number of active clustered atoms
        """
    def get_active_chains(self, tolerance: typing.SupportsFloat = 0.001) -> list[Molecule]:
        """
                  Get the crosslinker chains that are active.
        
                  :param tolerance: Springs under this length are considered inactive
                  :return: List of active crosslinker chains
        """
    def get_average_contour_length(self) -> float:
        """
                  Get the average contour length of all springs.
        
                  :return: The average contour length
        """
    def get_average_spring_length(self) -> float:
        """
                  Get the average length of the springs. Note that in contrast to :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor`,
                  this value is normalized by the number of springs rather than the number of chains.
        
                  :return: The average spring length
        """
    def get_crosslinker_universe(self) -> Universe:
        """
                  Returns the universe [of crosslinkers] with the positions of the current state of the simulation.
        
                  :return: A Universe object containing the crosslinkers with updated positions
        """
    def get_current_spring_distances(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Get the current spring distances.
        
                  :return: A vector containing the current spring distance vectors
        """
    def get_dangling_weight_fraction(self, tolerance: typing.SupportsFloat = 0.001) -> float:
        """
                  Compute the weight fraction of non-active springs.
        
                  Caution: ignores atom masses.
        
                  :param tolerance: Springs under this length are considered inactive
                  :return: The dangling weight fraction
        """
    def get_default_r0_square(self) -> float:
        """
                  Returns the value effectively used in :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor` for :math:`\\langle R_{0,\\eta}^2\\rangle`.
        
                  :return: The default R0 squared value used in gamma factor calculations
        """
    def get_effective_functionality_of_atoms(self, tolerance: typing.SupportsFloat = 0.001) -> dict[int, int]:
        """
                  Returns the number of active springs connected to each atom, atomId used as index.
        
                  :param tolerance: Springs under this length are considered inactive
                  :return: Vector with effective functionality for each atom (indexed by atom ID)
        """
    def get_exit_reason(self) -> ExitReason:
        """
                  Returns the reason for termination of the simulation.
        
                  :return: The exit reason enum value indicating why the simulation ended
        """
    def get_force(self) -> float:
        """
                  Returns the force at the current state of the simulation.
                  
                  :return: The current force value
        """
    def get_gamma_factor(self, b0_squared: typing.SupportsFloat = -1.0, nr_of_chains: typing.SupportsInt = -1) -> float:
        """
                  Computes the gamma factor as part of the ANT/MEHP formulism, i.e.:
        
                  :math:`\\Gamma = \\langle\\gamma_{\\eta}\\rangle`, with :math:`\\gamma_{\\eta} = \\\\frac{\\bar{r_{\\eta}}^2}{R_{0,\\eta}^2}`,
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
                  :param nr_of_chains: The value to normalize the sum of square distances by. Usually (and default if :math:`< 0`) the nr of chains.
        """
    def get_gamma_factors(self, b0_squared: typing.SupportsFloat = -1.0) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Computes the gamma factor for each spring as part of the ANT/MEHP formulism.
        
                  :math:`\\gamma_{\\eta} = \\\\frac{\\bar{r_{\\eta}}^2}{R_{0,\\eta}^2}`, with (here)
                  :math:`R_{0,\\eta}^2 = N_\\eta \\cdot ` the parameter `b0_squared`.
                  You can obtain this parameter e.g. by doing melt simulations at different lengths,
                  it's the slope you obtain.
        
                  :param b0_squared: Part of the denominator in the equation of :math:`\\Gamma`.
                       If :math:`-1.0` (default), the network is used for determination (which is not accurate), the system is assumed to be phantom.
                       For real systems, the value could be determined by :func:`~pylimer_tools_cpp.Universe.compute_mean_squared_end_to_end_distance()`
                       on the melt system, with subsequent division by the nr of bonds in the chain.
        
                  See also :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.get_gamma_factor` for the mean of these.
        """
    def get_ids_of_active_nodes(self, tolerance: typing.SupportsFloat = 0.001, minimum_nr_of_active_connections: typing.SupportsInt = 2, maximum_nr_of_active_connections: typing.SupportsInt = -1) -> list[int]:
        """
                  Get the atom ids of the nodes that are considered active.
        
                  :param tolerance: Springs under this length are considered inactive. A node is active if it has > 2 active springs.
                  :param minimum_nr_of_active_connections: Minimum number of active connections required for a node to be considered active
                  :param maximum_nr_of_active_connections: Maximum number of active connections allowed for a node to be considered active
                  :return: List of atom IDs for active nodes
        """
    def get_nr_of_active_nodes(self, tolerance: typing.SupportsFloat = 0.001, minimumNrOfActiveConnections: typing.SupportsInt = 2, maximumNrOfActiveConnections: typing.SupportsInt = -1) -> int:
        """
                  Get the number of active nodes remaining after running the simulation.
        
                  :param tolerance: Springs under this length are considered inactive
                  :param minimumNrOfActiveConnections: A node is active if it has equal or more than this number of active springs
                  :param maximumNrOfActiveConnections: A node is active if it has equal or less than this number of active springs.
                       Use a value < 0 to indicate that there is no maximum number of active connections
                  :return: The number of active nodes
        """
    def get_nr_of_active_springs(self, tolerance: typing.SupportsFloat = 0.001) -> int:
        """
                  Get the number of active springs remaining after running the simulation.
        
                  :param tolerance: Springs under this length are considered inactive
                  :return: The number of active springs
        """
    def get_nr_of_active_springs_connected(self, tolerance: typing.SupportsFloat = 0.05) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        """
                  Returns the number of active springs connected to each node.
        
                  :param tolerance: Springs under this length are considered inactive
                  :return: Vector with the number of active springs for each node
        """
    def get_nr_of_iterations(self) -> int:
        """
                  Returns the number of iterations used for force relaxation.
        
                  :return: The number of iterations performed during force relaxation
        """
    def get_nr_of_nodes(self) -> int:
        """
                  Get the number of nodes considered in this simulation.
        
                  :return: The number of nodes in the simulation
        """
    def get_nr_of_springs(self) -> int:
        """
                  Get the number of springs considered in this simulation.
        
                  :param tolerance: Springs under this length are considered inactive
                  :return: The number of springs in the simulation
        """
    def get_pressure(self) -> float:
        """
                  Returns the pressure at the current state of the simulation.
        
                  :return: The current pressure value
        """
    def get_residual_norm(self) -> float:
        """
                  Returns the residual norm at the current state of the simulation.
        
                  :return: The current residual norm value
        """
    def get_residuals(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Returns the residuals at the current state of the simulation.
        
                  :return: The current residual vector
        """
    def get_soluble_weight_fraction(self, tolerance: typing.SupportsFloat = 0.001) -> float:
        """
                  Compute the weight fraction of springs connected to active
                  springs (any depth).
        
                  Caution: ignores atom masses.
        
                  :param tolerance: Springs under this length are considered inactive
                  :return: The soluble weight fraction
        """
    def get_spring_contour_length(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Get the spring contour lengths.
        
                  :return: A vector containing the contour lengths of all springs
        """
    def get_spring_distances(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Get the current coordinate differences for all the springs.
        
                  :return: A vector of size 3*nrOfSprings, with each x, y, z values of the springs
        """
    def get_spring_lengths(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                  Get the current lengths for all the springs.
        
                  :return: A vector of size nrOfSprings, with each the norm of the distances
        """
    def get_stress_tensor(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 3]"]:
        """
                  Returns the stress tensor at the current state of the simulation.
        
                  :return: The current stress tensor matrix
        """
    def requires_another_run(self) -> bool:
        """
                  For performance reasons, the objective is only minimised within the distances of one box.
                  This means, that there is a possibility, e.g. for a single strand longer than two boxes,
                  that it would not be globally minimised.
        
                  If the final displacement of one of the atoms is close
                  (1e-3, configurable via :meth:`~pylimer_tools_cpp.MEHPForceRelaxation.config_rerun_epsilon`)
                  to the imposed min/max, after minimizing,
                  this method would return true.
        
                  :return: True if another run is suggested, False otherwise
        """
    def run_force_relaxation(self, algorithm: str = 'LD_MMA', max_nr_of_steps: typing.SupportsInt = 250000, x_tolerance: typing.SupportsFloat = 1e-12, f_tolerance: typing.SupportsFloat = 1e-09) -> None:
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
                  
                  :param force_evaluator: The new force evaluator to use
        """
    @property
    def network(self) -> SimplifiedNetwork:
        """
                  Get the network structure.
        
                  :return: The network structure used in the simulation
        """
class MaxDistanceProvider:
    """
    
         A generic implementation of a class, that shall provide a maximum distance for the MC sampling.
         
    """
    def get_max_distance(self, N: typing.SupportsFloat) -> float:
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
    def __getitem__(self, arg0: typing.SupportsInt) -> Atom:
        """
               Access an atom by its vertex index.
        """
    def __init__(self, arg0: Box, arg1: igraph_s, arg2: MoleculeType, arg3: collections.abc.Mapping[typing.SupportsInt, typing.SupportsFloat]) -> None:
        """
                 Construct a molecule from a graph structure.
                 
                 :param box: The simulation box containing this molecule
                 :param graph: Pointer to the igraph structure representing connectivity
                 :param type: The type of molecule (see MoleculeType enum)
                 :param mass_map: Map of atom types to their masses
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
    def compute_end_to_end_distance_with_derived_image_flags(self) -> float:
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
    def compute_end_to_end_vector_with_derived_image_flags(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
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
    def compute_total_vector(self, crosslinker_type: typing.SupportsInt = 2, close_loop: bool = True) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                       Compute the sum of all bond vectors.
                       
                       :param crosslinker_type: The type of crosslinker atoms
                       :param close_loop: Whether to close the loop for calculations
                       :return: The vector sum of all bonds
        """
    def compute_vector_from_to(self, atom_id_from: typing.SupportsInt, atom_id_to: typing.SupportsInt, crosslinker_type: typing.SupportsInt = 2, require_order: bool = True) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]:
        """
                       Compute the sum of all bond vectors between two specified atoms.
                       
                       :param atom_id_from: ID of the starting atom
                       :param atom_id_to: ID of the ending atom
                       :param crosslinker_type: The type of crosslinker atoms
                       :param require_order: Whether to require specific ordering
                       :return: The vector sum from start to end atom
        """
    def get_atom_by_id(self, atom_id: typing.SupportsInt) -> Atom:
        """
                    Get an atom by its id.
                    
                    :param atom_id: The atom ID to search for
                    :return: The atom with the specified ID
        """
    def get_atom_by_vertex_idx(self, vertex_idx: typing.SupportsInt) -> Atom:
        """
                    Get an atom for a specific vertex.
                    
                    :param vertex_idx: The vertex index to query
                    :return: The atom at the specified vertex
        """
    def get_atom_id_by_vertex_idx(self, vertex_id: typing.SupportsInt) -> int:
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
    def get_atoms_by_degree(self, degree: typing.SupportsInt) -> list[Atom]:
        """
                    Get the atoms that have the specified number of bonds.
                    
                    :param degree: The number of bonds (degree/functionality)
                    :return: List of atoms with the specified degree
        """
    def get_atoms_by_type(self, type: typing.SupportsInt) -> list[Atom]:
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
    def get_atoms_connected_to_vertex(self, vertex_idx: typing.SupportsInt) -> list[Atom]:
        """
                    Get the atoms connected to a specified vertex id.
                    
                    :param vertex_idx: The vertex index to query
                    :return: List of connected atoms
        """
    def get_atoms_lined_up(self, crosslinker_type: typing.SupportsInt = 2, assumed_coordinates: bool = False, close_loop: bool = False) -> list[Atom]:
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
    def get_edge_ids_from(self, vertex_id: typing.SupportsInt) -> list[int]:
        """
                 Get the edge IDs incident to a specific vertex.
                 
                 :param vertex_id: The vertex to query
                 :return: A vector of edge IDs connected to the vertex
        """
    def get_edge_ids_from_to(self, vertex_id_from: typing.SupportsInt, vertex_id_to: typing.SupportsInt) -> list[int]:
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
    def get_nr_of_edges_from_to(self, vertex_id_from: typing.SupportsInt, vertex_id_to: typing.SupportsInt, max_length: typing.SupportsInt = -1) -> int:
        """
                 Get the number of edges in the shortest path between two specific vertices.
        
                 If `max_length` is provided and positive, it will only consider paths up to that length.
                 
                 :param vertex_id_from: Starting vertex ID
                 :param vertex_id_to: Ending vertex ID  
                 :param max_length: Maximum path length to consider (-1 for no limit)
                 :return: Number of edges in the shortest path, or -1 if no path exists
        """
    def get_strand_ends(self, crosslinker_type: typing.SupportsInt = 2, close_loop: bool = False) -> list[Atom]:
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
    def get_vertex_idx_by_atom_id(self, atom_id: typing.SupportsInt) -> int:
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
class MoleculeType:
    """
    An enum representing the type of molecule/chain/strand.
    
    Members:
    
      UNDEFINED : This value indicates that either the property was not set or not discovered.
    
      NETWORK_STRAND : 
               A network strand is a strand in a network.
          
    
      PRIMARY_LOOP : 
               A primary loop is a network strand looping from and to the same crosslinker.
          
    
      DANGLING_CHAIN : 
               A dangling chain is a network strand where only one end is attached to a crosslinker.
          
    
      FREE_CHAIN : 
               A free chain is a strand not connected to any crosslinker.
          
    """
    DANGLING_CHAIN: typing.ClassVar[MoleculeType]  # value = <MoleculeType.DANGLING_CHAIN: 3>
    FREE_CHAIN: typing.ClassVar[MoleculeType]  # value = <MoleculeType.FREE_CHAIN: 4>
    NETWORK_STRAND: typing.ClassVar[MoleculeType]  # value = <MoleculeType.NETWORK_STRAND: 1>
    PRIMARY_LOOP: typing.ClassVar[MoleculeType]  # value = <MoleculeType.PRIMARY_LOOP: 2>
    UNDEFINED: typing.ClassVar[MoleculeType]  # value = <MoleculeType.UNDEFINED: 0>
    __members__: typing.ClassVar[dict[str, MoleculeType]]  # value = {'UNDEFINED': <MoleculeType.UNDEFINED: 0>, 'NETWORK_STRAND': <MoleculeType.NETWORK_STRAND: 1>, 'PRIMARY_LOOP': <MoleculeType.PRIMARY_LOOP: 2>, 'DANGLING_CHAIN': <MoleculeType.DANGLING_CHAIN: 3>, 'FREE_CHAIN': <MoleculeType.FREE_CHAIN: 4>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
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
    
        Gives access to somewhat fast queries on the neighbourhood of atoms.
        
        This class provides efficient spatial queries for finding atoms within
        a specified distance of each other.
        
    """
    def __init__(self, atoms: collections.abc.Sequence[Atom], box: Box, cutoff: typing.SupportsFloat) -> None:
        """
                 Instantiate a new neighbour list.
                 
                 :param atoms: Vector of atoms to include in the neighbour list
                 :param box: The simulation box
                 :param cutoff: Maximum distance for neighbour searches
        """
    def get_atoms_close_to(self, atom: Atom, upper_cutoff: typing.SupportsFloat = 1.0, lower_cutoff: typing.SupportsFloat = 0.0, unwrapped: bool = False, expect_self: bool = False) -> list[Atom]:
        """
                  List all atoms that are close to a given one. 
        
                  It is possible to request it within a new cutoff, 
                  though the underlying neighbour list will not be regenerated.
                  For performance reasons, it is recommended to initialize a 
                  new NeighbourList if you require a different cutoff, depending on your use case.
        
                  You can use a negative value for the upper_cutoff to use the cutoff used for 
                  filling the neighbour list buckets.
                  
                  :param atom: The reference atom
                  :param upper_cutoff: Maximum distance for neighbours
                  :param lower_cutoff: Minimum distance for neighbours
                  :param unwrapped: Whether to use unwrapped coordinates
                  :param expect_self: Whether to expect the atom itself in results
                  :return: List of neighbouring atoms
        """
    def remove_atom(self, atom: Atom, debug_hint: str = '') -> None:
        """
                 Remove an atom from this neighbour list.
                 It will not show up when querying for neighbours, 
                 but its neighbours cannot be queried either.
                 
                 :param atom: The atom to remove
                 :param debug_hint: Optional debug information
        """
class NoMaxDistanceProvider(MaxDistanceProvider):
    """
    
        For MC generation, to disable the neighbour list usage.
        
    """
    def __init__(self) -> None:
        ...
    def get_max_distance(self, N: typing.SupportsFloat) -> float:
        """
                 Get the maximum distance for a given N (always returns -1 to disable).
        
                 :param N: Number of segments (ignored).
                 :return: Always returns -1 to disable maximum distance checks.
        """
class NonGaussianSpringForceEvaluator(MEHPForceEvaluator):
    """
    
         This is equal to a spring evaluator for Langevin chains.
    
         The force for a certain spring is given by:
         :math:`f = 0.5 \\cdot \\\\frac{1}{l} \\scriptL^{-1}(\\frac{r}{N\\cdot l})`,
         where :math:`r` is the spring [between crosslinkers] length
         and :math:`\\scriptL^{-1}` the inverse langevin function.
    
         Please note that the inverse langevin is only approximated.
    
         Recommended optimization algorithm: "LD_MMA"
    
         :param kappa: The spring constant :math:`\\kappa`
         :param N: The number of links in a spring
         :param l: The  the length of a spring in the chain
        
    """
    def __init__(self, kappa: typing.SupportsFloat = 1.0, N: typing.SupportsFloat = 1.0, l: typing.SupportsFloat = 1.0) -> None:
        """
        Initialize this ForceEvaluator
        """
class NormalModeAnalyzer:
    """
    
        Compute the normal modes and predict the loss/storage moduli.
    
        Please cite :cite:t:`gusev_molecular_2024` if you use this method in your work.
        
    """
    def __getstate__(self) -> tuple:
        ...
    def __init__(self, spring_from: collections.abc.Sequence[typing.SupportsInt], spring_to: collections.abc.Sequence[typing.SupportsInt]) -> None:
        """
                 Initialize the NormalModeAnalyzer with the bonds (edges).
        
                 Constructs the connectivity matrix from the given edges.
                 
                 :param spring_from: Vector of starting node indices for springs/bonds
                 :param spring_to: Vector of ending node indices for springs/bonds
        """
    def __setstate__(self, arg0: tuple) -> None:
        ...
    def evaluate_loss_modulus(self, omega: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                 Evaluate the loss modulus :math:`G''(\\omega)`. Yet misses the conversion factor.
                 
                 :param omega: Angular frequencies
                 :return: Loss modulus values
        """
    def evaluate_storage_modulus(self, omega: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                 Evaluate the storage modulus :math:`G'(\\omega)`. Yet misses the conversion factor.
                 
                 :param omega: Angular frequencies
                 :return: Storage modulus values
        """
    def evaluate_stress_autocorrelation(self, t: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                 Evaluate stress autocorrelation :math:`C(t)`.
                 
                 :param t: The time at which to evaluate the stress autocorrelation
                 :return: Stress autocorrelation values
        """
    def find_all_eigenvalues(self, compute_eigenvectors: bool = False) -> None:
        """
                 Find all eigenvalues using a dense solver.
                 
                 :param compute_eigenvectors: Whether to also compute eigenvectors (default: False)
                 :return: True if computation was successful
        """
    def find_sparse_eigenvalues(self, nr_of_eigenvalues: typing.SupportsInt, compute_eigenvectors: bool = False) -> None:
        """
                 Find the k smallest eigenvalues using a sparse solver.
                 
                 :param nr_of_eigenvalues: Number of smallest eigenvalues to find
                 :param compute_eigenvectors: Whether to also compute eigenvectors (default: False)
                 :return: True if computation was successful
        """
    def get_eigenvalues(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        """
                 Get the eigenvalues.
                 
                 :return: Vector of eigenvalues
        """
    def get_eigenvectors(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, n]"]:
        """
                 Get eigenvectors.
                 
                 :return: Matrix of eigenvectors
        """
    def get_matrix(self) -> scipy.sparse.csc_matrix:
        """
                 Get the assembled connectivity matrix.
                 
                 :return: The connectivity matrix
        """
    def get_matrix_size(self) -> int:
        """
                 Get the size of the matrix (the maximum number of eigenvalues that could be queried).
                 
                 :return: Size of the connectivity matrix
        """
    def get_nr_of_soluble_clusters(self) -> int:
        """
                 Get the number of soluble clusters (eigenvalues = 0).
                 
                 :return: Number of soluble clusters
        """
    def set_eigenvalues(self, eigenvalues: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, 1]"]) -> None:
        """
                 Set the eigenvalues, e.g. if you use an external solver.
                 
                 :param eigenvalues: Vector of eigenvalues to set
        """
    def set_eigenvectors(self, eigenvectors: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[m, n]"]) -> None:
        """
                 Set eigenvectors, e.g. if you use an external solver.
                 
                 :param eigenvectors: Matrix of eigenvectors to set
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
    def output_every(self, arg0: typing.SupportsInt) -> None:
        ...
    @property
    def use_every(self) -> int:
        """
             For autocorrelation and averaging, how often to include values.
        
             Use a value of 1 to take average of or autocorrelate, respectively,
             all values encountered during the simulation or optimization procedure.
        """
    @use_every.setter
    def use_every(self, arg0: typing.SupportsInt) -> None:
        ...
class SLESolver:
    """
    Solver for sparse linear equation systems
    
    Members:
    
      DEFAULT : default
    
      SIMPLICIAL_LLT : SimplicialLLT
    
      SIMPLICIAL_LDLT : SimplicialLDLT
    
      SPARSE_LU : SparseLU
    
      SPARSE_QR : SparseQR
    
      CONJUGATE_GRADIENT : ConjugateGradient
    
      CONJUGATE_GRADIENT_DIAGONALIZED : ConjugateGradient, DiagonalPreconditioner
    
      CONJUGATE_GRADIENT_IDENTITY : ConjugateGradient, IdentityPreconditioner
    
      CONJUGATE_GRADIENT_INCOMPLETE_CHOLESKY : ConjugateGradient, IncompleteCholeskyPreconditioner
    
      LEAST_SQUARES_CONJUGATE_GRADIENT : LeastSquaresConjugateGradient
    
      LEAST_SQUARES_CONJUGATE_GRADIENT_DIAGONALIZED : LeastSquaresConjugateGradient, DiagonalPreconditioner
    
      LEAST_SQUARES_CONJUGATE_GRADIENT_IDENTITY : LeastSquaresConjugateGradient, IdentityPreconditioner
    
      BICGSTAB : BiCGSTAB
    
      BICGSTAB_DIAGONALIZED : BiCGSTAB, DiagonalPreconditioner
    
      BICGSTAB_IDENTITY : BiCGSTAB, IdentityPreconditioner
    
      BICGSTAB_INCOMPLETE_LU : BiCGSTAB, IncompleteLUTPreconditioner
    
      GRADIENT_DESCENT : GradientDescent
    
      GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT : GradientDescent (Barzilai-Borwein method, short time-step)
    
      GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG : GradientDescent (Barzilai-Borwein method, long time-step)
    
      GRADIENT_DESCENT_BARZILAI_BORWEIN_MOMENTUM : GradientDescent (Barzilai-Borwein & heavy ball method, selective time-step)
    """
    BICGSTAB: typing.ClassVar[SLESolver]  # value = <SLESolver.BICGSTAB: 12>
    BICGSTAB_DIAGONALIZED: typing.ClassVar[SLESolver]  # value = <SLESolver.BICGSTAB_DIAGONALIZED: 13>
    BICGSTAB_IDENTITY: typing.ClassVar[SLESolver]  # value = <SLESolver.BICGSTAB_IDENTITY: 14>
    BICGSTAB_INCOMPLETE_LU: typing.ClassVar[SLESolver]  # value = <SLESolver.BICGSTAB_INCOMPLETE_LU: 15>
    CONJUGATE_GRADIENT: typing.ClassVar[SLESolver]  # value = <SLESolver.CONJUGATE_GRADIENT: 5>
    CONJUGATE_GRADIENT_DIAGONALIZED: typing.ClassVar[SLESolver]  # value = <SLESolver.CONJUGATE_GRADIENT_DIAGONALIZED: 6>
    CONJUGATE_GRADIENT_IDENTITY: typing.ClassVar[SLESolver]  # value = <SLESolver.CONJUGATE_GRADIENT_IDENTITY: 7>
    CONJUGATE_GRADIENT_INCOMPLETE_CHOLESKY: typing.ClassVar[SLESolver]  # value = <SLESolver.CONJUGATE_GRADIENT_INCOMPLETE_CHOLESKY: 8>
    DEFAULT: typing.ClassVar[SLESolver]  # value = <SLESolver.DEFAULT: 0>
    GRADIENT_DESCENT: typing.ClassVar[SLESolver]  # value = <SLESolver.GRADIENT_DESCENT: 16>
    GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG: typing.ClassVar[SLESolver]  # value = <SLESolver.GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG: 18>
    GRADIENT_DESCENT_BARZILAI_BORWEIN_MOMENTUM: typing.ClassVar[SLESolver]  # value = <SLESolver.GRADIENT_DESCENT_BARZILAI_BORWEIN_MOMENTUM: 19>
    GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT: typing.ClassVar[SLESolver]  # value = <SLESolver.GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT: 17>
    LEAST_SQUARES_CONJUGATE_GRADIENT: typing.ClassVar[SLESolver]  # value = <SLESolver.LEAST_SQUARES_CONJUGATE_GRADIENT: 9>
    LEAST_SQUARES_CONJUGATE_GRADIENT_DIAGONALIZED: typing.ClassVar[SLESolver]  # value = <SLESolver.LEAST_SQUARES_CONJUGATE_GRADIENT_DIAGONALIZED: 10>
    LEAST_SQUARES_CONJUGATE_GRADIENT_IDENTITY: typing.ClassVar[SLESolver]  # value = <SLESolver.LEAST_SQUARES_CONJUGATE_GRADIENT_IDENTITY: 11>
    SIMPLICIAL_LDLT: typing.ClassVar[SLESolver]  # value = <SLESolver.SIMPLICIAL_LDLT: 2>
    SIMPLICIAL_LLT: typing.ClassVar[SLESolver]  # value = <SLESolver.SIMPLICIAL_LLT: 1>
    SPARSE_LU: typing.ClassVar[SLESolver]  # value = <SLESolver.SPARSE_LU: 3>
    SPARSE_QR: typing.ClassVar[SLESolver]  # value = <SLESolver.SPARSE_QR: 4>
    __members__: typing.ClassVar[dict[str, SLESolver]]  # value = {'DEFAULT': <SLESolver.DEFAULT: 0>, 'SIMPLICIAL_LLT': <SLESolver.SIMPLICIAL_LLT: 1>, 'SIMPLICIAL_LDLT': <SLESolver.SIMPLICIAL_LDLT: 2>, 'SPARSE_LU': <SLESolver.SPARSE_LU: 3>, 'SPARSE_QR': <SLESolver.SPARSE_QR: 4>, 'CONJUGATE_GRADIENT': <SLESolver.CONJUGATE_GRADIENT: 5>, 'CONJUGATE_GRADIENT_DIAGONALIZED': <SLESolver.CONJUGATE_GRADIENT_DIAGONALIZED: 6>, 'CONJUGATE_GRADIENT_IDENTITY': <SLESolver.CONJUGATE_GRADIENT_IDENTITY: 7>, 'CONJUGATE_GRADIENT_INCOMPLETE_CHOLESKY': <SLESolver.CONJUGATE_GRADIENT_INCOMPLETE_CHOLESKY: 8>, 'LEAST_SQUARES_CONJUGATE_GRADIENT': <SLESolver.LEAST_SQUARES_CONJUGATE_GRADIENT: 9>, 'LEAST_SQUARES_CONJUGATE_GRADIENT_DIAGONALIZED': <SLESolver.LEAST_SQUARES_CONJUGATE_GRADIENT_DIAGONALIZED: 10>, 'LEAST_SQUARES_CONJUGATE_GRADIENT_IDENTITY': <SLESolver.LEAST_SQUARES_CONJUGATE_GRADIENT_IDENTITY: 11>, 'BICGSTAB': <SLESolver.BICGSTAB: 12>, 'BICGSTAB_DIAGONALIZED': <SLESolver.BICGSTAB_DIAGONALIZED: 13>, 'BICGSTAB_IDENTITY': <SLESolver.BICGSTAB_IDENTITY: 14>, 'BICGSTAB_INCOMPLETE_LU': <SLESolver.BICGSTAB_INCOMPLETE_LU: 15>, 'GRADIENT_DESCENT': <SLESolver.GRADIENT_DESCENT: 16>, 'GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT': <SLESolver.GRADIENT_DESCENT_BARZILAI_BORWEIN_SHORT: 17>, 'GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG': <SLESolver.GRADIENT_DESCENT_BARZILAI_BORWEIN_LONG: 18>, 'GRADIENT_DESCENT_BARZILAI_BORWEIN_MOMENTUM': <SLESolver.GRADIENT_DESCENT_BARZILAI_BORWEIN_MOMENTUM: 19>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class SimpleSpringMEHPForceEvaluator(MEHPForceEvaluator):
    """
    
         This is equal to a spring evaluator for Gaussian chains.
    
         The force for a certain spring is given by:
         :math:`f = 0.5 \\cdot \\kappa r`,
         where :math:`r` is the spring [between crosslinkers] length.
    
         Recommended optimization algorithm: "LD_LBFGS"
    
         :param kappa: The spring constant :math:`\\kappa`
        
    """
    def __init__(self, kappa: typing.SupportsFloat = 1.0) -> None:
        ...
class SimplifiedBalance2Network:
    """
    
    A more efficient structure of the network for use in
    :obj:`~pylimer_tools_cpp.MEHPForceBalance2`.
    Consists usually only of the cross- and slip-links (and their connectivity),
    i.e., no "normal strand beads" in between, in order to reduce the degrees of freedom
    and therewith improve performance of the solver.
    
    A note on the terminology: a spring is the connection between two links (crosslink, entanglement-link/slip-link).
    A strand is a chain of connected links between two crosslinks.
    """
    @property
    def box_lengths(self) -> typing.Annotated[list[float], "FixedSize(3)"]:
        ...
    @property
    def coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    @property
    def link_is_entanglement(self) -> typing.Annotated[numpy.typing.NDArray[numpy.bool], "[m, 1]"]:
        ...
    @property
    def links_of_strand(self) -> list[list[int]]:
        ...
    @property
    def nr_of_crosslinks(self) -> int:
        ...
    @property
    def nr_of_links(self) -> int:
        ...
    @property
    def nr_of_springs(self) -> int:
        ...
    @property
    def nr_of_strands(self) -> int:
        ...
    @property
    def old_atom_ids(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def old_atom_types(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_box_offset(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    @property
    def spring_contour_length(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    @property
    def spring_coordinate_index_a(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_coordinate_index_b(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_index_a(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_index_b(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_is_entanglement(self) -> typing.Annotated[numpy.typing.NDArray[numpy.bool], "[m, 1]"]:
        ...
    @property
    def springs_of_strand(self) -> list[list[int]]:
        ...
    @property
    def strand_of_spring(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def strands_of_link(self) -> list[list[int]]:
        ...
class SimplifiedBalanceNetwork:
    """
    
         A more efficient structure of the network for use in MEHP force balance,
         namely :obj:`~pylimer_tools_cpp.MEHPForceBalance`, though also passable to
         namely :obj:`~pylimer_tools_cpp.MEHPForceBalance2`.
         Consists usually only of the cross- and slip-links (and their connectivity),
         i.e., no "normal strand beads" in between, in order to reduce the degrees of freedom
         and therewith improve performance of the solver.
    
         A note on the terminology: a spring is the connection between two links (crosslink, entanglement-link/slip-link).
         A strand is a chain of connected links between two crosslinks.
     
    """
    @property
    def box_lengths(self) -> typing.Annotated[list[float], "FixedSize(3)"]:
        ...
    @property
    def coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    @property
    def link_is_sliplink(self) -> typing.Annotated[numpy.typing.NDArray[numpy.bool], "[m, 1]"]:
        ...
    @property
    def links_of_strand(self) -> list[list[int]]:
        ...
    @property
    def nr_of_crosslink_swaps_endured(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def nr_of_crosslinks(self) -> int:
        ...
    @property
    def nr_of_links(self) -> int:
        ...
    @property
    def nr_of_springs(self) -> int:
        ...
    @property
    def nr_of_strands(self) -> int:
        ...
    @property
    def old_atom_ids(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_box_offset(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    @property
    def spring_coordinate_index_a(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_coordinate_index_b(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_index_a(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_index_b(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def springs_of_strand(self) -> list[list[int]]:
        ...
    @property
    def strand_contour_length(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
        ...
    @property
    def strand_coordinate_index_a(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def strand_coordinate_index_b(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def strand_index_a(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def strand_index_b(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def strand_of_spring(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def strands_of_link(self) -> list[list[int]]:
        ...
    @property
    def volume(self) -> float:
        ...
class SimplifiedNetwork:
    """
    
         A more efficient structure of the network for use in MEHP,
         namely :obj:`~pylimer_tools_cpp.MEHPForceRelaxation`.
         Consists usually only of the crosslinkers.
     
    """
    @property
    def assume_box_large_enough(self) -> bool:
        ...
    @property
    def assume_complete(self) -> bool:
        ...
    @property
    def box_lengths(self) -> typing.Annotated[list[float], "FixedSize(3)"]:
        ...
    @property
    def coordinates(self) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
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
    def old_atom_ids(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_coordinate_index_a(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_coordinate_index_b(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_index_a(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def spring_index_b(self) -> typing.Annotated[numpy.typing.NDArray[numpy.int32], "[m, 1]"]:
        ...
    @property
    def volume(self) -> float:
        ...
class StructureSimplificationMode:
    """
    How the structure shall be simplified during the optimization in order to remove non-trapped entanglement links.
    
    Members:
    
      NO_SIMPLIFICATION : No Simplification
    
      X2F_ONLY : Removal of two-functional crosslinks only
    
      INACTIVE_ONLY : Removal of inactive entanglement- and crosslinks only
    
      INACTIVE_THEN_X2F : Removal of inactive, and then two-functional entanglement- and crosslinks, one after the other
    
      X1F_X2F_THEN_INACTIVE : Removal of one- and twofunctional crosslinks, and then inactive entanglements and crosslinks. Deprecated, use INACTIVE_THEN_X2F instead
    """
    INACTIVE_ONLY: typing.ClassVar[StructureSimplificationMode]  # value = <StructureSimplificationMode.INACTIVE_ONLY: 2>
    INACTIVE_THEN_X2F: typing.ClassVar[StructureSimplificationMode]  # value = <StructureSimplificationMode.INACTIVE_THEN_X2F: 3>
    NO_SIMPLIFICATION: typing.ClassVar[StructureSimplificationMode]  # value = <StructureSimplificationMode.NO_SIMPLIFICATION: 0>
    X1F_X2F_THEN_INACTIVE: typing.ClassVar[StructureSimplificationMode]  # value = <StructureSimplificationMode.X1F_X2F_THEN_INACTIVE: 4>
    X2F_ONLY: typing.ClassVar[StructureSimplificationMode]  # value = <StructureSimplificationMode.X2F_ONLY: 1>
    __members__: typing.ClassVar[dict[str, StructureSimplificationMode]]  # value = {'NO_SIMPLIFICATION': <StructureSimplificationMode.NO_SIMPLIFICATION: 0>, 'X2F_ONLY': <StructureSimplificationMode.X2F_ONLY: 1>, 'INACTIVE_ONLY': <StructureSimplificationMode.INACTIVE_ONLY: 2>, 'INACTIVE_THEN_X2F': <StructureSimplificationMode.INACTIVE_THEN_X2F: 3>, 'X1F_X2F_THEN_INACTIVE': <StructureSimplificationMode.X1F_X2F_THEN_INACTIVE: 4>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt) -> Atom:
        """
               Access an atom by its vertex index.
        """
    def __getstate__(self) -> tuple:
        ...
    def __init__(self, Lx: typing.SupportsFloat, Ly: typing.SupportsFloat, Lz: typing.SupportsFloat) -> None:
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
    def add_angles(self, angles_from: collections.abc.Sequence[typing.SupportsInt], angles_via: collections.abc.Sequence[typing.SupportsInt], angles_to: collections.abc.Sequence[typing.SupportsInt], angle_types: collections.abc.Sequence[typing.SupportsInt]) -> None:
        """
                 Add angles to the Universe. No relation to the underlying graph, 
                 just a method to preserve read & write capabilities.
                 
                 :param angles_from: Vector of atom IDs for angle start points
                 :param angles_via: Vector of atom IDs for angle middle points
                 :param angles_to: Vector of atom IDs for angle end points
                 :param angle_types: Vector of angle types
        """
    def add_atoms(self, ids: collections.abc.Sequence[typing.SupportsInt], types: collections.abc.Sequence[typing.SupportsInt], x: collections.abc.Sequence[typing.SupportsFloat], y: collections.abc.Sequence[typing.SupportsFloat], z: collections.abc.Sequence[typing.SupportsFloat], nx: collections.abc.Sequence[typing.SupportsInt], ny: collections.abc.Sequence[typing.SupportsInt], nz: collections.abc.Sequence[typing.SupportsInt]) -> None:
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
    def add_bonds(self, bonds_from: collections.abc.Sequence[typing.SupportsInt], bonds_to: collections.abc.Sequence[typing.SupportsInt]) -> None:
        """
                 Add bonds to the underlying atoms, edges to the underlying graph. 
                 If the connected atoms are not found, the bonds are silently skipped.
                 
                 :param bonds_from: Vector of atom IDs for bond start points
                 :param bonds_to: Vector of atom IDs for bond end points
        """
    @typing.overload
    def add_bonds(self, nr_of_bonds: typing.SupportsInt, bonds_from: collections.abc.Sequence[typing.SupportsInt], bonds_to: collections.abc.Sequence[typing.SupportsInt], bond_types: collections.abc.Sequence[typing.SupportsInt], ignore_non_existent_atoms: bool = False, simplify_universe: bool = True) -> None:
        """
                 Add bonds to the underlying atoms, edges to the underlying graph.
                 
                 :param nr_of_bonds: Number of bonds to add
                 :param bonds_from: Vector of atom IDs for bond start points
                 :param bonds_to: Vector of atom IDs for bond end points
                 :param bond_types: Vector of bond types
                 :param ignore_non_existent_atoms: Whether to skip bonds to non-existent atoms
                 :param simplify_universe: Whether to simplify the universe after adding bonds
        """
    def add_bonds_with_types(self, bonds_from: collections.abc.Sequence[typing.SupportsInt], bonds_to: collections.abc.Sequence[typing.SupportsInt], bond_types: collections.abc.Sequence[typing.SupportsInt]) -> None:
        """
                 Add bonds to the underlying atoms, edges to the underlying graph. 
                 If the connected atoms are not found, the bonds are silently skipped.
                 
                 :param bonds_from: Vector of atom IDs for bond start points
                 :param bonds_to: Vector of atom IDs for bond end points
                 :param bond_types: Vector of bond types
        """
    def add_dihedral_angles(self, angles_from: collections.abc.Sequence[typing.SupportsInt], angles_via1: collections.abc.Sequence[typing.SupportsInt], angles_via2: collections.abc.Sequence[typing.SupportsInt], angles_to: collections.abc.Sequence[typing.SupportsInt], angle_types: collections.abc.Sequence[typing.SupportsInt]) -> None:
        """
                 Add dihedral angles to the Universe. No relation to the underlying graph, 
                 just a method to preserve read & write capabilities.
                 
                 :param angles_from: Vector of atom IDs for dihedral start points
                 :param angles_via1: Vector of atom IDs for first middle points
                 :param angles_via2: Vector of atom IDs for second middle points
                 :param angles_to: Vector of atom IDs for dihedral end points
                 :param angle_types: Vector of dihedral angle types
        """
    def compute_angles(self) -> list[float]:
        """
                 Compute the angle of each angle in the molecule, respecting periodic boundaries.
                 
                 :return: A list of angle values in radians
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
    def compute_dxs(self, atom_ids_to: collections.abc.Sequence[typing.SupportsInt], atom_ids_from: collections.abc.Sequence[typing.SupportsInt]) -> list[float]:
        """
                 Compute the dx distance for certain bonds (length in x direction).
                 
                 :param atom_ids_to: Vector of destination atom IDs
                 :param atom_ids_from: Vector of source atom IDs
                 :return: Vector of x-direction distances
        """
    def compute_dys(self, atom_ids_to: collections.abc.Sequence[typing.SupportsInt], atom_ids_from: collections.abc.Sequence[typing.SupportsInt]) -> list[float]:
        """
                 Compute the dy distance for certain bonds (length in y direction).
                 
                 :param atom_ids_to: Vector of destination atom IDs
                 :param atom_ids_from: Vector of source atom IDs
                 :return: Vector of y-direction distances
        """
    def compute_dzs(self, atom_ids_to: collections.abc.Sequence[typing.SupportsInt], atom_ids_from: collections.abc.Sequence[typing.SupportsInt]) -> list[float]:
        """
                 Compute the dz distance for certain bonds (length in z direction).
                 
                 :param atom_ids_to: Vector of destination atom IDs
                 :param atom_ids_from: Vector of source atom IDs
                 :return: Vector of z-direction distances
        """
    def compute_end_to_end_distances(self, crosslinker_type: typing.SupportsInt, derive_image_flags: bool = False) -> list[float]:
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
    def compute_mean_end_to_end_distance(self, crosslinker_type: typing.SupportsInt, derive_image_flags: bool = False) -> float:
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
    def compute_mean_squared_end_to_end_distance(self, crosslinker_type: typing.SupportsInt, only_those_with_two_crosslinkers: bool = False, derive_image_flags: bool = False) -> float:
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
    def compute_mean_strand_length(self, crosslinker_type: typing.SupportsInt) -> float:
        """
                      Compute the mean number of beads per strand.
                      
                      :param crosslinker_type: The type of crosslinker atoms
                      :return: Mean strand length
        """
    def compute_number_average_molecular_weight(self, crosslinker_type: typing.SupportsInt) -> float:
        """
                      Compute the number average molecular weight.
        
                      .. note:: 
                            Crosslinkers are ignored completely.
                            
                      :param crosslinker_type: The type of crosslinker atoms
                      :return: Number average molecular weight
        """
    def compute_polydispersity_index(self, crosslinker_type: typing.SupportsInt) -> float:
        """
                      Compute the polydispersity index: 
                      the weight average molecular weight over the number average molecular weight.
        """
    def compute_temperature(self, dimensions: typing.SupportsInt = 3, k_b: typing.SupportsFloat = 1.0) -> float:
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
    def compute_weight_average_molecular_weight(self, crosslinker_type: typing.SupportsInt) -> float:
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
    def contract_vertices_along_bond_type(self, bond_type: typing.SupportsInt) -> Universe:
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
    def count_atoms_in_skin_distance(self, distances: collections.abc.Sequence[typing.SupportsFloat], unwrapped: bool = False) -> list[int]:
        """
                  This is a function that may help you to compute the radial distribution function.
                  It loops through all atoms and counts neighbors within specified distance ranges.
        
                  :param distances: The edges of the bins for distance counting
                  :param unwrapped: Whether to measure the distance in unwrapped coordinates or as PBC-corrected distance
                  :return: Array of counts for each distance bin
        """
    def count_loop_lengths(self, max_length: typing.SupportsInt = -1) -> dict[int, int]:
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
    def find_loops(self, crosslinker_type: typing.SupportsInt, max_length: typing.SupportsInt = -1, skip_self_loops: bool = False) -> dict[int, list[list[Atom]]]:
        """
                    Decompose the Universe into loops.
                    The primary index specifies the degree of the loop.
        
                    CAUTION:
                       There are exponentially many paths between two crosslinkers of a network,
                       and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
                       You can use the `max_length` parameter to restrict the algorithm to only search for loops up to a certain length.
                       Use a negative value to find all loops and paths.
        """
    def find_minimal_order_loop_from(self, loop_start: typing.SupportsInt, loop_step1: typing.SupportsInt, max_length: typing.SupportsInt = -1, skip_self_loops: bool = False) -> list[Atom]:
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
    def get_atom(self, atom_id: typing.SupportsInt) -> Atom:
        """
                 Find an atom by its ID.
                 
                 :param atom_id: The ID of the atom to find
                 :return: The atom with the specified ID
        """
    def get_atom_by_vertex_id(self, vertex_id: typing.SupportsInt) -> Atom:
        """
              Find an atom by the ID of the vertex of the underlying graph.
              
              :param vertex_id: The vertex ID to query
              :return: The atom at the specified vertex
        """
    def get_atom_id_by_vertex_idx(self, vertex_id: typing.SupportsInt) -> int:
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
    def get_atoms_by_degree(self, functionality: typing.SupportsInt) -> list[Atom]:
        """
                    Get the atoms that have the specified number of bonds.
        """
    def get_atoms_by_type(self, atom_type: typing.SupportsInt) -> list[Atom]:
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
    def get_atoms_connected_to_vertex(self, vertex_idx: typing.SupportsInt) -> list[Atom]:
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
    def get_chains_with_crosslinker(self, crosslinker_type: typing.SupportsInt) -> list[Molecule]:
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
    def get_edge_ids_from(self, vertex_id: typing.SupportsInt) -> list[int]:
        """
                  Get the edge IDs incident to a specific vertex.
        
                  :param vertex_id: The ID of the vertex to query
                  :return: List of edge IDs connected to the vertex
        """
    def get_edge_ids_from_to(self, vertex_id_from: typing.SupportsInt, vertex_id_to: typing.SupportsInt) -> list[int]:
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
    def get_molecules(self, atom_type_to_omit: typing.SupportsInt) -> list[Molecule]:
        """
                  Decompose the Universe into molecules, which could be either chains, networks, or even lonely atoms.
                  
                  Reduces the Universe to a list of molecules. 
                  Specify the crosslinker_type to an existing type ID, 
                  then those atoms will be omitted, and this function returns chains instead.
        
                  :param atom_type_to_omit: The type of atom to omit from the universe to end up with the desired molecules (e.g., the type of the crosslinkers).
                  :return: List of molecules
        """
    def get_network_of_crosslinker(self, crosslinker_type: typing.SupportsInt) -> Universe:
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
    def get_nr_of_bonds_of_atom(self, arg0: typing.SupportsInt) -> int:
        """
               Count the number of immediate neighbors of an atom, specified by its ID.
        
               :param atom_id: The ID of the atom to query
               :return: Number of bonds connected to the atom
        """
    def get_nr_of_bonds_of_vertex(self, arg0: typing.SupportsInt) -> int:
        """
               Count the number of immediate neighbors of an atom, specified by its vertex ID.
        
               :param vertex_id: The vertex ID of the atom to query
               :return: Number of bonds connected to the vertex
        """
    def get_nr_of_dihedral_angles(self) -> int:
        """
                    Query the number of dihedral angles that have been added to this universe.
        """
    def get_nr_of_edges_from_to(self, vertex_id_from: typing.SupportsInt, vertex_id_to: typing.SupportsInt, max_length: typing.SupportsInt = -1) -> int:
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
    def get_vertex_idx_by_atom_id(self, atom_id: typing.SupportsInt) -> int:
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
    def has_atom_with_id(self, atom_id: typing.SupportsInt) -> bool:
        """
                 Check whether this universe contains an atom with the specified ID.
        
                 :param atom_id: The atom ID to check
                 :return: True if the atom exists in this universe, False otherwise
        """
    def has_infinite_strand(self, arg0: typing.SupportsInt, arg1: typing.SupportsInt) -> bool:
        """
                   Check whether there is a strand (with crosslinker) in the universe that loops through periodic images without coming back.
                   
                    .. warning::
                       There are exponentially many paths between two crosslinkers of a network,
                       and you may run out of memory when using this function, if your Universe/Network is lattice-like. 
                       
                   :return: True if infinite strands are detected, False otherwise
        """
    def hash_angle_type(self, angle_from: typing.SupportsInt, angle_via: typing.SupportsInt, angle_to: typing.SupportsInt) -> int:
        """
                  Convert the three integers to one long number/hash.
                  Used internally for duplicate detection.
                  
                  :param angle_from: First atom ID in the angle
                  :param angle_via: Middle atom ID in the angle
                  :param angle_to: Last atom ID in the angle
                  :return: Hash value for the angle type
        """
    def hash_dihedral_angle_type(self, angle_from: typing.SupportsInt, angle_via1: typing.SupportsInt, angle_via2: typing.SupportsInt, angle_to: typing.SupportsInt) -> int:
        """
                  Convert the four integers to one long number/hash.
                  Used internally for duplicate detection.
                  
                  :param angle_from: First atom ID in the dihedral
                  :param angle_via1: Second atom ID in the dihedral
                  :param angle_via2: Third atom ID in the dihedral
                  :param angle_to: Fourth atom ID in the dihedral
                  :return: Hash value for the dihedral angle type
        """
    def interpolate_edges(self, crosslinker_type: typing.SupportsInt, interpolation_factor: typing.SupportsFloat) -> list[tuple[int, int]]:
        """
                  Get more or less edges than currently present, interpolating between junctions.
                  
                  :param crosslinker_type: The type of crosslinker atoms
                  :param interpolation_factor: Factor for edge interpolation
                  :return: Interpolated edge structure
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
    def remove_atoms(self, atom_ids: collections.abc.Sequence[typing.SupportsInt]) -> None:
        """
                  Remove atoms and all associated bonds by their atom IDs. 
                  
                  :param atom_ids: Vector of atom IDs to remove
        """
    def remove_bonds(self, bonds_from: collections.abc.Sequence[typing.SupportsInt], bonds_to: collections.abc.Sequence[typing.SupportsInt]) -> None:
        """
                  Remove bonds by their connected atom IDs. 
                  
                  :param bonds_from: Vector of starting atom IDs
                  :param bonds_to: Vector of ending atom IDs
        """
    def remove_bonds_by_type(self, bond_type: typing.SupportsInt) -> None:
        """
                  Remove bonds with a specific type. 
                  
                  :param bond_type: The bond type to remove
        """
    def replace_atom(self, atom_id: typing.SupportsInt, replacement_atom: Atom) -> None:
        """
                  Replace the properties of an atom with the properties of another given atom.
                  
                  :param atom_id: ID of the atom to replace
                  :param replacement_atom: The atom with new properties
        """
    def replace_atom_type(self, atom_id: typing.SupportsInt, new_type: typing.SupportsInt) -> None:
        """
                  Replace the type of an atom with another type.
                  
                  :param atom_id: ID of the atom to modify
                  :param new_type: The new atom type
        """
    def resample_velocities(self, mean: typing.SupportsFloat, variance: typing.SupportsFloat, seed: str = '', is_2d: bool = False) -> None:
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
    def set_box_lengths(self, lx: typing.SupportsFloat, ly: typing.SupportsFloat, lz: typing.SupportsFloat, rescale_atoms: bool = False) -> None:
        """
                  Override the currently assigned box with one with the side lengths specified.
                  
                  :param lx: Length in x direction
                  :param ly: Length in y direction
                  :param lz: Length in z direction
                  :param rescale_atoms: Whether to rescale atom positions
        """
    def set_mass(self, atom_type: typing.SupportsInt, mass: typing.SupportsFloat) -> None:
        """
                 Set the mass for a specific atom type.
                 
                 :param atom_type: The atom type to set mass for
                 :param mass: The mass value to assign
        """
    def set_masses(self, mass_per_type: collections.abc.Mapping[typing.SupportsInt, typing.SupportsFloat]) -> None:
        """
                 Set the mass per type of atom.
                 
                 :param mass_per_type: Map of atom types to their masses
        """
    def set_timestep(self, timestep: typing.SupportsInt) -> None:
        """
                 Set the timestep when this Universe was captured.
                 
                 :param timestep: The timestep value
        """
    def set_vertex_property(self, vertex_id: typing.SupportsInt, property_name: str, value: typing.SupportsFloat) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt) -> Universe:
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
    def at_index(self, index: typing.SupportsInt) -> Universe:
        """
                 Get the Universe at the given index (as of in the sequence given 
                 by the dump file).
                 
                 :param index: The index of the universe to retrieve
                 :return: The Universe at the specified index
        """
    def compute_distance_autocorrelation_from_to(self, atom_ids_from: collections.abc.Sequence[typing.SupportsInt], atom_ids_to: collections.abc.Sequence[typing.SupportsInt], nr_of_origins: typing.SupportsInt = 25, reduce_memory: bool = False) -> dict[int, float]:
        """
                  Compute the autocorrelation of the dot product of the distance vector from certain to other atoms.
        
                  For example, this can be used to compute Eq. 4.51 from Masao Doi, Introduction to Polymer Physics, p. 74.
        """
    def compute_distance_from_to_atoms(self, atom_ids_from: collections.abc.Sequence[typing.SupportsInt], atom_ids_to: collections.abc.Sequence[typing.SupportsInt], reduce_memory: bool = False) -> list[float]:
        """
                  Compute the root square norm of all the (unwrapped!) distances for the given pair of atoms.
                  
                  Can be used to somewhat faster compute e.g. all the end-to-end or bond distances.
                  Pay attention that the image flags are correct, otherwise, this data may not be useable.
        """
    def compute_msd_for_atom_properties(self, atom_ids: collections.abc.Sequence[typing.SupportsInt], x_property: str, y_property: str, z_property: str, nr_of_origins: typing.SupportsInt = 25, reduce_memory: bool = False, max_tau: typing.SupportsInt = -1) -> dict[int, float]:
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
    def compute_msd_for_atoms(self, atom_ids: collections.abc.Sequence[typing.SupportsInt], nr_of_origins: typing.SupportsInt = 25, reduce_memory: bool = False, max_tau: typing.SupportsInt = -1) -> dict[int, float]:
        """
                  Compute the mean square displacement for atoms with the specified IDs.
                  
                  :param atom_ids: List of atom IDs for which to compute the MSD
                  :param nr_of_origins: Number of time origins to use for averaging. Higher values provide better statistics but increase computation time (default: 25)
                  :param reduce_memory: If True, reduces memory usage by forgetting universes after processing them (default: False)
                  :param max_tau: Maximum time lag (tau) to compute. If -1, computes for all possible tau values. For better statistics, consider setting this to approximately half the sequence length (default: -1)
                  :return: Dictionary mapping time lag (tau) to mean square displacement values
        """
    def compute_vector_from_to_atoms(self, atom_ids_from: collections.abc.Sequence[typing.SupportsInt], atom_ids_to: collections.abc.Sequence[typing.SupportsInt], reduce_memory: bool = False) -> list[typing.Annotated[numpy.typing.NDArray[numpy.float64], "[3, 1]"]]:
        """
                  Compute the (unwrapped!) distances for the given pair of atoms.
                  
                  Can be used to somewhat faster compute e.g. all the end-to-end or bond vectors.
                  Pay attention that the image flags are correct, otherwise, this data may not be useable.
        """
    def forget_at_index(self, index: typing.SupportsInt) -> None:
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
    def __init__(self, std_multiplier: typing.SupportsFloat, in_sqrt_multiplier: typing.SupportsFloat) -> None:
        ...
    def get_max_distance(self, N: typing.SupportsFloat) -> float:
        """
                 Get the maximum distance for a given N.
        
                 :param N: Number of segments.
        """
def do_linear_walk_chain_from_to(box: Box, from_coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], to_coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], chain_len: typing.SupportsInt, include_ends: bool = False) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
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
def do_random_walk(chain_len: typing.SupportsInt, bead_distance: typing.SupportsFloat = 1.0, mean_squared_bead_distance: typing.SupportsFloat = 1.0, seed: str = '') -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
    """
                Do a random walk, return the coordinates of each point visited.
    
                :param chain_len: Length of the chain to generate.
                :param bead_distance: Mean distance between consecutive beads (default: 1.0).
                :param mean_squared_bead_distance: Mean squared distance between consecutive beads (default: 1.0).
                :param seed: Random seed for reproducibility (default: empty string for random seed).
                :return: Coordinates of each point as a flat array (i.e., [x1, y1, z1, x2, y2, z2, ...]).
    """
def do_random_walk_chain_from_to(box: Box, from_coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], to_coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], chain_len: typing.SupportsInt, bead_distance: typing.SupportsFloat = 1.0, mean_squared_bead_distance: typing.SupportsFloat = 1.0, seed: str = '') -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
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
def do_random_walk_chain_from_to_mc(box: Box, from_coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], to_coordinates: typing.Annotated[numpy.typing.ArrayLike, numpy.float64, "[3, 1]"], chain_len: typing.SupportsInt, bead_distance: typing.SupportsFloat = 1.0, mean_squared_bead_distance: typing.SupportsFloat = 1.0, seed: str = '', n_iterations: typing.SupportsInt = 10000) -> typing.Annotated[numpy.typing.NDArray[numpy.float64], "[m, 1]"]:
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
def inverse_langevin(x: typing.SupportsFloat) -> float:
    """
         A somewhat accurate (for :math:`x \\in (-1, 1)`) implementation of the inverse Langevin function.
    
         Source: https://scicomp.stackexchange.com/a/30251
    
         :param x: Input value in the range (-1, 1)
         :return: Inverse Langevin function value
    """
def randomly_sample_entanglements(universe: Universe, nr_of_samples: typing.SupportsInt, upper_cutoff: typing.SupportsFloat, lower_cutoff: typing.SupportsFloat = 0, minimum_nr_of_samples: typing.SupportsInt = 0, same_strand_cutoff: typing.SupportsFloat = 3.0, seed: str = '', crosslinker_type: typing.SupportsInt = 2, ignore_crosslinks: bool = True, filter_dangling_and_soluble: bool = False) -> AtomPairEntanglements:
    """
        Randomly find pairs of atoms that are close together and could be
        entanglements
    
        :param universe: The universe of atoms from which to sample entanglements from.
        :param nr_of_samples: The number of pairs of atoms to randomly sample.
        :param upper_cutoff: The maximum distance between atoms for a pair to be considered a potential entanglement.
        :param lower_cutoff: The minimum distance between atoms for a pair to be considered a potential entanglement.
        :param minimum_nr_of_samples: The minimum number of entanglements to be found.
        :param same_strand_cutoff: The maximum distance between atoms on the same strand for a pair to be considered a potential entanglement.
        :param seed: A seed for the random number generator.
        :param crosslinker_type: The type of crosslinker to consider when finding entanglements. Used for the splitting into strands.
        :param ignore_crosslinks: Whether to ignore crosslinks when finding entanglements. Careful: if you don't ignore them, the same-strand policy might not work correctly, since each crosslink should actually be associated with more than one strand.
        :param filter_dangling_and_soluble: Whether to filter out dangling chains and soluble crosslinks when finding entanglements.
          This means, entanglements involving an obviously (1st order) dangling or soluble chain are 
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
__version__: str = '0.3.9'
