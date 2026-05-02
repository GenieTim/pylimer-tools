<a id="sphx-glr-auto-examples-readers-writers"></a>

# File I/O: Readers & Writers

pylimer-tools provides comprehensive support for reading and writing various LAMMPS file formats, enabling seamless integration with molecular dynamics workflows.

The library supports four main file I/O operations:

1. **Reading LAMMPS data files** - Complete system configurations
2. **Reading LAMMPS dump files** - Trajectory and snapshot data
3. **Reading LAMMPS output files** - Measurements and simulation results
4. **Writing LAMMPS data files** - Modified or generated systems

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="Read output from pylimer-tools&#x27; own simulators.">  <div class="sphx-glr-thumbnail-title">pylimer-tools Output Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read correlation functions from fix ave/correlate and fix ave/correlate/long commands.">  <div class="sphx-glr-thumbnail-title">Correlated Averages Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Process histogram data from fix ave/hist commands.">  <div class="sphx-glr-thumbnail-title">Histogram Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="All readers support automatic caching for improved performance.">  <div class="sphx-glr-thumbnail-title">Performance and Caching</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Convert LAMMPS units to SI or other unit systems:">  <div class="sphx-glr-thumbnail-title">Unit Conversion</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read LAMMPS log files containing thermodynamic output from the thermo command.">  <div class="sphx-glr-thumbnail-title">Log & Thermo File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Take a later entry from a LAMMPS trajectory and save it as a data file.">  <div class="sphx-glr-thumbnail-title">Extract Trajectory Frame to Data File</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The DumpFileReader processes LAMMPS trajectory files. Is has a practical interface to read frames from a trajectory file and access them as Universe objects in the UniverseSequence class.">  <div class="sphx-glr-thumbnail-title">Dump File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read output from LAMMPS fix ave/time commands.">  <div class="sphx-glr-thumbnail-title">Time-Averaged Data Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="LAMMPS bond swapping can significantly speed up equilibration for polymer networks. However, to preserve the strand lengths, you need to ensure that the molecule indices are set correctly. Here&#x27;s how you could re-set the molecule indices for bond swapping with constant chain lengths.">  <div class="sphx-glr-thumbnail-title">Convert for Bond Swapping</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The DataFileParser handles LAMMPS data files containing complete system definitions. The UniverseSequence provides a convenient way to have these data files read into the Universe objects.">  <div class="sphx-glr-thumbnail-title">Data File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The UniverseSequence class provides automatic memory management.">  <div class="sphx-glr-thumbnail-title">Memory-Efficient Dump File Reading</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The DataFileWriter generates LAMMPS-compatible data files.">  <div class="sphx-glr-thumbnail-title">Data File Writer</div>
</div>
<!-- thumbnail-parent-div-close --></div>
