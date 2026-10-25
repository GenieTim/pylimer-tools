<a id="sphx-glr-auto-examples-readers-writers"></a>

# File I/O: Readers & Writers

pylimer-tools provides comprehensive support for reading and writing various LAMMPS file formats, enabling seamless integration with molecular dynamics workflows.

The library supports four main file I/O operations:

1. **Reading LAMMPS data files** - Complete system configurations
2. **Reading LAMMPS dump files** - Trajectory and snapshot data
3. **Reading LAMMPS output files** - Measurements and simulation results
4. **Writing LAMMPS data files** - Modified or generated systems

<div id='sg-tag-list' class='sphx-glr-tag-list'></div><div class="sphx-glr-thumbnails">
<!-- thumbnail-parent-div-open --><div class="sphx-glr-thumbcontainer" tooltip="Convert LAMMPS units to SI or other unit systems:">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_unit_conversion_thumb.png)

[Unit Conversion](plot_unit_conversion.md)

  <div class="sphx-glr-thumbnail-title">Unit Conversion</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read output from pylimer-tools&#x27; own simulators.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_pylimer_tools_output_reader_thumb.png)

[pylimer-tools Output Reader](plot_pylimer_tools_output_reader.md)

  <div class="sphx-glr-thumbnail-title">pylimer-tools Output Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read correlation functions from fix ave/correlate and fix ave/correlate/long commands.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_correlated_averages_reader_thumb.png)

[Correlated Averages Reader](plot_correlated_averages_reader.md)

  <div class="sphx-glr-thumbnail-title">Correlated Averages Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Process histogram data from fix ave/hist commands.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_histogram_reader_thumb.png)

[Histogram Reader](plot_histogram_reader.md)

  <div class="sphx-glr-thumbnail-title">Histogram Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="All readers support automatic caching for improved performance.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_performance_caching_thumb.png)

[Performance and Caching](plot_performance_caching.md)

  <div class="sphx-glr-thumbnail-title">Performance and Caching</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Take a later entry from a LAMMPS trajectory and save it as a data file.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_convert_trajectory_to_frame_thumb.png)

[Extract Trajectory Frame to Data File](plot_convert_trajectory_to_frame.md)

  <div class="sphx-glr-thumbnail-title">Extract Trajectory Frame to Data File</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read LAMMPS log files containing thermodynamic output from the thermo command.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_log_thermo_file_reader_thumb.png)

[Log & Thermo File Reader](plot_log_thermo_file_reader.md)

  <div class="sphx-glr-thumbnail-title">Log & Thermo File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The DumpFileReader processes LAMMPS trajectory files. Is has a practical interface to read frames from a trajectory file and access them as Universe objects in the UniverseSequence class.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_dump_file_reader_thumb.png)

[Dump File Reader](plot_dump_file_reader.md)

  <div class="sphx-glr-thumbnail-title">Dump File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="Read output from LAMMPS fix ave/time commands.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_time_averaged_data_reader_thumb.png)

[Time-Averaged Data Reader](plot_time_averaged_data_reader.md)

  <div class="sphx-glr-thumbnail-title">Time-Averaged Data Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="LAMMPS bond swapping can significantly speed up equilibration for polymer networks. However, to preserve the strand lengths, you need to ensure that the molecule indices are set correctly. Here&#x27;s how you could re-set the molecule indices for bond swapping with constant chain lengths.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_convert_for_bond_swapping_thumb.png)

[Convert for Bond Swapping](plot_convert_for_bond_swapping.md)

  <div class="sphx-glr-thumbnail-title">Convert for Bond Swapping</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The DataFileParser handles LAMMPS data files containing complete system definitions. The UniverseSequence provides a convenient way to have these data files read into the Universe objects.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_data_file_reader_thumb.png)

[Data File Reader](plot_data_file_reader.md)

  <div class="sphx-glr-thumbnail-title">Data File Reader</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The UniverseSequence class provides automatic memory management.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_memory_efficient_dump_file_reading_thumb.png)

[Memory-Efficient Dump File Reading](plot_memory_efficient_dump_file_reading.md)

  <div class="sphx-glr-thumbnail-title">Memory-Efficient Dump File Reading</div>
</div><div class="sphx-glr-thumbcontainer" tooltip="The DataFileWriter generates LAMMPS-compatible data files.">![](auto_examples/readers_writers/images/thumb/sphx_glr_plot_data_file_writer_thumb.png)

[Data File Writer](plot_data_file_writer.md)

  <div class="sphx-glr-thumbnail-title">Data File Writer</div>
</div>
<!-- thumbnail-parent-div-close --></div>
