# Command Line Interface

Apart from the Python API, pylimer-tools provides a command line interface (CLI) for various functionalities.
This allows you to perform tasks directly from the terminal without writing Python scripts.

## Network Generation

Generate crosslinked polymer networks using Monte Carlo procedures.

## Network Analysis

Analyze existing polymer networks to compute structural properties.

### pylimer-analyse-networks

Basic CLI application reading all passed files, outputting some stats on the structures therein

Arguments:
: - files: list of files to read

### Usage

```shell
pylimer-analyse-networks [OPTIONS] [FILES]...
```

### Options

### --crosslinker-type <crosslinker_type>

### --polymer-name <polymer_name>

Name of the polymer to use for parameter retrieval.

* **Options:**
  pi7 | pib413 | pi75 | pe | aphma | app463 | apams | apma | app | app348 | sbr | pi50 | pc | apva | peo | pi34 | ipp | aps | spp | apoa | ips | p6n | cispbd | pdms | pet | pib | pbd20 | pe413 | apmma | pi20 | pbd98 | apea | pom | cispi | ptfe | app413

### Arguments

### FILES

Optional argument(s)

## LAMMPS Data Statistics

Compute basic statistics from LAMMPS data files.

### pylimer-basic-lammps-stats

Basic CLI application reading all passed files, outputting some stats on the structures therein

Arguments:
: - files: list of files to read

### Usage

```shell
pylimer-basic-lammps-stats [OPTIONS] [FILES]...
```

### Arguments

### FILES

Optional argument(s)

## Random Displacement

Randomly displace atoms in a structure file.

### pylimer-displace-randomly

Basic CLI application iterating all atoms in a file, displacing them by a bit.

Arguments:
: - file: The file to read (and write, with prefix “random-displaced-“)
  - max_displacement: The maximum displacement

### Usage

```shell
pylimer-displace-randomly [OPTIONS] FILE [MAX_DISPLACEMENT]
```

### Arguments

### FILE

Required argument

### MAX_DISPLACEMENT

Optional argument
