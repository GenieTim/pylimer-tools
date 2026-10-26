# Command Line Interface

Apart from the Python API, pylimer-tools provides a command line interface (CLI) for various functionalities.
This allows you to perform tasks directly from the terminal without writing Python scripts.

## Network Generation

Generate crosslinked polymer networks using Monte Carlo procedures.

### pylimer-generate-network

Generate a crosslinked polymer network using Monte Carlo procedure.

This command generates a crosslinked polymer network with specified parameters.
You must specify either –target-p (crosslinker conversion) or –target-wsol
(soluble fraction), but not both.

Example:
: pylimer-generate-network –polymer-name PDMS –parameter-type GAUSSIAN –target-p 0.8 –n-chains-crosslinks 100 –target-file network.dat

### Usage

```shell
pylimer-generate-network [OPTIONS]
```

### Options

### --polymer-name <polymer_name>

Name of the polymer to generate the network for

* **Options:**
  apoa | cispbd | pbd98 | pi20 | pe413 | apea | pi75 | app348 | apams | pc | peo | pe | app463 | apma | apmma | pet | pib | apva | cispi | ips | pi50 | aps | ipp | pi7 | pib413 | pbd20 | aphma | pom | sbr | spp | p6n | app | ptfe | app413 | pi34 | pdms

### --parameter-type <parameter_type>

Type of parameters to use for the polymer model

* **Options:**
  GAUSSIAN | KG_LJ | KUHN

### --target-p <target_p>

Target crosslinker conversion (p). Mutually exclusive with –target-wsol.

### --target-wsol <target_wsol>

Target soluble fraction (w_sol). Mutually exclusive with –target-p.

### --n-chains-crosslinkers <n_chains_crosslinkers>

Number of crosslinkers

### --target-f <target_f>

Crosslinker functionality

### --n-solvent-chains <n_solvent_chains>

Number of solvent chains

### --n-beads-per-solvent-chain <n_beads_per_solvent_chain>

Beads per solvent chain

### --n-chains-1 <n_chains_1>

Number of chains

### --n-beads-per-chain-1 <n_beads_per_chain_1>

Beads per chain

### --n-chains-2 <n_chains_2>

Number of other chains

### --n-beads-per-chain-2 <n_beads_per_chain_2>

Beads per other chain

### --n-mono-chains <n_mono_chains>

Number of monofunctional chains

### --n-mono-beads-per-chain <n_mono_beads_per_chain>

Beads per monofunctional chain

### --n-beads-per-xlink <n_beads_per_xlink>

Beads per crosslinker (xlink)

### --remove-wsol

Remove soluble fraction after generation

### --functionalize-discrete, --no-functionalize-discrete

Functionalize crosslinks discrete. Applicable to crosslinkers longer than one bead. If enabled, the crosslinker functionality is distributed in discrete chunks. If disabled, the crosslinker functionality is distributed randomly across the crosslinker beads.

### --disable-primary-loops

Disable primary loops in the network

### --disable-secondary-loops

Disable secondary loops in the network

### --z-score-std-mult <z_score_std_mult>

Z-score standard deviation multiplier for max distance

### --target-file <target_file>

**Required** Output file for the generated network

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
  apoa | cispbd | pbd98 | pi20 | pe413 | apea | pi75 | app348 | apams | pc | peo | pe | app463 | apma | apmma | pet | pib | apva | cispi | ips | pi50 | aps | ipp | pi7 | pib413 | pbd20 | aphma | pom | sbr | spp | p6n | app | ptfe | app413 | pi34 | pdms

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
