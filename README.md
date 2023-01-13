smalldisco
==========

<img src="misc/logo.png" width="60">

smalldisco discovers putative siRNA regions in the genome, as well as quantification of small RNA tails.

Check out our preprint at https://www.biorxiv.org/content/10.1101/2022.07.15.500275v1.

### Installation

Installation steps:

1. Clone the smalldisco repository
2. Install Tailor
3. Create a conda environment for running smalldisco

First, clone this repository to your machine:

```console
$ git clone https://github.com/ianvcaldas/smalldisco.git
$ cd smalldisco
```

Smalldisco depends on the [Tailor](https://github.com/jhhung/Tailor) program for detecting read tails. By default, the Tailor repository is set up as a [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) of smalldisco. To initialize it, run this command from within the smalldisco folder:

```console
$ git submodule update --init
```

This will clone the Tailor repository as a folder. If you're running on Linux, an executable binary is available in `Tailor/bin/tailor_v1.1_linux_static`. If you're running on a different operating system, you will have to compile Tailor from source yourself; see its documentation for instructions.

To run smalldisco, you will need to install its software dependencies. We provide the file `environment.yaml` with all the dependencies, and we recommend using the [conda package manager](https://docs.conda.io/en/latest/) to install them in an isolated environment. You can then activate the environment with `conda activate smalldisco`, run the program, and deactivate the environment with `conda deactivate` when you're done. Check the [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) for more information on conda environments. To create a new conda environment called `smalldisco`, run:

```console
$ conda env create -f environment.yaml
```

To test that the installation worked, run the `--help` command in the conda environment. If you see the smalldisco help, starting with its usage string, the software has been installed successfully.

```console
$ conda activate smalldisco
$ python smalldisco.py --help

Usage: smalldisco.py [OPTIONS] COMMAND [ARGS]...
...
```

### Usage

Use `python smalldisco.py --help` for main usage instructions. Smalldisco has two commands, `sirna` and `tail`, whose usage can be checked with `python smalldisco.py sirna --help` and `python smalldisco.py tail --help`, respectively.

Behind the scenes, `smalldisco.py` is a wrapper script. Smalldisco is actually implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. Workflow definition and helper scripts are in the `workflow` folder, and the program will not work if `smalldisco.py` and `workflow` are not in the same directory.

**Tailor integration:** When running the `tail` command, we assume by default that the path to the Tailor executable is `Tailor/bin/tailor_v1.1_linux_static`. This assumes that you are on Linux, are in the smalldisco repository, (e.g. `cd smalldisco`), and the Tailor submodule has been initialized as described in the installation section. If you are running smalldisco from a different folder, or have a custom Tailor installation, you must specify a path to a valid Tailor executable, for instance:

```console
$ python smalldisco.py tail --tailor_command /usr/local/bin/tailor
```
