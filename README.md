smalldisco
==========

<img src="misc/logo.png" width="60">

smalldisco discovers putative siRNA regions in the genome, as well as quantification of small RNA tails.


### Installation

smalldisco is implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow. All of its dependencies are listed in the `environment.yaml` file.

We recommend using the [conda package manager](https://docs.conda.io/en/latest/) to install requirements and run the program. To create a new conda environment called `smalldisco` with all the dependencies, run:

    conda env create -f environment.yaml

You can then activate the environment with `conda activate smalldisco`, run the program, and deactivate the environment with `conda deactivate` when you're done.

Check the [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) for more information on conda environments.


### Running

To run, first you need to set up your analysis. This is done in the `config.yaml`, and we recommend following the instructions in that file. After the configuration is done, run smalldisco with:

    python smalldisco.py <mode> --configfile=<config> [optional params]
	
Where `<mode>` can be one of:

* `align`: Use HISAT2 to align FASTQ input to a reference genome.
* `sirna`: Produce a list of putative siRNA regions.
* `tail`: Quantify read tails for a given list of sRNA regions.
* `run_all`: Run all the above steps.

Specifying a configuration file is mandatory. If you've filled in the provided configuration with the parameters of your analysis, run with `--configfile=config.yaml`.

The optional parameters are passed along to Snakemake; see below.


### Optional Snakemake parameters

Since smalldisco is implemented as a Python script that invokes a Snakemake workflow, it benefits from all the advantages that Snakemake provides. It will continue a run from an intermediate set out outputs if those are available, and will pass along any extra parameters to the Snakemake command line.

All Snakemake command line options are described in its [documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html). The most useful ones for running smalldisco are:

* `-cX`: Run with `X` threads.
* `-n`: Perform a dry run, indicating what steps of the pipeline would be run, without running them.
* `--notemp`: Keep all intermediate output, including large BAM and SAM files made by the pipeline.
* `--config KEY=VALUE`: Overwrite the value of the `KEY` parameter in the configuration file with `VALUE`.
