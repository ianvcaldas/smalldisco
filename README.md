smalldisco
==========

<p float="left">
  <img src="misc/logo.png" width="200" />
  <img src="misc/logo2.jpg" width="300" /> 
</p>


smalldisco discovers putative siRNA regions in the genome and quantifies tails of small RNA reads.

Check out our preprint at https://www.biorxiv.org/content/10.1101/2022.07.15.500275v1 for details on how it works and examples of use cases.

You can also check out the smalldisco vingette located in the Wiki tab.

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

This will clone the Tailor repository as a folder. If you're running on Linux, you should be all set: an executable binary is available in `Tailor/bin/tailor_v1.1_linux_static`. If you're running on a different operating system, however, you will have to compile it from source yourself. See the Tailor documentation for instructions on how to do that.

To run smalldisco, you will need to install its software dependencies. They are listed in the file `environment.yaml`, which specifies a [conda](https://docs.conda.io/en/latest/) environment. We recommend using the  [mamba package mananger](https://mamba.readthedocs.io/) to install the environment instead of conda, as conda can sometimes fail when there are packages from conda-forge. To create a new conda environment called `smalldisco` using mamba, run:

```console
$ mamba env create -f environment.yaml
```

You can then activate the environment with `conda activate smalldisco`, run the program, and deactivate the environment with `conda deactivate` when you're done. You can also run a single command in the environment with `conda run -n smalldisco <command>`. Check the [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) for more information on conda environments.

### Usage

Use `python smalldisco.py --help` for main usage instructions. Smalldisco has two commands, `sirna` and `tail`, whose usage can be checked with `python smalldisco.py sirna --help` and `python smalldisco.py tail --help`, respectively.

Behind the scenes, smalldisco is actually implemented as a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline and `smalldisco.py` is just a wrapper script. The pipeline code and helper scripts are in the `workflow` folder, and the program will not work if `smalldisco.py` and `workflow` are not in the same directory.

Overall, to run both modes of smalldisco for a full siRNA identification and tailing analysis - `sirna` and `tai`l - the user needs to provide their small RNA-seq alignment BAM files of interest, a GTF/GFF reference annotation file, and a FASTA genome reference file. Below, we outline the details of these modes.

<img src="misc/Workflow.png">

#### `sirna` mode

`sirna` identifies genomic regions to which user-inputted reads map antisense. To this end, sirna requires:
* GTF/GFF reference annotation file
  * Please be sure that your GTF/GFF file follows the standard tab-separated, nine-column format (https://useast.ensembl.org/info/website/upload/gff.html). 
* Your small RNA-seq alignment alignment BAM files of interest in one folder
  * You can use any workflow or mapping program to create your sRNA-seq BAM files for analysis. However, be sure that your small RNA reads are trimmed of 3’ adapters (see vignette for example).

The required and optional arguements can be viewed with `python smalldisco.py sirna --help`

```
Usage: smalldisco.py sirna [OPTIONS] BAMFOLDER

  Find siRNA regions from antisense reads.

  This command generates putative siRNA regions based on .bam-formatted read
  alignment files in the folder BAMFOLDER.

Options:
  -o, --out PATH                  Name of output file in BED format.
                                  [default: sirna.bed]
  -a, --annotation PATH           Genome annotation file in either GTF or GFF
                                  format.
  -k, --annotation_kind [GTF|GFF]
                                  Format of genome annotation file.  [default:
                                  GTF]
  -f, --feature TEXT              Feature type in the annotation file assumed
                                  to contain siRNA regions.  [default: CDS]
  -r X                            Minimum amount of overlapping reads to
                                  create a putative siRNA region.  [default:
                                  10; x>=1]
  -s X                            Minimum size, in base pairs, of a putative
                                  siRNA.  [default: 10; x>=1]
  --help                          Show this message and exit.
  ```

The user must define the feature type (`-f`) from which they wish to map antisense reads to in their GTF/GFF file. For example, if one is interested in canonical siRNAs, the user could define their feature type as “CDS”. 

#### `tail` mode

`tail` uses Tailor to identify non-templated nucleotides (i.e., tails) on the 3’ end of small RNA reads. This mode requires: 
* FASTA genome reference file (-g)
* Your small RNA-seq alignment alignment BAM files of interest in one folder
* List of small RNA regions in BED format
  * The user can use the output BED file of putative siRNAs from sirna or a predefined BED file of another small RNA type.

For users interested in non-siRNA small RNAs (such as miRNAs or piRNAs): one way to create a bed file for a specific small RNA type is to obtain a GTF/GFF file for only that small RNA type or to filter a genomic GTF/GFF file for your small RNA type of interest. Then, a tool such as `gtf2bed` (part of BEDOPs https://bedops.readthedocs.io/en/latest/index.html) can be used to convert the GTF/GFF to a BED file.


#### Tailor integration

When running the `tail` command, we assume by default that the path to the Tailor executable is `Tailor/bin/tailor_v1.1_linux_static`. This will work if you are on Linux, are in the smalldisco repository, (e.g. `cd smalldisco`), and the Tailor submodule has been initialized as described in the installation section. If you are running smalldisco from a different folder, or have a custom Tailor installation, you must specify a path to a valid Tailor executable, for instance:

```console
$ python smalldisco.py tail --tailor_command /usr/local/bin/tailor
```

### Example

This repository comes with example data in the `example` folder. To test that smalldisco is installed correctly, you can run a toy analysis with:

```console
$ bash example-run.sh
```
This should create two outputs: `example/sirna.bed` with locations of putative siRNA regions and `example/tails.tsv` with quantified read tails.
