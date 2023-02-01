import subprocess
import shutil
import json
from pathlib import Path

import click

TEMPDIR = Path(".smalldisco")


@click.group()
@click.option(
    "-c",
    "--cores",
    default=1,
    show_default=True,
    help="Number of parallel cores to use.",
)
@click.option(
    "--rmtemp/--no-rmtemp",
    default=True,
    show_default=True,
    help="If enabled, delete the folder .smalldisco with intermediate files after running.",
)
def smalldisco(cores, rmtemp):
    pass


@smalldisco.command()
@click.argument("bamfolder", type=click.Path(exists=True))
@click.option(
    "-o",
    "--out",
    help="Name of output file in BED format.",
    type=click.Path(),
    default="sirna.bed",
    show_default=True,
)
@click.option(
    "-a",
    "--annotation",
    help="Genome annotation file in either GTF or GFF format.",
    type=click.Path(exists=True),
    prompt="Please provide a path to an annotation file in GTF or GFF format",
)
@click.option(
    "-k",
    "--annotation_kind",
    type=click.Choice(("GTF", "GFF"), case_sensitive=False),
    default="GTF",
    show_default=True,
    help="Format of genome annotation file.",
)
@click.option(
    "-f",
    "--feature",
    default="CDS",
    show_default=True,
    help="Feature type in the annotation file assumed to contain siRNA regions.",
)
@click.option(
    "-r",
    "sirna_min_reads",
    metavar="X",
    default=10,
    type=click.IntRange(min=1),
    help="Minimum amount of overlapping reads to create a putative siRNA region.",
    show_default=True,
)
@click.option(
    "-s",
    "sirna_min_size",
    metavar="X",
    default=10,
    type=click.IntRange(min=1),
    help="Minimum size, in base pairs, of a putative siRNA.",
    show_default=True,
)
@click.pass_context
def sirna(
    context,
    bamfolder,
    out,
    annotation,
    annotation_kind,
    feature,
    sirna_min_reads,
    sirna_min_size,
):
    """Find siRNA regions from antisense reads.

    This command generates putative siRNA regions based on .bam-formatted read alignment files in the folder BAMFOLDER."""
    params = context.parent.params.copy()
    params.update(context.params)
    run_smalldisco_pipeline(params, "sirna.smk", "siRNA discovery")


@smalldisco.command()
@click.argument("bedfile", type=click.Path(exists=True))
@click.argument("bamfolder", type=click.Path(exists=True))
@click.option(
    "-o",
    "--out",
    help="Name of output file, in TSV (tab-separated values) format.",
    type=click.Path(),
    default="tails.tsv",
    show_default=True,
)
@click.option(
    "-g",
    "--genome",
    help="Reference genome in FASTA format.",
    type=click.Path(exists=True),
    prompt="Please provide a reference genome in FASTA format",
)
@click.option(
    "--tails-antisense/--tails-all",
    default=True,
    help="Whether to quantify tails for antisense reads only or for all reads.",
    show_default=True,
)
@click.option(
    "--tailor_command",
    type=click.Path(exists=True),
    default="Tailor/bin/tailor_v1.1_linux_static",
    help="Path to Tailor binary executable.",
    show_default=True,
)
@click.option(
    "--tailor-min-prefix",
    metavar="X",
    type=click.IntRange(min=1),
    default=18,
    help=(
        "Minimum number of base pairs matching exactly in a read alignment before a tail can start. "
        "Equivalent to Tailor's '-l' command-line parameter."
    ),
    show_default=True,
)
@click.pass_context
def tail(
    context,
    bedfile,
    bamfolder,
    out,
    genome,
    tails_antisense,
    tailor_command,
    tailor_min_prefix,
):
    """Quantify tails of reads aligning to specified genome regions.

    This command quantifies tails from read alignments in .bam format found in BAMFOLDER. Only reads that overlap with certain genome regions, specified in .bed format in BEDFILE, are considered."""
    params = context.parent.params.copy()
    params.update(context.params)
    run_smalldisco_pipeline(params, "tail.smk", "read tail quantification")


def run_smalldisco_pipeline(params, snakefile, description):
    command = [
        "snakemake",
        "-c" + str(params["cores"]),
        "--configfile",
        str(TEMPDIR / "config.json"),
        "--snakefile",
        str(Path(__file__).parent / "workflow" / snakefile),
    ]
    click.echo(
        f"Running smalldisco {description} pipeline with command:\n" + " ".join(command)
    )
    Path.mkdir(TEMPDIR, exist_ok=True)
    create_configfile(params)
    subprocess.run(command)
    if params["rmtemp"]:
        shutil.rmtree(TEMPDIR)


def create_configfile(params):
    with open(TEMPDIR / "config.json", "w") as f:
        json.dump(params, f)


if __name__ == "__main__":
    smalldisco()
