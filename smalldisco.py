import subprocess
import shutil
import json
from pathlib import Path

import click

TEMPDIR = Path(".smalldisco")
SIRNAFILE = "sirna.bed"
TAILSFILE = "tails.tsv"

@click.group()
@click.option(
    "-c",
    "--cores",
    default=1,
    show_default=True,
    help="Number of parallel cores to use.",
)
@click.option(
    "--snakefile",
    type=click.Path(exists=True),
    default="./workflow/Snakefile",
    show_default=True,
    help="Snakefile with the smalldisco pipeline.",
)
@click.option(
    "--rmtemp/--no-rmtemp",
    default=True,
    show_default=True,
    help="After running, delete the folder .smalldisco with intermediate files.",
)
def smalldisco():
    pass

@smalldisco.command()
@click.argument("bamfolder", type=click.Path(exists=True))
@click.option(
    "-a",
    "--annotation",
    help="Genome annotation file in GTF or GFF format.",
    type=click.Path(exists=True),
    prompt="Genome annotation file in GTF or GFF format",
)
@click.option(
    "-k",
    "--annotation_kind",
    type=click.Choice(("GTF", "GFF"), case_sensitive=False),
    default="GTF",
    show_default=True,
    help="Type of genome annotation.",
)
@click.option(
    "-f",
    "--feature",
    default="CDS",
    show_default=True,
    help="Annotation feature to search for siRNAs in.",
)
@click.option(
    "-r",
    "sirna_min_reads",
    default=10,
    type=click.IntRange(min=1),
    help="Minimum amount of overlapping reads to create a putative siRNA region.",
    show_default=True,
)
@click.option(
    "-s",
    "sirna_min_size",
    default=10,
    type=click.IntRange(min=1),
    help="Minimum size in base pairs of putative siRNA regions.",
    show_default=True,
)
def sirna(
    bamfolder,
    genome,
    annotation,
    annotation_kind,
    feature,
    sirna_min_reads,
    sirna_min_size,
    tails_antisense,
    tailor_command,
    tailor_min_prefix,
    cores,
    snakefile,
    rmtemp,
):
    """Finds siRNA based on the alignments in BAMFOLDER."""
    command = [
        "snakemake",
        "-c" + str(cores),
        "--configfile",
        str(TEMPDIR / "config.json"),
        "--snakefile",
        snakefile,
        "locate_sirna_regions",
    ]
    click.echo("Running smalldisco pipeline with command:\n" + " ".join(command))
    Path.mkdir(TEMPDIR, exist_ok=True)
    create_configfile(
        bamfolder,
        genome,
        annotation,
        annotation_kind,
        feature,
        sirna_min_reads,
        sirna_min_size,
        tails_antisense,
        tailor_command,
        tailor_min_prefix,
    )
    subprocess.run(command)
    if rmtemp:
        shutil.rmtree(TEMPDIR)

@click.argument("bamfolder", type=click.Path(exists=True))
@click.option(
    "-g",
    "--genome",
    help="Reference genome in FASTA format.",
    type=click.Path(exists=True),
    prompt="Reference genome in FASTA format",
)
@click.option(
    "--tails-antisense/--tails-all",
    default=True,
    help="Quantify tails for antisense reads only, or for all tails.",
    show_default=True,
)
@click.option(
    "--tailor_command",
    type=click.Path(exists=True),
    default="Tailor/bin/tailor_v1.1_linux_static",
    help="Path to Tailor executable.",
    show_default=True,
)
@click.option(
    "--tailor-min-prefix",
    type=click.IntRange(min=1),
    default=18,
    help="Minimum number of bp that must align to genome before a tail is detected (Tailor parameter).",
    show_default=True,
)
@click.option(
    "-c",
    "--cores",
    default=1,
    show_default=True,
    help="Number of parallel cores to use.",
)
@click.option(
    "--snakefile",
    type=click.Path(exists=True),
    default="./workflow/Snakefile",
    show_default=True,
    help="Snakefile with the smalldisco pipeline.",
)
@click.option(
    "--rmtemp/--no-rmtemp",
    default=True,
    show_default=True,
    help="After running, delete the folder .smalldisco with intermediate files.",
)
@smalldisco.command()
def tail():
    """Get tails of reads intersecting with specified genome regions"""
    click.echo("Tails!")

def create_configfile(
    bamfolder,
    genome,
    annotation,
    annotation_kind,
    feature,
    sirna_min_reads,
    sirna_min_size,
    tails_antisense,
    tailor_command,
    tailor_min_prefix,
):
    data = {
        "bam_dir": bamfolder,
        "genome": genome,
        "annotation": annotation,
        "annotation_kind": annotation_kind,
        "sirna_source": feature,
        "min_reads_per_sirna": sirna_min_reads,
        "min_sirna_size": sirna_min_size,
        "sirna_output": SIRNAFILE,
        "small_rna_reference": SIRNAFILE,
        "tails_antisense_only": tails_antisense,
        "tailor_command": tailor_command,
        "tailor_min_prefix_match": tailor_min_prefix,
        "tails_output": TAILSFILE,
    }
    print(data)
    with open(TEMPDIR / "config.json", "w") as f:
        json.dump(data, f)


if __name__ == "__main__":
    smalldisco()
