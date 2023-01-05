from pathlib import Path

OUTDIR = Path(".smalldisco")
BAMDIR = Path(config["bamfolder"])
FASTQDIR = OUTDIR/"tail"/"fastq"
TAILORINDEXPREFIX = OUTDIR/"tail"/"genomeindex"
TMPDIR = str(OUTDIR/"tmp")

def sample_names():
    return [f.stem for f in Path(config["bamfolder"]).iterdir() if f.suffix == ".bam"]

rule get_tailing:
  input: expand(str(OUTDIR/"tail/tail-distributions/{sample}.tsv"), sample=sample_names())
  output: config["out"]
  script: "scripts/merge-tables-of-tails.py"


rule get_sample_tails:
  input: OUTDIR/"tail/uniquely-mapped-tagged-reads/{sample}.sam"
  output: OUTDIR/"tail/tail-distributions/{sample}.tsv"
  script: "scripts/get-table-of-tails.py"


rule get_uniquely_mapped_reads:
  input: OUTDIR/"tail/tagged/{sample}_tagged.bam"
  output: OUTDIR/"tail/uniquely-mapped-tagged-reads/{sample}.sam"
  shell:
    "samtools view {input} | awk '$12 == \"NH:i:1\"' | grep \"YB:Z\" > {output}"

def get_map_antisense(wildcards):
  if config["tails_antisense"]:
    return "-S"
  else:
    return ""

rule tag_reads_with_mapped_loci:
  """Produce a BAM alignment where each read is tagged with the region it maps to."""
  input:
    database = config["bedfile"],
    sample = OUTDIR/"tail/tailor/{sample}.bam"
  output:
    OUTDIR/"tail/tagged/{sample}_tagged.bam"
  params:
    map_antisense = get_map_antisense
  shell:
    "bedtools tag -i {input.sample} {params.map_antisense} -files {input.database} -names > {output} ; "


rule tailor_one_sample:
  input:
    fastq = FASTQDIR/"{sample}.fastq",
    index_ran = TAILORINDEXPREFIX.parent/"index.log"
  output:
    sam = OUTDIR/"tail/tailor/{sample}.sam",
    bam = OUTDIR/"tail/tailor/{sample}.bam"
  threads: 8
  log: OUTDIR/"tail/tailor/{sample}.log"
  params:
    index = TAILORINDEXPREFIX,
    tailor_command = config["tailor_command"],
    tailor_min = int(config["tailor_min_prefix"])
  shell:
    "{params.tailor_command} map -l {params.tailor_min} -n {threads} -i {input.fastq} -p {params.index} -o {output.sam} &> {log} ; "
    "samtools sort -@ {threads} -o {output.bam} {output.sam} ;"


rule fastq_from_bam:
  input: Path(config["bamfolder"])/"{sample}.bam"
  output: FASTQDIR/"{sample}.fastq"
  shell:
    "samtools fastq {input} > {output}"


rule tailor_index:
  input: config["genome"]
  output: TAILORINDEXPREFIX.parent/"index.log"
  threads: 4
  params:
    index = TAILORINDEXPREFIX,
    tailor_command = config["tailor_command"]
  shell:
    "{params.tailor_command} build -i {input} -p {params.index} &> {output} ; "
