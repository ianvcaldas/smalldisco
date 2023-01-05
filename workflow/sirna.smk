from pathlib import Path

OUTDIR = Path(".smalldisco")
BAMDIR = Path(config["bamfolder"])
TMPDIR = str(OUTDIR/"tmp")

def sample_names():
    return [f.stem for f in Path(config["bamfolder"]).iterdir() if f.suffix == ".bam"]

rule locate_sirna_regions:
    """Create a database of sirna regions with unique names."""
    input: OUTDIR/"sirna/sirna-regions-annotated.bed"
    output: config["out"]
    script: "scripts/clean-sirna-regions.py"


rule intersect_sirna_and_exons:
    input:
        sirna_regions = OUTDIR/"sirna/putative-sirna-regions.bed",
        exon_list = OUTDIR/"sirna/exon-list.bed"
    output:
        OUTDIR/"sirna/sirna-regions-annotated.bed"
    shell:
        "bedtools intersect -a {input.sirna_regions} -b {input.exon_list} -wo > {output} ; "


rule find_putative_sirna_regions:
    """Produce putative sirna regions by merging overlapping reads that map antisense to exons."""
    input: OUTDIR/"sirna/all-samples_antisense-reads.bam"
    output: OUTDIR/"sirna/putative-sirna-regions.bed"
    params:
        filter_threshold = config["sirna_min_reads"]
    resources:
        tmpdir = TMPDIR
    shell:
        # Merge overlapping reads into putative regions.
        "bedtools merge -c 1 -o count -i {input} > ${{TMPDIR}}/merged.bed ; "
        # Filter merged BED file to regions with a minimum number of reads.
        "awk '$4 >= {params.filter_threshold}' ${{TMPDIR}}/merged.bed > {output} ; "


rule merge_antisense_reads_from_samples:
    """Merge all BAMs with reads that map antisense to exons into one multi-sample BAM."""
    input: expand(str(OUTDIR/"sirna/antisense-reads-samples/{sample}_antisense-reads.bam"), sample=sample_names())
    output: temp(OUTDIR/"sirna/all-samples_antisense-reads.bam")
    shell:
        "samtools merge {output} {input} ; "


rule map_antisense_reads:
    """Obtain reads from one sample mapped antisense to exons."""
    input:
        exons = OUTDIR/"sirna/exon-list.bed",
        alignment = BAMDIR/"{sample}.bam"
    resources:
        tmpdir = TMPDIR
    output: temp(OUTDIR/"sirna/antisense-reads-samples/{sample}_antisense-reads.bam")
    shell:
        "samtools sort {input.alignment} > ${{TMPDIR}}/sorted.bam ; "
        "bedtools intersect -S -a ${{TMPDIR}}/sorted.bam -b {input.exons} > {output} ; "


rule remove_exon_overlaps:
    """Produce exon list where regions that overlap on the same strand are
    merged together, and regions of overlap across different strands are removed."""
    input: OUTDIR/"sirna/exon-list-with-overlaps.bed"
    output: OUTDIR/"sirna/exon-list.bed"
    params:
        min_region_size = config["sirna_min_size"]
    resources:
        tmpdir = TMPDIR
    shell:
        # Sort input by chromosome and start position
        "sort -k1,1 -k2,2n {input} > ${{TMPDIR}}/sorted.bed ; "
        # Merge all overlapping features on the same strand first.
        "bedtools merge -i ${{TMPDIR}}/sorted.bed -s -c 4,5,6,1 -o distinct,distinct,distinct,count > ${{TMPDIR}}/merged.bed ; "
        # After this step, every possible overlap can only be across different strands.
        # So now we remove the regions of overlap across strands.
        # First, separate the merged features into the ones in the plus strand and the ones in the minus strand.
        "awk '$6 == \"-\"' ${{TMPDIR}}/merged.bed > ${{TMPDIR}}/minus.bed ; "
        "awk '$6 == \"+\"' ${{TMPDIR}}/merged.bed > ${{TMPDIR}}/plus.bed ; "
        # The intersection of those are the areas where there is some overlap across strands.
        "bedtools intersect -a ${{TMPDIR}}/minus.bed -b ${{TMPDIR}}/plus.bed > ${{TMPDIR}}/overlaps.bed ; "
        # Finally, take those areas out.
        "bedtools subtract -a ${{TMPDIR}}/merged.bed -b ${{TMPDIR}}/overlaps.bed > ${{TMPDIR}}/subtracted.bed ; "
        # Filter database to remove tiny regions, then replace "," characters in the
        # database by "+" characters such that read tagging preserves multi-region
        # entries in the database.
        "cat ${{TMPDIR}}/subtracted.bed | awk '($3 - $2) > {params.min_region_size}' |  tr ',' '+' > {output} ; "


rule make_exon_list:
    """Make a list of exon locations in the genome that could have a siRNA region."""
    input: config["annotation"]
    output: OUTDIR/"sirna/exon-list-with-overlaps.bed"
    script: "scripts/make-exon-list.py"
