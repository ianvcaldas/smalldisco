#!/usr/bin/env python
# coding: utf-8

from collections import Counter
import pandas as pd

raw = pd.read_table(
    snakemake.input[0],
    header=None,
    names=[
        "chrom",
        "start",
        "end",
        "num_reads",
        "feat_chrom",
        "feat_start",
        "feat_end",
        "feat_name",
        "score",
        "strand",
        "num_exons_in_feat",
        "bp_overlap",
    ],
)

# If an overlapping group of reads (with coordinates 'start' and 'end') is so
# large it covers more than one exon (with coordinates 'feat_start' and
# 'feat_end'), it will show up as two rows in the raw file with the same start
# and end coordinates. Since both of those cases represent one single putative
# siRNA region, we remove the duplicated rows by keeping only the first one,
# which will give the putative region its name.
# This happens to only a tiny handful of regions.
dedupped = raw.drop_duplicates(["chrom", "start", "end"], keep="first")

# From previous analyses, we know that there might be two rows with different
# exon coordinates (feat_start and feat_end will be different) but the same
# exon name (feat_name will be the same). This is okay, since in the next step
# those different coordinates will be assigned different siRNA IDs. So no
# worries about that situation.

base_ids = [s.split("+")[0] for s in dedupped.feat_name]
c = Counter()
sirna_ids = []
for bid in base_ids:
    c[bid] += 1
    sirna_ids.append(f"sirna_{bid}_{str(c[bid]).zfill(2)}")

bed_cols = ["chrom", "start", "end", "sirna_name", "score", "strand", "num_reads"]

df = dedupped.assign(
    sirna_name=sirna_ids,
    frame=".",
)[bed_cols]

# Add comment line to the top to indicate column names
with open(snakemake.output[0], "w") as f:
    f.write("# COLS:" + " ".join(bed_cols) + "\n")
df.to_csv(snakemake.output[0], mode="a", sep="\t", index=False, header=False)
