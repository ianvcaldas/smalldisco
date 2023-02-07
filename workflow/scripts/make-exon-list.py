#!/usr/bin/env python
# coding: utf-8

import re
import csv
import pandas as pd


# Extract info from attribute field and turn into columns
if snakemake.config["annotation_kind"].lower() == "gtf":
    pattern = re.compile(r'(\w+) "(.*?)";')
    id_name = "gene_id"
elif snakemake.config["annotation_kind"].lower() == "gff":
    pattern = re.compile(r"(\w+)=(.*?);")
    id_name = "ID"


def extract(ix, s, patt):
    return [(ix, *attributes) for attributes in patt.findall(s)]


annot = pd.read_table(snakemake.input[0], comment="#", header=None)
annot.columns = [
    "chrom",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]

parsed_attributes = []
for ix, value in annot.attribute.iteritems():
    parsed_attributes.extend(extract(ix, value, pattern))

attributes_df = pd.DataFrame.from_records(
    parsed_attributes, columns=["index", "attribute", "value"]
)
attributes_df = attributes_df.loc[attributes_df.attribute == id_name]
attributes_df = attributes_df.drop(attributes_df.loc[attributes_df.value == ""].index)
attributes_df = attributes_df.pivot_table(
    index="index", columns="attribute", values="value", aggfunc=lambda x: ";".join(x)
)

data = annot.loc[annot.feature == snakemake.config["feature"]].join(attributes_df)

# Backup IDs in case getting them from the annotation fails
backup_ids = pd.Series(
    [f"{i}-{j}-{k}" for i, j, k in zip(data.chrom, data.start, data.end)],
    index=data.index,
    dtype='object'
)
data[id_name].fillna(backup_ids, inplace=True)

# Save in BED format
# We don't use quotes in the output to stay compatible with the GTF and BED formats.
bed_cols = ["chrom", "start", "end", id_name, "score", "strand", "frame"]
data[bed_cols].to_csv(
    snakemake.output[0], sep="\t", header=False, index=False, quoting=csv.QUOTE_NONE
)
