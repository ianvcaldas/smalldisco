from collections import Counter
import pandas as pd

tags = []

with open(snakemake.input[0], "r") as f:
    for line in f:
        # Skip SAM header
        if line.startswith("@"):
            continue
        # Get read info
        elements = line.strip().split()
        # 1st element is the read name
        read_name = elements[0]
        # 11th element onward in the SAM are the custom tags
        read_tags = elements[11:]
        read_tags = {s.split(":")[0]: s.split(":")[2] for s in read_tags}
        # We don't care about multiple-mapped reads or reads that didn't map to a RNA locus
        if ("YB" not in read_tags) or (int(read_tags["NH"]) != 1):
            print("IMPOSSIBLE!", args.samplename, read_name, read_tags)
            continue
        # The following reads are uniquely mapped to a locus of interest
        else:
            # We remove the NH tag counting how many loci the reads mapped to
            # because we want all read_tag dictionaries to look the same from now
            # on, no matter what NH was.
            read_tags.pop("NH")
            if "TL" not in read_tags:
                read_tags["TL"] = "untailed"
            # If a read maps to multiple overlapping loci, take only the first one.
            if "," in read_tags["YB"]:
                these_tags = [
                    {"YB": mapped_locus, "TL": read_tags["TL"]}
                    for mapped_locus in read_tags["YB"].split(",")
                ]
                read_tags = these_tags[0]
            # Determine read type and log it if it's a read type of interest
            tags.append(read_tags)


tag_counter = Counter(tuple(sorted(t.items())) for t in tags)

records = []
for tag in tag_counter:
    records.append((tag[0][1], tag[1][1], tag_counter[tag]))

counts = pd.DataFrame.from_records(
    records, columns=["read_tail", "region", "num_reads"]
)

all_tails = (
    pd.DataFrame(counts.loc[counts.read_tail != "untailed"])
    .groupby("region")
    .sum()
    .num_reads.reset_index()
    .assign(read_tail="all_tails")
)

counts = (
    pd.concat([counts, all_tails])
    .sort_values(by=["region", "read_tail"])[["region", "read_tail", "num_reads"]]
    .reset_index(drop=True)
    .assign(sample=snakemake.wildcards["sample"])
)

counts.to_csv(snakemake.output[0], sep="\t", index=False)
