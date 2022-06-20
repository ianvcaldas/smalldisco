import pandas as pd

df = (
    pd.concat([pd.read_table(f) for f in snakemake.input])
    .sort_values(by=["region", "sample", "read_tail"])
)
df.to_csv(snakemake.output[0], sep='\t', index=False)
