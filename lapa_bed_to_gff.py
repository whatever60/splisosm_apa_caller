import argparse

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input bed file")
parser.add_argument("-o", "--output", help="output gff file")
args = parser.parse_args()

df = pd.read_table(
    args.input,
    header=None,
    names=[
        "chr",
        "start",
        "end",
        "polya_site",
        "count",
        "strand",
        "feature",
        "gene_id",
        "tpm",
        "gene_count",
        "usage",
        "frac_a",
        "signal",
        "annotated_site",
    ],
)

df["score"] = df["usage"]
df["frame"] = "."
attributes = [
    "polya_site",
    "count",
    "gene_id",
    "tpm",
    "gene_count",
    "usage",
    "frac_a",
    "signal",
    "annotated_site",
]
df = df.assign(
    attribute=df.apply(
        lambda row: ";".join(f"{attr}={row[attr]}" for attr in attributes), axis=1
    )
)
df["attribute"] += df.apply(
    lambda row: ";gene_id_polya_site=" + row["gene_id"] + "_" + str(row["polya_site"]),
    axis=1,
)
df["source"] = "lapa_bed_to_gff"
# GFF is 1-based, BED is 0-based
df["start"] = df["start"] + 1

df = df[
    [
        "chr",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]
].sort_values(["chr", "start"])

df.to_csv(args.output, sep="\t", index=False, header=False)
