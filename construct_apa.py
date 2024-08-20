"""
Script to process APA data for a single sample. The script maps peaks to genes, filters ambiguous peaks, 
and processes the data to prepare it for downstream analysis, such as in spatial transcriptomics.
"""

import os
import argparse

import pandas as pd
import polars as pl
import numpy as np
import scanpy as sc
import scipy.sparse as ss
from tqdm.auto import tqdm

from utils_annot import get_tx2gene
from utils import load_annotation


def construct_apa_dataset(
    df_count: pl.DataFrame,
    gene2peaks: dict[str, list[str]],
    cells: list[str] | np.ndarray | None,
) -> dict[int, tuple[np.ndarray, pd.DataFrame]]:
    # for each num_peaks, for example num_peaks == 2, construct a np matrix of shape
    # [num_genes_with_num_peaks, num_peaks, num_cells], save to npz format
    if cells is None:
        cells = df_count["cell"].unique()
    gene_num_peaks = pd.Series({gene: len(peaks) for gene, peaks in gene2peaks.items()})
    ret = {}
    for num_peaks in gene_num_peaks.unique():
        if num_peaks == 1:
            continue
        num_cells = len(cells)
        genes_with_num_peaks = gene_num_peaks[
            gene_num_peaks == num_peaks
        ].index.to_numpy()
        count = df_count.filter(pl.col("gene").is_in(genes_with_num_peaks))
        count_arrs = []
        for gene in tqdm(genes_with_num_peaks):
            count_g = count.filter(pl.col("gene") == gene)
            peaks = gene2peaks[gene]
            apaxcell = (
                pl.Series(peaks)
                .to_frame("apa_site_name")
                .join(pl.Series(cells).to_frame("cell"), how="cross")
            )
            count_gene_arr = (
                apaxcell.join(count_g, on=["cell", "apa_site_name"], how="left")[
                    "counts"
                ]
                .to_numpy()
                .reshape(num_peaks, num_cells)
            )
            count_arrs.append(count_gene_arr)
        count_arr = np.nan_to_num(np.stack(count_arrs), 0).astype(int)
        ret[num_peaks] = (
            count_arr,
            pd.DataFrame.from_dict(
                {g: gene2peaks[g] for g in genes_with_num_peaks},
                orient="index",
                columns=[f"apa_peak_{i}" for i in range(num_peaks)],
            ),
        )
    return ret


def read_spot_info(tissue_position_file: str) -> pd.DataFrame:
    """Read the tissue position file from space ranger output.

    The format of this file changed when space ranger was updated to 2.0, when it added
        header and changed file name from tissue_positions_list.csv to
        tissue_positions.csv. Also for visium HD, the file changed from csv format to
        parquet. The function here deals with all these cases.

    Number of rows in this file is:
        - 4,992 for 6.5mm capture area
        - 14,336 for 11mm capture area
        - 11,222,500 for 2um Visium HD
    """
    if tissue_position_file.endswith(".parquet"):
        spot_info = pd.read_parquet(tissue_position_file)
    else:
        has_header = False
        with open(tissue_position_file, "r") as f:
            first_line = f.readline()
            if first_line.startswith("barcode"):
                has_header = True
        if has_header:
            names = None
        else:
            names = [
                "barcode",
                "in_tissue",
                "array_row",
                "array_col",
                "pxl_row_in_fullres",
                "pxl_col_in_fullres",
            ]
        spot_info = pd.read_csv(tissue_position_file, names=names)
    assert spot_info.shape[0] in [
        4992,
        14336,
        11222500,
    ], f"Unexpected number of rows: {spot_info.shape[0]}"
    return spot_info


def main(args):
    df = load_annotation(args.annot)
    print(f"Loaded annotation file: {df.feature.value_counts()} features identified.")

    tx2gene = get_tx2gene(df, flavor=args.annot_flavor)
    # gene2name = get_gene2name(df, flavor=args.annot_flavor)

    df_peak = pd.read_table(args.peak, index_col=0)
    print(f"Peak data loaded: {df_peak.shape[0]} rows, {df_peak.shape[1]} columns.")

    # first extract all parent transcript IDs of 3' UTRs in the peak df
    # then map these transcript to gene
    if "gff" in args.annot:
        parent_transcript_ids = [
            j.split("=")[1]
            for attr in df_peak.utr_info.str.split(";")
            for j in attr
            if j.startswith("Parent=")
        ]
    elif "gtf" in args.annot:
        print(
            "WARNING: GTF format not yet well supported, unexpected results may occur."
        )
        parent_transcript_ids = [
            j.split()[1].strip('"')
            for attrs in df_peak.utr_info.str.split(";;")
            for attr in attrs
            for j in attr.split("; ")
            if j.startswith("transcript_id ")
        ]
    else:
        raise NotImplementedError
    parent_gene_ids = [tx2gene[i] for i in parent_transcript_ids if i in tx2gene]

    print(f"{len(parent_transcript_ids)} transcripts found in the peak data.")
    print(
        f"{sum([i in tx2gene for i in parent_transcript_ids])} transcripts mapped to genes."
    )

    print(f"Parent gene ID counts:\n{pd.Series(parent_gene_ids).value_counts()}")

    df_peak["parent_transcripts"] = df_peak.utr_info.apply(
        lambda x: ";".join(
            [
                attr.split("=")[-1].split(":")[-1]
                for attr in x.split(";")
                if attr.startswith("Parent=")
            ]
        )
    )
    df_peak["parent_genes"] = df_peak.parent_transcripts.apply(
        lambda x: ";".join(list(set(tx2gene[i] for i in x.split(";"))))
    )

    # extract peaks that unanimously map to a single gene (persumably that's all peaks)
    peak2genes = df_peak.parent_genes.to_dict()
    peaks_unambiguous = df_peak.query("~parent_genes.str.contains(';')").index.tolist()
    print(
        f"Peak-to-gene mapping statistics:\n{(df_peak.parent_genes.str.count(';') + 1).value_counts()}"
    )

    gene2peaks: dict[str, list[str]] = {}
    for peak_name, row in df_peak.query("~parent_genes.str.contains(';')").iterrows():
        gene = row.parent_genes
        if gene in gene2peaks:
            gene2peaks[gene].append(peak_name)
        else:
            gene2peaks[gene] = [peak_name]
    print(
        f"Number of peaks per gene:\n{pd.Series(gene2peaks).apply(len).value_counts()}"
    )

    df_count = (
        pl.read_csv(
            args.count,
            separator="\t",
            new_columns=["apa_site_name", "cell", "counts"],
        )
        .filter(pl.col("apa_site_name").is_in(peaks_unambiguous))  # remove bad peaks
        .with_columns(gene=pl.col("apa_site_name").replace(peak2genes))
    )

    spot_info = read_spot_info(args.spot_info)
    cells = spot_info.barcode.map(lambda x:x[:-2])
    print(f"Number of cells in APA count table: {len(cells)}")
    apa_data = construct_apa_dataset(df_count, gene2peaks, cells)

    # save results
    sample = args.sample
    if sample is not None:
        # npz_path_template = os.path.join(
        #     args.output_dir,
        #     f"{args.sample}_mtx_apa_{{}}.npz" if args.sample else "mtx_apa_{}.npz",
        # )
        # genes_apa_path_template = os.path.join(
        #     args.output_dir,
        #     f"{args.sample}_genes_apa_{{}}.tsv" if args.sample else "genes_apa_{}.tsv",
        # )
        # spot_info_save_path = os.path.join(
        #     args.output_dir, f"{args.sample}_cells.csv" if args.sample else "cells.csv"
        # )
        # gene_ids_save_path = os.path.join(
        #     args.output_dir,
        #     f"{args.sample}_gene_ids.csv" if args.sample else "gene_ids.csv",
        # )
        # gene_count_save_path = os.path.join(
        #     args.output_dir,
        #     f"{args.sample}_gene_count.csv" if args.sample else "gene_count.csv",
        # )
        npz_path_template = os.path.join(args.output_dir, f"{sample}_mtx_apa_{{}}.npz")
        genes_apa_path_template = os.path.join(
            args.output_dir, f"{sample}_genes_apa_{{}}.tsv"
        )
        spot_info_save_path = os.path.join(args.output_dir, f"{sample}_cells.csv")
        gene_ids_save_path = os.path.join(args.output_dir, f"{sample}_gene_ids.csv")
        gene_count_save_path = os.path.join(args.output_dir, f"{sample}_gene_count.csv")
    else:
        npz_path_template = os.path.join(args.output_dir, "mtx_apa_{}.npz")
        genes_apa_path_template = os.path.join(args.output_dir, "genes_apa_{}.tsv")
        spot_info_save_path = os.path.join(args.output_dir, "cells.csv")
        gene_ids_save_path = os.path.join(args.output_dir, "gene_ids.csv")
        gene_count_save_path = os.path.join(args.output_dir, "gene_count.csv")

    os.makedirs(args.output_dir, exist_ok=True)
    for num_peaks, (count_arr, genes_with_num_peaks) in apa_data.items():
        np.savez_compressed(npz_path_template.format(num_peaks), count_arr)
        genes_with_num_peaks.to_csv(genes_apa_path_template.format(num_peaks), sep="\t")

    adata = sc.read_10x_h5(args.adata)
    count = adata.X
    gene_ids = adata.var.gene_ids

    # Construct output file paths
    spot_info.to_csv(spot_info_save_path)
    gene_ids.to_csv(gene_ids_save_path, header=False)
    ss.save_npz(gene_count_save_path, count)

    print(f"Spot information loaded: {spot_info.shape[0]} spots.")
    print(f"adata shape: {adata.shape}")
    print(f"Spot information saved to {spot_info_save_path}")
    print(f"Gene IDs saved to {gene_ids_save_path}")
    print(f"Gene count matrix saved to {gene_count_save_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process APA data for a single sample."
    )

    parser.add_argument(
        "-a",
        "--annot",
        type=str,
        required=True,
        help="Path to the annotation file (GFF/GTF).",
    )
    parser.add_argument(
        "-f",
        "--annot-flavor",
        type=str,
        default="gencode-gff3",
        help="Flavor/type of annotation used (e.g., gencode-gff3).",
    )
    parser.add_argument(
        "-p",
        "--peak",
        type=str,
        required=True,
        help="Path to the peak file for the sample.",
    )
    parser.add_argument(
        "-c",
        "--count",
        type=str,
        required=True,
        help="Path to the UMI tools count file.",
    )
    parser.add_argument(
        "-s",
        "--spot-info",
        type=str,
        required=True,
        help="Path to the spot information file.",
    )
    parser.add_argument(
        "-d",
        "--adata",
        type=str,
        required=True,
        help="Path to the 10x Genomics filtered feature matrix (HDF5 format).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        required=True,
        help="Directory to save the processed output files.",
    )
    parser.add_argument(
        "-S",
        "--sample",
        type=str,
        default=None,
        help="Sample identifier (e.g., s1). Optional.",
    )

    args = parser.parse_args()
    main(args)
