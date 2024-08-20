from pathlib import Path
from io import StringIO
import warnings

import numpy as np
import pandas as pd
import polars as pl
import pyranges as pr


def load_region(file: str | StringIO, file_type: str | None = None) -> pd.DataFrame:
    """
    Load a BED or narrowPeak file into a DataFrame.

    Args:
    file (str | StringIO): Path to the file or StringIO object.
    file_type (str): Type of the file, either 'bed' or 'narrowPeak'.

    Returns:
    pd.DataFrame: DataFrame containing the file data.
    """
    if file_type is None:
        if isinstance(file, str):
            file_type = deduce_file_type(file)
        else:
            raise ValueError(
                "Please provide the file type, otherwise provide the file path as a string."
            )
    if file_type == "bed":
        columns = ["chrom", "start", "end", "name", "score", "strand", "peak"]
    elif file_type == "narrowPeak":
        columns = [
            "chrom",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "signalValue",
            "pValue",
            "qValue",
            "peak",
        ]
    else:
        raise ValueError(
            f"Unsupported file type {file_type}. Please use 'bed' or 'narrowPeak'."
        )

    df = pd.read_csv(file, sep="\t", header=None, names=columns)
    return df


def load_annotation(
    file: str | StringIO, annotation_type: str | None = None, feature: str | None = None
) -> pd.DataFrame:
    """
    Load a GFF3 or GTF file into a DataFrame and filter for 3' UTR features.

    Args:
    file (str | StringIO): Path to the GFF3 or GTF file or StringIO object.
    annotation_type (str): Type of the annotation file, either 'gff3', 'gtf', 'gff3.gz', or 'gtf.gz'.

    Returns:
    pd.DataFrame: DataFrame containing the 3' UTR data.
    """
    if annotation_type is None:
        if isinstance(file, str):
            annotation_type = deduce_file_type(file)
        else:
            raise ValueError(
                "Please provide the annotation file type, otherwise provide the file path as a string."
            )
    if not annotation_type in [
        "gff3",
        "gff3.gz",
        "gtf",
        "gtf.gz",
        "gff",
        "gff.gz",
    ]:
        raise ValueError(
            "Unsupported annotation file type. "
            "Please use 'gff3', 'gtf', 'gff3.gz', 'gtf.gz', 'gff', or 'gff.gz'."
        )
    annot_df = pl.read_csv(
        file,
        separator="\t",
        comment_prefix="#",
        # has_header=None,
        new_columns=[
            "chrom",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ],
    )
    if feature is not None:
        annot_df = annot_df.filter(pl.col("feature") == feature).to_pandas()
    else:
        annot_df = annot_df.to_pandas()
    return annot_df


def deduce_file_type(file: str) -> str:
    """
    Deduce the file type based on the file extension.

    Args:
    file (str): Path to the file.

    Returns:
    str: Deduced file type.
    """
    suffixes = Path(file).suffixes
    if suffixes[-1] == ".gz":
        suffix = suffixes[-2] + ".gz"
    else:
        suffix = suffixes[-1]
    return suffix.lstrip(".")


def find_overlaps(
    bed_df: pd.DataFrame,
    utr_df: pd.DataFrame,
    strand_specific: bool,
    feat: str,
    overlap_thres: float = 0.5,
) -> pd.DataFrame:
    """
    Find overlaps between BED regions and 3' UTRs, retaining those with at least 50% overlap.

    Args:
    bed_df (pd.DataFrame): DataFrame containing BED regions.
    utr_df (pd.DataFrame): DataFrame containing 3' UTRs.
    strand_specific (bool): Whether to consider strand specificity.

    Returns:
    pd.DataFrame: DataFrame with overlapping regions that meet the 50% overlap criterion.
    """
    # Convert DataFrames to PyRanges objects
    bed_df = bed_df.rename(
        columns={
            "chrom": "Chromosome",
            "start": "Start",
            "end": "End",
            "strand": "Strand",
            "name": "Name",
            "score": "Score",
        }
    )

    utr_df = utr_df.rename(
        columns={
            "chrom": "Chromosome",
            "start": "Start",
            "end": "End",
            "strand": "Strand",
            "attribute": "Attribute",
        }
    )

    if strand_specific:
        strandedness = "same"
    else:
        strandedness = False

    # Find overlaps
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bed_pr = pr.PyRanges(bed_df)
        utr_pr = pr.PyRanges(utr_df)
        overlaps = bed_pr.join(utr_pr, strandedness=strandedness, suffix=f"_{feat}")

    # Calculate overlap lengths and filter by 50% criterion
    # overlaps = overlaps.insert(
    #     (
    #         np.minimum(overlaps.End, overlaps[f"End_{feat}"])
    #         - np.maximum(overlaps.Start, overlaps[f"Start_{feat}"])
    #     ).rename("overlap_len")
    # ).insert((overlaps.End - overlaps.Start).rename("bed_len"))
    overlaps = overlaps.df
    overlaps["overlap_len"] = np.minimum(
        overlaps.End, overlaps[f"End_{feat}"]
    ) - np.maximum(overlaps.Start, overlaps[f"Start_{feat}"])
    overlaps["bed_len"] = overlaps.End - overlaps.Start
    overlaps[f"{feat}_len"] = overlaps[f"End_{feat}"] - overlaps[f"Start_{feat}"]
    overlaps = overlaps[
        ((overlaps.overlap_len / overlaps.bed_len) >= overlap_thres)
        | ((overlaps.overlap_len / overlaps[f"{feat}_len"]) >= overlap_thres)
    ]

    # Prepare the final DataFrame
    result_df = overlaps[
        [
            "Chromosome",
            "Start",
            "End",
            "Name",
            "Score",
            "Strand",
            "peak",
            f"Start_{feat}",
            f"End_{feat}",
            f"Strand_{feat}",
            "Attribute",
        ]
    ]
    result_df.columns = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "peak",
        f"{feat}_start",
        f"{feat}_end",
        f"{feat}_strand",
        f"{feat}_info",
    ]

    return result_df


def read_genomic_region(
    file: str | StringIO,
    annotation_file: str | StringIO,
    *,
    strand_specific: bool = True,
    feature: str,
) -> pd.DataFrame:
    """
    Process BED or narrowPeak and GFF3 or GTF files to find BED regions with significant overlaps with 3' UTRs.

    Args:
    file (str | StringIO): Path to the BED or narrowPeak file or StringIO object.
    annotation_file (str | StringIO): Path to the GFF3 or GTF file or StringIO object.
    file_type (str, optional): Type of the BED or narrowPeak file. If not provided, deduced from file name.
    annotation_file_type (str, optional): Type of the GFF3 or GTF file. If not provided, deduced from file name.
    strand_specific (bool): Whether to consider strand specificity. Default is True.

    Returns:
    pd.DataFrame: DataFrame with BED or narrowPeak regions having significant overlaps with 3' UTRs and relevant 3' UTR information.
    """
    bed_df = load_region(file)
    annotation_df = load_annotation(annotation_file, feature=feature)
    print(annotation_df.query("attribute.str.contains('Exosc5')"))
    result_df = find_overlaps(
        bed_df, annotation_df, strand_specific=strand_specific, feat=feature
    )
    return result_df


# Simulate dummy data for testing
def generate_dummy_bed_gff3() -> tuple[StringIO, StringIO]:
    """
    Generate dummy BED and GFF3 data files for testing.

    Returns:
    tuple[StringIO, StringIO]: StringIO objects for the dummy BED and GFF3 data files.
    """
    bed_data = """chr1\t100\t200\tregion1\t0\t+
chr1\t150\t250\tregion2\t0\t+
chr1\t300\t400\tregion3\t0\t-
chr2\t100\t200\tregion4\t0\t+
"""
    gff3_data = """chr1\tsource\tthree_prime_UTR\t180\t220\t.\t+\t.\tID=utr1
chr1\tsource\tthree_prime_UTR\t170\t210\t.\t+\t.\tID=utr2
chr1\tsource\tthree_prime_UTR\t390\t450\t.\t-\t.\tID=utr3
chr2\tsource\tthree_prime_UTR\t50\t150\t.\t+\t.\tID=utr4
"""
    bed_file = StringIO(bed_data)
    gff3_file = StringIO(gff3_data)
    return bed_file, gff3_file


# Generate dummy data for narrowPeak format
def generate_dummy_narrowPeak() -> tuple[StringIO, StringIO]:
    """
    Generate dummy narrowPeak and GFF3 data files for testing.

    Returns:
    tuple[StringIO, StringIO]: StringIO objects for the dummy narrowPeak and GFF3 data files.
    """
    narrowPeak_data = """chr1\t100\t200\tpeak1\t0\t+\t10.0\t5.0\t2.0\t150
chr1\t150\t250\tpeak2\t0\t+\t20.0\t6.0\t3.0\t200
chr1\t300\t400\tpeak3\t0\t-\t15.0\t4.5\t1.5\t350
chr2\t100\t200\tpeak4\t0\t+\t18.0\t5.5\t2.5\t125
"""
    gff3_data = """chr1\tsource\tthree_prime_UTR\t180\t220\t.\t+\t.\tID=utr1
chr1\tsource\tthree_prime_UTR\t170\t210\t.\t+\t.\tID=utr2
chr1\tsource\tthree_prime_UTR\t390\t450\t.\t-\t.\tID=utr3
chr2\tsource\tthree_prime_UTR\t50\t150\t.\t+\t.\tID=utr4
"""
    narrowPeak_file = StringIO(narrowPeak_data)
    gff3_file = StringIO(gff3_data)
    return narrowPeak_file, gff3_file


# Generate dummy data for GTF format
def generate_dummy_gtf() -> tuple[StringIO, StringIO]:
    """
    Generate dummy BED and GTF data files for testing.

    Returns:
    tuple[StringIO, StringIO]: StringIO objects for the dummy BED and GTF data files.
    """
    bed_data = """chr1\t100\t200\tregion1\t0\t+
chr1\t150\t250\tregion2\t0\t+
chr1\t300\t400\tregion3\t0\t-
chr2\t100\t200\tregion4\t0\t+
"""
    gtf_data = """chr1\tsource\tthree_prime_UTR\t180\t220\t.\t+\t.\tgene_id "utr1"; transcript_id "utr1";
chr1\tsource\tthree_prime_UTR\t170\t210\t.\t+\t.\tgene_id "utr2"; transcript_id "utr2";
chr1\tsource\tthree_prime_UTR\t390\t450\t.\t-\t.\tgene_id "utr3"; transcript_id "utr3";
chr2\tsource\tthree_prime_UTR\t50\t150\t.\t+\t.\tgene_id "utr4"; transcript_id "utr4";
"""
    bed_file = StringIO(bed_data)
    gtf_file = StringIO(gtf_data)
    return bed_file, gtf_file


if __name__ == "__main__":
    import argparse

    # 4 arguments for input narrowpeak file, input gff3 file, output tsv file, output narrowpeak file
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input narrowpeak file")
    parser.add_argument("-a", "--annotation", help="Input gff3 file")
    parser.add_argument("-o", "--output", help="Output tsv file")
    parser.add_argument("-n", "--narrowpeak", help="Output narrowpeak file")
    args = parser.parse_args()

    np_path = args.input
    annot_path = args.annotation
    output_tsv = args.output
    output_np = args.narrowpeak

    feat = "three_prime_UTR" if "gff" in annot_path else "UTR"
    df_peak = read_genomic_region(
        np_path,
        annotation_file=annot_path,
        strand_specific=True,
        feature=feat,
    )
    # aggregate utr_info into list by grouping entries with same chrom, start, end,
    # strand and peak. Drop utr_start, utr_end, utr_strand
    df_peak_g = []
    for name, group in df_peak.groupby(["chrom", "start", "end", "strand", "peak"]):
        df_peak_g.append(
            {
                "chrom": name[0],
                "start": name[1],
                "end": name[2],
                "strand": name[3],
                "peak": name[4],
                "utr_info": ";".join(group[f"{feat}_info"]),
            }
        )
    df_peak_g = pd.DataFrame(df_peak_g)
    df_peak_g
    df_peak_g.index = (
        df_peak_g.chrom
        + "_"
        + df_peak_g.start.astype(str)
        + "_"
        + df_peak_g.end.astype(str)
        + "_"
        + df_peak_g.peak.astype(str)
    )
    print(
        f"Found {df_peak_g.shape[0]} RNA-seq coverage peaks across "
        f"{df_peak[['chrom', 'start', 'end', 'strand']].drop_duplicates().shape[0]} 3' UTR regions"
    )
    # save as tsv
    df_peak_g.to_csv(output_tsv, sep="\t")
    # save as narrowPeak (add dummy columns and remove columns not needed)
    df_peak_g.reset_index(names=["name"]).assign(
        score=999, signalValue=".", pValue=".", qValue="."
    )[
        [
            "chrom",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "signalValue",
            "pValue",
            "qValue",
            "peak",
        ]
    ].sort_values(
        ["chrom", "start", "end", "strand", "peak"]
    ).to_csv(
        output_np, sep="\t", index=False, header=False
    )

    # check whether there are overlapping regions (after dropping duplicates)
    for s in ["+", "-"]:
        chr2arr = {
            c: np.zeros(df_peak_g.query(f"chrom == '{c}'").end.max())
            for c in df_peak_g.chrom.unique()
        }
        for name, row in (
            df_peak_g.query(f"strand == '{s}'")[["chrom", "start", "end", "strand"]]
            .drop_duplicates()
            .iterrows()
        ):
            chrom = row["chrom"]
            start = row["start"]
            end = row["end"]
            chr2arr[chrom][start:end] += 1
        for chr_ in chr2arr:
            max_count = pd.Series(chr2arr[chr_]).value_counts().index.max()
            if max_count > 1:
                print(f"WARNING: {chr_} has overlapping regions")

    # # Generate dummy data
    # bed_file, gff3_file = generate_dummy_bed_gff3()

    # # Process dummy data with strand specificity (BED file)
    # result_df_strand_specific = read_genomic_region(
    #     bed_file,
    #     gff3_file,
    #     file_type="bed",
    #     annotation_file_type="gff3",
    #     strand_specific=True,
    # )
    # print("With Strand Specificity (BED file):")
    # print(result_df_strand_specific)

    # # Generate dummy narrowPeak data
    # narrowPeak_file, gff3_file = generate_dummy_narrowPeak()

    # # Process dummy data with strand specificity (narrowPeak file)
    # result_df_narrowPeak_strand_specific = read_genomic_region(
    #     narrowPeak_file,
    #     gff3_file,
    #     file_type="narrowPeak",
    #     annotation_file_type="gff3",
    #     strand_specific=True,
    # )
    # print("\nWith Strand Specificity (narrowPeak file):")
    # print(result_df_narrowPeak_strand_specific)

    # # Generate dummy GTF data
    # bed_file, gtf_file = generate_dummy_gtf()

    # # Process dummy data with strand specificity (GTF file)
    # result_df_gtf_strand_specific = read_genomic_region(
    #     bed_file,
    #     gtf_file,
    #     file_type="bed",
    #     annotation_file_type="gtf",
    #     strand_specific=True,
    # )
    # print("\nWith Strand Specificity (GTF file):")
    # print(result_df_gtf_strand_specific)

    # # Process dummy data without strand specificity (BED file)
    # bed_file, gff3_file = (
    #     generate_dummy_bed_gff3()
    # )  # Regenerate because StringIO object has been read
    # result_df_non_strand_specific = read_genomic_region(
    #     bed_file,
    #     gff3_file,
    #     file_type="bed",
    #     annotation_file_type="gff3",
    #     strand_specific=False,
    # )
    # print("\nWithout Strand Specificity (BED file):")
    # print(result_df_non_strand_specific)
