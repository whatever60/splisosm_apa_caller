from io import StringIO
import os
import subprocess

import pysam
import pandas as pd
import numpy as np
from tqdm import tqdm
from intervaltree import IntervalTree


def get_trees(df: pd.DataFrame) -> dict[str, IntervalTree]:
    trees = {}
    for (chr_, start, end, strand), group in df.groupby(["chrom", "start", "end", "strand"]):
        peaks = group["peak"].to_numpy()
        strands = group["strand"].to_numpy()
        if chr_ not in trees:
            trees[chr_] = IntervalTree()
        trees[chr_][start:end] = {"peaks": peaks, "strands": strands}
    return trees


def get_bam_length(bam: pysam.AlignmentFile | str) -> int:
    """
    Get the total number of reads in an indexed BAM file by summing the number of reads on each contig.

    This function takes a BAM file (either as a pysam.AlignmentFile object or a file 
        path) and uses the BAM index to retrieve the number of reads for each contig. 
        It then sums these values to return the total number of reads.

    Args:
        bam (pysam.AlignmentFile | str): The BAM file, either as a pysam.AlignmentFile 
            object or a file path.

    Returns:
        int: The total number of reads in the BAM file.
    """
    if isinstance(bam, pysam.AlignmentFile):
        bam = bam.filename
    return (
        pd.read_table(StringIO(pysam.idxstats(bam)), header=None)
        .iloc[:, 2:4]
        .fillna(0)
        .sum()
        .sum()
    )


def add_tag(
    bam_in: pysam.AlignmentFile,
    bam_out: pysam.AlignmentFile,
    trees: dict[str, IntervalTree],
) -> None:
    for read in tqdm(
        bam_in.fetch(until_eof=True),
        desc="Processing reads",
        total=get_bam_length(bam_in),
    ):
        if (not read.is_unmapped) and read.reference_name in trees:
            # get 0-based position of the 3' end of the read
            # if read.is_reverse:
            #     read_end = read.reference_start
            # else:
            #     read_end = read.reference_end
            # get 0-based position of center of the read
            read_strand = "-" if read.is_reverse else "+"
            read_end = (read.reference_start + read.reference_end) // 2
            intervals = trees[read.reference_name][read_end]
            if intervals:
                interval = next(iter(intervals))
                peaks = interval.data["peaks"]
                strands = interval.data["strands"]
                # take only the same strand peaks
                peaks = peaks[strands == read_strand]
                if not peaks.size:
                    continue
                closest_peak = peaks[np.argmin(np.abs(peaks - read_end))]
                tag_value = f"{read.reference_name}_{interval.begin}_{interval.end}_{closest_peak}"
                read.set_tag("UP", tag_value, value_type="Z")
                bam_out.write(read)


if __name__ == "__main__":
    # Read the TSV file into a pandas DataFrame
    # data_dir = "/mnt/c/aws_data/data/10x/visium_fresh_frozen_adult_mouse_brain"
    # file_base = (
    #     "Visium_Fresh_Frozen_Adult_Mouse_Brain_possorted_genome_nuclear_bam.sorted"
    # )
    # bam_in = f"{data_dir}/{file_base}.bam"  # same input to macs3 callpeak
    # peak_file = f"{data_dir}/{file_base}.macs3_peaks.tsv"  # output from macs3 callpeak
    # bam_out = f"{data_dir}/{file_base}.macs3_peaks.bam"
    import argparse
    
    parser = argparse.ArgumentParser()
    # parser.add_argument("-s", "--sample", help="sample name")
    parser.add_argument("-i", "--input_bam", help="input bam file")
    parser.add_argument("-p", "--peak_file", help="peak file")
    parser.add_argument("-o", "--output_bam", help="output bam file")
    parser.add_argument("-f", "--force", action="store_true", help="force overwrite")
    args = parser.parse_args()

    peak_file = args.peak_file
    bam_in = args.input_bam
    bam_out = args.output_bam
    force = args.force

    if os.path.isfile(bam_out) and not force:
        raise FileExistsError(f"Output file {bam_out} already exists. Use -f to overwrite.")

    # peak_file = f"/mnt/c/aws_data/data/arora_nc_2022/macs3_out/{sample}_peaks.tsv"
    # bam_in = f"/mnt/c/aws_data/data/arora_nc_2022/bam/{sample}_bam.bam"
    # bam_out = f"/mnt/c/aws_data/data/arora_nc_2022/bam/{sample}_bam.peaks.bam"

    regions_df = pd.read_csv(
        peak_file,
        sep="\t",
        header=0,
        names=["region", "chrom", "start", "end", "strand", "peak", "utr_info"],
    )

    # Create the IntervalTree for each chromosome
    trees = get_trees(regions_df)

    # Add the UP tag to the reads
    input_bam = pysam.AlignmentFile(bam_in, "rb")
    output_bam = pysam.AlignmentFile(bam_out, "wb", template=input_bam)
    add_tag(input_bam, output_bam, trees)
    input_bam.close()
    output_bam.close()
    subprocess.run(["sambamba", "index",  bam_out])
