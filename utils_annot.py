"""
Utilities for parsing genome annotation files.
"""

import re
import pandas as pd
from joblib import Parallel, delayed
from tqdm.auto import tqdm


def _parse_attributes_gff3(attribute: str) -> tuple[str, str]:
    """
    Parse the attribute string from a GFF3 annotation file to extract the parent gene ID and the element ID.

    Args:
        attribute (str): Attribute string from a GFF3 file, containing multiple key-value pairs separated by semicolons.

    Returns:
        tuple: A tuple containing the parent gene ID and the element ID.
    """
    parent_gene_id, element_id = None, None
    for attr in attribute.split(";"):
        if attr.startswith("Parent=gene:"):
            parent_gene_id = attr.split(":")[1]
        elif attr.startswith("ID=transcript:"):
            element_id = attr.split(":")[1]
    return parent_gene_id, element_id


def _parse_attributes_gtf(attribute: str) -> tuple[str, str]:
    """
    Parse the attribute string from a GTF annotation file to extract the parent gene ID and the element ID.

    Args:
        attribute (str): Attribute string from a GTF file, containing multiple key-value pairs separated by semicolons.

    Returns:
        tuple: A tuple containing the parent gene ID and the element ID.
    """
    parent_gene_id, element_id = None, None
    for attr in attribute.split(";"):
        key_value = attr.strip().split(" ")
        if len(key_value) == 2:
            key, value = key_value
            value = value.strip('"')
            if key == "gene_id":
                parent_gene_id = value
            elif key == "transcript_id":
                element_id = value
    return parent_gene_id, element_id


def _parse_attributes_gencode_gff3(attributes: str) -> tuple[str, str]:
    """
    Parse the attributes column for gencode-gff3 to extract the parent gene ID and transcript ID.

    Args:
        attributes (str): The attributes string from a GFF3 entry.

    Returns:
        tuple: A tuple containing the parent gene ID and transcript ID.
    """
    parent_gene_id, element_id, gene_id = None, None, None
    for attr in attributes.split(";"):
        if attr.startswith("Parent="):
            parent_gene_id = attr.split("=")[1]
        elif attr.startswith("ID="):
            element_id = attr.split("=")[1]
        elif attr.startswith("gene_id="):
            gene_id = attr.split("=")[1]
    if gene_id != parent_gene_id:
        parent_gene_id = None
    return parent_gene_id, element_id

def _check_format(df: pd.DataFrame, flavor: str, num_rows: int = 100) -> bool:
    """
    Check if the annotation format of the DataFrame matches the expected flavor.

    Args:
        df (pd.DataFrame): DataFrame containing annotation data.
        flavor (str): Flavor of the annotation ("ensembl-gff3", "space_ranger-gtf", or "gencode-gff3").
        num_rows (int): Number of rows to check. Default is 100.

    Returns:
        bool: True if the format matches the expected flavor, False otherwise.

    Raises:
        ValueError: If an invalid flavor is provided.
    """
    if flavor == "ensembl-gff3":
        expected_pattern = {
            "transcript": re.compile(r"Parent=gene:[^;]+.*ID=transcript:[^;]+"),
            "gene": re.compile(r"ID=gene:[^;]+.*Name=[^;]+"),
        }
    elif flavor == "space_ranger-gtf":
        expected_pattern = {
            "transcript": re.compile(r"gene_id \"[^\"]+\";.*transcript_id \"[^\"]+\";"),
            "gene": re.compile(r"gene_id \"[^\"]+\";.*gene_name \"[^\"]+\";"),
        }
    elif flavor == "gencode-gff3":
        expected_pattern = {
            "transcript": re.compile(r"ID=[^;]+.*Parent=[^;]+"),
            "gene": re.compile(r"ID=[^;]+.*gene_name=[^;]+"),
        }
    else:
        raise ValueError(
            "Invalid flavor argument. Must be 'ensembl-gff3', 'space_ranger-gtf', or 'gencode-gff3'."
        )

    right = False
    for i in range(min(num_rows, len(df))):
        attribute = df.iloc[i]["attribute"]
        if expected_pattern["transcript"].search(attribute) or expected_pattern[
            "gene"
        ].search(attribute):
            right = True
    return right


def get_tx2gene(
    df: pd.DataFrame, flavor: str = "ensembl-gff3", n_jobs: int = 4
) -> dict[str, str]:
    """
    Extract a mapping from transcript IDs to gene IDs from the annotation DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing annotation data.
        flavor (str): Flavor of the annotation ("ensembl-gff3" or "space_ranger-gtf").
        n_jobs (int): Number of parallel jobs to run. Default is 4.

    Returns:
        dict: A dictionary mapping transcript IDs to gene IDs.

    Raises:
        ValueError: If an invalid flavor is provided.

    This function checks the format of the DataFrame to ensure it matches the expected
    pattern for the given flavor. It then processes each row to extract the parent gene ID
    and the element ID, creating a dictionary that maps each transcript ID to its corresponding
    gene ID. Inconsistent gene names for the same transcript ID are reported.
    """
    if flavor not in {"ensembl-gff3", "space_ranger-gtf", "gencode-gff3"}:
        raise ValueError(
            "Invalid flavor argument. Must be 'ensembl-gff3' or 'space_ranger-gtf' or 'gencode-gff3'."
        )

    format_right = _check_format(df, flavor)
    if not format_right:
        print(
            f"Warning: The attributes in annotation do not match the expected format for {flavor}."
        )

    if flavor == "ensembl-gff3":
        parse_attributes = _parse_attributes_gff3
    elif flavor == "gencode-gff3":
        parse_attributes = _parse_attributes_gencode_gff3
    elif flavor == "space_ranger-gtf":
        # since gtf file does not support hierarchical structure, and gtf does have a
        # type of feature called "transcript" (in comparison gff3 decompose it into
        # mRNA, lncRNA, etc.), we filter to transcript
        df = df.query("feature == 'transcript'")
        parse_attributes = _parse_attributes_gtf

    else:
        raise ValueError(
            "Invalid flavor argument. Must be 'ensembl-gff3' or 'space_ranger-gtf'."
        )

    def process_row(attributes: str) -> tuple[str, str]:
        parent_gene_id, element_id = parse_attributes(attributes)
        return parent_gene_id, element_id

    results = Parallel(n_jobs=n_jobs)(
        delayed(process_row)(attr) for attr in tqdm(df["attribute"].tolist())
    )

    to_gene = {}
    for parent_gene_id, element_id in results:
        if parent_gene_id and element_id:
            if element_id in to_gene:
                if not to_gene[element_id] == parent_gene_id:
                    print(
                        f"Inconsistent gene names for {element_id}: {to_gene[element_id]} vs {parent_gene_id}"
                    )
            else:
                to_gene[element_id] = parent_gene_id

    return to_gene


def _parse_gene_attributes_gff3(attribute: str) -> tuple[str, str]:
    """
    Parse the attribute string from a GFF3 annotation file to extract the gene ID and gene name.

    Args:
        attribute (str): Attribute string from a GFF3 file, containing multiple key-value pairs separated by semicolons.

    Returns:
        tuple: A tuple containing the gene ID and gene name.
    """
    gene_id, gene_name = None, None
    for attr in attribute.split(";"):
        if attr.startswith("ID=gene:"):
            gene_id = attr.split(":")[1]
        elif attr.startswith("Name="):
            gene_name = attr.split("=")[1]
    return gene_id, gene_name


def _parse_gene_attributes_gtf(attribute: str) -> tuple[str, str]:
    """
    Parse the attribute string from a GTF annotation file to extract the gene ID and gene name.

    Args:
        attribute (str): Attribute string from a GTF file, containing multiple key-value pairs separated by semicolons.

    Returns:
        tuple: A tuple containing the gene ID and gene name.
    """
    gene_id, gene_name = None, None
    for attr in attribute.split(";"):
        key_value = attr.strip().split(" ")
        if len(key_value) == 2:
            key, value = key_value
            value = value.strip('"')
            if key == "gene_id":
                gene_id = value
            elif key == "gene_name":
                gene_name = value
    return gene_id, gene_name


def _parse_gene_attributes_gencode_gff3(attributes: str) -> tuple[str, str]:
    """
    Parse the attributes column for gencode-gff3 to extract the gene ID and gene name.

    Args:
        attributes (str): The attributes string from a GFF3 entry.

    Returns:
        tuple: A tuple containing the gene ID and gene name.
    """
    gene_id, gene_name = None, None
    for attr in attributes.split(";"):
        if attr.startswith("gene_id="):
            gene_id = attr.split("=")[1]
        elif attr.startswith("gene_name="):
            gene_name = attr.split("=")[1]
    return gene_id, gene_name
    


def get_gene2name(
    df: pd.DataFrame, flavor: str = "ensembl-gff3", n_jobs: int = 4
) -> dict[str, str]:
    """
    Extract a mapping from gene IDs to gene names from the annotation DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing annotation data.
        flavor (str): Flavor of the annotation ("ensembl-gff3" or "space_ranger-gtf").
        n_jobs (int): Number of parallel jobs to run. Default is 4.

    Returns:
        dict: A dictionary mapping gene IDs to gene names.

    Raises:
        ValueError: If an invalid flavor is provided.

    This function checks the format of the DataFrame to ensure it matches the expected
    pattern for the given flavor. It then processes each row to extract the gene ID
    and gene name, creating a dictionary that maps each gene ID to its corresponding
    gene name. Inconsistent gene names for the same gene ID are reported.
    """
    if flavor not in {"ensembl-gff3", "space_ranger-gtf", "gencode-gff3"}:
        raise ValueError(
            "Invalid flavor argument. Must be 'ensembl-gff3' or 'space_ranger-gtf' or 'gencode-gff3'."
        )

    format_right = _check_format(df, flavor)
    if not format_right:
        print(
            f"Warning: The attributes in annotation do not match the expected format for {flavor}."
        )

    if flavor == "ensembl-gff3":
        parse_gene_attributes = _parse_gene_attributes_gff3
    elif flavor == "gencode-gff3":
        parse_gene_attributes = _parse_gene_attributes_gencode_gff3
    elif flavor == "space_ranger-gtf":
        parse_gene_attributes = _parse_gene_attributes_gtf
    else:
        raise ValueError(
            "Invalid flavor argument. Must be 'ensembl-gff3' or 'space_ranger-gtf'."
        )

    def process_gene_row(attributes: str) -> tuple[str, str]:
        gene_id, gene_name = parse_gene_attributes(attributes)
        return gene_id, gene_name

    gene_rows = df.query("feature == 'gene'")
    results = Parallel(n_jobs=n_jobs)(
        delayed(process_gene_row)(attr)
        for attr in tqdm(gene_rows["attribute"].tolist())
    )

    gene2name = {}
    for gene_id, gene_name in results:
        if gene_id is not None:
            if gene_id in gene2name:
                if not gene2name[gene_id] == gene_name:
                    print(
                        f"Inconsistent gene names for {gene_id}: {gene2name[gene_id]} vs {gene_name}"
                    )
            else:
                gene2name[gene_id] = gene_name

    return gene2name


if __name__ == "__main__":
    from utils import load_annotation

    # a sample of 10000 lines
    gtf_file = "/mnt/c/aws_data/data/10x/space_ranger/reference/refdata-gex-GRCh38-2020-A/genes/genes.sorted.sample.gtf"
    gff3_file = "/mnt/c/aws_data/data/ensembl/pub/release-98/gff3/mus_musculus/Mus_musculus.GRCm38.98.chr_prefix.sample.gff3.gz"
    gff3_gencode_file = "/mnt/c/aws_data/data/gencode/Gencode_human/release_32/gencode.v32.annotation.sample.gff3"

    df_gtf = load_annotation(gtf_file)
    to_gene_gtf = get_tx2gene(df_gtf, flavor="space_ranger-gtf")
    gene2name_gtf = get_gene2name(df_gtf, flavor="space_ranger-gtf")
    print(
        f"Number of transcript IDs: {len(to_gene_gtf)}, Number of gene IDs: {len(gene2name_gtf)}"
    )

    df_gff3 = load_annotation(gff3_file)
    to_gene_gff3 = get_tx2gene(df_gff3, flavor="ensembl-gff3")
    gene2name_gff3 = get_gene2name(df_gff3, flavor="ensembl-gff3")
    print(
        f"Number of transcript IDs: {len(to_gene_gff3)}, Number of gene IDs: {len(gene2name_gff3)}"
    )

    df_gff3_gencode = load_annotation(gff3_gencode_file)
    to_gene_gff3_gencode = get_tx2gene(df_gff3_gencode, flavor="gencode-gff3")
    gene2name_gff3_gencode = get_gene2name(df_gff3_gencode, flavor="gencode-gff3")
    print(
        f"Number of transcript IDs: {len(to_gene_gff3_gencode)}, Number of gene IDs: {len(gene2name_gff3_gencode)}"
    )
