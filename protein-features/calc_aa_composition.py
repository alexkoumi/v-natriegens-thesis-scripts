"""
Purpose:
    Calculate amino acid composition for unique and non-unique Vibrio natriegens genes.

Description:
    This script reads the PPanGGOLiN matrix.csv file to classify genes as unique
    or non-unique based on the "No. isolates" column. It then parses the
    corresponding RefSeq GenBank files for each strain, extracts CDS protein
    translations by locus tag, and calculates amino acid counts and percentages
    for each gene.

Input:
    - PPanGGOLiN matrix.csv file.
    - RefSeq GenBank folders containing genomic.gbff files for each strain.

Output:
    - Excel workbook containing:
        - unique gene amino acid composition
        - non-unique gene amino acid composition
        - amino acid percentage summary statistics
        - run summary
    - Log file reporting missing files, missing genes, and duplicate entries.

Required Python packages:
    - pandas
    - biopython
    - tqdm
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update PANGENOME_DIR, MATRIX_PATH, OUT_DIR, EXCEL_PATH, and LOG_PATH
    before running the script.
"""

import os
import re
from collections import Counter
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# Paths
PANGENOME_DIR = "path/to/refseq_genbank_folders"
MATRIX_PATH = "path/to/ppanggolin_output/matrix.csv"

OUT_DIR = "path/to/output/amino_acid_composition"
os.makedirs(OUT_DIR, exist_ok=True)

EXCEL_PATH = os.path.join(OUT_DIR, "name-of-file.xlsx")

# Amino acids
AA_ORDER = [
    ("A", "Alanine"), ("R", "Arginine"), ("N", "Asparagine"), ("D", "Aspartic acid"),
    ("C", "Cysteine"), ("Q", "Glutamine"), ("E", "Glutamic acid"), ("G", "Glycine"),
    ("H", "Histidine"), ("I", "Isoleucine"), ("L", "Leucine"), ("K", "Lysine"),
    ("M", "Methionine"), ("F", "Phenylalanine"), ("P", "Proline"), ("S", "Serine"),
    ("T", "Threonine"), ("W", "Tryptophan"), ("Y", "Tyrosine"), ("V", "Valine"),
]

AA_CODES = [aa for aa, name in AA_ORDER]
AA_SET = set(AA_CODES)


# Helper functions
def split_locus_tags(cell_value):
    if pd.isna(cell_value):
        return []

    text = str(cell_value).strip()

    if not text:
        return []

    return [tag.strip().strip('"') for tag in re.findall(r'[^\s"]+', text)]


def read_gene_list_from_matrix(matrix_path):
    matrix = pd.read_csv(matrix_path)

    strain_columns = list(matrix.columns[14:])
    gene_lookup_by_strain = {strain: {} for strain in strain_columns}

    for _, row in matrix.iterrows():
        if int(row["No. isolates"]) == 1:
            gene_type = "unique"
        else:
            gene_type = "non-unique"

        for strain in strain_columns:
            locus_tags = split_locus_tags(row[strain])

            for locus_tag in locus_tags:
                gene_lookup_by_strain[strain][locus_tag] = gene_type

    return strain_columns, gene_lookup_by_strain


def parse_genbank_proteins(gbff_path):
    proteins = {}

    for record in SeqIO.parse(gbff_path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue

            qualifiers = feature.qualifiers
            locus_tags = qualifiers.get("locus_tag", [])
            translations = qualifiers.get("translation", [])

            if not locus_tags or not translations:
                continue

            locus_tag = locus_tags[0].strip()
            sequence = translations[0].replace(" ", "").replace("\n", "").upper()

            proteins[locus_tag] = sequence

    return proteins


def compute_aa_counts_and_percentages(sequence):
    sequence = sequence.upper()
    valid_residues = [aa for aa in sequence if aa in AA_SET]
    valid_length = len(valid_residues)

    counts = Counter(valid_residues)

    count_dict = {}
    percentage_dict = {}

    for aa in AA_CODES:
        count = counts.get(aa, 0)
        count_dict[aa] = count

        if valid_length == 0:
            percentage_dict[aa] = 0.0
        else:
            percentage_dict[aa] = (count / valid_length) * 100

    return count_dict, percentage_dict, valid_length


def build_percentage_summary(df_detail):
    summary_rows = []

    for gene_type in ["unique", "non-unique"]:
        subset = df_detail[df_detail["Gene occurrence type"] == gene_type]

        for aa_code, aa_name in AA_ORDER:
            column = f"{aa_name} (%)"
            values = subset[column].dropna()

            if values.empty:
                mean_value = float("nan")
                sd_value = float("nan")
                median_value = float("nan")
                q1 = float("nan")
                q3 = float("nan")
                iqr = float("nan")
            else:
                mean_value = values.mean()
                sd_value = values.std(ddof=1) if len(values) > 1 else 0.0
                median_value = values.median()
                q1 = values.quantile(0.25)
                q3 = values.quantile(0.75)
                iqr = q3 - q1

            summary_rows.append({
                "Gene occurrence type": gene_type,
                "Amino acid": aa_name,
                "Amino acid code": aa_code,
                "Mean (%)": round(mean_value, 6) if pd.notna(mean_value) else mean_value,
                "SD (%)": round(sd_value, 6) if pd.notna(sd_value) else sd_value,
                "Median (%)": round(median_value, 6) if pd.notna(median_value) else median_value,
                "Q1 (%)": round(q1, 6) if pd.notna(q1) else q1,
                "Q3 (%)": round(q3, 6) if pd.notna(q3) else q3,
                "IQR (%)": round(iqr, 6) if pd.notna(iqr) else iqr,
                "Gene count": len(subset),
            })

    return pd.DataFrame(summary_rows)

# Main
def main():
    strain_columns, gene_lookup_by_strain = read_gene_list_from_matrix(MATRIX_PATH)

    rows = []

    for strain in tqdm(strain_columns, desc="Processing strains"):
        gbff_path = Path(PANGENOME_DIR) / strain / "genomic.gbff"

        if not gbff_path.exists():
            print(f"Missing GenBank file for strain: {strain}")
            continue

        proteins = parse_genbank_proteins(str(gbff_path))
        genes_for_strain = gene_lookup_by_strain[strain]

        for locus_tag, gene_type in genes_for_strain.items():
            sequence = proteins.get(locus_tag)

            if sequence is None:
                continue

            counts, percentages, sequence_length = compute_aa_counts_and_percentages(sequence)

            row = {
                "Strain accession": strain,
                "Locus tag": locus_tag,
                "Sequence length (aa)": sequence_length,
                "Gene occurrence type": gene_type,
            }

            for aa_code, aa_name in AA_ORDER:
                row[f"{aa_name} count"] = counts[aa_code]
                row[f"{aa_name} (%)"] = round(percentages[aa_code], 6)

            rows.append(row)

    df_detail = pd.DataFrame(rows)

    df_unique = df_detail[df_detail["Gene occurrence type"] == "unique"].copy()
    df_non_unique = df_detail[df_detail["Gene occurrence type"] == "non-unique"].copy()
    df_summary = build_percentage_summary(df_detail)

    with pd.ExcelWriter(EXCEL_PATH, engine="openpyxl") as writer:
        df_unique.to_excel(writer, sheet_name="Unique_genes", index=False)
        df_non_unique.to_excel(writer, sheet_name="Non_unique_genes", index=False)
        df_summary.to_excel(writer, sheet_name="AA_percentage_summary", index=False)

    print(f"Excel written to: {EXCEL_PATH}")

if __name__ == "__main__":
    main()
