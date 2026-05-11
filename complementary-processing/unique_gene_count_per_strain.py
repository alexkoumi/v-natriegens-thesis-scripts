"""
Purpose:
    Identify strain-specific genes from the PPanGGOLiN matrix.csv output file.

Description:
    This script reads the PPanGGOLiN gene presence/absence matrix and selects genes
    where the "No. isolates" column is equal to 1. These genes are interpreted as
    strain-specific genes because they are present in only one strain.

    For each strain, the script extracts the corresponding locus tags, counts the
    number of strain-specific genes, and records duplicated genes where applicable.

Input:
    - PPanGGOLiN matrix.csv file.

Outputs:
    - A fasta-formatted file containing strain names as headers, followed by the gene locus tag.
    - An Excel file summarising the number of strain-specific genes per strain.

Required Python packages:
    - pandas
    - tqdm
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update the INPUT_MATRIX_CSV, OUTPUT_GENE_LIST, and OUTPUT_EXCEL_SUMMARY
    variables before running the script.
"""

import pandas as pd
from tqdm import tqdm

# File paths
matrix_path = "path/to/ppanggolin_output/matrix.csv"
output_fasta = "path/to/output/unique_gene_list_per_strain.fasta"
output_excel = "path/to/output/unique_gene_count_per_strain.xlsx"


# Read PPanGGOLiN matrix and identify strain-specific genes
print("reading matrix file")
matrix_file = pd.read_csv(matrix_path)
print("identifying unique genes")
unique_genes = matrix_file[matrix_file["No. isolates"] == 1]
strain_columns = matrix_file.columns[14:]


# Extract strain-specific locus tags for each strain
strain_with_locus_tags = {}
duplicate_counts = {}
print("extracting unique genes for each strain")
for strain in tqdm(strain_columns, desc="processing strains"):
    strain_specific_tags = unique_genes[strain].dropna().tolist()
    cleaned_tags = []
    duplicates = 0
    for entry in strain_specific_tags:
        for raw_tag in str(entry).split():
            tag = raw_tag.strip('"')
            if tag:
                if raw_tag.startswith('"') and raw_tag.endswith('"'):
                    cleaned_tags.append(f"{tag} (duplicate)")
                    duplicates += 1
                else:
                    cleaned_tags.append(tag)
    strain_with_locus_tags[strain] = cleaned_tags
    duplicate_counts[strain] = duplicates


# Write per-strain gene list output
print("writing FASTA file")
with open(output_fasta, "w") as fasta_file:
    for strain, locus_tags in tqdm(strain_with_locus_tags.items(), desc="Writing FASTA"):
        fasta_file.write(f">{strain} (count = {len(locus_tags)})\n")
        for locus_tag in locus_tags:
            fasta_file.write(f"{locus_tag}\n")
        fasta_file.write("\n")

# Write Excel summary table
print("Writing Excel file")
data = {
    "Strain": list(strain_with_locus_tags.keys()),
    "Unique gene count": [len(v) for v in strain_with_locus_tags.values()],
    "Duplicate unique genes": [duplicate_counts[s] for s in strain_with_locus_tags.keys()],
}
df_counts = pd.DataFrame(data)
df_counts.to_excel(output_excel, index=False)

# Print summary
total_unique_locus_tags = sum(len(v) for v in strain_with_locus_tags.values())
print(f"Total unique genes including duplicates: {total_unique_locus_tags}")
print(f"Number of strains processed: {len(strain_with_locus_tags)}")
print(f"Outputs generated:\n - {output_fasta}\n - {output_excel}")