"""
Purpose:
    Count the number of genes assigned to each PPanGGOLiN pangenome partition.

Description:
    This script reads the PPanGGOLiN matrix.csv file and counts the number of gene
    locus tags assigned to each pangenome partition category. It then writes a text
    summary reporting gene counts per partition and the total number of counted genes.

Input:
    - PPanGGOLiN presence/absence matrix file.

Output:
    - Text file containing gene counts per pangenome partition.

Required Python packages:
    - pandas
"""

import pandas as pd
import re

# File paths
matrix_path = "path/to/ppanggolin_output/matrix.csv"
output_txt = "path/to/output/ppanggolin_gene_counts_per_partition.txt"

# User-defined columns
partition_col = "Non-unique Gene name"
first_genome_column_index = 14

print("Reading PPanGGOLiN matrix file")
df = pd.read_csv(matrix_path)

# Genome columns contain locus tags
genome_cols = list(df.columns[first_genome_column_index:])

# Function to count gene IDs in a cell
def count_genes(cell):
    if pd.isna(cell):
        return 0

    s = str(cell).strip().strip('"')

    if not s:
        return 0

    parts = re.split(r"[;,\s]+", s)

    return len([p for p in parts if p])


# Count genes per family
print("Counting gene occurrences per family")

df["gene_count"] = df[genome_cols].applymap(count_genes).sum(axis=1)

# Sum by partition
print("Summarising gene counts per partition")

partition_counts = (
    df.groupby(partition_col)["gene_count"]
    .sum()
    .sort_values(ascending=False)
)

total_genes = partition_counts.sum()

# Write output
print("Writing output file")

with open(output_txt, "w") as out:
    out.write("PPanGGOLiN gene counts per partition\n\n")

    for partition, count in partition_counts.items():
        out.write(f"{partition}: {count}\n")

    out.write("-" * 25 + "\n")
    out.write(f"Total genes: {total_genes}\n")

print("Done.")
print(partition_counts)
print(f"Total genes: {total_genes}")
print(f"Output written to: {output_txt}")
