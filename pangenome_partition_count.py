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

# User-defined inputs
domain = "domain name"
input_file = f"{domain}_reference_proteomes.xlsx"
df = pd.read_excel(input_file)

#Clean column names for consistency
df.columns = df.columns.str.strip()

#Extract Phylum, Family, and Species from 'Taxonomic lineage'
def extract_taxon(lineage, position):
    try:
        parts = [p.strip() for p in str(lineage).split(",")]
        if len(parts) >= position:
            return parts[position - 1]
        else:
            return "Unknown"
    except Exception:
        return "Unknown"

def extract_last_taxon(lineage):
    try:
        parts = [p.strip() for p in str(lineage).split(",")]
        return parts[-1] if parts else "Unknown"
    except Exception:
        return "Unknown"

df["Phylum"] = df["Taxonomic lineage"].apply(lambda x: extract_taxon(x, 4))
df["Family"] = df["Taxonomic lineage"].apply(lambda x: extract_taxon(x, 7))
df["Genus"] = df["Taxonomic lineage"].apply(lambda x: extract_taxon(x, 8))
df["Species"] = df["Taxonomic lineage"].apply(extract_last_taxon)

#Create summaries
def summarize_by(column_name):
    return (
        df.groupby(column_name)
        .agg(
            Proteome_count=("Proteome Id", "count"),
            Total_protein_count=("Protein count", "sum")
        )
        .reset_index()
        .sort_values(by="Proteome_count", ascending=False)
    )

species_summary = summarize_by("Species")
family_summary = summarize_by("Family")
phylum_summary = summarize_by("Phylum")
genus_summary = summarize_by("Genus")

#Save results to Excel
output_file = f"{domain}_reference_proteomes_counts.xlsx"
with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
    species_summary.to_excel(writer, index=False, sheet_name="Species_summary")
    family_summary.to_excel(writer, index=False, sheet_name="Family_summary")
    phylum_summary.to_excel(writer, index=False, sheet_name="Phylum_summary")
    genus_summary.to_excel(writer, index=False, sheet_name="Genus_summary")
print("Analysis complete, results saved to:", output_file)
