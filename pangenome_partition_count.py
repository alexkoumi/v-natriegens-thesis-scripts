"""
Purpose:
    Summarise the taxonomic composition of UniProt reference proteomes.

Description:
    This script reads an Excel file containing reference proteomes for a
    selected domain and calculates the number of proteomes and total number of
    proteins represented at species, genus, family, and phylum levels.

Input:
    - Excel file containing reference proteomes.
    - Required columns:
        - "Proteome Id"
        - "Protein count"
        - "Taxonomic lineage"

Output:
    - Excel workbook containing separate summary sheets for:
        - Species-level counts
        - Genus-level counts
        - Family-level counts
        - Phylum-level counts

Required Python packages:
    - pandas
    - openpyxl
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
