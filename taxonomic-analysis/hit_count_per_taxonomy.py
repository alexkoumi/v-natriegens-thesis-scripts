"""
Script name:
    hit_count_per_taxonomy.py

Purpose:
    Generate taxonomic hit-count summaries from MMseqs2 taxonomy TSV files.

Description:
    This script reads a taxonomy-annotated MMseqs2 TSV file in chunks and counts
    the number of hits assigned to each species, genus, and family. It writes the
    resulting summaries to an Excel workbook with separate sheets for each
    taxonomic rank. It also generates a report listing entries with missing
    species, genus, or family assignments.

Input:
    - Taxonomy-annotated MMseqs2 TSV file containing taxonomy information.
    - Expected columns include:
        - query
        - target
        - Species
        - Genus
        - Family

Output:
    - Excel file containing species-, genus-, and family-level hit-count summaries.
    - Text report listing entries with missing taxonomic ranks.

Required Python packages:
    - pandas
    - tqdm
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update the domain, dataset, and base_dir values before running the script.
"""

import pandas as pd
import os
from tqdm import tqdm

# file paths
domain = "domain name"
dataset = "dataset name (all/top)"
base_dir = f"/path/to/mmseqs_hits_taxonomy/{dataset}_taxonomy/{domain}"
input_file = os.path.join(base_dir, f"{domain}_{dataset}_taxonomy.tsv")
output_excel = os.path.join(base_dir, f"{domain}_{dataset}_taxonomy_counts.xlsx")
output_report = os.path.join(base_dir, f"{domain}_{dataset}_taxonomy_missing_ranks.txt")

# parameters
CHUNK_SIZE = 100000   # process in chunks for memory

# main workflow
def process_taxonomy_tsv():
    print(f"Reading taxonomy TSV: {input_file}")
    if not os.path.exists(input_file):
        raise SystemExit(f"File not found: {input_file}")

    species_counts, genus_counts, family_counts = {}, {}, {}
    missing_species, missing_genus, missing_family = [], [], []

    for chunk in tqdm(pd.read_csv(input_file, sep="\t", chunksize=CHUNK_SIZE), desc="Processing chunks", unit="chunk"):
        for _, row in chunk.iterrows():
            q, t = str(row.get("query", "")).strip(), str(row.get("target", "")).strip()
            s = (str(row.get("Species", "")).strip() if pd.notna(row.get("Species")) else "")
            g = (str(row.get("Genus", "")).strip() if pd.notna(row.get("Genus")) else "")
            f = (str(row.get("Family", "")).strip() if pd.notna(row.get("Family")) else "")

            # Count only non-empty ranks
            if s:
                species_counts[s] = species_counts.get(s, 0) + 1
            else:
                missing_species.append((q, t))

            if g:
                genus_counts[g] = genus_counts.get(g, 0) + 1
            else:
                missing_genus.append((q, t))

            if f:
                family_counts[f] = family_counts.get(f, 0) + 1
            else:
                missing_family.append((q, t))

    #save results
    df_species = pd.DataFrame(sorted(species_counts.items()), columns=["Species", "Hit_Count"])
    df_genus = pd.DataFrame(sorted(genus_counts.items()), columns=["Genus", "Hit_Count"])
    df_family = pd.DataFrame(sorted(family_counts.items()), columns=["Family", "Hit_Count"])

    with pd.ExcelWriter(output_excel) as writer:
        df_species.to_excel(writer, sheet_name="Species", index=False)
        df_genus.to_excel(writer, sheet_name="Genus", index=False)
        df_family.to_excel(writer, sheet_name="Family", index=False)

    print(f"Excel summary written: {output_excel}")

    # write the report file
    with open(output_report, "w") as f:
        f.write(f"=== {domain.upper()} TAXONOMY SUMMARY ===\n")
        f.write(f"Total distinct Species: {len(species_counts)}\n")
        f.write(f"Total distinct Genus: {len(genus_counts)}\n")
        f.write(f"Total distinct Family: {len(family_counts)}\n\n")

        f.write(f"--- Missing Species ({len(missing_species)} entries) ---\n")
        for q, t in missing_species:
            f.write(f"{q}\t{t}\n")

        f.write(f"\n--- Missing Genus ({len(missing_genus)} entries) ---\n")
        for q, t in missing_genus:
            f.write(f"{q}\t{t}\n")

        f.write(f"\n--- Missing Family ({len(missing_family)} entries) ---\n")
        for q, t in missing_family:
            f.write(f"{q}\t{t}\n")

    print(f"Missing-rank report written: {output_report}")

# run
if __name__ == "__main__":
    process_taxonomy_tsv()