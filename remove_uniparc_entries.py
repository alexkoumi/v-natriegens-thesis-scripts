"""
Purpose:
    Remove UniParc-associated or taxonomically unresolved MMseqs2 hits from results file.

Description:
    This script reads a taxonomy-annotated MMseqs2 hit table and removes entries whose
    target accessions are present in a pre-existing removal cache. In the original
    analysis, this cache was used to identify accessions corresponding to UniParc-only
    entries or entries without usable UniProt taxonomy information.

    The script also removes rows with blank TaxID values or completely blank taxonomy
    fields across the selected taxonomy columns.

Input:
    - Taxonomy-annotated MMseqs2 hit table in TSV format.
    - JSON cache containing target accessions to remove.

Output:
    - Cleaned taxonomy TSV file with unresolved/UniParc-associated entries removed.
    - Updated JSON cache containing the removed accessions.

Required Python packages:
    - pandas

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update TAXONOMY_TSV, EXISTING_REMOVAL_CACHE, CLEANED_OUTPUT_TSV, and
    UPDATED_REMOVAL_CACHE before running the script.
"""

import pandas as pd
import json
import os

# file paths
taxonomy_file = "path/to/input/taxonomy.tsv"
existing_cache = "path/to/input/taxonomy_cache.json"
output_file = "path/to/output/top_hits_taxonomy_cleaned.tsv"
new_cache = "path/to/output/removed_total_accessions.json"

# Read taxonomy file
print(f"Reading taxonomy file: {taxonomy_file}")
if not os.path.exists(taxonomy_file):
    raise SystemExit(f"File not found: {taxonomy_file}")

df = pd.read_csv(taxonomy_file, sep="\t", dtype=str).fillna("")
print(f"Total rows before cleaning: {len(df)}")

# Load existing removal cache
if os.path.exists(existing_cache):
    with open(existing_cache, "r") as f:
        cached_accessions = set(json.load(f))
    print(f"Loaded {len(cached_accessions)} accessions from cache.")
else:
    cached_accessions = set()
    print("No existing cache found.")

# Identify rows to remove
taxonomy_cols = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
remove_mask = (
    df["target"].isin(cached_accessions)
    | (df["TaxID"].str.strip() == "")
    | df[taxonomy_cols].apply(lambda row: all(x.strip() == "" for x in row), axis=1)
)

removed_df = df[remove_mask]
df_clean = df[~remove_mask]

print(f"Removed {len(removed_df)} rows (cached + blank taxonomy/TaxID).")
print(f"Remaining: {len(df_clean)}")

# Save cleaned taxonomy file
df_clean.to_csv(output_file, sep="\t", index=False)
print(f"Cleaned taxonomy saved: {output_file}")

# Save updated removal cache
removed_accessions = sorted(set(removed_df["target"].dropna().tolist()))
combined_cache = sorted(cached_accessions.union(removed_accessions))

with open(new_cache, "w") as f:
    json.dump(combined_cache, f, indent=2)

print(f"Combined cache saved: {new_cache}")
print("\nProcessing complete.")