"""
Purpose:
    Split bacterial MMseqs2 taxonomy results into Vibrio and non-Vibrio entries.

Description:
    This script reads a taxonomy-annotated MMseqs2 TSV file in chunks and separates
    rows where the Genus column is equal to "Vibrio" from rows assigned to other
    genera. The outputs are written as two separate TSV files.

Input:
    - Taxonomy-annotated MMseqs2 TSV file.
    - Expected column:
        - Genus

Output:
    - TSV file containing rows assigned to the genus Vibrio.
    - TSV file containing rows not assigned to the genus Vibrio.

Required Python packages:
    - pandas

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update input_file, output_vibrio, and output_nonvibrio before running.
"""

import pandas as pd

# file paths
input_file = "path/to/bacterial_taxonomy.tsv"
output_vibrio = "path/to/output/vibrio_taxonomy.tsv"
output_nonvibrio = "path/to/output/non_vibrio_taxonomy.tsv"

# Parameters
chunksize = 200_000

print(f"Processing single TSV file: {input_file}")

first_chunk = True

for chunk_index, chunk in enumerate(pd.read_csv(input_file, sep="\t", chunksize=chunksize), start=1):
    print(f"  Processing chunk {chunk_index} with {len(chunk)} rows")

    # Split by genus
    vibrio_chunk = chunk[chunk["Genus"] == "Vibrio"]
    nonvibrio_chunk = chunk[chunk["Genus"] != "Vibrio"]

    # Save chunks to separate files
    if first_chunk:
        vibrio_chunk.to_csv(output_vibrio, sep="\t", index=False, mode="w")
        nonvibrio_chunk.to_csv(output_nonvibrio, sep="\t", index=False, mode="w")
        first_chunk = False
    else:
        vibrio_chunk.to_csv(output_vibrio, sep="\t", index=False, mode="a", header=False)
        nonvibrio_chunk.to_csv(output_nonvibrio, sep="\t", index=False, mode="a", header=False)

    print(f" Wrote {len(vibrio_chunk)} Vibrio rows and {len(nonvibrio_chunk)} non-Vibrio rows")

print("\nFinished creating:")
print(f" - {output_vibrio}")
print(f" - {output_nonvibrio}")