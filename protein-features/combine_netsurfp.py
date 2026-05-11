"""
Purpose:
    Combine per-residue NetSurfP output files into one gene-level summary table.

Description:
    This script reads NetSurfP result CSV files from two folders:
        - one folder containing unique gene results
        - one folder containing non-unique gene results

    Each NetSurfP file contains per-residue predictions for one protein. For each
    protein, the script calculates the average RSA, ASA, phi, psi, and disorder
    values across residues. It also calculates the fraction of residues assigned
    to each Q3 and Q8 secondary structure state.

Input:
    - Folder containing NetSurfP CSV files for unique genes.
    - Folder containing NetSurfP CSV files for non-unique genes.

Output:
    - Combined CSV file with one row per gene.

Required Python packages:
    - pandas
    - tqdm

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update UNIQUE_ROOT, NONUNIQUE_ROOT, and OUTPUT_FILE before running.
"""

from pathlib import Path
import pandas as pd
from tqdm import tqdm

# Paths
UNIQUE_ROOT = "path/to/netsurfp_results/Vnat_unique_all"
NONUNIQUE_ROOT = "path/to/netsurfp_results/Vnat_nonunique_all"
OUTPUT_FILE = "path/to/output/netsurfp_results_combined.csv"

# Helper functions
def find_csv_files(root_folder):
    """
    Find all CSV files inside a folder and its subfolders.
    """
    return sorted(root_folder.rglob("*.csv"))

def clean_column_names(df):
    """
    Remove accidental spaces from column names.
    """
    df.columns = df.columns.str.strip()
    return df


def extract_locus_tag(csv_file, df):
    """
    Extract locus tag from the NetSurfP file.

    First, try the id column.
    If that is not available, use the folder or file name.
    """
    if "id" in df.columns and len(df) > 0:
        gene_id = str(df["id"].iloc[0]).strip().replace(">", "")

        if gene_id:
            return gene_id

    name = csv_file.stem

    if "_" in name:
        parts = name.split("_", 1)
        return parts[1]

    parent_name = csv_file.parent.name

    if "_" in parent_name:
        parts = parent_name.split("_", 1)
        return parts[1]

    return name

def calculate_fraction(df, column, state):
    """
    Calculate the fraction of residues assigned to a specific Q3 or Q8 state.
    """
    if column not in df.columns or len(df) == 0:
        return 0.0

    values = df[column].astype(str).str.strip()
    return (values == state).mean()


def calculate_mean(df, column):
    """
    Calculate the mean value of a numeric column.
    """
    if column not in df.columns:
        return float("nan")

    values = pd.to_numeric(df[column], errors="coerce")
    return values.mean()


def summarise_netsurfp_file(csv_file, gene_occurrence_type):
    """
    Read one NetSurfP CSV file and return one gene-level summary row.
    """
    df = pd.read_csv(csv_file)
    df = clean_column_names(df)

    locus_tag = extract_locus_tag(csv_file, df)

    row = {
        "locus_tag": locus_tag,
        "gene_occurrence_type": gene_occurrence_type,
        "rsa": calculate_mean(df, "rsa"),
        "asa": calculate_mean(df, "asa"),
        "phi": calculate_mean(df, "phi"),
        "psi": calculate_mean(df, "psi"),
        "disorder": calculate_mean(df, "disorder"),
        "frac_q3_H": calculate_fraction(df, "q3", "H"),
        "frac_q3_E": calculate_fraction(df, "q3", "E"),
        "frac_q3_C": calculate_fraction(df, "q3", "C"),
        "frac_q8_G": calculate_fraction(df, "q8", "G"),
        "frac_q8_H": calculate_fraction(df, "q8", "H"),
        "frac_q8_I": calculate_fraction(df, "q8", "I"),
        "frac_q8_B": calculate_fraction(df, "q8", "B"),
        "frac_q8_E": calculate_fraction(df, "q8", "E"),
        "frac_q8_S": calculate_fraction(df, "q8", "S"),
        "frac_q8_T": calculate_fraction(df, "q8", "T"),
        "frac_q8_C": calculate_fraction(df, "q8", "C"),
    }

    return row


def process_folder(root_folder, gene_occurrence_type):
    """
    Process all NetSurfP CSV files in one root folder.
    """
    csv_files = find_csv_files(root_folder)
    rows = []

    print(f"Found {len(csv_files)} CSV files for {gene_occurrence_type} genes.")

    for csv_file in tqdm(csv_files, desc=f"Processing {gene_occurrence_type} genes"):
        try:
            row = summarise_netsurfp_file(csv_file, gene_occurrence_type)
            rows.append(row)
        except Exception as error:
            print(f"Could not process file: {csv_file}")
            print(f"Reason: {error}")

    return rows


# Main
def main():
    all_rows = []

    all_rows.extend(process_folder(UNIQUE_ROOT, "unique"))
    all_rows.extend(process_folder(NONUNIQUE_ROOT, "non-unique"))

    combined_df = pd.DataFrame(all_rows)

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    combined_df.to_csv(OUTPUT_FILE, index=False)

    print("Done.")
    print(f"Combined NetSurfP file written to: {OUTPUT_FILE}")
    print(f"Total genes processed: {len(combined_df)}")

if __name__ == "__main__":
    main()