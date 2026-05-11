"""
Purpose:
    Add locus tags to GO annotation Excel sheets using GenBank protein IDs.

Description:
    This script reads selected sheets from an Excel workbook, looks for the
    "genbank protein id" column, and uses a reference GenBank file to map each
    GenBank protein ID to its corresponding locus tag. The script then adds or
    overwrites a column named "locus tag" in each processed sheet.

Input:
    - Excel workbook containing GO annotation sheets.
    - Reference genome GenBank file in .gbff format.

Output:
    - Updated Excel workbook containing locus tags.

Required Python packages:
    - pandas
    - biopython
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update EXCEL_IN, EXCEL_OUT, and GBFF_PATH before running the script.
"""

import os
import re
import pandas as pd
from Bio import SeqIO

# Paths
EXCEL_IN  = "path/to/excel/with/reference genome annotations.xlsx"
EXCEL_OUT = "path/to/the new/excel/with/reference genome annotations.xlsx"
GBFF_PATH = "path/to/reference_genome/genomic.gbff"

SHEETS_TO_PROCESS = ["GenBank_GO", "UniProt_GO"]

PROTEIN_COL = "genbank protein id"
LOCUS_COL   = "locus tag"


# Helper functions
def clean_cell(x):
    if x is None:
        return ""
    if isinstance(x, float) and pd.isna(x):
        return ""
    s = str(x).strip()
    return "" if s.lower() in {"nan", "none"} else s


def extract_first_protein_id(cell):
    cell = clean_cell(cell)
    if not cell:
        return ""

    parts = re.split(r"[;\s,]+", cell)

    for p in parts:
        p = p.strip()
        if p:
            return p

    return ""


def build_protein_to_locus_map(gbff_path):
    mapping = {}

    for rec in SeqIO.parse(gbff_path, "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue

            q = feat.qualifiers
            protein_ids = q.get("protein_id", [])
            locus_tags  = q.get("locus_tag", [])

            if not protein_ids or not locus_tags:
                continue

            pid = protein_ids[0].strip()
            lt  = locus_tags[0].strip()

            if pid and lt and pid not in mapping:
                mapping[pid] = lt

    return mapping

# Main
def main():
    if not os.path.isfile(GBFF_PATH):
        raise SystemExit(f"ERROR: GenBank file not found:\n  {GBFF_PATH}")

    if not os.path.isfile(EXCEL_IN):
        raise SystemExit(f"ERROR: Excel file not found:\n  {EXCEL_IN}")

    print("Parsing GenBank file...")
    pid_to_locus = build_protein_to_locus_map(GBFF_PATH)
    print(f"Indexed protein_id to locus_tag mappings: {len(pid_to_locus):,}")

    xls = pd.ExcelFile(EXCEL_IN)
    sheet_names = xls.sheet_names if SHEETS_TO_PROCESS is None else SHEETS_TO_PROCESS

    updated = {}
    total_rows = 0
    filled = 0
    missing_pid_col_sheets = []

    for sh in sheet_names:
        df = pd.read_excel(EXCEL_IN, sheet_name=sh)

        if PROTEIN_COL not in df.columns:
            missing_pid_col_sheets.append(sh)
            updated[sh] = df
            continue

        locus_vals = []

        for v in df[PROTEIN_COL].tolist():
            pid = extract_first_protein_id(v)
            lt = pid_to_locus.get(pid, "")
            locus_vals.append(lt)

        df[LOCUS_COL] = locus_vals

        total_rows += len(df)
        filled += sum(1 for x in locus_vals if clean_cell(x) != "")

        updated[sh] = df

    out_dir = os.path.dirname(EXCEL_OUT)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with pd.ExcelWriter(EXCEL_OUT, engine="openpyxl") as writer:
        for sh, df in updated.items():
            df.to_excel(writer, sheet_name=sh, index=False)

    print("\nDone.")
    print(f"Input:  {EXCEL_IN}")
    print(f"Output: {EXCEL_OUT}")
    print(f"Processed sheets: {len(sheet_names)}")
    print(f"Total rows scanned: {total_rows:,}")
    print(f"Rows with locus tag filled: {filled:,}")

    if missing_pid_col_sheets:
        print("\nWarning: These sheets did not have the column 'genbank protein id' and were left unchanged:")
        for sh in missing_pid_col_sheets:
            print(f"  - {sh}")

if __name__ == "__main__":
    main()