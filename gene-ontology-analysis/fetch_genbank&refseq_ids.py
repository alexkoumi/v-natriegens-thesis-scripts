"""
Purpose:
    Add GenBank protein IDs and RefSeq protein IDs to a GO annotation workbook.

Description:
    This script reads an Excel workbook containing strain names and locus tags.
    For each row, it opens the corresponding strain GenBank file, finds the CDS
    feature with the matching locus tag, and extracts:
        - GenBank protein ID from the protein_id qualifier
        - RefSeq protein ID from the inference qualifier

Input:
    - Excel workbook containing strain and locus tag columns.
    - GenBank files arranged in strain-specific folders.

Output:
    - Updated Excel workbook containing GenBank protein IDs and RefSeq protein IDs.
    - Log file listing missing GenBank files, missing locus tags, or missing RefSeq IDs.

Required Python packages:
    - pandas
    - biopython
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update EXCEL_FILE, OUTPUT_FILE, BASE_DIR, and LOG_FILE before running.
"""

import os
import re
import pandas as pd
from Bio import SeqIO

# Paths
EXCEL_FILE = "path/to/Excel/workbook/with/strain/locus tag.xlsx"
OUTPUT_FILE = "path/to/Excel/workbook/with/strain/locus tagupdated.xlsx"
BASE_DIR = "path/to/refseq_records_folders"
LOG_FILE =  "path/to/output/missing_refseq_ids.log"

# If None: script will auto-pick GO_annotations if present, otherwise first sheet.
SHEET_NAME = None  # or "GO_annotations"
REFSEQ_REGEX = re.compile(r"RefSeq:([A-Za-z0-9_.]+)")

#helper functions
def pick_sheet(xlsx_path: str, preferred: str | None) -> str:
    xls = pd.ExcelFile(xlsx_path)
    sheets = xls.sheet_names

    if preferred and preferred in sheets:
        return preferred

    if "GO_annotations" in sheets:
        return "GO_annotations"

    # fallback: first sheet
    return sheets[0]

def build_locus_map(gbff_path: str) -> dict:
    locus_map = {}
    for record in SeqIO.parse(gbff_path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue

            q = feature.qualifiers
            locus_tags = q.get("locus_tag", [])
            if not locus_tags:
                continue

            locus_tag = locus_tags[0]
            protein_id = q.get("protein_id", [None])[0]

            refseq_id = None
            for inf in q.get("inference", []):
                m = REFSEQ_REGEX.search(inf)
                if m:
                    refseq_id = m.group(1)
                    break

            locus_map[locus_tag] = {"protein_id": protein_id, "refseq_id": refseq_id}

    return locus_map

# Main
def main():
    if not os.path.isfile(EXCEL_FILE):
        raise SystemExit(f"ERROR: Excel file not found: {EXCEL_FILE}")

    sheet = pick_sheet(EXCEL_FILE, SHEET_NAME)
    print(f"Using sheet: {sheet}")

    df = pd.read_excel(EXCEL_FILE, sheet_name=sheet)

    # Make sure these columns are treated as strings (prevents dtype warnings)
    for col in ["genbank protein id", "refseq protein id"]:
        if col in df.columns:
            df[col] = df[col].astype("string")

    required = ["strain", "locus tag", "genbank protein id", "refseq protein id"]
    for c in required:
        if c not in df.columns:
            raise SystemExit(f"ERROR: Missing column '{c}' in sheet '{sheet}'")

    strain_cache = {}
    missing_refseq = []

    for i, row in df.iterrows():
        if pd.isna(row["strain"]) or pd.isna(row["locus tag"]):
            continue

        strain = str(row["strain"]).strip()
        locus_tag = str(row["locus tag"]).strip()

        gbff_path = os.path.join(BASE_DIR, strain, "genomic.gbff")
        if not os.path.isfile(gbff_path):
            missing_refseq.append((strain, locus_tag, "", "missing_genomic.gbff"))
            continue

        if strain not in strain_cache:
            strain_cache[strain] = build_locus_map(gbff_path)

        hits = strain_cache[strain].get(locus_tag)
        if not hits:
            missing_refseq.append((strain, locus_tag, "", "locus_not_found"))
            continue

        # update genbank protein id
        if hits["protein_id"]:
            cur = df.at[i, "genbank protein id"]
            if pd.isna(cur) or str(cur).strip() != hits["protein_id"]:
                df.at[i, "genbank protein id"] = hits["protein_id"]

        # update refseq protein id
        if hits["refseq_id"]:
            cur = df.at[i, "refseq protein id"]
            if pd.isna(cur) or str(cur).strip() != hits["refseq_id"]:
                df.at[i, "refseq protein id"] = hits["refseq_id"]
        else:
            missing_refseq.append((strain, locus_tag, hits["protein_id"] or "", "no_refseq_found"))

    # Save output
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    df.to_excel(OUTPUT_FILE, sheet_name=sheet, index=False)

    # Write log
    os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)
    with open(LOG_FILE, "w", encoding="utf-8") as f:
        f.write("strain\tlocus_tag\tgenbank_protein_id\treason\n")
        for s, lt, pid, reason in missing_refseq:
            f.write(f"{s}\t{lt}\t{pid}\t{reason}\n")

    print("Done.")
    print("Output:", OUTPUT_FILE)
    print("Log:", LOG_FILE)


if __name__ == "__main__":
    main()
