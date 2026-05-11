"""
Purpose:
    Extract locally available Gene Ontology (GO) annotations from the RefSeq
    GenBank file of the Vibrio natriegens reference genome.

Description:
    This script parses the reference genome GenBank file, identifies CDS features,
    and extracts locally available GO annotations from the following GenBank
    qualifiers:
        - GO_process
        - GO_function
        - GO_component

    The extracted GO annotations are written to an Excel workbook sheet.

Input:
    - Reference genome GenBank file in .gbff format.

Output:
    - Excel workbook containing a GenBank-derived GO annotation sheet.

Required Python packages:
    - pandas
    - biopython
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update INPUT_EXCEL, OUTPUT_EXCEL, and GBFF_PATH before running the script.
"""

import os
import pandas as pd
from Bio import SeqIO
import re

# PATHS
INPUT_EXCEL = "path/to/existing/excel file.xlsx"
OUTPUT_EXCEL = "path/to/outputv/excel.xlsx"
GBFF_PATH = "path/to/reference_genome/genomic.gbff"
SHEET_NAME = "genbank"

# ----------------------------
# Helper functions
# ----------------------------
GO_ID_PATTERN = re.compile(r"(GO:\d{7})")

def extract_go_id(text):
    match = GO_ID_PATTERN.search(text or "")
    return match.group(1) if match else ""

def extract_label(text):
    text = text or ""
    if "[" in text:
        return text.split("[", 1)[0].strip()
    if " - " in text:
        return text.split(" - ", 1)[1].strip()
    return text.strip()

def format_term(label, goid):
    return f"{label} [{goid}]" if label else f"[{goid}]"

def deduplicate_preserve_order(items):
    seen = set()
    result = []
    for x in items:
        if x not in seen:
            seen.add(x)
            result.append(x)
    return result

# ----------------------------
# Parse GenBank file
# ----------------------------
def parse_genbank_structured_go(gbff_path):
    if not os.path.exists(gbff_path):
        raise FileNotFoundError(f"GenBank file not found: {gbff_path}")

    print(f"Parsing GenBank file: {gbff_path}")

    records = []
    cds_counter = 0

    with open(gbff_path, "r") as handle:
        for rec in SeqIO.parse(handle, "genbank"):
            for feat in rec.features:
                if feat.type != "CDS":
                    continue

                cds_counter += 1

                locus_tags = feat.qualifiers.get("locus_tag", [])
                if not locus_tags:
                    continue

                locus = locus_tags[0]
                protein_id = feat.qualifiers.get("protein_id", [""])[0]

                bp_raw = feat.qualifiers.get("GO_process", [])
                mf_raw = feat.qualifiers.get("GO_function", [])
                cc_raw = feat.qualifiers.get("GO_component", [])

                def process_terms(raw_terms):
                    formatted = []
                    for t in raw_terms:
                        goid = extract_go_id(t)
                        label = extract_label(t)
                        if goid:
                            formatted.append(format_term(label, goid))
                    return deduplicate_preserve_order(formatted)

                bp = process_terms(bp_raw)
                mf = process_terms(mf_raw)
                cc = process_terms(cc_raw)

                row = {
                    "locus_tag": locus,
                    "GenBank protein ID": protein_id,
                    "Gene Ontology (biological process)": "; ".join(bp),
                    "source (bp)": "; ".join(["GenBank"] * len(bp)) if bp else "",
                    "evidence (bp)": "",
                    "Gene Ontology (molecular function)": "; ".join(mf),
                    "source (mf)": "; ".join(["GenBank"] * len(mf)) if mf else "",
                    "evidence (mf)": "",
                    "Gene Ontology (cellular component)": "; ".join(cc),
                    "source (cc)": "; ".join(["GenBank"] * len(cc)) if cc else "",
                    "evidence (cc)": "",
                }

                records.append(row)

    print(f"[INFO] CDS features processed: {cds_counter}")
    print(f"[INFO] Rows created: {len(records)}")

    return pd.DataFrame(records).sort_values("locus_tag")

# ----------------------------
# Write to Excel
# ----------------------------
def write_genbank_sheet(df, input_excel, output_excel, sheet_name):
    print(f"Writing sheet '{sheet_name}' to Excel...")

    if os.path.exists(input_excel):
        with pd.ExcelWriter(output_excel, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        with pd.ExcelWriter(output_excel, engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    print("Excel writing complete.")

# ----------------------------
# Main
# ----------------------------
def main():
    print("\n=== Structured GenBank GO Extraction (Reference Genome) ===\n")
    df = parse_genbank_structured_go(GBFF_PATH)
    write_genbank_sheet(df, INPUT_EXCEL, OUTPUT_EXCEL, SHEET_NAME)
    print("\nCompleted successfully.\n")

if __name__ == "__main__":
    main()