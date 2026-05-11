"""
Purpose:
    Insert RefSeq protein IDs into a UniProtKB GO annotation workbook.

Description:
    This script reads a GO annotation Excel workbook containing UniProt identifiers
    and a separate UniProt-to-RefSeq mapping file. It matches UniProt identifiers
    between the two files and adds a RefSeq protein ID column to the GO annotation
    sheet.

Input:
    - GO annotation Excel workbook containing UniProt identifiers.
    - UniProt-to-RefSeq mapping Excel file.

Output:
    - Updated GO annotation Excel workbook containing a RefSeq protein ID column.
    - Text file listing UniProt entries that could not be mapped.

Required Python packages:
    - pandas
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update GO_XLSX, MAP_XLSX, OUTPUT_XLSX, and UNMAPPED_TXT before running.
"""

import pandas as pd

# Paths
GO_XLSX = "path/to/uniprotkb_reference_vnat_genome_GO.xlsx"
MAP_XLSX = "path/to/uniprot_to_refseq_mapping.xlsx"

GO_SHEET_NAME = "GO_annotations"

GO_ID_COL = "Entry Name"
MAP_TO_COL = "To"
OUT_REFSEQ_COL = "RefSeq Protein ID"

OUTPUT_XLSX = "path/to/output/uniprotkb_reference_vnat_genome_GO_updated.xlsx"
UNMAPPED_TXT = "path/to/output/unmapped_uniprot_entry_names.txt"


# Helper functions
def extract_uniprot_key(entry_name):
    """
    Convert A0AAN1CUF2_VIBNA to A0AAN1CUF2.
    If there is no underscore, return the full value.
    """
    if entry_name is None or pd.isna(entry_name):
        return ""

    s = str(entry_name).strip()

    if not s:
        return ""

    return s.split("_", 1)[0]

def pick_from_column(map_df):
    """
    Select the column containing the UniProt identifiers in the mapping file.
    """
    candidates = ["From", "from", "Entry", "Entry Name", "UniProt", "UniProtKB"]

    for c in candidates:
        if c in map_df.columns:
            return c

    raise ValueError(f"Could not find a UniProt 'From' column. Columns: {list(map_df.columns)}")

# Main
def main():
    print("Loading GO workbook...")
    go_sheets = pd.read_excel(GO_XLSX, sheet_name=None)
    print(f"Sheets found: {list(go_sheets.keys())}")

    if GO_SHEET_NAME not in go_sheets:
        raise ValueError(f"Sheet '{GO_SHEET_NAME}' not found in GO workbook.")

    go_df = go_sheets[GO_SHEET_NAME]
    print(f"GO rows loaded from '{GO_SHEET_NAME}': {len(go_df)}")

    if GO_ID_COL not in go_df.columns:
        raise ValueError(f"GO sheet missing '{GO_ID_COL}'. Columns: {list(go_df.columns)}")

    print("Loading mapping workbook...")
    map_df = pd.read_excel(MAP_XLSX)
    print(f"Mapping rows loaded: {len(map_df)}")

    if MAP_TO_COL not in map_df.columns:
        raise ValueError(f"Mapping file missing '{MAP_TO_COL}'. Columns: {list(map_df.columns)}")

    map_from_col = pick_from_column(map_df)
    print(f"Using mapping source column: {map_from_col}")

    go_df["_uniprot_key"] = go_df[GO_ID_COL].apply(extract_uniprot_key)
    map_df["_uniprot_key"] = map_df[map_from_col].apply(extract_uniprot_key)

    map_df = map_df.dropna(subset=["_uniprot_key", MAP_TO_COL])
    map_df = map_df[map_df["_uniprot_key"].astype(str).str.len() > 0]

    mapping = (
        map_df.drop_duplicates(subset=["_uniprot_key"], keep="first")
              .set_index("_uniprot_key")[MAP_TO_COL]
              .to_dict()
    )

    print(f"Unique UniProt keys in mapping: {len(mapping)}")

    go_df[OUT_REFSEQ_COL] = go_df["_uniprot_key"].map(mapping).fillna("")

    matched_mask = go_df[OUT_REFSEQ_COL].astype(str).str.len() > 0
    matched = int(matched_mask.sum())
    total = len(go_df)
    unmapped = total - matched

    print(f"Matched RefSeq IDs for {matched}/{total} rows.")
    print(f"Unmapped rows: {unmapped}")

    unmapped_entries = go_df.loc[~matched_mask, GO_ID_COL].astype(str).tolist()

    with open(UNMAPPED_TXT, "w", encoding="utf-8") as f:
        for x in unmapped_entries:
            f.write(x + "\n")

    print(f"Wrote unmapped Entry Names to: {UNMAPPED_TXT}")

    go_df = go_df.drop(columns=["_uniprot_key"])

    go_sheets[GO_SHEET_NAME] = go_df

    print(f"Writing updated workbook: {OUTPUT_XLSX}")

    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
        for sheet_name, df in go_sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    print("Done")

if __name__ == "__main__":
    main()