"""
Purpose:
    Fetch and merge Gene Ontology (GO) annotations for genes using local GenBank
    annotations and UniProtKB-derived annotations.

Description:
    This script extracts locally available GO terms from RefSeq GenBank files using
    strain names and locus tags. It also maps RefSeq or GenBank protein identifiers
    to UniProtKB accessions using the UniProt ID-mapping REST API, retrieves
    UniProtKB GO annotations, and merges GenBank-derived and UniProtKB-derived GO
    terms per gene.

    The merged output removes duplicate GO terms based on GO ID and records the
    source of each term.

Input:
    - Excel workbook containing gene identifiers.
    - RefSeq GenBank folders containing genomic.gbff files for each strain.

Output:
    - Excel workbook with three sheets:
        1. GenBank_GO
        2. UniProt_GO
        3. Merged_GO
    - Missing/mapping log file.
    - Local JSON cache files for UniProt mapping and GO retrieval.

Required Python packages:
    - pandas
    - requests
    - biopython
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update EXCEL_FILE, OUTPUT_FILE, GENBANK_ROOT, MISSING_LOG_FILE, and
    CACHE_DIR before running the script.
"""

import os
import re
import json
import time
from io import StringIO
import pandas as pd
import requests
from Bio import SeqIO

#file paths
EXCEL_FILE  = "path/to/input_GO_annotation_workbook.xlsx"
INPUT_SHEET = "Uniprot_GO"
OUTPUT_FILE = "path/to/output_GO_annotation_workbook_FETCHED.xlsx"
GENBANK_ROOT = "path/to/refseq_genbank_folders"
GBFF_NAME    = "genomic.gbff"
MISSING_LOG_FILE = "path/to/output/missing_GO_or_uniprot_mapping.log"
CACHE_DIR = "path/to/output/.cache"
MAP_CACHE_JSON = os.path.join(CACHE_DIR, "filename.json")
GO_CACHE_JSON  = os.path.join(CACHE_DIR, "filename.json")

UNIPROT_IDMAP_URL    = "https://rest.uniprot.org/idmapping/run"
UNIPROT_STATUS_URL   = "https://rest.uniprot.org/idmapping/status/"
UNIPROT_RESULT_URL   = "https://rest.uniprot.org/idmapping/results/"
UNIPROT_STREAM_URL   = "https://rest.uniprot.org/uniprotkb/stream"

BATCH_SIZE  = 200
SLEEP_SEC   = 0.35

COL_STRAIN = "strain"
COL_LOCUS  = "locus tag"
COL_REFSEQ = "refseq protein id"
COL_GB_PID = "genbank protein id"
COL_UNI_ACC = "uniprotkb id"

GO_COLS = {
    "bp": "Gene Ontology (biological process)",
    "mf": "Gene Ontology (molecular function)",
    "cc": "Gene Ontology (cellular component)",
}
SRC_COLS = {"bp": "source (bp)", "mf": "source (mf)", "cc": "source (cc)"}
EV_COLS  = {"bp": "evidence (bp)", "mf": "evidence (mf)", "cc": "evidence (cc)"}

GO_ID_RE = re.compile(r"(GO:\d{7})")
GO_ID_FROM_NOTE_RE = re.compile(r"GO:(\d{7})")
UNIPROT_GO_ID_IN_BRACKETS = re.compile(r"\[(GO:\d{7})\]")

def ensure_dir(p: str):
    if p:
        os.makedirs(p, exist_ok=True)

def load_json(path: str) -> dict:
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return {}
    return {}

def save_json(path: str, obj: dict):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)
    os.replace(tmp, path)

def clean(x) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and pd.isna(x):
        return ""
    s = str(x).strip()
    return "" if s.lower() in {"nan", "none"} else s

def pick_sheet(xlsx: str, preferred: str | None) -> str:
    xls = pd.ExcelFile(xlsx)
    if preferred and preferred in xls.sheet_names:
        return preferred
    if "GO_annotations" in xls.sheet_names:
        return "GO_annotations"
    return xls.sheet_names[0]

def as_string_cols(df: pd.DataFrame, cols: list[str]):
    for c in cols:
        if c in df.columns:
            df[c] = df[c].astype("string")
        else:
            df[c] = pd.Series([""] * len(df), dtype="string")

def require_input_cols(df: pd.DataFrame):
    missing = [c for c in [COL_STRAIN, COL_LOCUS, COL_REFSEQ, COL_GB_PID] if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: missing required columns: {missing}")

def clear_go_columns(dfx: pd.DataFrame):
    for aspect in ("bp", "mf", "cc"):
        dfx[GO_COLS[aspect]]  = ""
        dfx[SRC_COLS[aspect]] = ""
        dfx[EV_COLS[aspect]]  = ""


# RefSeq record parsing
def extract_go_from_note(note: str) -> list[str]:
    return [f"GO:{m}" for m in GO_ID_FROM_NOTE_RE.findall(note or "")]

def parse_gbff_locus_index(gbff_path: str) -> dict:
    locus_index = {}
    for rec in SeqIO.parse(gbff_path, "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            if "locus_tag" not in q:
                continue
            locus_tag = q["locus_tag"][0].strip()

            bp = list(q.get("GO_process", []))
            mf = list(q.get("GO_function", []))
            cc = list(q.get("GO_component", []))

            for n in q.get("note", []):
                bp.extend(extract_go_from_note(n))

            locus_index[locus_tag] = {"bp": bp, "mf": mf, "cc": cc}
    return locus_index


# UniProt mapping
def uniprot_map_refseq_to_uniprot(refseq_ids: list[str], map_cache: dict) -> dict:
    out = {}
    to_map = [r for r in refseq_ids if r and r not in map_cache]
    if not to_map:
        return out

    for i in range(0, len(to_map), BATCH_SIZE):
        batch = to_map[i:i + BATCH_SIZE]
        try:
            job = requests.post(
                UNIPROT_IDMAP_URL,
                data={"from": "RefSeq_Protein", "to": "UniProtKB", "ids": ",".join(batch)},
                timeout=30,
            )
            job.raise_for_status()
            job_id = job.json().get("jobId")
        except Exception:
            continue

        ok = False
        for _ in range(40):
            time.sleep(SLEEP_SEC)
            try:
                r = requests.get(UNIPROT_STATUS_URL + job_id, timeout=30)
                if r.status_code != 200:
                    continue
                st = r.json().get("jobStatus")
                if st == "FINISHED":
                    ok = True
                    break
                if st == "FAILED":
                    break
            except Exception:
                continue
        if not ok:
            continue

        try:
            tsv = requests.get(UNIPROT_RESULT_URL + job_id + "?format=tsv", timeout=30).text
            lines = tsv.strip().split("\n")
            if len(lines) <= 1:
                continue
            for line in lines[1:]:
                cols = line.split("\t")
                if len(cols) >= 2:
                    ref, uni = cols[0].strip(), cols[1].strip()
                    if uni and not uni.startswith("UPI"):
                        out[ref] = uni
        except Exception:
            continue

    map_cache.update(out)
    return out


# Fetching GO terms from UniProt
def _split_semicolon_list(series) -> list[str]:
    out = []
    for v in series.dropna().astype(str).tolist():
        out.extend([x.strip() for x in v.split(";") if x.strip()])
    return out

def fetch_uniprot_go(uniprot_acc: str, go_cache: dict):
    if not uniprot_acc:
        return {"bp": ([], []), "mf": ([], []), "cc": ([], [])}

    if uniprot_acc in go_cache:
        return go_cache[uniprot_acc]

    url = (
        f"{UNIPROT_STREAM_URL}?query=accession:{uniprot_acc}"
        f"&fields=go_p,go_f,go_c,go_evidence&format=tsv"
    )

    try:
        r = requests.get(url, timeout=30)
        if r.status_code != 200 or not r.text.strip():
            out = {"bp": ([], []), "mf": ([], []), "cc": ([], [])}
            go_cache[uniprot_acc] = out
            return out

        df = pd.read_csv(StringIO(r.text), sep="\t")
        col_bp = GO_COLS["bp"]
        col_mf = GO_COLS["mf"]
        col_cc = GO_COLS["cc"]
        col_ev = "Gene Ontology (evidence)"

        bp = _split_semicolon_list(df[col_bp]) if col_bp in df.columns else []
        mf = _split_semicolon_list(df[col_mf]) if col_mf in df.columns else []
        cc = _split_semicolon_list(df[col_cc]) if col_cc in df.columns else []
        ev = _split_semicolon_list(df[col_ev]) if col_ev in df.columns else []

        out = {"bp": (bp, ev), "mf": (mf, ev), "cc": (cc, ev)}
        go_cache[uniprot_acc] = out
        return out

    except Exception:
        out = {"bp": ([], []), "mf": ([], []), "cc": ([], [])}
        go_cache[uniprot_acc] = out
        return out

# Merge and deduplicate by GO ID
def go_id_from_any(term: str) -> str:
    term = term or ""
    m = UNIPROT_GO_ID_IN_BRACKETS.search(term)
    if m:
        return m.group(1)
    m = GO_ID_RE.search(term)
    return m.group(1) if m else ""

def label_from_any(term: str) -> str:
    t = (term or "").strip()
    if "[" in t and "]" in t and "GO:" in t:
        return t.split("[", 1)[0].strip()
    if t.startswith("GO:") and " - " in t:
        return t.split(" - ", 1)[1].strip()
    if t.startswith("GO:"):
        return ""
    return t

def merge_by_go_id(gb_terms: list[str], uni_terms: list[str]):
    merged = {}  # goid -> {"label": str, "sources": set()}

    for t in gb_terms or []:
        goid = go_id_from_any(t)
        if not goid:
            continue
        merged.setdefault(goid, {"label": "", "sources": set()})
        merged[goid]["sources"].add("GenBank")
        if not merged[goid]["label"]:
            merged[goid]["label"] = label_from_any(t)

    for t in uni_terms or []:
        goid = go_id_from_any(t)
        if not goid:
            continue
        merged.setdefault(goid, {"label": "", "sources": set()})
        merged[goid]["sources"].add("UniProtKB")
        ulab = label_from_any(t)
        if ulab:
            merged[goid]["label"] = ulab

    goids = sorted(merged.keys())
    items = []
    src_items = []
    for goid in goids:
        lab = merged[goid]["label"]
        items.append(f"{lab} [{goid}]" if lab else f"[{goid}]")
        srcs = sorted(merged[goid]["sources"], key=lambda x: 0 if x == "GenBank" else 1)
        src_items.append(", ".join(srcs))

    return "; ".join(items), "; ".join(src_items), merged

def main():
    ensure_dir(CACHE_DIR)
    map_cache = load_json(MAP_CACHE_JSON)
    go_cache  = load_json(GO_CACHE_JSON)

    sheet = pick_sheet(EXCEL_FILE, INPUT_SHEET)
    df_in = pd.read_excel(EXCEL_FILE, sheet_name=sheet)
    require_input_cols(df_in)

    as_string_cols(df_in, [COL_STRAIN, COL_LOCUS, COL_REFSEQ, COL_GB_PID, COL_UNI_ACC])
    as_string_cols(df_in, [*GO_COLS.values(), *SRC_COLS.values(), *EV_COLS.values()])

    df_gb  = df_in.copy()
    df_uni = df_in.copy()
    df_mrg = df_in.copy()

    clear_go_columns(df_gb)
    clear_go_columns(df_uni)
    clear_go_columns(df_mrg)
    df_uni[COL_UNI_ACC] = ""

    strain_index_cache = {}

    query_ids = []
    row_query_id = []
    for i, row in df_in.iterrows():
        refseq = clean(row[COL_REFSEQ])
        gbpid  = clean(row[COL_GB_PID])
        qid = refseq if refseq else gbpid
        row_query_id.append((i, qid))
        if qid:
            query_ids.append(qid)

    uniprot_map_refseq_to_uniprot(sorted(set(query_ids)), map_cache)
    save_json(MAP_CACHE_JSON, map_cache)

    missing_rows = []

    for i, qid in row_query_id:
        strain = clean(df_in.at[i, COL_STRAIN])
        locus  = clean(df_in.at[i, COL_LOCUS])

        gbff_path = os.path.join(GENBANK_ROOT, strain, GBFF_NAME)
        if strain not in strain_index_cache:
            if os.path.isfile(gbff_path):
                try:
                    strain_index_cache[strain] = parse_gbff_locus_index(gbff_path)
                except Exception:
                    strain_index_cache[strain] = None
            else:
                strain_index_cache[strain] = None

        locus_index = strain_index_cache.get(strain)
        if locus_index is None:
            gb = {"bp": [], "mf": [], "cc": []}
            missing_rows.append((strain, locus, qid or "", "", "missing_or_unreadable_genomic.gbff"))
        else:
            gb = locus_index.get(locus)
            if gb is None:
                gb = {"bp": [], "mf": [], "cc": []}
                missing_rows.append((strain, locus, qid or "", "", "locus_tag_not_found_in_gbff"))

        for aspect in ("bp", "mf", "cc"):
            gb_terms = gb.get(aspect, [])
            if gb_terms:
                df_gb.at[i, GO_COLS[aspect]] = "; ".join(gb_terms)
                df_gb.at[i, SRC_COLS[aspect]] = "; ".join(["GenBank"] * len(gb_terms))
                df_gb.at[i, EV_COLS[aspect]] = ""

        uni_acc = map_cache.get(qid, "") if qid else ""
        if uni_acc:
            uni = fetch_uniprot_go(uni_acc, go_cache)
        else:
            uni = {"bp": ([], []), "mf": ([], []), "cc": ([], [])}
            if qid:
                missing_rows.append((strain, locus, qid, "", "no_uniprot_mapping_for_query_id"))

        if uni_acc:
            df_uni.at[i, COL_UNI_ACC] = uni_acc

        for aspect in ("bp", "mf", "cc"):
            terms, ev = uni[aspect]
            if terms:
                df_uni.at[i, GO_COLS[aspect]] = "; ".join(terms)
                df_uni.at[i, SRC_COLS[aspect]] = "; ".join(["UniProtKB"] * len(terms))
                ev_unique = sorted({clean(x) for x in ev if clean(x)})
                df_uni.at[i, EV_COLS[aspect]] = "; ".join(ev_unique)

        for aspect in ("bp", "mf", "cc"):
            gb_terms = gb.get(aspect, [])
            uni_terms = uni[aspect][0]
            merged_terms_str, merged_sources_str, merged_dict = merge_by_go_id(gb_terms, uni_terms)

            if merged_dict:
                df_mrg.at[i, GO_COLS[aspect]] = merged_terms_str
                df_mrg.at[i, SRC_COLS[aspect]] = merged_sources_str

                if any("UniProtKB" in v["sources"] for v in merged_dict.values()):
                    ev_unique = sorted({clean(x) for x in uni[aspect][1] if clean(x)})
                    df_mrg.at[i, EV_COLS[aspect]] = "; ".join(ev_unique)

    save_json(GO_CACHE_JSON, go_cache)

    ensure_dir(os.path.dirname(OUTPUT_FILE))
    with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
        df_gb.to_excel(writer, sheet_name="GenBank_GO", index=False)
        df_uni.to_excel(writer, sheet_name="UniProt_GO", index=False)
        df_mrg.to_excel(writer, sheet_name="Merged_GO", index=False)

    ensure_dir(os.path.dirname(MISSING_LOG_FILE))
    with open(MISSING_LOG_FILE, "w", encoding="utf-8") as f:
        f.write("strain\tlocus_tag\tquery_protein_id\tuniprotkb_id\treason\n")
        for s, lt, qid, uni, reason in missing_rows:
            f.write(f"{clean(s)}\t{clean(lt)}\t{clean(qid)}\t{clean(uni)}\t{clean(reason)}\n")

    print("Done (GO fetching only).")
    print(f"Read:  {EXCEL_FILE} (sheet: {sheet})")
    print(f"Wrote: {OUTPUT_FILE} (GenBank_GO, UniProt_GO, Merged_GO)")
    print(f"Missing log: {MISSING_LOG_FILE}")

if __name__ == "__main__":
    main()