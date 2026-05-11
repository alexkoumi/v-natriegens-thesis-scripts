"""
Purpose:
    Summarise Gene Ontology (GO) annotation coverage, depth, source contribution,
    aspect composition, and most frequent GO terms.

Description:
    This script reads an Excel workbook containing GenBank-derived and UniProtKB-
    derived GO annotation sheets. It calculates summary statistics for mappable
    genes with locus tags, including:
        - number of genes with at least one GO term
        - number of genes with no GO terms
        - annotation depth per gene
        - GenBank and UniProtKB source contribution
        - BP/MF/CC aspect composition
        - top GO terms per aspect

    It also creates an annotation-depth Excel file and can summarise UniProtKB
    entries without locus tags if present.

Input:
    - Excel workbook containing GO annotation sheets.

Output:
    - Text summary report.
    - Excel file containing per-gene GO annotation depth.

Required Python packages:
    - pandas
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update EXCEL_FILE, SUMMARY_LOG_FILE, and DEPTH_XLSX_OUT before running.
"""

import os
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Paths
EXCEL_FILE = "path/to/GO_annotation_workbook.xlsx"
SUMMARY_LOG_FILE = "path/to/output/GO_annotation_summary.txt"
DEPTH_XLSX_OUT = "path/to/output/GO_annotation_depth.xlsx"
TOP_N_TERMS = 10

SHEET_GENBANK = "GenBank_GO"
SHEET_UNIPROT = "UniProt_GO"

COL_STRAIN = "strain"
COL_LOCUS = "locus tag"

GO_COLS = {
    "bp": "gene ontology (biological process)",
    "mf": "gene ontology (molecular function)",
    "cc": "gene ontology (cellular component)",
}

UNIPROT_ID_CANDIDATES = [
    "uniprotkb id",
    "uniprotkb_id",
    "uniprotkb accession",
    "uniprot id",
    "uniprot accession",
]


# Helper functions
def clean(x):
    if x is None:
        return ""
    if isinstance(x, float) and pd.isna(x):
        return ""
    s = str(x).strip()
    return "" if s.lower() in {"nan", "none"} else s


def split_semicolon(s):
    s = clean(s)
    if not s:
        return []
    return [p.strip() for p in s.split(";") if p.strip()]


def safe_percent(n, d):
    return "0.00%" if d == 0 else f"{(n / d) * 100:.2f}%"


def compute_depth_stats(depths):
    if not depths:
        return {"mean": 0.0, "median": 0.0, "sd": 0.0, "q1": 0.0, "q3": 0.0, "iqr": 0.0}

    s = pd.Series(depths, dtype=float)
    q1 = float(s.quantile(0.25))
    q3 = float(s.quantile(0.75))
    sd = float(s.std(ddof=1)) if len(s) > 1 else 0.0

    return {
        "mean": float(s.mean()),
        "median": float(s.median()),
        "sd": sd,
        "q1": q1,
        "q3": q3,
        "iqr": float(q3 - q1),
    }


def compute_aspect_combination_counts(all_genes, combined_aspect_terms):
    counts = Counter({
        "bp_only": 0,
        "mf_only": 0,
        "cc_only": 0,
        "bp_mf_only": 0,
        "bp_cc_only": 0,
        "mf_cc_only": 0,
        "bp_mf_cc": 0,
    })

    for gk in all_genes:
        has_bp = len(combined_aspect_terms["bp"].get(gk, set())) > 0
        has_mf = len(combined_aspect_terms["mf"].get(gk, set())) > 0
        has_cc = len(combined_aspect_terms["cc"].get(gk, set())) > 0

        if has_bp and not has_mf and not has_cc:
            counts["bp_only"] += 1
        elif has_mf and not has_bp and not has_cc:
            counts["mf_only"] += 1
        elif has_cc and not has_bp and not has_mf:
            counts["cc_only"] += 1
        elif has_bp and has_mf and not has_cc:
            counts["bp_mf_only"] += 1
        elif has_bp and has_cc and not has_mf:
            counts["bp_cc_only"] += 1
        elif has_mf and has_cc and not has_bp:
            counts["mf_cc_only"] += 1
        elif has_bp and has_mf and has_cc:
            counts["bp_mf_cc"] += 1

    return counts


def make_column_names_lowercase(df):
    if df is not None:
        df.columns = [str(c).lower() for c in df.columns]


def ensure_cols(df):
    if df is None:
        return

    for c in [COL_STRAIN, COL_LOCUS, *GO_COLS.values()]:
        if c not in df.columns:
            df[c] = ""


def gene_key_from_row(row):
    return (clean(row.get(COL_STRAIN, "")), clean(row.get(COL_LOCUS, "")))


def pick_uniprot_id_column(df):
    cols_lower = {c.lower(): c for c in df.columns}

    for cand in UNIPROT_ID_CANDIDATES:
        if cand in cols_lower:
            return cols_lower[cand]

    return None


def write_summary_log(out_path, summary, top_n):
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("=== GO ANNOTATION SUMMARY ===\n\n")

        f.write("1a) Annotation coverage\n")
        f.write(f"Total genes processed: {summary['total_genes']}\n")
        f.write(
            f"Genes with ≥1 GO term (combined): {summary['genes_with_any_go']} "
            f"({safe_percent(summary['genes_with_any_go'], summary['total_genes'])})\n"
        )
        f.write(
            f"Genes with 0 GO terms (combined): {summary['genes_with_no_go']} "
            f"({safe_percent(summary['genes_with_no_go'], summary['total_genes'])})\n\n"
        )

        f.write("1b) Annotation depth (GO terms per annotated gene)\n")
        d = summary["depth"]
        f.write(f"  Mean: {d['mean']:.3f}\n")
        f.write(f"  Median: {d['median']:.3f}\n")
        f.write(f"  SD: {d['sd']:.3f}\n")
        f.write(f"  Q1: {d['q1']:.3f}\n")
        f.write(f"  Q3: {d['q3']:.3f}\n")
        f.write(f"  IQR: {d['iqr']:.3f}\n\n")

        f.write("2a) Source contribution (gene-level; any aspect)\n")
        f.write("  Total\n")
        for src, cnt in summary["gene_src_total"].most_common():
            f.write(f"     - {src}: {cnt}\n")

        for aspect in ("bp", "mf", "cc"):
            f.write(f"     {aspect.upper()}:\n")
            for src, cnt in summary["gene_src_by_aspect"][aspect].most_common():
                f.write(f"             - {src}: {cnt} ({safe_percent(cnt, summary['total_genes'])})\n")
        f.write("\n")

        f.write("2a) Overlap of sources (gene-level; any aspect)\n")
        f.write(f"  Genes with ≥1 UniProtKB term: {summary['genes_with_uniprot_any']}\n")
        f.write(f"  Genes with ≥1 GenBank term: {summary['genes_with_genbank_any']}\n")
        f.write(f"  Genes with both sources: {summary['genes_with_both_sources_any']}\n")
        f.write(f"  UniProtKB-only genes: {summary['genes_with_uniprot_only_any']}\n")
        f.write(f"  GenBank-only genes: {summary['genes_with_genbank_only_any']}\n\n")

        f.write("2b) Source contribution (term-level; GO term entries from source)\n")
        f.write("  Total (all aspects):\n")
        for src, cnt in summary["term_src_total"].most_common():
            f.write(f"    - {src}: {cnt}\n")

        for aspect in ("bp", "mf", "cc"):
            f.write(f"  {aspect.upper()}:\n")
            for src, cnt in summary["term_src_by_aspect"][aspect].most_common():
                f.write(f"    - {src}: {cnt}\n")
        f.write("\n\n")

        f.write("2c) GO aspect composition (gene-level; combined)\n")
        acc = summary["aspect_combo_counts"]
        f.write(f"  BP only: {acc['bp_only']} ({safe_percent(acc['bp_only'], summary['total_genes'])})\n")
        f.write(f"  MF only: {acc['mf_only']} ({safe_percent(acc['mf_only'], summary['total_genes'])})\n")
        f.write(f"  CC only: {acc['cc_only']} ({safe_percent(acc['cc_only'], summary['total_genes'])})\n")
        f.write(f"  BP + MF only: {acc['bp_mf_only']} ({safe_percent(acc['bp_mf_only'], summary['total_genes'])})\n")
        f.write(f"  BP + CC only: {acc['bp_cc_only']} ({safe_percent(acc['bp_cc_only'], summary['total_genes'])})\n")
        f.write(f"  MF + CC only: {acc['mf_cc_only']} ({safe_percent(acc['mf_cc_only'], summary['total_genes'])})\n")
        f.write(f"  BP + MF + CC: {acc['bp_mf_cc']} ({safe_percent(acc['bp_mf_cc'], summary['total_genes'])})\n\n")

        f.write(f"Top {top_n} GO terms per aspect (combined; count = genes with term)\n")
        for aspect in ("bp", "mf", "cc"):
            f.write(f"  {aspect.upper()}:\n")
            top = summary["top_terms"][aspect].most_common(top_n)

            if not top:
                f.write("    (none)\n")
            else:
                for term, cnt in top:
                    f.write(f"    - {term}: {cnt}\n")

        extra = summary.get("unmappable_uniprot", None)

        if extra is None:
            return

        f.write("\n\n")
        f.write("UniProt entries without locus tags (unmappable)\n")
        f.write(f"Total UniProt entries (unmappable): {extra['total_entries']}\n")
        f.write(
            f"Entries with ≥1 GO term (any aspect): {extra['entries_with_any_go']} "
            f"({safe_percent(extra['entries_with_any_go'], extra['total_entries'])})\n"
        )
        f.write(
            f"Entries with 0 GO terms (all aspects): {extra['entries_with_no_go']} "
            f"({safe_percent(extra['entries_with_no_go'], extra['total_entries'])})\n\n"
        )

        f.write("Entries with ≥1 GO term per aspect\n")
        for aspect in ("bp", "mf", "cc"):
            f.write(
                f"  {aspect.upper()}: {extra['entries_with_aspect'][aspect]} "
                f"({safe_percent(extra['entries_with_aspect'][aspect], extra['total_entries'])})\n"
            )
        f.write("\n")

        f.write("GO term entries per aspect (term-level; deduplicated within entry)\n")
        for aspect in ("bp", "mf", "cc"):
            f.write(f"  {aspect.upper()}: {extra['term_entries_by_aspect'][aspect]}\n")
        f.write("\n")

        f.write(f"Top {top_n} GO terms per aspect (unmappable UniProt; count = entries with term)\n")
        for aspect in ("bp", "mf", "cc"):
            f.write(f"  {aspect.upper()}:\n")
            top = extra["top_terms"][aspect].most_common(top_n)

            if not top:
                f.write("    (none)\n")
            else:
                for term, cnt in top:
                    f.write(f"    - {term}: {cnt}\n")


# Main
def main():
    xls = pd.ExcelFile(EXCEL_FILE)

    df_gb = pd.read_excel(EXCEL_FILE, sheet_name=SHEET_GENBANK) if SHEET_GENBANK in xls.sheet_names else None
    df_up = pd.read_excel(EXCEL_FILE, sheet_name=SHEET_UNIPROT) if SHEET_UNIPROT in xls.sheet_names else None

    if df_gb is None and df_up is None:
        df_gb = pd.read_excel(EXCEL_FILE, sheet_name=xls.sheet_names[0])

    make_column_names_lowercase(df_gb)
    make_column_names_lowercase(df_up)

    ensure_cols(df_gb)
    ensure_cols(df_up)

    gb_gene_any = set()
    up_gene_any = set()
    gb_gene_aspect = {"bp": set(), "mf": set(), "cc": set()}
    up_gene_aspect = {"bp": set(), "mf": set(), "cc": set()}

    term_src_total = Counter()
    term_src_by_aspect = {"bp": Counter(), "mf": Counter(), "cc": Counter()}

    top_terms = {"bp": Counter(), "mf": Counter(), "cc": Counter()}

    gb_terms_union = {}
    up_terms_union = {}
    combined_union = {}
    combined_aspect_terms = {"bp": {}, "mf": {}, "cc": {}}

    def add_gene_terms(term_dict, gk, terms):
        if not terms:
            return
        term_dict.setdefault(gk, set()).update(terms)

    def process_mappable(df, source):
        if df is None:
            return

        for _, row in df.iterrows():
            gk = gene_key_from_row(row)

            if gk[1] == "":
                continue

            has_any = False

            for aspect in ("bp", "mf", "cc"):
                terms = split_semicolon(row.get(GO_COLS[aspect], ""))

                if not terms:
                    continue

                has_any = True

                if source == "GenBank":
                    gb_gene_aspect[aspect].add(gk)
                else:
                    up_gene_aspect[aspect].add(gk)

                combined_aspect_terms[aspect].setdefault(gk, set()).update(terms)

                for t in terms:
                    term_src_by_aspect[aspect][source] += 1
                    term_src_total[source] += 1

                for t in set(terms):
                    top_terms[aspect][t] += 1

                if source == "GenBank":
                    add_gene_terms(gb_terms_union, gk, terms)
                else:
                    add_gene_terms(up_terms_union, gk, terms)

            if has_any:
                if source == "GenBank":
                    gb_gene_any.add(gk)
                else:
                    up_gene_any.add(gk)

    process_mappable(df_gb, "GenBank")
    process_mappable(df_up, "UniProtKB")

    all_genes = set()

    if df_gb is not None:
        for _, r in df_gb[[COL_STRAIN, COL_LOCUS]].iterrows():
            if clean(r[COL_LOCUS]) != "":
                all_genes.add((clean(r[COL_STRAIN]), clean(r[COL_LOCUS])))

    if df_up is not None:
        for _, r in df_up[[COL_STRAIN, COL_LOCUS]].iterrows():
            if clean(r[COL_LOCUS]) != "":
                all_genes.add((clean(r[COL_STRAIN]), clean(r[COL_LOCUS])))

    for gk in all_genes:
        gbset = gb_terms_union.get(gk, set())
        upset = up_terms_union.get(gk, set())
        combined_union[gk] = set(gbset) | set(upset)

    aspect_combo_counts = compute_aspect_combination_counts(all_genes, combined_aspect_terms)

    total_genes = len(all_genes)
    depth_list = []
    genes_with_any_go = 0

    for gk in all_genes:
        n = len(combined_union.get(gk, set()))

        if n > 0:
            genes_with_any_go += 1
            depth_list.append(n)

    gene_src_total = Counter({
        "GenBank": len(gb_gene_any),
        "UniProtKB": len(up_gene_any)
    })

    gene_src_by_aspect = {
        "bp": Counter({"GenBank": len(gb_gene_aspect["bp"]), "UniProtKB": len(up_gene_aspect["bp"])}),
        "mf": Counter({"GenBank": len(gb_gene_aspect["mf"]), "UniProtKB": len(up_gene_aspect["mf"])}),
        "cc": Counter({"GenBank": len(gb_gene_aspect["cc"]), "UniProtKB": len(up_gene_aspect["cc"])}),
    }

    genes_with_uniprot_any = len(up_gene_any)
    genes_with_genbank_any = len(gb_gene_any)
    genes_with_both_sources_any = len(gb_gene_any & up_gene_any)
    genes_with_uniprot_only_any = len(up_gene_any - gb_gene_any)
    genes_with_genbank_only_any = len(gb_gene_any - up_gene_any)

    unmappable = {
        "total_entries": 0,
        "entries_with_any_go": 0,
        "entries_with_no_go": 0,
        "entries_with_aspect": {"bp": 0, "mf": 0, "cc": 0},
        "term_entries_by_aspect": {"bp": 0, "mf": 0, "cc": 0},
        "top_terms": {"bp": Counter(), "mf": Counter(), "cc": Counter()},
    }

    if df_up is not None:
        uni_id_col = pick_uniprot_id_column(df_up)

        entry_terms_any = {}
        entry_terms_by_aspect = {"bp": {}, "mf": {}, "cc": {}}

        for idx, row in df_up.iterrows():
            locus = clean(row.get(COL_LOCUS, ""))

            if locus != "":
                continue

            if uni_id_col:
                entry_key = clean(row.get(uni_id_col, ""))
            else:
                entry_key = ""

            if entry_key == "":
                entry_key = f"ROW_{idx}"

            any_terms_here = set()

            for aspect in ("bp", "mf", "cc"):
                terms = set(split_semicolon(row.get(GO_COLS[aspect], "")))

                if terms:
                    entry_terms_by_aspect[aspect].setdefault(entry_key, set()).update(terms)
                    any_terms_here |= terms

            entry_terms_any.setdefault(entry_key, set()).update(any_terms_here)

        unmappable["total_entries"] = len(entry_terms_any)

        for entry_key, all_terms in entry_terms_any.items():
            if len(all_terms) > 0:
                unmappable["entries_with_any_go"] += 1

        unmappable["entries_with_no_go"] = unmappable["total_entries"] - unmappable["entries_with_any_go"]

        for aspect in ("bp", "mf", "cc"):
            entries_in_aspect = 0

            for entry_key, terms in entry_terms_by_aspect[aspect].items():
                if terms:
                    entries_in_aspect += 1
                    unmappable["term_entries_by_aspect"][aspect] += len(terms)

                    for t in terms:
                        unmappable["top_terms"][aspect][t] += 1

            unmappable["entries_with_aspect"][aspect] = entries_in_aspect

    summary = {
        "total_genes": total_genes,
        "genes_with_any_go": genes_with_any_go,
        "genes_with_no_go": total_genes - genes_with_any_go,
        "depth": compute_depth_stats(depth_list),
        "aspect_combo_counts": aspect_combo_counts,
        "gene_src_total": gene_src_total,
        "gene_src_by_aspect": gene_src_by_aspect,
        "term_src_total": term_src_total,
        "term_src_by_aspect": term_src_by_aspect,
        "genes_with_uniprot_any": genes_with_uniprot_any,
        "genes_with_genbank_any": genes_with_genbank_any,
        "genes_with_both_sources_any": genes_with_both_sources_any,
        "genes_with_uniprot_only_any": genes_with_uniprot_only_any,
        "genes_with_genbank_only_any": genes_with_genbank_only_any,
        "top_terms": top_terms,
        "unmappable_uniprot": unmappable,
    }

    out_dir = os.path.dirname(SUMMARY_LOG_FILE)

    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    write_summary_log(SUMMARY_LOG_FILE, summary, TOP_N_TERMS)

    depth_rows = []

    for gk in sorted(all_genes):
        gb_n = len(gb_terms_union.get(gk, set()))
        up_n = len(up_terms_union.get(gk, set()))
        comb_n = len(combined_union.get(gk, set()))

        depth_rows.append({
            "strain": gk[0],
            "locus tag": gk[1],
            "genbank_terms": gb_n,
            "uniprot_terms": up_n,
            "combined_terms": comb_n,
        })

    depth_df = pd.DataFrame(depth_rows)

    depth_dir = os.path.dirname(DEPTH_XLSX_OUT)

    if depth_dir:
        os.makedirs(depth_dir, exist_ok=True)

    with pd.ExcelWriter(DEPTH_XLSX_OUT, engine="openpyxl") as writer:
        depth_df.to_excel(writer, sheet_name="annotation_depth", index=False)

    print("Done.")
    print(f"Read:  {EXCEL_FILE}")
    print(f"Wrote: {SUMMARY_LOG_FILE}")
    print(f"Wrote: {DEPTH_XLSX_OUT}")


if __name__ == "__main__":
    main()