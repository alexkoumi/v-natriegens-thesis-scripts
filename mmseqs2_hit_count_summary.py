"""
Purpose:
    Count the number of MMseqs2 hits detected for each strain-specific
    Vibrio natriegens gene across major reference proteome domains.

Description:
    This script uses a master list of strain-specific genes to ensure that genes
    with zero MMseqs2 hits are retained in the final output. It reads taxonomy-
    annotated MMseqs2 result files, counts the number of hits per query gene for
    each domain, and writes a combined Excel summary table.

    The script also generates a log file summarising, for each domain:
        - number of genes with at least one hit
        - number of genes with zero hits
        - number of distinct target proteins

Input:
    - Master strain-specific gene list file.
    - Domain-specific MMseqs2 taxonomy TSV files with at least the columns:
        - query
        - target

Output:
    - Excel file containing per-domain hit counts for each query gene.
    - Log file summarising hit detection across domains.
    - Text files listing genes with zero hits for each domain.

Required Python packages:
    - pandas
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update base_dir, output_excel, log_file, master_file, and the paths in
    files before running the script.
"""

import os
import pandas as pd

# File paths
base_dir = "path/to/mmseqs_hits_taxonomy"
output_excel = os.path.join(base_dir, "hit_count_per_query.xlsx")
log_file = os.path.join(base_dir, "domain_summary_log.txt")
MASTER_FILE = "path/to/unique_gene_list_per_strain.fasta"
files = {
    "Archaea": os.path.join(base_dir, "archaea", "archaea_all_hits_taxonomy.tsv"),
    "Eukaryota": os.path.join(base_dir, "eukaryotic", "eukaryotic_all_hits_taxonomy.tsv"),
    "Viral": os.path.join(base_dir, "viral", "viral_all_hits_taxonomy.tsv"),
    "Bacterial": os.path.join(base_dir, "bacterial", "bacterial_all_hits_taxonomy.tsv"),
}

# Normalisation functions
def normalize_master_gene(gene_id: str):
    gene_id = str(gene_id).strip()
    if "(duplicate)" in gene_id:
        return None
    return gene_id

def normalize_tsv_query(gene_id: str):
    gene_id = str(gene_id).strip()
    if gene_id.endswith("_1"):
        gene_id = gene_id[:-2]
    return gene_id

# Read master strain-specific gene list
def read_master_gene_list(path: str) -> pd.DataFrame:
    rows = []
    strain = None

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                strain = line[1:].strip()
                continue

            gene = normalize_master_gene(line)
            if gene:
                rows.append({"strain": strain, "query": gene})

    return pd.DataFrame(rows)

# Analyse each domain-specific MMseqs2 result file
def analyze_domain(tsv_path, domain, master_queries, log_handle):
    if not os.path.exists(tsv_path):
        log_handle.write(f"{domain}: FILE MISSING\n")
        return None

    df = pd.read_csv(tsv_path, sep="\t", usecols=["query", "target"], dtype=str)
    df["query"] = df["query"].map(normalize_tsv_query)

    # distinct proteins
    distinct_targets = df["target"].nunique()

    # genes with ≥1 hit
    genes_with_hits = set(df["query"].unique())

    all_genes = set(master_queries)
    zero_hit_genes = sorted(all_genes - genes_with_hits)

    # write orphan file
    orphan_file = os.path.join(base_dir, f"{domain}_zero_hit_genes.txt")
    with open(orphan_file, "w") as f:
        for g in zero_hit_genes:
            f.write(g + "\n")

    log_handle.write(
        f"{domain}:\n"
        f"  Genes with ≥1 hit: {len(genes_with_hits)}\n"
        f"  Genes with 0 hits: {len(zero_hit_genes)}\n"
        f"  Distinct target proteins: {distinct_targets}\n\n"
    )

    return df["query"].value_counts().rename(f"{domain}_hits")

# Main workflow
def main():
    print("\n--- Domain homology summary ---\n")

    master = read_master_gene_list(MASTER_FILE)
    master_queries = master["query"].tolist()

    merged = master.copy()

    with open(log_file, "w") as log:
        log.write("Domain Summary\n\n")
        log.write(f"Total non-duplicated strain-specific genes: {len(master_queries)}\n\n")

        for domain, path in files.items():
            print(f"Processing {domain}...")
            counts = analyze_domain(path, domain, master_queries, log)

            if counts is not None:
                counts = counts.reset_index()
                counts.columns = ["query", f"{domain}_hits"]
                merged = merged.merge(counts, on="query", how="left")

        # fill missing with zero
        for c in merged.columns:
            if c.endswith("_hits"):
                merged[c] = merged[c].fillna(0).astype(int)
    merged.to_excel(output_excel, index=False)
    print(f"\nExcel written: {output_excel}")
    print(f"Summary log: {log_file}")
    print("\n--- Done ---\n")
if __name__ == "__main__":
    main()