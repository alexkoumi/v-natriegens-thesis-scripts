"""
Purpose:
    Fetch taxonomy for all target proteins in MMseqs2 result files.

Description:
    This script reads MMseqs2 result files in `.m8` format and extracts the query
    and target protein accession columns. It then retrieves TaxID information for
    target proteins using the UniProtKB REST API and retrieves taxonomic lineage
    information using NCBI Entrez Taxonomy, with an NCBI Datasets fallback.

    The final output is a taxonomy-annotated TSV file containing the query protein,
    target protein, TaxID, and major taxonomic ranks.

Input:
    - MMseqs2 `.m8` result files.

Output:
    - Taxonomy-annotated TSV file containing:
        query, target, TaxID, Kingdom, Phylum, Class, Order, Family, Genus, Species
    - Summary report listing retrieved and failed taxonomy entries.
    - JSON cache storing retrieved TaxID-to-lineage results.

Required Python packages:
    - pandas
    - requests
    - tqdm
    - biopython

Notes:
    Local file paths and personal email address used during the original analysis
    were replaced with placeholder values. Update dataset_name, results_dir,
    out_dir, shared_cache_path, and Entrez.email before running the script.
"""

import os
import time
import json
import re
import requests
import pandas as pd
from io import StringIO
from tqdm import tqdm
from Bio import Entrez

# Paths
dataset_name = "domain name"
results_dir = f"/path/to/mmseqs2_{dataset_name}_results"
out_dir = f"/path/to/output/all_hits_taxonomy/{dataset_name}"

# Shared cache for domains
shared_cache_path = "/path/to/output/taxonomy_cache.json"
os.makedirs(os.path.dirname(shared_cache_path), exist_ok=True)

# Parameters
Entrez.email = "your.email@example.com"
UNIPROT_BATCH = 250
NCBI_BATCH = 100
REQUEST_TIMEOUT = 60
PAUSE = 0.2
RETRIES = 3

# Helper functions
def canonicalize(acc: str) -> str:
    return acc.split("-")[0] if acc else acc

def clean_species(scientific_name: str) -> str:
    if not scientific_name:
        return ""

    toks = scientific_name.strip().split()

    for i in range(len(toks) - 2, -1, -1):
        if re.match(r"^[A-Z][a-z-]+$", toks[i]) and re.match(r"^[a-z][a-z0-9.-]*$", toks[i + 1]):
            return f"{toks[i]} {toks[i + 1]}"

    return scientific_name.strip()

def _is_blank_lineage(v):
    if not isinstance(v, list) or len(v) != 7:
        return True

    return all(not str(x).strip() for x in v)

# Read files
def read_all_m8(folder: str) -> pd.DataFrame:
    files = [f for f in os.listdir(folder) if f.endswith(".m8")]

    if not files:
        raise SystemExit(f"No .m8 files found in {folder}")

    frames = []

    for f in tqdm(files, desc="Reading .m8 files", unit="file"):
        df = pd.read_csv(
            os.path.join(folder, f),
            sep="\t",
            header=None,
            usecols=[0, 1],
            names=["query", "target"],
            dtype=str
        )
        frames.append(df)

    return pd.concat(frames, ignore_index=True)

def load_cache(path):
    if os.path.exists(path):
        try:
            with open(path, "r") as f:
                return json.load(f)
        except Exception:
            return {}

    return {}

def save_cache(path, data):
    tmp = path + ".tmp"

    with open(tmp, "w") as f:
        json.dump(data, f)

    os.replace(tmp, path)

# UniProt TaxIDs
def uniprot_stream_taxids(accessions):
    result = {}
    accs = sorted(set(canonicalize(a) for a in accessions))
    for i in tqdm(range(0, len(accs), UNIPROT_BATCH), desc="Fetching TaxIDs (UniProt)", unit="batch"):
        batch = accs[i:i + UNIPROT_BATCH]
        params = {"accessions": ",".join(batch), "fields": "accession,organism_id", "format": "tsv"}
        text = None
        for attempt in range(RETRIES):
            try:
                r = requests.get("https://rest.uniprot.org/uniprotkb/accessions", params=params, timeout=REQUEST_TIMEOUT)
                if r.status_code == 200:
                    text = r.text
                    break
                time.sleep(PAUSE * (attempt + 1))
            except requests.RequestException:
                time.sleep(PAUSE * (attempt + 1))
        if text is None:
            continue
        df = pd.read_csv(StringIO(text), sep="\t")
        if "Entry" in df.columns and "Organism (ID)" in df.columns:
            for _, row in df.iterrows():
                result[str(row["Entry"]).strip()] = str(row["Organism (ID)"]).strip()
        time.sleep(PAUSE)
    return result

# NCBI Taxonomy fallback
def fetch_lineage_entrez_record(rec):
    lin_map = {r.get("Rank", "").lower(): r.get("ScientificName", "") for r in rec.get("LineageEx", [])}
    if not lin_map:
        return None
    species = lin_map.get("species") or clean_species(rec.get("ScientificName", ""))
    return [lin_map.get(r, "") for r in ["kingdom", "phylum", "class", "order", "family", "genus"]] + [species]

def fetch_lineage_datasets_single(taxid: str):
    try:
        r = requests.get(f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{taxid}", timeout=REQUEST_TIMEOUT)
        if not r.ok:
            return None
        data = r.json()
        lineage = data.get("taxonomy_nodes", [{}])[0].get("lineage", [])
        ranks = {x.get("rank", "").lower(): x.get("organism_name", "") for x in lineage}
        return [ranks.get(k, "") for k in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]]
    except Exception:
        return None

def fetch_ncbi_taxonomy_batched(taxids, existing_cache):
    new_items, failed = {}, []
    entrez_ok, datasets_ok = 0, 0
    ids = [t for t in taxids if t and (t not in existing_cache or _is_blank_lineage(existing_cache.get(t)))]

    for i in tqdm(range(0, len(ids), NCBI_BATCH), desc="Fetching taxonomy (Entrez+fallback)", unit="batch"):
        batch = ids[i:i + NCBI_BATCH]
        try:
            h = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
            records = Entrez.read(h)
            h.close()
        except Exception:
            records = []

        rec_by_taxid = {r.get("TaxId", ""): r for r in records}
        for tx in batch:
            rec = rec_by_taxid.get(tx)
            used = False
            if rec:
                ent_lineage = fetch_lineage_entrez_record(rec)
                if ent_lineage:
                    new_items[tx] = ent_lineage
                    entrez_ok += 1
                    used = True
            if not used:
                ds_lineage = fetch_lineage_datasets_single(tx)
                if ds_lineage:
                    new_items[tx] = ds_lineage
                    datasets_ok += 1
                    used = True
            if not used:
                failed.append(tx)
            time.sleep(PAUSE)
    return new_items, failed, entrez_ok, datasets_ok

# main workflow
def main():
    df = read_all_m8(results_dir)
    df["target"] = df["target"].astype(str).map(canonicalize)
    uniq_targets = sorted(set(df["target"]))
    acc2taxid = uniprot_stream_taxids(uniq_targets)
    df["TaxID"] = df["target"].map(acc2taxid).fillna("")
    failed_taxid_accs = [a for a in uniq_targets if a not in acc2taxid]
    cache = load_cache(shared_cache_path)

    taxids_needed = [t for t in df["TaxID"].unique() if t]
    new_lin, failed_ncbi, entrez_ok, datasets_ok = fetch_ncbi_taxonomy_batched(taxids_needed, cache)
    if new_lin:
        cache.update(new_lin)
        save_cache(shared_cache_path, cache)

    cols = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    for i, c in enumerate(cols):
        df[c] = df["TaxID"].map(lambda x: cache.get(x, [""] * 7)[i])

    os.makedirs(out_dir, exist_ok=True)
    out_tsv = os.path.join(out_dir, f"{dataset_name}_all_hits_taxonomy.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)

    report = os.path.join(out_dir, f"{dataset_name}_all_hits_taxonomy_report.txt")
    with open(report, "w") as f:
        f.write(f"Dataset: {dataset_name}\n")
        f.write(f"Distinct targets: {len(uniq_targets)}\n")
        f.write(f"Retrieved TaxIDs: {len(acc2taxid)}\n")
        f.write(f"Failed TaxID retrievals: {len(failed_taxid_accs)}\n")
        if failed_taxid_accs:
            f.write("Proteins with missing TaxIDs:\n")
            f.write("\n".join(sorted(failed_taxid_accs)) + "\n\n")
        f.write(f"TaxIDs with lineage via Entrez: {entrez_ok}\n")
        f.write(f"TaxIDs with lineage via Datasets fallback: {datasets_ok}\n")
        f.write(f"Failed taxonomy retrievals: {len(failed_ncbi)}\n")
        if failed_ncbi:
            f.write("List of TaxIDs with missing taxonomy:\n")
            f.write("\n".join(sorted(failed_ncbi)) + "\n\n")
    print(f"\nTSV: {out_tsv}\nReport: {report}\nCache: {shared_cache_path}\n")
if __name__ == "__main__":
    main()
