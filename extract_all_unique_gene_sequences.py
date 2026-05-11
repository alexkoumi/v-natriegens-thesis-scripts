"""
Purpose:
    Extract amino acid sequences for strain-specific genes and create the FASTA file
    used as input for MMseqs2 searches.

Description:
    This script reads a per-strain list of strain-specific locus tags, then searches
    the corresponding genome annotation folder for each strain. For each locus tag,
    the script reads the strain's GenBank-format annotation file and extracts the
    amino acid sequence from the CDS feature's "translation" qualifier.

    All extracted amino acid sequences are written to a single FASTA file:
    all_unique_genes.fasta.

Input:
    - unique_gene_list_per_strain.fasta

    - Genome annotation folders:
        A directory containing one folder per strain. Each strain folder should be
        named using the same strain name as the headers in unique_gene_list_per_strain.fasta
        and should contain at least one GenBank-format annotation file.

Output:
    - all_unique_genes.fasta:
        FASTA file containing amino acid sequences for all extracted strain-specific genes.

    - missing_unique_gene_sequences.tsv:
        Report listing locus tags that could not be found

Required Python packages:
    - biopython
    - tqdm

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update "UNIQUE_GENE_LIST", "GENOME_ROOT_DIR", "OUTPUT_FASTA", and
    "MISSING_REPORT" before running the script.
"""

from collections import defaultdict
from pathlib import Path
import csv
from Bio import SeqIO
from tqdm import tqdm

# File paths
UNIQUE_GENE_LIST = "path/to/unique_gene_list_per_strain.fasta"
GENOME_ROOT_DIR = "path/to/genome_folders"
OUTPUT_FASTA = "path/to/output/all_unique_genes.fasta"
MISSING_REPORT = "path/to/output/missing_unique_gene_sequences.tsv"

# Helper functions
def clean_locus_tag(locus_tag_line):
    return (
        locus_tag_line
        .replace(" (duplicate)", "")
        .strip()
    )

def parse_unique_gene_list(unique_gene_list_file):
    strain_to_locus_tags = defaultdict(list)
    current_strain = None
    with open(unique_gene_list_file, "r") as input_file:
        for line in input_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                current_strain = header.split(" (count")[0].strip()
                strain_to_locus_tags[current_strain] = []
            else:
                if current_strain is None:
                    continue
                locus_tag = clean_locus_tag(line)
                if locus_tag:
                    strain_to_locus_tags[current_strain].append(locus_tag)
    return dict(strain_to_locus_tags)

def find_genbank_files(strain_folder):
    genbank_files = []
    for pattern in ("*.gbff"):
        genbank_files.extend(sorted(strain_folder.rglob(pattern)))
    return genbank_files

def load_cds_translations(genbank_files):
    locus_tag_to_sequence = {}
    for genbank_file in genbank_files:
        for record in SeqIO.parse(str(genbank_file), "genbank"):
            for feature in record.features:
                if feature.type != "CDS":
                    continue

                locus_tags = feature.qualifiers.get("locus_tag", [])
                translations = feature.qualifiers.get("translation", [])

                if not locus_tags:
                    continue

                locus_tag = locus_tags[0]

                if translations:
                    locus_tag_to_sequence[locus_tag] = translations[0]

    return locus_tag_to_sequence


def wrap_sequence(sequence, line_length=60):
    return "\n".join(
        sequence[i:i + line_length] for i in range(0, len(sequence), line_length)
    )

def make_unique_fasta_header(strain_name, locus_tag, header_counts):
    base_header = f"{strain_name}|{locus_tag}"
    header_counts[base_header] += 1

    if header_counts[base_header] == 1:
        return base_header

    return f"{base_header}|copy{header_counts[base_header]}"

# Main
print("Reading per-strain unique gene list...")
strain_to_locus_tags = parse_unique_gene_list(UNIQUE_GENE_LIST)

OUTPUT_FASTA.parent.mkdir(parents=True, exist_ok=True)
MISSING_REPORT.parent.mkdir(parents=True, exist_ok=True)

missing_entries = []
fasta_header_counts = defaultdict(int)
total_extracted_sequences = 0

print("Extracting amino acid sequences")

with open(OUTPUT_FASTA, "w") as fasta_output:
    for strain_name, locus_tags in tqdm(strain_to_locus_tags.items(), desc="Processing strains"):
        strain_folder = GENOME_ROOT_DIR / strain_name
        if not strain_folder.exists():
            for locus_tag in locus_tags:
                missing_entries.append([strain_name, locus_tag, "strain folder not found"])
            continue

        genbank_files = find_genbank_files(strain_folder)

        if not genbank_files:
            for locus_tag in locus_tags:
                missing_entries.append([strain_name, locus_tag, "GenBank file not found"])
            continue

        locus_tag_to_sequence = load_cds_translations(genbank_files)

        for locus_tag in locus_tags:
            amino_acid_sequence = locus_tag_to_sequence.get(locus_tag)

            if amino_acid_sequence is None:
                missing_entries.append([strain_name, locus_tag, "locus tag or translation not found"])
                continue

            fasta_header = make_unique_fasta_header(
                strain_name,
                locus_tag,
                fasta_header_counts
            )

            fasta_output.write(f">{fasta_header}\n")
            fasta_output.write(f"{wrap_sequence(amino_acid_sequence)}\n")

            total_extracted_sequences += 1
# Write missing sequence report
with open(MISSING_REPORT, "w", newline="") as report_file:
    writer = csv.writer(report_file, delimiter="\t")
    writer.writerow(["Strain", "Locus_tag", "Reason"])
    writer.writerows(missing_entries)

# Print summary
print("Sequence extraction complete.")
print(f"Number of strains processed: {len(strain_to_locus_tags)}")
print(f"Total amino acid sequences extracted: {total_extracted_sequences}")
print(f"Missing entries reported: {len(missing_entries)}")
print("Outputs generated:")
print(f" - {OUTPUT_FASTA}")
print(f" - {MISSING_REPORT}")