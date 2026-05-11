"""
Purpose:
    Summarise SignalP 6.0 predictions for strain-specific Vibrio natriegens genes.

Description:
    This script reads zipped SignalP 6.0 output files for strain-specific genes.
    For each strain, it locates the prediction_results.txt file inside the ZIP
    archive, parses the SignalP predictions, and generates one summary row per
    strain containing:
        - total number of proteins analysed
        - total number of proteins predicted to contain signal peptides
        - counts for each SignalP signal peptide class
        - mean cleavage position
        - standard deviation of cleavage position

Input:
    - ZIP files containing SignalP 6.0 prediction_results.txt files.

Output:
    - CSV file containing one SignalP summary row per strain.

Required Python packages:
    - None. This script uses only the Python standard library.

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update INPUT_DIR and OUTPUT_CSV before running the script.
"""

import os
import zipfile
import csv
import re
from statistics import mean, stdev

# Paths
INPUT_DIR = "path/to/signalp_results/unique_genes"
OUTPUT_CSV = "path/to/output/signalp_results_summary_unique.csv"

PREDICTION_RESULTS_FILENAME = "prediction_results.txt"

SP_CLASS_COLUMNS = [
    "SP(Sec/SPI)",
    "LIPO(Sec/SPII)",
    "TAT(Tat/SPI)",
    "TATLIPO(Tat/SPII)",
    "PILIN(Sec/SPIII)",
]

CS_REGEX = re.compile(r"CS pos:\s*(\d+)-(\d+)")


# Helper functions
def find_prediction_results_member(zip_file):
    candidates = []

    for name in zip_file.namelist():
        if name.endswith(PREDICTION_RESULTS_FILENAME):
            candidates.append(name)

    if not candidates:
        return None

    return sorted(candidates, key=len)[0]


def parse_prediction_results(file_obj):
    rows = []
    header = None
    column_index = {}

    for raw_line in file_obj:
        line = raw_line.decode("utf-8", errors="replace").rstrip("\n")

        if not line.strip():
            continue

        if line.lstrip().startswith("#"):
            cleaned = line.lstrip("#").strip()
            parts = cleaned.split("\t")

            if "ID" in parts and "Prediction" in parts:
                header = [p.strip() for p in parts]
                column_index = {name: i for i, name in enumerate(header)}

            continue

        parts = line.split("\t")

        if header is None and ("ID" in parts and "Prediction" in parts):
            header = [p.strip() for p in parts]
            column_index = {name: i for i, name in enumerate(header)}
            continue

        if header is None:
            continue

        row = {}

        for column, index in column_index.items():
            if index < len(parts):
                row[column] = parts[index].strip()
            else:
                row[column] = ""

        rows.append(row)

    return rows


def extract_cleavage_position(cs_field):
    if not cs_field:
        return None

    match = CS_REGEX.search(cs_field)

    if not match:
        return None

    return int(match.group(1))


def summarize_strain(rows, strain_name):
    total_genes = 0
    total_sp_count = 0
    per_class_counts = {column: 0 for column in SP_CLASS_COLUMNS}
    cleavage_positions = []

    for row in rows:
        total_genes += 1
        prediction = row.get("Prediction", "OTHER")

        if prediction != "OTHER":
            total_sp_count += 1

            max_class = None
            max_probability = -1.0

            for column in SP_CLASS_COLUMNS:
                try:
                    probability = float(row.get(column, "0") or "0")
                except ValueError:
                    probability = 0.0

                if probability > max_probability:
                    max_probability = probability
                    max_class = column

            if max_class:
                per_class_counts[max_class] += 1

            cleavage_position = extract_cleavage_position(row.get("CS Position", ""))

            if cleavage_position is not None:
                cleavage_positions.append(cleavage_position)

    if cleavage_positions:
        mean_cleavage_position = f"{mean(cleavage_positions):.2f}"

        if len(cleavage_positions) > 1:
            sd_cleavage_position = f"{stdev(cleavage_positions):.2f}"
        else:
            sd_cleavage_position = "not_available"
    else:
        mean_cleavage_position = "not_available"
        sd_cleavage_position = "not_available"

    return {
        "strain": strain_name,
        "total_genes": total_genes,
        "total_SP_count": total_sp_count,
        "SP(Sec/SPI)_count": per_class_counts["SP(Sec/SPI)"],
        "LIPO(Sec/SPII)_count": per_class_counts["LIPO(Sec/SPII)"],
        "TAT(Tat/SPI)_count": per_class_counts["TAT(Tat/SPI)"],
        "TATLIPO(Tat/SPII)_count": per_class_counts["TATLIPO(Tat/SPII)"],
        "PILIN(Sec/SPIII)_count": per_class_counts["PILIN(Sec/SPIII)"],
        "mean_cleavage_position": mean_cleavage_position,
        "sd_cleavage_position": sd_cleavage_position,
    }


# Main
def main():
    summaries = []

    for filename in sorted(os.listdir(INPUT_DIR)):
        if not filename.lower().endswith(".zip"):
            continue

        strain_name = filename[:-4]
        zip_path = os.path.join(INPUT_DIR, filename)

        try:
            with zipfile.ZipFile(zip_path, "r") as zip_file:
                member = find_prediction_results_member(zip_file)

                if member is None:
                    print(f"Warning: no prediction_results.txt in {filename}")
                    continue

                with zip_file.open(member) as input_file:
                    rows = parse_prediction_results(input_file)

            print(f"{strain_name}: parsed {len(rows)} genes")

            summaries.append(summarize_strain(rows, strain_name))

        except Exception as error:
            print(f"Error: failed to process {filename}: {error}")

    fieldnames = [
        "strain",
        "total_genes",
        "total_SP_count",
        "SP(Sec/SPI)_count",
        "LIPO(Sec/SPII)_count",
        "TAT(Tat/SPI)_count",
        "TATLIPO(Tat/SPII)_count",
        "PILIN(Sec/SPIII)_count",
        "mean_cleavage_position",
        "sd_cleavage_position",
    ]

    with open(OUTPUT_CSV, "w", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summaries)

    print(f"Done. Wrote {len(summaries)} strain summaries to: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()