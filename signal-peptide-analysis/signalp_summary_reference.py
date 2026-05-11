"""
Purpose:
    Summarise SignalP 6.0 predictions for the Vibrio natriegens reference genome.

Description:
    This script reads SignalP 6.0 prediction_results.txt files from one or more
    output folders for the reference genome. It combines the predictions and
    produces one strain-level summary row containing:
        - total number of proteins analysed
        - total number of proteins predicted to contain signal peptides
        - counts for each SignalP signal peptide class
        - mean cleavage position
        - standard deviation of cleavage position

Input:
    - SignalP 6.0 output folders containing prediction_results.txt files.

Output:
    - CSV file containing one summary row for the reference genome.

Required Python packages:
    - None. This script uses only the Python standard library.

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update INPUT_DIR, OUTPUT_CSV, and STRAIN_NAME before running the script.
"""

import os
import csv
import re
from statistics import mean, stdev

# Paths
INPUT_DIR = "path/to/signalp_results/reference_genome"
OUTPUT_CSV = "path/to/output/signalp_results_reference_genome_summary.csv"

PREDICTION_RESULTS_FILENAME = "prediction_results.txt"
STRAIN_NAME = "vibrio_natriegens_reference_genome"

SP_CLASS_COLUMNS = [
    "SP(Sec/SPI)",
    "LIPO(Sec/SPII)",
    "TAT(Tat/SPI)",
    "TATLIPO(Tat/SPII)",
    "PILIN(Sec/SPIII)",
]

CS_REGEX = re.compile(r"CS pos:\s*(\d+)-(\d+)")


# Helper functions
def parse_prediction_results(file_path):
    rows = []
    header = None
    column_index = {}

    with open(file_path, "r", encoding="utf-8", errors="replace") as input_file:
        for line in input_file:
            line = line.rstrip("\n")

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


# Main
def main():
    all_rows = []

    for entry in sorted(os.listdir(INPUT_DIR)):
        run_path = os.path.join(INPUT_DIR, entry)

        if not os.path.isdir(run_path):
            continue

        prediction_file = os.path.join(run_path, PREDICTION_RESULTS_FILENAME)

        if not os.path.isfile(prediction_file):
            print(f"Warning: missing prediction_results.txt in {entry}")
            continue

        rows = parse_prediction_results(prediction_file)
        print(f"{entry}: added {len(rows)} proteins")

        all_rows.extend(rows)

    total_genes = 0
    total_sp_count = 0
    per_class_counts = {column: 0 for column in SP_CLASS_COLUMNS}
    cleavage_positions = []

    for row in all_rows:
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

    summary = {
        "strain": STRAIN_NAME,
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

    fieldnames = list(summary.keys())

    with open(OUTPUT_CSV, "w", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(summary)

    print(f"Wrote reference-genome summary to: {OUTPUT_CSV}")
if __name__ == "__main__":
    main()
