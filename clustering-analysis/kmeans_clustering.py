"""
Purpose:
    Perform k-means clustering on different feature sets for unique and non-unique
    Vibrio natriegens genes.

Description:
    This script reads a combined feature table and runs MiniBatchKMeans clustering
    separately on predefined feature sets, including amino acid composition,
    RSA/disorder, Q3 secondary structure, Q8 secondary structure, and all combined
    features. Features are standardized before clustering.

Input:
    - Excel or CSV file containing gene identifiers, gene occurrence type, and
      numerical feature columns.

Output:
    - One Excel file per feature set containing:
        - strain accession
        - locus tag
        - gene occurrence type
        - assigned cluster

Required Python packages:
    - pandas
    - scikit-learn
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update INPUT_FILE, SHEET_NAME, and OUTPUT_DIR before running the script.
"""

from pathlib import Path
import os
import pandas as pd
from sklearn.cluster import MiniBatchKMeans
from sklearn.preprocessing import StandardScaler

# Paths
INPUT_FILE = "path/to/combined_feature_table.xlsx"
SHEET_NAME = "Sheet1"

OUTPUT_DIR = "path/to/output/kmeans_clustering"

OUTPUT_COLUMNS = [
    "Strain accession",
    "Locus tag",
    "Gene occurrence type",
]

# Feature sets and selected k values
FEATURE_CONFIG = {
    "aa": {
        "k": "user-defined k value",
        "features": [
            "Alanine (%)", "Arginine (%)", "Asparagine (%)",
            "Cysteine (%)", "Glutamic acid (%)", "Glycine (%)",
            "Lysine (%)", "Methionine (%)", "Serine (%)",
            "Threonine (%)", "Tryptophan (%)", "Tyrosine (%)", "Valine (%)"
        ],
    },
    "disorder_rsa": {
        "k": "user-defined k value",
        "features": ["rsa", "disorder"],
    },
    "q3": {
        "k": "user-defined k value",
        "features": ["frac_q3_H", "frac_q3_E", "frac_q3_C"],
    },
    "q8": {
        "k": "user-defined k value",
        "features": [
            "frac_q8_H", "frac_q8_E", "frac_q8_C", "frac_q8_T",
            "frac_q8_G", "frac_q8_S", "frac_q8_B", "frac_q8_I"
        ],
    },
    "all_features": {
        "k": "user-defined k value",
        "features": [
            "Alanine (%)", "Arginine (%)", "Asparagine (%)",
            "Cysteine (%)", "Glutamic acid (%)", "Glycine (%)",
            "Lysine (%)", "Methionine (%)", "Serine (%)",
            "Threonine (%)", "Tryptophan (%)", "Tyrosine (%)", "Valine (%)",
            "rsa", "disorder", "asa",
            "frac_q3_H", "frac_q3_E", "frac_q3_C",
            "frac_q8_H", "frac_q8_E", "frac_q8_C", "frac_q8_T",
            "frac_q8_G", "frac_q8_S", "frac_q8_B", "frac_q8_I",
            "phi", "psi"
        ],
    },
}


# Helper functions
def load_input_file(file_path, sheet_name=None):
    file_path = Path(file_path)
    suffix = file_path.suffix.lower()

    if suffix in [".xlsx", ".xls"]:
        df = pd.read_excel(file_path, sheet_name=sheet_name, engine="openpyxl")
    elif suffix == ".csv":
        df = pd.read_csv(file_path)
    else:
        raise ValueError(f"Unsupported file type: {suffix}. Use .xlsx, .xls, or .csv")

    df.columns = df.columns.str.strip()
    return df


def validate_required_columns(df, required_columns):
    missing = []

    for col in required_columns:
        if col not in df.columns:
            missing.append(col)

    if missing:
        raise ValueError(
            "These required columns were not found in the input file:\n"
            + "\n".join(missing)
            + "\n\nAvailable columns are:\n"
            + "\n".join(df.columns.astype(str))
        )


def run_single_feature_set(df, feature_set_name, feature_columns, k, output_dir):
    print("\nRunning feature set:", feature_set_name)
    print("Selected k:", k)

    validate_required_columns(df, feature_columns + OUTPUT_COLUMNS)

    working_df = df.copy()

    working_df[feature_columns] = working_df[feature_columns].apply(pd.to_numeric, errors="coerce")

    before_drop = len(working_df)
    working_df = working_df.dropna(subset=feature_columns).copy()
    dropped = before_drop - len(working_df)

    if dropped > 0:
        print(f"Dropped {dropped:,} rows with missing values in selected features")

    if working_df.empty:
        raise ValueError(f"No rows remain after dropping missing values for feature set: {feature_set_name}")

    scaler = StandardScaler()
    X = scaler.fit_transform(working_df[feature_columns])

    kmeans = MiniBatchKMeans(
        n_clusters=k,
        init="k-means++",
        n_init=10,
        max_iter=300,
        tol=0.0001,
        random_state=0,
    )

    labels = kmeans.fit_predict(X)
    working_df["Cluster"] = labels

    output_df = working_df[OUTPUT_COLUMNS + ["Cluster"]].copy()

    os.makedirs(output_dir, exist_ok=True)

    output_file = Path(output_dir) / f"kmeans_clustering_results_{feature_set_name}_k={k}.xlsx"
    output_df.to_excel(output_file, index=False)

    print(f"Rows written: {len(output_df):,}")
    print(f"Output saved to: {output_file}")


# Main
def main():
    print(f"Loading input file: {INPUT_FILE}")

    if SHEET_NAME is not None:
        print(f"Reading sheet: {SHEET_NAME}")

    df = load_input_file(INPUT_FILE, SHEET_NAME)

    print(f"Rows loaded: {len(df):,}")

    validate_required_columns(df, OUTPUT_COLUMNS)

    for feature_set_name, config in FEATURE_CONFIG.items():
        run_single_feature_set(
            df=df,
            feature_set_name=feature_set_name,
            feature_columns=config["features"],
            k=config["k"],
            output_dir=OUTPUT_DIR,
        )

    print("\nAll clustering analyses completed.")

if __name__ == "__main__":
    main()