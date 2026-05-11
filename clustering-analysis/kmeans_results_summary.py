"""
Purpose:
    Summarise k-means clustering results and generate PCA plots for each feature set.

Description:
    This script reads a combined feature table and previously generated k-means
    cluster assignment files. For each user-defined feature set, it merges the
    cluster assignments with the original feature values, calculates cluster
    centroids and cluster sizes, performs PCA on standardised feature values, and
    writes a summary workbook.

Input:
    - Combined feature table in Excel or CSV format.
    - K-means cluster assignment files for each feature set.

Output:
    - Excel summary workbook for each feature set containing:
        - cluster centroids
        - cluster sizes
        - PCA coordinates
        - PCA loadings
    - PCA plot PNG file for each feature set.

Required Python packages:
    - pandas
    - matplotlib
    - scikit-learn
    - openpyxl

Notes:
    Local file paths and dataset-specific feature names were replaced with placeholder
    values. Update FEATURE_FILE, FEATURE_SHEET, FEATURE_CONFIG, and feature column
    names before running the script.
"""

import os
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Paths
FEATURE_FILE = "path/to/combined_feature_table.xlsx"
FEATURE_SHEET = "Sheet1"

# Feature sets and cluster assignment files
# Replace placeholder feature names and cluster_file paths with the exact values used in your analysis.
FEATURE_CONFIG = {
    "feature_set_1": {
        "k": "user-defined k value",
        "cluster_file": "path/to/kmeans_clustering_results_feature_set_1.xlsx",
        "features": [
            "feature_column_1",
            "feature_column_2",
            "feature_column_3",
        ],
    },
    "feature_set_2": {
        "k": "user-defined k value",
        "cluster_file": "path/to/kmeans_clustering_results_feature_set_2.xlsx",
        "features": [
            "feature_column_4",
            "feature_column_5",
            "feature_column_6",
        ],
    },
}

# Constants
MERGE_COLUMNS = ["Strain accession", "Locus tag", "Gene occurrence type"]
CLUSTER_COLUMN = "Cluster"

# Cluster colours
CLUSTER_COLORS = {
    0: "red",
    1: "yellow",
    2: "blue",
    3: "green",
}

# Marker by gene type
GENE_TYPE_MARKERS = {
    "unique": "o",
    "non-unique": "^",
}

GENE_TYPE_LABELS = {
    "unique": "Strain-specific genes",
    "non-unique": "Non-strain-specific genes",
}


# Helper functions
def load_input_table(file_path, sheet_name=0):
    path = Path(file_path)
    suffix = path.suffix.lower()

    if suffix in [".xlsx", ".xls"]:
        df = pd.read_excel(path, sheet_name=sheet_name, engine="openpyxl")
    elif suffix == ".csv":
        df = pd.read_csv(path)
    else:
        raise ValueError(f"Unsupported input format: {suffix}")

    df.columns = df.columns.str.strip()

    return df


def validate_columns(df, required_columns, df_name):
    missing = []

    for col in required_columns:
        if col not in df.columns:
            missing.append(col)

    if missing:
        raise KeyError(
            f"The following required columns are missing from {df_name}:\n- "
            + "\n- ".join(missing)
        )


def normalize_gene_occurrence_values(df):
    df = df.copy()

    df["Gene occurrence type"] = (
        df["Gene occurrence type"]
        .astype(str)
        .str.strip()
        .str.lower()
        .replace({
            "nonunique": "non-unique",
            "non_unique": "non-unique",
            "non unique": "non-unique",
        })
    )

    return df


def build_cluster_sizes(df):
    cluster_sizes_total = (
        df[CLUSTER_COLUMN]
        .value_counts()
        .rename("total_count")
        .sort_index()
    )

    cluster_sizes_unique = (
        df[df["Gene occurrence type"] == "unique"][CLUSTER_COLUMN]
        .value_counts()
        .rename("strain-specific_count")
        .sort_index()
    )

    cluster_sizes_nonunique = (
        df[df["Gene occurrence type"] == "non-unique"][CLUSTER_COLUMN]
        .value_counts()
        .rename("non-strain-specific_count")
        .sort_index()
    )

    cluster_sizes = pd.concat(
        [cluster_sizes_total, cluster_sizes_unique, cluster_sizes_nonunique],
        axis=1
    ).fillna(0).astype(int)

    return cluster_sizes


def make_pca_plot(df, feature_columns, pca_plot_file):
    scaler = StandardScaler()
    X = scaler.fit_transform(df[feature_columns])

    pca = PCA(n_components=2)
    coordinates = pca.fit_transform(X)

    plot_df = df.copy()
    plot_df["PC1"] = coordinates[:, 0]
    plot_df["PC2"] = coordinates[:, 1]

    pca_loadings = pd.DataFrame(
        pca.components_,
        columns=feature_columns,
        index=["PC1", "PC2"]
    )

    pca_coordinates = plot_df[
        ["Strain accession", "Locus tag", "Gene occurrence type", CLUSTER_COLUMN, "PC1", "PC2"]
    ].copy()

    plt.figure(figsize=(9, 7))

    gene_types_present = sorted(plot_df["Gene occurrence type"].dropna().unique())
    clusters_present = sorted(plot_df[CLUSTER_COLUMN].dropna().unique())

    for gene_type in gene_types_present:
        marker = GENE_TYPE_MARKERS.get(gene_type, "o")
        subset = plot_df[plot_df["Gene occurrence type"] == gene_type]

        for cluster_id in clusters_present:
            cluster_subset = subset[subset[CLUSTER_COLUMN] == cluster_id]

            if cluster_subset.empty:
                continue

            plt.scatter(
                cluster_subset["PC1"],
                cluster_subset["PC2"],
                c=CLUSTER_COLORS.get(cluster_id, "grey"),
                marker=marker,
                s=14,
                alpha=0.8,
                edgecolors="black",
                linewidths=0.2,
            )

    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()

    legend_elements = []

    if "unique" in gene_types_present:
        legend_elements.append(
            Line2D(
                [0], [0],
                marker="o",
                color="w",
                label="Strain-specific genes",
                markerfacecolor="black",
                markeredgecolor="black",
                markersize=7
            )
        )

    if "non-unique" in gene_types_present:
        legend_elements.append(
            Line2D(
                [0], [0],
                marker="^",
                color="w",
                label="Non-strain-specific genes",
                markerfacecolor="black",
                markeredgecolor="black",
                markersize=7
            )
        )

    for cluster_id in clusters_present:
        legend_elements.append(
            Line2D(
                [0], [0],
                marker="s",
                color=CLUSTER_COLORS.get(cluster_id, "grey"),
                label=f"Cluster {cluster_id}",
                markersize=10,
                linewidth=0
            )
        )

    plt.legend(handles=legend_elements, loc="upper right", frameon=True)
    plt.savefig(pca_plot_file, dpi=300)
    plt.close()

    return pca_coordinates, pca_loadings


def run_single_summary(feature_df, feature_set_name, config):
    k = config["k"]
    feature_columns = config["features"]
    cluster_file = config["cluster_file"]

    output_dir = os.path.dirname(cluster_file)
    os.makedirs(output_dir, exist_ok=True)

    summary_output_file = os.path.join(
        output_dir,
        f"kmeans_clustering_results_summary_{feature_set_name}_k={k}.xlsx"
    )

    pca_plot_file = os.path.join(
        output_dir,
        f"pca_plot_{feature_set_name}_k={k}.png"
    )

    print(f"\nProcessing feature set: {feature_set_name}")
    print(f"Reading cluster file: {cluster_file}")

    cluster_df = load_input_table(cluster_file)

    cluster_df = normalize_gene_occurrence_values(cluster_df)
    feature_df = normalize_gene_occurrence_values(feature_df)

    validate_columns(
        cluster_df,
        MERGE_COLUMNS + [CLUSTER_COLUMN],
        f"cluster file for {feature_set_name}"
    )

    validate_columns(
        feature_df,
        MERGE_COLUMNS + feature_columns,
        f"feature file for {feature_set_name}"
    )

    merged_df = pd.merge(
        cluster_df,
        feature_df[MERGE_COLUMNS + feature_columns],
        on=MERGE_COLUMNS,
        how="left",
        validate="one_to_one"
    )

    missing_feature_rows = merged_df[feature_columns].isna().any(axis=1).sum()

    if missing_feature_rows > 0:
        raise ValueError(
            f"{missing_feature_rows} rows in {feature_set_name} could not be matched back to the feature table."
        )

    print("Computing cluster centroids")
    cluster_centroids = merged_df.groupby(CLUSTER_COLUMN)[feature_columns].mean()

    print("Calculating cluster sizes")
    cluster_sizes = build_cluster_sizes(merged_df)

    print("Performing PCA")
    pca_coordinates, pca_loadings = make_pca_plot(
        df=merged_df,
        feature_columns=feature_columns,
        pca_plot_file=pca_plot_file
    )

    print(f"PCA plot saved to: {pca_plot_file}")

    print(f"Writing summary Excel file: {summary_output_file}")

    with pd.ExcelWriter(summary_output_file, engine="openpyxl") as writer:
        cluster_centroids.to_excel(writer, sheet_name="cluster_centroids")
        cluster_sizes.to_excel(writer, sheet_name="cluster_sizes")
        pca_coordinates.to_excel(writer, sheet_name="pca_coordinates", index=False)
        pca_loadings.to_excel(writer, sheet_name="pca_loadings")

    print("Summary file written successfully")


# Main
def main():
    print("Loading feature file")

    feature_df = load_input_table(FEATURE_FILE, FEATURE_SHEET)

    validate_columns(feature_df, MERGE_COLUMNS, "feature file")

    for feature_set_name, config in FEATURE_CONFIG.items():
        run_single_summary(feature_df, feature_set_name, config)

    print("\nAll summary files and PCA plots completed.")

if __name__ == "__main__":
    main()