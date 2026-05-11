"""
Purpose:
    Evaluate k-means clustering stability using Jaccard similarity across repeated
    clustering runs.

Description:
    This script reads a combined feature table and performs repeated k-means
    clustering for user-defined feature sets and k-value ranges. For each k value,
    clustering is repeated several times using different random seeds. Cluster
    similarity between runs is then estimated using mean Jaccard similarity after
    optimal cluster matching.

Input:
    - Excel or CSV file containing gene-level feature values.

Output:
    - Excel files containing Jaccard stability scores for each tested k value.
    - PNG plots showing Jaccard stability across k values.

Required Python packages:
    - pandas
    - numpy
    - matplotlib
    - scipy
    - scikit-learn
    - openpyxl

Notes:
    Local file paths and dataset-specific feature names were replaced with placeholder
    values. Update INPUT_FILE, SHEET_NAME, OUTPUT_DIR, and FEATURE_CONFIG before
    running the script.
"""

import os
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# Paths
INPUT_FILE = "path/to/combined_feature_table.xlsx"
SHEET_NAME = "Sheet1"
OUTPUT_DIR = "path/to/output/jaccard_stability"

RUNS_PER_K = 10

# Plot settings
LOWER_AXIS_MIN = 0.0
LOWER_AXIS_MAX = 0.01
UPPER_AXIS_MIN = 0.60
UPPER_AXIS_MAX = 1.00

UPPER_MAJOR_TICK_STEP = 0.05
UPPER_MINOR_TICK_STEP = 0.01

PLOT_DPI = 300
PLOT_FIGSIZE = (8, 6)
LINE_COLOR = "0.5"

# Feature sets and k ranges
# Replace these placeholder feature names with the exact column names in your input file.
FEATURE_CONFIG = {
    "feature_set_1": {
        "features": [
            "feature_column_1",
            "feature_column_2",
            "feature_column_3",
        ],
        "k_min": "user-defined k value",
        "k_max": "used-defined k value",
    },
    "feature_set_2": {
        "features": [
            "feature_column_4",
            "feature_column_5",
            "feature_column_6",
        ],
        "k_min": 3,
        "k_max": 5,
    },
}

# Helper functions
def load_input_table(input_file, sheet_name=None):
    suffix = Path(input_file).suffix.lower()

    if suffix in [".xlsx", ".xls"]:
        return pd.read_excel(input_file, sheet_name=sheet_name, engine="openpyxl")

    if suffix == ".csv":
        return pd.read_csv(input_file)

    raise ValueError(f"Unsupported input format: {suffix}")

def cluster_jaccard(labels_a, labels_b):
    clusters_a = np.unique(labels_a)
    clusters_b = np.unique(labels_b)

    jaccard_matrix = np.zeros((len(clusters_a), len(clusters_b)), dtype=float)

    for i, cluster_a in enumerate(clusters_a):
        set_a = set(np.where(labels_a == cluster_a)[0])

        for j, cluster_b in enumerate(clusters_b):
            set_b = set(np.where(labels_b == cluster_b)[0])
            union = len(set_a | set_b)

            if union > 0:
                score = len(set_a & set_b) / union
            else:
                score = 0.0

            jaccard_matrix[i, j] = score

    row_ind, col_ind = linear_sum_assignment(1 - jaccard_matrix)
    matched_scores = jaccard_matrix[row_ind, col_ind]

    return float(np.mean(matched_scores))

def make_broken_axis_plot(k_values, stability_scores, output_path, runs_per_k):
    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=PLOT_FIGSIZE,
        gridspec_kw={"height_ratios": [30, 0.25], "hspace": 0.08}
    )

    ax_top.plot(k_values, stability_scores, marker="o", color=LINE_COLOR, linewidth=1.5)

    ax_top.set_ylim(UPPER_AXIS_MIN, UPPER_AXIS_MAX)
    ax_bottom.set_ylim(0.0, 0.001)

    ax_top.set_ylabel(f"Mean Jaccard similarity ({runs_per_k} runs)")
    ax_bottom.set_xlabel("k value")

    ax_bottom.set_xticks(k_values)
    ax_bottom.set_yticks([0.0])
    ax_bottom.set_yticklabels(["0"])

    ax_top.yaxis.set_major_locator(MultipleLocator(UPPER_MAJOR_TICK_STEP))
    ax_top.yaxis.set_minor_locator(MultipleLocator(UPPER_MINOR_TICK_STEP))

    ax_top.grid(True, which="major", linewidth=0.8)
    ax_top.grid(True, which="minor", linewidth=0.4, alpha=0.5)
    ax_bottom.grid(False)

    ax_top.spines["bottom"].set_visible(False)
    ax_bottom.spines["top"].set_visible(False)
    ax_top.tick_params(labeltop=False, bottom=False)
    ax_bottom.xaxis.tick_bottom()
    ax_bottom.spines["right"].set_visible(False)

    d = 0.015

    kwargs = dict(transform=ax_top.transAxes, color="k", clip_on=False, linewidth=1.0)
    ax_top.plot((-d, +d), (-d, +d), **kwargs)
    ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)

    kwargs.update(transform=ax_bottom.transAxes)
    ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    plt.savefig(output_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close()

# Load data
print("Loading dataset")

df = load_input_table(INPUT_FILE, SHEET_NAME)
df.columns = df.columns.str.strip()

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Run Jaccard stability analysis
for feature_set_name, config in FEATURE_CONFIG.items():
    print(f"\nRunning feature set: {feature_set_name}")

    selected_features = config["features"]
    k_min = config["k_min"]
    k_max = config["k_max"]
    missing_features = []

    for feature in selected_features:
        if feature not in df.columns:
            missing_features.append(feature)

    if missing_features:
        raise KeyError(
            f"Feature set '{feature_set_name}' is missing these columns:\n- "
            + "\n- ".join(missing_features)
        )

    print("Extracting features")

    X_raw = df[selected_features].apply(pd.to_numeric, errors="coerce")

    rows_before = len(X_raw)
    X_raw = X_raw.dropna()
    rows_after = len(X_raw)

    print(f"Rows before missing-value removal: {rows_before}")
    print(f"Rows after missing-value removal: {rows_after}")

    if X_raw.empty:
        raise ValueError(
            f"No rows remain after removing missing values for feature set: {feature_set_name}"
        )

    print("Standardising data")
    scaler = StandardScaler()
    X = scaler.fit_transform(X_raw)
    k_values = list(range(k_min, k_max + 1))
    stability_scores = []

    print(f"Running Jaccard stability analysis ({RUNS_PER_K} runs per k)")

    for k in k_values:
        print(f"k = {k}")

        all_labels = []

        for run_number in range(RUNS_PER_K):
            model = KMeans(
                n_clusters=k,
                init="k-means++",
                n_init=10,
                max_iter=300,
                tol=1e-4,
                random_state=run_number,
            )

            model.fit(X)
            all_labels.append(model.labels_)

        pair_scores = []

        for i, j in combinations(range(RUNS_PER_K), 2):
            score = cluster_jaccard(all_labels[i], all_labels[j])
            pair_scores.append(score)

        mean_score = float(np.mean(pair_scores))

        print(f"Mean stability: {mean_score:.4f}")

        stability_scores.append(mean_score)

    stability_table_path = os.path.join(
        OUTPUT_DIR,
        f"jaccard_stability_table_{feature_set_name}.xlsx"
    )

    print(f"Saving Jaccard table: {stability_table_path}")

    pd.DataFrame({
        "k": k_values,
        "jaccard_stability": stability_scores
    }).to_excel(stability_table_path, index=False)

    stability_plot_path = os.path.join(
        OUTPUT_DIR,
        f"jaccard_stability_plot_{feature_set_name}.png"
    )

    print(f"Generating plot: {stability_plot_path}")

    make_broken_axis_plot(
        k_values=k_values,
        stability_scores=stability_scores,
        output_path=stability_plot_path,
        runs_per_k=RUNS_PER_K
    )
print("\nAll Jaccard analyses completed.")