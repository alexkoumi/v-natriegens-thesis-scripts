"""
Purpose:
    Perform elbow analysis for k-means clustering across multiple feature sets.

Description:
    This script reads a combined feature table, extracts predefined feature sets,
    standardises the selected features, and runs MiniBatchKMeans clustering across
    a range of k values. For each feature set, it records the inertia values and
    generates an elbow plot.

Input:
    - Excel or CSV file containing gene-level feature values.

Output:
    - Excel files containing inertia values for each tested k value.
    - Elbow plot PNG files for each feature set.

Required Python packages:
    - pandas
    - numpy
    - matplotlib
    - scikit-learn
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update INPUT_FILE, SHEET_NAME, and OUTPUT_DIR before running the script.
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, AutoMinorLocator
from sklearn.cluster import MiniBatchKMeans
from sklearn.preprocessing import StandardScaler

# Paths
INPUT_FILE = "path/to/vnat_combined_features.xlsx"
SHEET_NAME = "Sheet1"
OUTPUT_DIR = "path/to/output/elbow_analysis"

K_MIN = "user-defined minimum k value"
K_MAX = "user-defined minimum k value"
PLOT_DPI = 300
PLOT_FIGSIZE = (8, 6)
LINE_COLOR = "0.5"

# Feature sets
FEATURE_CONFIG = {
    "aa": {
        "features": [
            "Alanine (%)", "Arginine (%)", "Asparagine (%)",
            "Cysteine (%)", "Glutamic acid (%)", "Glycine (%)",
            "Lysine (%)", "Methionine (%)", "Serine (%)",
            "Threonine (%)", "Tryptophan (%)", "Tyrosine (%)", "Valine (%)"
        ]
    },
    "disorder_rsa": {
        "features": ["rsa", "disorder"]
    },
    "q3": {
        "features": ["frac_q3_H", "frac_q3_E", "frac_q3_C"]
    },
    "q8": {
        "features": [
            "frac_q8_H", "frac_q8_E", "frac_q8_C", "frac_q8_T",
            "frac_q8_G", "frac_q8_S", "frac_q8_B", "frac_q8_I"
        ]
    },
    "all_features": {
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
        ]
    }
}


# Helper functions
def load_input_table(input_file, sheet_name=None):
    suffix = Path(input_file).suffix.lower()

    if suffix in [".xlsx", ".xls"]:
        return pd.read_excel(input_file, sheet_name=sheet_name, engine="openpyxl")

    if suffix == ".csv":
        return pd.read_csv(input_file)

    raise ValueError(f"Unsupported input format: {suffix}")


def make_dynamic_broken_axis_elbow_plot(k_values, inertia_scaled, plot_path):
    y = np.array(inertia_scaled, dtype=float)

    y_min = float(y.min())
    y_max = float(y.max())
    y_range = y_max - y_min

    if y_range == 0:
        y_range = max(abs(y_max) * 0.05, 1e-6)

    lower_axis_max = max(y_min * 0.02, 1e-4)
    upper_axis_min = y_min - (0.06 * y_range)
    upper_axis_max = y_max + (0.05 * y_range)

    if upper_axis_min <= lower_axis_max:
        upper_axis_min = y_min * 0.95

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=PLOT_FIGSIZE,
        gridspec_kw={"height_ratios": [30, 0.25], "hspace": 0.08}
    )

    ax_top.plot(k_values, inertia_scaled, marker="o", color=LINE_COLOR, linewidth=1.5)

    ax_top.set_ylim(upper_axis_min, upper_axis_max)
    ax_bottom.set_ylim(0.0, lower_axis_max)

    ax_top.set_ylabel("Inertia (× 10⁶)")
    ax_bottom.set_xlabel("k value")

    ax_bottom.set_xticks(k_values)
    ax_bottom.set_yticks([0.0])
    ax_bottom.set_yticklabels(["0"])

    ax_top.yaxis.set_major_locator(MaxNLocator(nbins=6))
    ax_top.yaxis.set_minor_locator(AutoMinorLocator(2))

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

    plt.savefig(plot_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close()


# Load data
print("Loading dataset")

df = load_input_table(INPUT_FILE, SHEET_NAME)
df.columns = df.columns.str.strip()

os.makedirs(OUTPUT_DIR, exist_ok=True)


# Run elbow analysis
for feature_set_name, config in FEATURE_CONFIG.items():
    print(f"\nRunning feature set: {feature_set_name}")

    selected_features = config["features"]

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

    k_values = list(range(K_MIN, K_MAX + 1))
    inertia = []

    print("Running elbow analysis")

    for k in k_values:
        print(f"k = {k}")

        model = MiniBatchKMeans(
            n_clusters=k,
            init="k-means++",
            n_init=10,
            max_iter=300,
            tol=1e-4,
            batch_size=2000,
            random_state=0,
        )

        model.fit(X)
        inertia.append(model.inertia_)

    inertia_scaled = [value / 1e6 for value in inertia]

    inertia_table_path = os.path.join(
        OUTPUT_DIR,
        f"elbow_inertia_{feature_set_name}.xlsx"
    )

    print(f"Saving inertia table: {inertia_table_path}")

    pd.DataFrame({
        "k": k_values,
        "inertia_scaled_millions": inertia_scaled
    }).to_excel(inertia_table_path, index=False)

    plot_path = os.path.join(
        OUTPUT_DIR,
        f"elbow_plot_{feature_set_name}.png"
    )

    print(f"Generating plot: {plot_path}")

    make_dynamic_broken_axis_elbow_plot(
        k_values=k_values,
        inertia_scaled=inertia_scaled,
        plot_path=plot_path
    )

print("\nAll elbow analyses completed.")