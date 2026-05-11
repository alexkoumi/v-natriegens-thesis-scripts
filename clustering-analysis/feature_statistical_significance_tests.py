"""
Purpose:
    Test whether selected gene-level features differ between unique and non-unique
    Vibrio natriegens genes.

Description:
    This script reads a combined feature table and compares unique genes against
    non-unique genes for amino acid composition, RSA, ASA, disorder probability,
    Q3/Q8 secondary structure fractions, and phi/psi angle features.

    For each feature, the script performs:
        - Welch's t-test
        - Mann-Whitney U test
        - Kolmogorov-Smirnov test

    False-discovery-rate correction is applied to each set of p-values. The script
    also writes summary statistics and generates violin and box plots for each
    feature.

Input:
    - Excel or CSV file containing gene occurrence labels and numerical feature columns.

Output:
    - CSV file containing statistical test results and FDR-adjusted p-values.
    - Text file containing summary statistics for each feature.
    - Violin plot PNG files.
    - Box plot PNG files.

Required Python packages:
    - pandas
    - numpy
    - matplotlib
    - seaborn
    - scipy
    - statsmodels
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update INPUT_FILE, SHEET_NAME, GROUP_COLUMN, and OUTPUT_DIR before running.
"""

import re
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind, mannwhitneyu, ks_2samp
from statsmodels.stats.multitest import multipletests

# Paths
INPUT_FILE = "path/to/combined_feature_table.xlsx"
SHEET_NAME = "Sheet1"
GROUP_COLUMN = "Gene occurrence type"

OUTPUT_DIR = "path/to/output/statistical_tests"
OUTPUT_STATS_FILE = "feature_significance_results.csv"
SUMMARY_FILE = "feature_summary_statistics.txt"
VIOLIN_PLOT_DIR_NAME = "violin_plots"
BOX_PLOT_DIR_NAME = "box_plots"

UNIQUE_LABEL = "unique"
NON_UNIQUE_LABELS = ["non-unique", "nonunique", "non unique", "non_unique"]

# Feature columns to test
# Replace these placeholders with the exact feature column names in your input file.
FEATURES = [
    "feature_column_1",
    "feature_column_2",
    "feature_column_3",
]

# Y-axis labels for plots
# Replace these placeholders with readable labels for the features listed above.
FEATURE_LABELS = {
    "feature_column_1": "Feature 1",
    "feature_column_2": "Feature 2",
    "feature_column_3": "Feature 3",
}

PLOT_DPI = 300
PLOT_FIGSIZE = (6, 4)
GREY_COLOR = "0.6"
X_AXIS_LABEL = "Gene occurrence type"
X_TICK_LABELS = ["Strain-specific", "Non-strain-specific"]


# Helper functions
def safe_name(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("_")


def load_input_table(input_file, sheet_name=None):
    suffix = Path(input_file).suffix.lower()

    if suffix in [".xlsx", ".xls"]:
        return pd.read_excel(input_file, sheet_name=sheet_name, engine="openpyxl")

    if suffix == ".csv":
        return pd.read_csv(input_file)

    raise ValueError(f"Unsupported input format: {suffix}")


def standardize_group_value(value):
    if pd.isna(value):
        return None

    value_str = str(value).strip().lower()
    normalized = re.sub(r"[\s_]+", "-", value_str)
    normalized = re.sub(r"-+", "-", normalized)

    if normalized == UNIQUE_LABEL:
        return "unique"

    valid_non_unique = set()

    for label in NON_UNIQUE_LABELS:
        cleaned_label = re.sub(r"[\s_]+", "-", label.strip().lower())
        valid_non_unique.add(cleaned_label)

    if normalized in valid_non_unique:
        return "non-unique"

    return None


def summarize_series(series):
    series = pd.to_numeric(series, errors="coerce").dropna()

    if len(series) == 0:
        return {
            "n": 0,
            "mean": np.nan,
            "std": np.nan,
            "median": np.nan,
            "q1": np.nan,
            "q3": np.nan,
        }

    if len(series) > 1:
        standard_deviation = float(series.std(ddof=1))
    else:
        standard_deviation = np.nan

    return {
        "n": int(series.shape[0]),
        "mean": float(series.mean()),
        "std": standard_deviation,
        "median": float(series.median()),
        "q1": float(series.quantile(0.25)),
        "q3": float(series.quantile(0.75)),
    }


def safe_t_test(v1, v2):
    if len(v1) == 0 or len(v2) == 0:
        return np.nan, np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        stat, p_value = ttest_ind(v1, v2, equal_var=False, nan_policy="omit")

    return float(stat), float(p_value)


def safe_mann_whitney(v1, v2):
    if len(v1) == 0 or len(v2) == 0:
        return np.nan, np.nan

    try:
        stat, p_value = mannwhitneyu(v1, v2, alternative="two-sided")
        return float(stat), float(p_value)
    except ValueError:
        return np.nan, np.nan


def safe_ks(v1, v2):
    if len(v1) == 0 or len(v2) == 0:
        return np.nan, np.nan

    stat, p_value = ks_2samp(v1, v2)

    return float(stat), float(p_value)


def apply_fdr(df, p_value_column, output_column):
    pvals = pd.to_numeric(df[p_value_column], errors="coerce")
    valid_mask = pvals.notna()
    adjusted = np.full(len(df), np.nan)

    if valid_mask.any():
        adjusted[valid_mask] = multipletests(pvals[valid_mask], method="fdr_bh")[1]

    df[output_column] = adjusted

    return df


# Output directories
output_dir = Path(OUTPUT_DIR)
violin_dir = output_dir / VIOLIN_PLOT_DIR_NAME
box_dir = output_dir / BOX_PLOT_DIR_NAME

output_dir.mkdir(parents=True, exist_ok=True)
violin_dir.mkdir(parents=True, exist_ok=True)
box_dir.mkdir(parents=True, exist_ok=True)

output_stats_path = output_dir / OUTPUT_STATS_FILE
summary_path = output_dir / SUMMARY_FILE


# Load data
print(f"Loading file: {INPUT_FILE}")

df = load_input_table(INPUT_FILE, SHEET_NAME)
df.columns = df.columns.str.strip()

if GROUP_COLUMN not in df.columns:
    raise KeyError(f"Group column not found: {GROUP_COLUMN}")

print(
    "Raw group labels found in input:",
    sorted(df[GROUP_COLUMN].dropna().astype(str).str.strip().unique().tolist())
)

missing_features = []

for feature in FEATURES:
    if feature not in df.columns:
        missing_features.append(feature)

if missing_features:
    raise KeyError(
        "The following feature columns were not found in the input file:\n- "
        + "\n- ".join(missing_features)
    )

analysis_df = df.copy()
analysis_df["__group__"] = analysis_df[GROUP_COLUMN].apply(standardize_group_value)
analysis_df = analysis_df[analysis_df["__group__"].isin(["unique", "non-unique"])].copy()

if analysis_df.empty:
    raise ValueError("No rows remained after standardizing group labels.")

g1 = analysis_df[analysis_df["__group__"] == "unique"].copy()
g2 = analysis_df[analysis_df["__group__"] == "non-unique"].copy()

print(f"Unique genes: {len(g1)}")
print(f"Non-unique genes: {len(g2)}")


# Statistical tests
results = []

for feature in FEATURES:
    print(f"Testing feature: {feature}")

    v1 = pd.to_numeric(g1[feature], errors="coerce").dropna()
    v2 = pd.to_numeric(g2[feature], errors="coerce").dropna()

    t_stat, t_p = safe_t_test(v1, v2)
    mw_stat, mw_p = safe_mann_whitney(v1, v2)
    ks_stat, ks_p = safe_ks(v1, v2)

    results.append({
        "feature": feature,
        "t_test_statistic": t_stat,
        "t_test_p_value": t_p,
        "mann_whitney_u_statistic": mw_stat,
        "mann_whitney_p_value": mw_p,
        "kolmogorov_smirnov_statistic": ks_stat,
        "kolmogorov_smirnov_p_value": ks_p,
    })

stats_df = pd.DataFrame(results)

stats_df = apply_fdr(stats_df, "t_test_p_value", "t_test_fdr")
stats_df = apply_fdr(stats_df, "mann_whitney_p_value", "mann_whitney_fdr")
stats_df = apply_fdr(stats_df, "kolmogorov_smirnov_p_value", "kolmogorov_smirnov_fdr")

stats_df.to_csv(output_stats_path, index=False)

print(f"Saved significance results to: {output_stats_path}")


# Generate plots
sns.set_style("whitegrid")

plot_df = analysis_df.copy()
plot_df["__group__"] = pd.Categorical(
    plot_df["__group__"],
    categories=["unique", "non-unique"],
    ordered=True,
)

for feature in FEATURES:
    print(f"Plotting feature: {feature}")

    y_label = FEATURE_LABELS.get(feature, feature)
    file_stub = safe_name(feature)

    plt.figure(figsize=PLOT_FIGSIZE)

    ax = sns.violinplot(
        data=plot_df,
        x="__group__",
        y=feature,
        inner="quartile",
        color=GREY_COLOR,
        cut=0,
        linewidth=1,
    )

    ax.set_xlabel(X_AXIS_LABEL)
    ax.set_ylabel(y_label)
    ax.set_xticklabels(X_TICK_LABELS)

    plt.tight_layout()
    plt.savefig(violin_dir / f"{file_stub}_violin.png", dpi=PLOT_DPI, bbox_inches="tight")
    plt.close()

    plt.figure(figsize=PLOT_FIGSIZE)

    ax = sns.boxplot(
        data=plot_df,
        x="__group__",
        y=feature,
        color=GREY_COLOR,
        showfliers=False,
        width=0.5,
    )

    ax.set_xlabel(X_AXIS_LABEL)
    ax.set_ylabel(y_label)
    ax.set_xticklabels(X_TICK_LABELS)

    plt.tight_layout()
    plt.savefig(box_dir / f"{file_stub}_box.png", dpi=PLOT_DPI, bbox_inches="tight")
    plt.close()


# Summary text file
with open(summary_path, "w", encoding="utf-8") as handle:
    handle.write("Median, IQR, mean and standard deviation values for each feature\n")
    handle.write("Computed separately for unique and non-unique genes\n\n")

    for feature in FEATURES:
        unique_vals = pd.to_numeric(g1[feature], errors="coerce").dropna()
        non_unique_vals = pd.to_numeric(g2[feature], errors="coerce").dropna()

        unique_summary = summarize_series(unique_vals)
        non_unique_summary = summarize_series(non_unique_vals)

        handle.write(f"=== {feature} ===\n")

        handle.write("unique:\n")
        handle.write(f" n: {unique_summary['n']}\n")
        handle.write(f" mean: {unique_summary['mean']}\n")
        handle.write(f" standard deviation: {unique_summary['std']}\n")
        handle.write(f" median: {unique_summary['median']}\n")
        handle.write(f" Q1 (25%): {unique_summary['q1']}\n")
        handle.write(f" Q3 (75%): {unique_summary['q3']}\n")

        handle.write("non-unique:\n")
        handle.write(f" n: {non_unique_summary['n']}\n")
        handle.write(f" mean: {non_unique_summary['mean']}\n")
        handle.write(f" standard deviation: {non_unique_summary['std']}\n")
        handle.write(f" median: {non_unique_summary['median']}\n")
        handle.write(f" Q1 (25%): {non_unique_summary['q1']}\n")
        handle.write(f" Q3 (75%): {non_unique_summary['q3']}\n\n")

print("All processes completed")
