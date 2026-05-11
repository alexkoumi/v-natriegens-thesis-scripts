"""
Purpose:
    Summarise NetSurfP-derived structural features for unique and non-unique
    Vibrio natriegens genes.

Description:
    This script reads a cleaned combined NetSurfP results file and compares unique
    and non-unique genes for structural accessibility, disorder probability,
    secondary structure composition, and selected integrated structural features.

    The script generates three Excel workbooks:
        1. rsa_disorder
        2. secondary_structure
        3. combined_features

Input:
    - Cleaned combined NetSurfP results CSV file.

Output:
    - Excel workbook summarising RSA, ASA, and disorder probability.
    - Excel workbook summarising Q3 secondary structure fractions.
    - Excel workbook summarising combined structural features and standardized
      feature differences.

Required Python packages:
    - pandas
    - numpy
    - scipy
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update INPUT_FILE, rsa_disorder, secondary_structure, and
    combined_features before running the script.
"""

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

# File paths
INPUT_FILE = "path/to/combined_netsurfp_results_file.csv"

rsa_disorder = "path/to/output/netsurfp_structural_accessibility_disorder.xlsx"
secondary_structure = "path/to/output/netsurfp_secondary_structure_composition.xlsx"
combined_features = "path/to/output/netsurfp_integrated_feature_comparison.xlsx"


# Load data
df = pd.read_csv(INPUT_FILE, sep=None, engine="python", encoding="utf-8-sig")
df.columns = df.columns.str.strip()


# Helper functions
def find_column(df, candidates, required=True):
    cols_lower = {c.lower(): c for c in df.columns}

    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]

    if required:
        raise KeyError(
            f"Could not find any of these columns: {candidates}\n"
            f"Available columns: {list(df.columns)}"
        )

    return None


def summarize_by_group(df, feature_col, feature_name, group_col="gene_occu"):
    rows = []

    unique_vals = df.loc[df[group_col] == "unique", feature_col].dropna()
    nonunique_vals = df.loc[df[group_col] == "non-unique", feature_col].dropna()

    for group_name, values in [("unique", unique_vals), ("non-unique", nonunique_vals)]:
        q1 = values.quantile(0.25)
        q3 = values.quantile(0.75)
        iqr = q3 - q1

        rows.append({
            "Feature": feature_name,
            "Gene occurrence type": group_name,
            "n": len(values),
            "Median": values.median(),
            "Q1": q1,
            "Q3": q3,
            "IQR": iqr,
            "Minimum": values.min(),
            "Maximum": values.max(),
            "Mean": values.mean(),
            "Standard deviation": values.std(ddof=1)
        })

    if len(unique_vals) > 0 and len(nonunique_vals) > 0:
        u_stat, p_val = mannwhitneyu(unique_vals, nonunique_vals, alternative="two-sided")
    else:
        u_stat, p_val = np.nan, np.nan

    test_row = {
        "Feature": feature_name,
        "Gene occurrence type": "Mann-Whitney U test",
        "n": "",
        "Median": "",
        "Q1": "",
        "Q3": "",
        "IQR": "",
        "Minimum": "",
        "Maximum": "",
        "Mean": "",
        "Standard deviation": "",
        "U statistic": u_stat,
        "p-value": p_val
    }

    out = pd.DataFrame(rows)

    for col in ["U statistic", "p-value"]:
        if col not in out.columns:
            out[col] = ""

    out = pd.concat([out, pd.DataFrame([test_row])], ignore_index=True)

    col_order = [
        "Feature", "Gene occurrence type", "n",
        "Median", "Q1", "Q3", "IQR",
        "Minimum", "Maximum", "Mean", "Standard deviation",
        "U statistic", "p-value"
    ]

    return out[col_order]


def make_boxplot_table(df, feature_col, group_col="gene_occu"):
    unique_vals = df.loc[df[group_col] == "unique", feature_col].dropna().reset_index(drop=True)
    nonunique_vals = df.loc[df[group_col] == "non-unique", feature_col].dropna().reset_index(drop=True)

    return pd.DataFrame({
        "Unique genes": unique_vals,
        "Non-unique genes": nonunique_vals
    })


def add_run_summary(df, group_col="gene_occu"):
    counts = df[group_col].value_counts(dropna=False).to_dict()

    return pd.DataFrame({
        "Metric": [
            "Total rows",
            "Unique genes",
            "Non-unique genes"
        ],
        "Value": [
            len(df),
            counts.get("unique", 0),
            counts.get("non-unique", 0)
        ]
    })


# Detect and standardize core columns
group_col = "gene_occurrence_type"
id_col = "locus_tag"

if group_col not in df.columns:
    raise ValueError(f"Column '{group_col}' not found. Columns available: {list(df.columns)}")

df[group_col] = (
    df[group_col]
    .astype(str)
    .str.strip()
    .str.lower()
    .replace({
        "nonunique": "non-unique",
        "non_unique": "non-unique",
        "non unique": "non-unique"
    })
)

df = df[df[group_col].isin(["unique", "non-unique"])].copy()


# Detect feature columns
rsa_col = find_column(df, ["rsa"])
asa_col = find_column(df, ["asa"])
disorder_col = find_column(df, ["disorder"])

q3_h_col = find_column(df, ["frac_q3_h"])
q3_e_col = find_column(df, ["frac_q3_e"])
q3_c_col = find_column(df, ["frac_q3_c"])

phi_col = find_column(df, ["phi"])
psi_col = find_column(df, ["psi"])

feature_cols = [
    rsa_col, asa_col, disorder_col,
    q3_h_col, q3_e_col, q3_c_col,
    phi_col, psi_col
]

for col in feature_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")


# Rsa and disorder
summary_371 = pd.concat([
    summarize_by_group(df, rsa_col, "RSA", group_col),
    summarize_by_group(df, asa_col, "ASA", group_col),
    summarize_by_group(df, disorder_col, "Disorder probability", group_col)
], ignore_index=True)

rsa_box = make_boxplot_table(df, rsa_col, group_col)
asa_box = make_boxplot_table(df, asa_col, group_col)
disorder_box = make_boxplot_table(df, disorder_col, group_col)

raw_371_cols = [group_col, rsa_col, asa_col, disorder_col]

if id_col is not None:
    raw_371_cols = [id_col] + raw_371_cols

raw_371 = df[raw_371_cols].copy()

rename_371 = {
    group_col: "Gene occurrence type",
    rsa_col: "RSA",
    asa_col: "ASA",
    disorder_col: "Disorder probability"
}

if id_col is not None:
    rename_371[id_col] = "ID"

raw_371.rename(columns=rename_371, inplace=True)

run_summary = add_run_summary(df, group_col)

with pd.ExcelWriter(rsa_disorder, engine="openpyxl") as writer:
    summary_371.to_excel(writer, sheet_name="Summary_statistics", index=False)
    rsa_box.to_excel(writer, sheet_name="RSA_boxplot_data", index=False)
    asa_box.to_excel(writer, sheet_name="ASA_boxplot_data", index=False)
    disorder_box.to_excel(writer, sheet_name="Disorder_boxplot_data", index=False)
    raw_371.to_excel(writer, sheet_name="Raw_long_format", index=False)
    run_summary.to_excel(writer, sheet_name="Run_summary", index=False)


# Secondary structure
summary_372 = pd.concat([
    summarize_by_group(df, q3_h_col, "Helix fraction (Q3_H)", group_col),
    summarize_by_group(df, q3_e_col, "Strand fraction (Q3_E)", group_col),
    summarize_by_group(df, q3_c_col, "Coil fraction (Q3_C)", group_col)
], ignore_index=True)

helix_box = make_boxplot_table(df, q3_h_col, group_col)
strand_box = make_boxplot_table(df, q3_e_col, group_col)
coil_box = make_boxplot_table(df, q3_c_col, group_col)

raw_372_cols = [group_col, q3_h_col, q3_e_col, q3_c_col]

if id_col is not None:
    raw_372_cols = [id_col] + raw_372_cols

raw_372 = df[raw_372_cols].copy()

rename_372 = {
    group_col: "Gene occurrence type",
    q3_h_col: "Helix fraction (Q3_H)",
    q3_e_col: "Strand fraction (Q3_E)",
    q3_c_col: "Coil fraction (Q3_C)"
}

if id_col is not None:
    rename_372[id_col] = "ID"

raw_372.rename(columns=rename_372, inplace=True)

with pd.ExcelWriter(secondary_structure, engine="openpyxl") as writer:
    summary_372.to_excel(writer, sheet_name="Summary_statistics", index=False)
    helix_box.to_excel(writer, sheet_name="Helix_boxplot_data", index=False)
    strand_box.to_excel(writer, sheet_name="Strand_boxplot_data", index=False)
    coil_box.to_excel(writer, sheet_name="Coil_boxplot_data", index=False)
    raw_372.to_excel(writer, sheet_name="Raw_long_format", index=False)
    run_summary.to_excel(writer, sheet_name="Run_summary", index=False)


# Combined features
integrated_features = {
    "RSA": rsa_col,
    "ASA": asa_col,
    "Disorder": disorder_col,
    "Helix fraction (Q3_H)": q3_h_col,
    "Strand fraction (Q3_E)": q3_e_col,
    "Coil fraction (Q3_C)": q3_c_col,
    "Phi": phi_col,
    "Psi": psi_col
}

heatmap_rows = []

for feature_name, col in integrated_features.items():
    feature_values = df[col].dropna()

    mean_all = feature_values.mean()
    std_all = feature_values.std()

    z_values = (df[col] - mean_all) / std_all

    z_unique = z_values[df[group_col] == "unique"].mean()
    z_nonunique = z_values[df[group_col] == "non-unique"].mean()

    difference = z_unique - z_nonunique

    heatmap_rows.append({
        "Feature": feature_name,
        "Standardized difference": difference
    })

heatmap_df = pd.DataFrame(heatmap_rows).set_index("Feature")

summary_373 = pd.concat([
    summarize_by_group(df, rsa_col, "RSA", group_col),
    summarize_by_group(df, asa_col, "ASA", group_col),
    summarize_by_group(df, disorder_col, "Disorder", group_col),
    summarize_by_group(df, q3_h_col, "Helix fraction (Q3_H)", group_col),
    summarize_by_group(df, q3_e_col, "Strand fraction (Q3_E)", group_col),
    summarize_by_group(df, q3_c_col, "Coil fraction (Q3_C)", group_col),
    summarize_by_group(df, phi_col, "Phi", group_col),
    summarize_by_group(df, psi_col, "Psi", group_col)
], ignore_index=True)

raw_373_cols = [group_col] + feature_cols

if id_col is not None:
    raw_373_cols = [id_col] + raw_373_cols

raw_373 = df[raw_373_cols].copy()

rename_373 = {
    group_col: "Gene occurrence type",
    rsa_col: "RSA",
    asa_col: "ASA",
    disorder_col: "Disorder",
    q3_h_col: "Helix fraction (Q3_H)",
    q3_e_col: "Strand fraction (Q3_E)",
    q3_c_col: "Coil fraction (Q3_C)",
    phi_col: "Phi",
    psi_col: "Psi"
}

if id_col is not None:
    rename_373[id_col] = "ID"

raw_373.rename(columns=rename_373, inplace=True)

with pd.ExcelWriter(combined_features, engine="openpyxl") as writer:
    heatmap_df.to_excel(writer, sheet_name="Heatmap_matrix")
    summary_373.to_excel(writer, sheet_name="Supporting_statistics", index=False)
    raw_373.to_excel(writer, sheet_name="Raw_long_format", index=False)
    run_summary.to_excel(writer, sheet_name="Run_summary", index=False)


# Final report
print("Done.")
print(f"Created: {rsa_disorder}")
print(f"Created: {secondary_structure}")
print(f"Created: {combined_features}")
print("\nDetected columns:")
print(f"Group column: {group_col}")

if id_col is not None:
    print(f"ID column: {id_col}")

print(f"RSA: {rsa_col}")
print(f"ASA: {asa_col}")
print(f"Disorder: {disorder_col}")
print(f"Q3_H: {q3_h_col}")
print(f"Q3_E: {q3_e_col}")
print(f"Q3_C: {q3_c_col}")
print(f"Phi: {phi_col}")
print(f"Psi: {psi_col}")