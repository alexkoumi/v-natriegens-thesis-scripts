"""
Purpose:
    Run Ward hierarchical clustering on selected NetSurfP-derived features.

Description:
    This script reads a combined NetSurfP results file, standardizes selected
    numerical features, constructs a k-nearest-neighbour connectivity graph, and
    applies Ward hierarchical clustering using a selected number of clusters.

Input:
    - Combined NetSurfP results CSV file.

Output:
    - Excel file containing the original input table with an added Ward cluster
      assignment column.

Required Python packages:
    - pandas
    - scikit-learn
    - openpyxl

Notes:
    Local file paths used during the original analysis were replaced with placeholder
    paths. Update input_file, output_file, features, and k before running.
"""

import os
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import kneighbors_graph

k = "user-defined k value"
cluster_column = f"ward_k={k}"

# Input
input_file = "path/to/netsurfp_results_combined.csv"

df = pd.read_csv(input_file)
df.columns = df.columns.str.strip()

# Features
features = ["user-defined feature set"]

# Standardise
scaler = StandardScaler()
X = scaler.fit_transform(df[features])

# Connectivity graph
knn_graph = kneighbors_graph(
    X,
    n_neighbors=10,
    mode="connectivity",
    include_self=False
)

# Ward hierarchical clustering
model = AgglomerativeClustering(
    n_clusters=k,
    metric="euclidean",
    linkage="ward",
    connectivity=knn_graph
)

labels = model.fit_predict(X)
df[cluster_column] = labels

# Output
output_file = f"path/to/output/"
os.makedirs(os.path.dirname(output_file), exist_ok=True)

df.to_excel(output_file, index=False)

print(f"Ward clustering complete, saved to {output_file}")