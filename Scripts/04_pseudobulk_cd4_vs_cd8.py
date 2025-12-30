#!/usr/bin/env python3
"""
04_pseudobulk_cd4_vs_cd8.py

Create a pseudo-bulk expression matrix by summing raw counts
across CD4 and CD8 T cells.
"""

import os
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp

# -----------------------------
# Paths
# -----------------------------
IN_CD4 = "data_GSE149689/processed/adata_cd4.h5ad"
IN_CD8 = "data_GSE149689/processed/adata_cd8.h5ad"

OUT_DIR = "results/pseudobulk"
OUT_FILE = os.path.join(OUT_DIR, "pseudobulk_CD4_vs_CD8_counts.tsv")

# -----------------------------
# Load data
# -----------------------------
print("Loading CD4 and CD8 AnnData objects...")
adata_cd4 = sc.read_h5ad(IN_CD4)
adata_cd8 = sc.read_h5ad(IN_CD8)

# -----------------------------
# Sum raw counts
# -----------------------------
def sum_counts(adata):
    X = adata.X
    if sp.issparse(X):
        return np.asarray(X.sum(axis=0)).ravel()
    else:
        return X.sum(axis=0)

cd4_sum = sum_counts(adata_cd4)
cd8_sum = sum_counts(adata_cd8)

# -----------------------------
# Build output table
# -----------------------------
genes = (
    adata_cd4.var["gene_id"].astype(str)
    if "gene_id" in adata_cd4.var.columns
    else adata_cd4.var_names.astype(str)
)

gene_names = (
    adata_cd4.var["gene_name"].astype(str)
    if "gene_name" in adata_cd4.var.columns
    else adata_cd4.var_names.astype(str)
)

df = pd.DataFrame({
    "gene_id": genes,
    "gene_name": gene_names,
    "CD4": cd4_sum.astype(int),
    "CD8": cd8_sum.astype(int),
})

print("Total counts:")
print(df[["CD4", "CD8"]].sum())

# -----------------------------
# Save
# -----------------------------
os.makedirs(OUT_DIR, exist_ok=True)
df.to_csv(OUT_FILE, sep="\t", index=False)

print(f"Saved pseudo-bulk matrix to: {OUT_FILE}")

#python scripts/04_pseudobulk_cd4_vs_cd8.py
