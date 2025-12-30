#!/usr/bin/env python3
"""
02_qc_filtering.py

Compute quality control metrics and apply minimal QC filtering
to remove low-quality cells from the raw AnnData object.
"""

import os
import scanpy as sc

# -----------------------------
# Paths
# -----------------------------
IN_FILE = "data_GSE149689/processed/adata_raw.h5ad"
OUT_DIR = "data_GSE149689/processed"
OUT_FILE = os.path.join(OUT_DIR, "adata_qc.h5ad")

# -----------------------------
# Load data
# -----------------------------
print("Loading raw AnnData...")
adata = sc.read_h5ad(IN_FILE)
print("Initial shape:", adata.shape)

# -----------------------------
# QC metrics
# -----------------------------
print("Computing QC metrics...")

# Identify mitochondrial genes
adata.var["mt"] = adata.var["gene_name"].astype(str).str.startswith("MT-")

# Compute QC metrics
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    inplace=True
)

print("QC summary:")
print("n_genes_by_counts (min / median / max):",
      int(adata.obs["n_genes_by_counts"].min()),
      int(adata.obs["n_genes_by_counts"].median()),
      int(adata.obs["n_genes_by_counts"].max()))

print("pct_counts_mt (min / median / max):",
      round(float(adata.obs["pct_counts_mt"].min()), 2),
      round(float(adata.obs["pct_counts_mt"].median()), 2),
      round(float(adata.obs["pct_counts_mt"].max()), 2))

# -----------------------------
# Minimal filtering
# -----------------------------
print("Applying minimal QC filtering...")

adata_qc = adata[
    (adata.obs["n_genes_by_counts"] >= 200) &
    (adata.obs["pct_counts_mt"] <= 20)
].copy()

print("After QC shape:", adata_qc.shape)
print("Cells removed:", adata.n_obs - adata_qc.n_obs)

# -----------------------------
# Save
# -----------------------------
os.makedirs(OUT_DIR, exist_ok=True)
adata_qc.write_h5ad(OUT_FILE)

print(f"Saved QC-filtered AnnData to: {OUT_FILE}")

#python scripts/02_qc_filtering.py
