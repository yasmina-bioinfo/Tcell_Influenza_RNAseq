#!/usr/bin/env python3
"""
03_select_tcells_cd4_cd8.py

Select T cells and separate CD4-like and CD8-like populations
using marker-based gene signature scoring.
"""

import os
import scanpy as sc

# -----------------------------
# Paths
# -----------------------------
IN_FILE = "data_GSE149689/processed/adata_qc.h5ad"
OUT_DIR = "data_GSE149689/processed"

OUT_T  = os.path.join(OUT_DIR, "adata_tcells.h5ad")
OUT_CD4 = os.path.join(OUT_DIR, "adata_cd4.h5ad")
OUT_CD8 = os.path.join(OUT_DIR, "adata_cd8.h5ad")

# -----------------------------
# Load data
# -----------------------------
print("Loading QC-filtered AnnData...")
adata = sc.read_h5ad(IN_FILE)
print("Initial shape:", adata.shape)

# -----------------------------
# Prepare gene symbols for scoring
# -----------------------------
adata.var_names = adata.var["gene_name"].astype(str)
adata.var_names_make_unique()

# -----------------------------
# Normalization (for scoring only)
# -----------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# -----------------------------
# Marker genes
# -----------------------------
t_markers = ["TRAC", "CD3D", "CD3E"]
cd4_markers = ["CD4", "IL7R", "CCR7", "LTB", "MAL"]
cd8_markers = ["CD8A", "CD8B", "NKG7", "GZMB", "PRF1", "CTSW"]

genes_present = set(adata.var_names)

t_markers = [g for g in t_markers if g in genes_present]
cd4_markers = [g for g in cd4_markers if g in genes_present]
cd8_markers = [g for g in cd8_markers if g in genes_present]

print("T markers present:", t_markers)
print("CD4 markers present:", cd4_markers)
print("CD8 markers present:", cd8_markers)

# -----------------------------
# Scoring
# -----------------------------
sc.tl.score_genes(adata, t_markers, score_name="score_T")
sc.tl.score_genes(adata, cd4_markers, score_name="score_CD4")
sc.tl.score_genes(adata, cd8_markers, score_name="score_CD8")

# -----------------------------
# Selection
# -----------------------------
adata_tcells = adata[adata.obs["score_T"] > 0].copy()
adata_cd4 = adata_tcells[adata_tcells.obs["score_CD4"] > adata_tcells.obs["score_CD8"]].copy()
adata_cd8 = adata_tcells[adata_tcells.obs["score_CD8"] >= adata_tcells.obs["score_CD4"]].copy()

print("T cells shape:", adata_tcells.shape)
print("CD4-like shape:", adata_cd4.shape)
print("CD8-like shape:", adata_cd8.shape)

# -----------------------------
# Make objects safe for h5ad saving
# -----------------------------
def make_h5ad_safe(a):
    a.obs.index = a.obs.index.astype(str)
    a.var.index = a.var.index.astype(str)

    for col in a.obs.columns:
        if a.obs[col].dtype == "object":
            a.obs[col] = a.obs[col].astype(str)

    for col in a.var.columns:
        if a.var[col].dtype == "object":
            a.var[col] = a.var[col].astype(str)

    a.var.index.name = None
    return a

adata_tcells = make_h5ad_safe(adata_tcells)
adata_cd4 = make_h5ad_safe(adata_cd4)
adata_cd8 = make_h5ad_safe(adata_cd8)

# -----------------------------
# Save
# -----------------------------
os.makedirs(OUT_DIR, exist_ok=True)

adata_tcells.write_h5ad(OUT_T)
adata_cd4.write_h5ad(OUT_CD4)
adata_cd8.write_h5ad(OUT_CD8)

print("Saved:")
print(OUT_T)
print(OUT_CD4)
print(OUT_CD8)

#python scripts/03_select_tcells_cd4_cd8.py
