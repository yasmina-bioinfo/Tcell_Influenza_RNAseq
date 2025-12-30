#!/usr/bin/env python3
"""
01_load_data.py

Load raw 10X Genomics files (matrix.mtx, barcodes.tsv, features.tsv),
create an AnnData object, and save it as adata_raw.h5ad.

This script performs NO filtering and NO normalization.
"""

import os
import scanpy as sc
import pandas as pd
from scipy.io import mmread
import scipy.sparse as sp

# -----------------------------
# Paths
# -----------------------------
RAW_DIR = "data_GSE149689/raw/GSE149689"
OUT_DIR = "data_GSE149689/processed"
OUT_FILE = os.path.join(OUT_DIR, "adata_raw.h5ad")

# -----------------------------
# Load raw data
# -----------------------------
print("Loading raw 10X files...")

X = mmread(os.path.join(RAW_DIR, "matrix.mtx")).T
X = X.tocsr()  # required by AnnData

obs = pd.read_csv(
    os.path.join(RAW_DIR, "barcodes.tsv"),
    header=None,
    sep="\t"
)
obs.columns = ["barcode"]

var = pd.read_csv(
    os.path.join(RAW_DIR, "features.tsv"),
    header=None,
    sep="\t"
)
var.columns = ["gene_id", "gene_name", "feature_type"]

# -----------------------------
# Create AnnData object
# -----------------------------
adata = sc.AnnData(X=X, obs=obs, var=var)

print("AnnData object created:")
print(adata)

# -----------------------------
# Save
# -----------------------------
os.makedirs(OUT_DIR, exist_ok=True)
adata.write_h5ad(OUT_FILE)

print(f"Saved raw AnnData to: {OUT_FILE}")

#python scripts/01_load_data.py
