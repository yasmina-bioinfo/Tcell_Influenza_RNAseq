#!/usr/bin/env python3
"""
00_check_metadata.py

Check whether sample-level or patient-level metadata
are available in the AnnData object.
"""

import scanpy as sc

adata = sc.read_h5ad("data_GSE149689/processed/adata_qc.h5ad")

keywords = ["patient", "donor", "sample", "condition", "infect", "control", "batch", "orig"]

matching_cols = [
    col for col in adata.obs.columns
    if any(k in col.lower() for k in keywords)
]

print("Metadata columns related to samples / conditions:")
print(matching_cols if matching_cols else "None found")

# Exécution : python scripts/00_check_metadata.py
#Output : Metadata columns related to samples / conditions: None found
