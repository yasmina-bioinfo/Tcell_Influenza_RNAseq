# T cell transcriptomic profiling in Influenza infection (scRNA-seq)

**Dataset:** GEO accession **GSE149689**  
**Study title:** Single-cell RNA sequencing of PBMCs from patients with Influenza, COVID-19, and healthy controls  
**Organism:** *Homo sapiens*  
**Technology:** 10x Genomics single-cell RNA-seq  

---

## Project overview

This project implements a complete and reproducible single-cell RNA-seq analysis pipeline
to explore transcriptional programs in human T cells (CD4 and CD8) using the public dataset
**GSE149689**.

The initial objective was to investigate host immune responses to viral infection.
During the analysis, dataset-specific limitations related to metadata availability were
identified, which guided the scope and interpretation of downstream analyses.

This repository focuses on methodological rigor, reproducibility, and critical evaluation
of dataset suitability with respect to the biological question.

---

## Analysis pipeline

The analysis was performed using a stepwise and reproducible pipeline implemented in Python,
with each script corresponding to a clearly defined analytical stage.

### 0. Metadata verification (`00_check_metadata.py`)
Sample-level metadata were explicitly inspected to assess whether biological replicates
or condition labels (e.g. infected vs control) were available at the single-cell level.

No exploitable columns related to patient, donor, sample, or condition were found,
which directly informed the downstream analytical strategy.

---

### 1. Data loading (`01_load_data.py`)
Raw 10x Genomics files (`matrix.mtx`, `barcodes.tsv`, `features.tsv`) were loaded and
encapsulated into an AnnData object without any filtering or normalization.
This step ensures a faithful representation of the original data.

**Output:**
- `adata_raw.h5ad`

---

### 2. Quality control and filtering (`02_qc_filtering.py`)
Quality control metrics were computed at the single-cell level, including:
- number of detected genes per cell,
- total UMI counts,
- percentage of mitochondrial gene expression.

Cells with low gene detection (< 200 genes) or high mitochondrial content (> 20%)
were removed to exclude low-quality or dying cells.

**Output:**
- `adata_qc.h5ad`

---

### 3. T cell selection and CD4/CD8 separation (`03_select_tcells_cd4_cd8.py`)
After QC, cells were biologically annotated using canonical marker-based gene signatures.
T cells were identified first, then separated into CD4-like and CD8-like populations
based on relative expression scores.

Normalization and log-transformation were applied **only** for marker scoring purposes,
not for downstream count-based analyses.

**Outputs:**
- `adata_tcells.h5ad`
- `adata_cd4.h5ad`
- `adata_cd8.h5ad`

---

### 4. Pseudo-bulk construction (`04_pseudobulk_cd4_vs_cd8.py`)
Raw gene expression counts were aggregated across all CD4 and CD8 T cells separately
to generate bulk-like expression profiles.

This pseudo-bulk representation enables a direct comparison of transcriptional
programs between CD4 and CD8 T cell populations.

**Output:**
- `pseudobulk_CD4_vs_CD8_counts.tsv`

---

## Results summary

After quality control and biological selection, a large number of high-quality T cells
were retained and separated into CD4-like and CD8-like populations.

Marker-based scoring confirmed biologically coherent transcriptional profiles:
CD8-like T cells displayed higher expression of cytotoxic and effector-associated genes,
while CD4-like T cells were enriched for helper and memory-associated markers.
These results are consistent with established immunological knowledge and validate
the robustness of the cell selection strategy.

A pseudo-bulk representation was constructed by aggregating raw gene counts across
all CD4 and CD8 T cells separately. This produced bulk-like expression profiles
allowing descriptive comparison of transcriptional programs between CD4 and CD8 T cells.

Due to dataset constraints, results are interpreted descriptively rather than through
formal differential expression statistics.

---

## Limitations of the dataset

Despite the successful implementation of a complete single-cell RNA-seq analysis pipeline,
this study is subject to important dataset-specific limitations that constrain the scope
of biological interpretation.

First, the available 10x Genomics files do not contain exploitable sample-level metadata.
In particular, there is no information linking individual cells to a specific patient,
donor, or experimental condition (e.g. infected vs control). As a result, biological
replicates cannot be defined, and inter-individual variability cannot be assessed.

Second, although the original study includes both infected individuals and healthy controls,
this information is not accessible at the single-cell level in the processed data.
Consequently, a direct comparison of transcriptional responses between infected and
non-infected conditions cannot be performed using the current dataset.

Third, the absence of biological replicates limits the analysis to descriptive comparisons.
While pseudo-bulk profiles can be constructed to compare CD4 and CD8 T cell populations,
formal differential expression analysis with statistical inference is not supported.

These limitations reflect constraints inherent to the structure and metadata availability
of the public dataset rather than shortcomings of the analysis pipeline.

---

## Conclusion and lessons learned

This project demonstrates the implementation of a robust and reproducible scRNA-seq
analysis pipeline, from raw data loading to biologically meaningful cell population
selection and pseudo-bulk construction.

Importantly, it highlights the necessity of critically evaluating dataset design and
metadata availability before attempting statistical inference. When the data were found
to be insufficient to address the initial biological question (infected vs control),
the analysis scope was deliberately limited to avoid overinterpretation.

This work therefore emphasizes methodological rigor, transparency, and scientific
decision-making as essential components of computational biology research.

---

## Repository structure

```text
Tcell_Influenza_RNAseq/
├── data_GSE149689/
│   ├── raw/
│   └── processed/
├── Scripts/
│   ├── 00_check_metadata.py
│   ├── 01_load_data.py
│   ├── 02_qc_filtering.py
│   ├── 03_select_tcells_cd4_cd8.py
│   └── 04_pseudobulk_cd4_vs_cd8.py
├── results/
│   └── pseudobulk/
│       └── pseudobulk_CD4_vs_CD8_counts.tsv
└── README.md
