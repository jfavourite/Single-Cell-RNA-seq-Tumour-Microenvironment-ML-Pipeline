# Single-Cell RNA-seq Tumour Microenvironment Pipeline — Breast Cancer

A complete single-cell analysis pipeline applied to the Wu et al. (2021) breast cancer
single-cell atlas (GSE176078), combining batch-corrected representation learning (scVI),
unsupervised clustering, and a translational machine learning layer linking tumour
microenvironment (TME) composition to immune phenotype.

## Dataset

- **Source:** Wu et al. 2021, *Nature Genetics* — [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078)
- **26 breast cancer patients** spanning ER+, HER2+, and TNBC subtypes
- **~100,000 cells**, ~29,700 genes, with author-curated cell type annotations at three
  levels of granularity (`celltype_major`, `celltype_minor`, `celltype_subset`)

## Pipeline overview

1. **Data loading** — 26 per-patient count matrices loaded and concatenated into a single
   `AnnData` object (100,064 cells × 29,733 genes)
2. **Quality control** — filtered on gene counts (200–6000) and mitochondrial percentage
   (<20%), retaining 98,593 cells (98.5%)
3. **Normalisation & HVG selection** — library-size normalisation, log1p transform, and
   batch-aware selection of 3,000 highly variable genes
4. **PCA baseline** — linear dimensionality reduction; UMAP reveals strong
   patient-specific batch effects, particularly within the malignant cell population
5. **scVI** — a variational autoencoder trained on raw counts (30,000-cell subsample,
   `n_latent=30`, negative binomial likelihood, `batch_key='sample'`) to learn a
   batch-corrected latent representation
6. **Leiden clustering** — 17 clusters identified in the scVI latent space, achieving
   **>93% concordance** with author-curated cell type annotations across all major TME
   populations
7. **TME composition** — per-patient cell-type fraction matrix (9 major populations × 26
   patients)
8. **Subtype classification** — tested whether TME composition predicts ER+/HER2+/TNBC
   subtype (LOOCV)
9. **Immune hot/cold classification** — reframed toward predicting CD8+ T-cell
   infiltration status from the non-T-cell TME composition (LOOCV)
10. **SHAP explainability** — identified and biologically interpreted the dominant
    predictors of immune phenotype
11. **Benchmarking** — silhouette score comparison of PCA vs scVI embeddings, with
    discussion of metric limitations for batch-corrected data

## Key results

### Clustering validation

Unsupervised Leiden clustering on the scVI latent space recovered all 9 major TME
populations (T-cells, Cancer Epithelial, Myeloid, CAFs, Endothelial, PVL, B-cells,
Plasmablasts, Normal Epithelial) with >93% purity per cluster, including **four distinct
T-cell subclusters** discovered without supervision — consistent with known CD8+, CD4+,
regulatory, and exhausted T-cell states.

### Subtype classification — a negative result, correctly interpreted

Predicting molecular subtype (ER+/HER2+/TNBC) from TME composition alone achieved only
**50% LOOCV accuracy** (baseline 42%), with 0% recall for HER2+ (n=5). This is interpreted
as a biological limitation rather than a methodological failure: molecular subtype is
defined by receptor expression *within malignant cells*, not by the surrounding
immune/stromal composition. This null result motivated reframing toward immune phenotype.

### Immune hot vs cold classification — the primary finding

Using CD8+ T-cell fraction (median split) to define "immune-hot" vs "immune-cold" tumours,
and a feature matrix of the **remaining 8 major cell-type fractions** (T-cell categories
excluded to avoid circularity), an XGBoost classifier achieved **80.8% LOOCV accuracy**
on a balanced 13/13 split.

### SHAP — biological interpretation of the immune phenotype model

- **B-cell fraction** is the dominant predictor. Tumours with B-cell fraction above ~4%
  consistently push toward "hot"; near-zero B-cell fraction pushes toward "cold". This is
  consistent with **tertiary lymphoid structures (TLS)**, where B-cell and T-cell
  co-localisation reflects an organised anti-tumour immune response.
- **Cancer Epithelial fraction** is the second strongest predictor, with an inverse
  relationship — high malignant cell fraction (>40%) consistently predicts "cold". This
  matches the established **immune exclusion by tumour mass** mechanism in breast cancer.

### PCA vs scVI benchmark

Silhouette scores alone showed PCA with a *higher* raw cell-type silhouette (0.25 vs
0.14) and similar batch silhouettes (-0.04 vs -0.06) — at face value, inconclusive.
However, this is interpreted as a limitation of silhouette score for batch-corrected
embeddings: PCA's higher silhouette likely reflects residual patient-specific structure
(visible as per-patient islands in the UMAP) that inflates apparent separation without
reflecting generalisable biology. The >93% ground-truth cluster concordance achieved on
the scVI embedding — across all 26 patients — is the more meaningful validation.

## Tech stack

- **scanpy / anndata** — single-cell data structures and standard preprocessing
- **scvi-tools** — variational autoencoder for batch-corrected representation learning
- **leidenalg / python-igraph** — graph-based clustering
- **scikit-learn** — Leave-One-Out cross-validation, RandomForest, silhouette scoring
- **xgboost** — gradient-boosted classifiers
- **shap** — model explainability
- **matplotlib / seaborn** — visualisation

## Repository structure

```
scRNA_TME_pipeline_clean.ipynb   # full annotated pipeline notebook
results/
  adata_raw.h5ad                 # checkpoint: raw concatenated counts
  adata_pca_baseline.h5ad        # checkpoint: post-QC, normalised, HVG, PCA
  adata_scvi.h5ad                # checkpoint: scVI latent space + UMAP
  adata_scvi_clustered.h5ad      # checkpoint: + Leiden clusters + annotations
  scvi_model/                    # trained scVI model
  *.png                          # all figures referenced above
```

## Notes on methodology

- A 30,000-cell subsample was used for scVI training to keep CPU training time
  tractable (~10 min vs >30 min on the full 98,593 cells). All 9 major cell populations
  remain well represented (smallest population retains >900 cells).
- All classification results use **Leave-One-Out Cross-Validation** due to the small
  number of patients (n=26), which is the most reliable evaluation strategy at this
  sample size.
- The immune hot/cold feature matrix explicitly excludes all T-cell populations to avoid
  circularity between the label definition and the predictive features.
