# AROHI v2 — Neoadjuvant Breast Cancer pCR Prediction Pipeline
 
AROHI v2 is a two-stage (R + Python) pipeline that predicts pathologic complete
response (**pCR** vs. residual disease, **RD**) to neoadjuvant chemotherapy in
breast cancer from public gene-expression cohorts. It ingests and harmonizes
31 GEO datasets, derives pathway- and transcription-factor–level features via
differential expression and enrichment analysis, and trains a stacked
GLM + XGBoost ensemble that is evaluated with nested, leakage-safe
cross-validation.
 
The pipeline is split into two phases that hand off through CSV files on disk:
 
- **R (Phases 0–3):** data ingestion, harmonization, differential expression,
  pathway/TF signature discovery, and feature scoring.
- **Python (Phases 4–5):** feature fusion, feature selection, model training,
  and evaluation.
## Pipeline Overview
 
```
Phase 0  Ingest & merge 31 GEO datasets ─────────────── R/ingest_multi_dataset.r
Phase 1  Stratified train/test split (70/30, locked)    R/split_harmonize.r
Phase 2  limma differential expression                  R/limma_de.r
         fgsea pathway enrichment (Hallmark)             R/pathway_enrichment.r
         Signature/leading-edge gene extraction          R/leading_edge_extract.r
Phase 3  ssGSEA pathway scoring                          R/signature_scoring.r
         TF activity scoring (DoRothEA + VIPER)          R/tf_activity.r
         AROHI Signature Gene Panel construction          R/signature_gene_panel.r
──────────────────────────────────────────────────────────────────────────────
Phase 4  Feature fusion (clinical + pathway + TF scores)  feature_fusion.py
         ElasticNetCV feature selection                   feature_selection.py
Phase 5  GLM + XGBoost training, Optuna tuning,            model_training.py
         stacked ensemble, SHAP/coefficient explanation    evaluation.py
         End-to-end orchestration                          run_pipeline.py
```
 
## Repository Structure
 
```
sp4cefruity-arohi-v2/
├── run_pipeline.py          # Orchestrates Phases 4-5 end to end
├── feature_fusion.py        # Merges clinical, ssGSEA, leading-edge & TF features
├── feature_selection.py     # Train-only preprocessing + ElasticNetCV selection
├── model_training.py        # GLM + XGBoost training, Optuna tuning, stacking
├── evaluation.py            # Dual-model prediction, metrics, SHAP/GLM explanation
└── R/
    ├── dependencies.R           # Checks/installs required CRAN + Bioconductor packages
    ├── ingest_multi_dataset.r   # Downloads & merges the 31 GSE cohorts
    ├── normalize.r              # Expression normalization utilities
    ├── probe_map.r              # Microarray probe → gene mapping
    ├── ensembl_map.r            # Ensembl ID mapping
    ├── split_harmonize.r        # Stratified train/test split + cross-cohort harmonization
    ├── limma_de.r                # Batch-corrected differential expression (limma)
    ├── pathway_enrichment.r     # fgsea enrichment against Hallmark gene sets
    ├── leading_edge_extract.r   # Leading-edge gene matrix extraction
    ├── signature_scoring.r      # ssGSEA scoring of pathway signatures
    ├── tf_activity.r            # Transcription-factor activity (DoRothEA + VIPER)
    ├── signature_gene_panel.r   # Builds the final AROHI signature gene panel
    └── run_pipeline.r           # Orchestrates Phases 0-3 end to end
```
 
## Requirements
 
**R** (see `R/dependencies.R` for the authoritative list):
`here`, `BiocManager`, `Biobase`, `limma`, `edgeR`, `DESeq2`, `GEOquery`,
`affy`, `affyio`, `fgsea`, `msigdbr`, `GSVA`, `pheatmap`, `preprocessCore`,
`AnnotationDbi`, `biomaRt`, `dorothea`, `viper`, `ggplot2`, `R.utils`
 
**Python**: `pandas`, `numpy`, `scikit-learn`, `xgboost`, `optuna`, `shap`,
`matplotlib`, `joblib`
 
Install everything with:
 
```r
source("R/dependencies.R")
install_project_dependencies()
```
 
```bash
pip install pandas numpy scikit-learn xgboost optuna shap matplotlib joblib
```
 
## Usage
 
### 1. Run the R pipeline (Phases 0–3)
 
From the project root:
 
```r
source("R/run_pipeline.r")
```
 
This downloads and merges the 31 GEO cohorts, performs the locked 70/30
train/test split, runs differential expression and pathway enrichment, and
writes the ssGSEA, leading-edge, and TF activity score matrices to
`outputs/scores/` (consumed by the Python phase).
 
### 2. Run the Python pipeline (Phases 4–5)
 
```bash
python run_pipeline.py
```
 
Edit the `CONFIG` block at the top of `run_pipeline.py` to point at your
input paths (clinical metadata, score matrices, labels) before running.
 
Optional flags for cohort-agnostic modeling ablations:
 
```bash
# Compare all-features vs. gse_id-dropped vs. gse_id-dropped + top-150 MI features
python run_pipeline.py --ablate
 
# Run a single variant
python run_pipeline.py --drop-gse-id [--k 150] [--n-trials 60]
```
 
Outputs (fused feature matrices, selected features, trained model bundle,
evaluation plots and metrics) are written to `outputs/python/`.
 
## Modeling Notes
 
- **Leakage control:** imputation and scaling (`feature_selection.py`) are
  fit on the training fold only and applied transform-only elsewhere.
- **Honest evaluation:** hyperparameter tuning (Optuna, seed `7`) and the
  reported cross-validated metrics (seed `123`) use non-overlapping CV folds
  to avoid optimistic bias.
- **Ensemble:** predictions are produced from a GLM (ElasticNet-selected
  features), an XGBoost model (full feature set), and a stacked ensemble of
  the two.
- **Interpretability:** `evaluation.py` provides SHAP-based explanation for
  XGBoost and coefficient-based explanation for the GLM.
