if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here", repos = "https://cloud.r-project.org")
}
setwd(here::here())

source("R/normalize.r")
source("R/ingest_multi_dataset.r")
source("R/probe_map.r")
source("R/ensembl_map.r")
source("R/split_harmonize.r")
source("R/limma_de.r")
source("R/pathway_enrichment.r")
source("R/signature_scoring.r")
source("R/leading_edge_extract.r")
source("R/tf_activity.r")
source("R/signature_gene_panel.r")

dir.create("outputs/cache", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/qc",    recursive = TRUE, showWarnings = FALSE)

.save <- function(obj, filename, dir = "outputs") {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, filename)
  write.csv(obj, path)
  message(sprintf("[SAVE] %s", path))
}

message("[PHASE 0] Ingesting all GEO datasets...")
merged <- merge_datasets(DATASET_LIST)

MERGED_EXPR   <- merged$expr_mat
MERGED_LABELS <- merged$labels
MERGED_META   <- merged$metadata

plot_dataset_summary(MERGED_META)

message("[PHASE 1] Stratified train/test split across datasets and labels...")
splits <- split_train_test(MERGED_EXPR, MERGED_LABELS, MERGED_META, train_prop = 0.70, seed = 42L)

message("[PHASE 1] Test set locked (outputs/gse_test_locked.csv, written by split_train_test). Will not be touched until Phase 5.")

gse_train    <- splits$train
gse_test_raw <- splits$test
labels_train <- splits$labels_train

harmonized <- choose_harmonisation(gse_train, gse_test_raw, MERGED_META)

.save(harmonized$train, "gse_train.csv")
.save(harmonized$test,  "gse_test_harmonized_locked.csv")

message("[PHASE 2] Differential expression (limma, batch-corrected by gse_id)...")
de_result <- run_limma_de(
  expr_mat   = harmonized$train,
  labels     = labels_train,
  meta       = MERGED_META,
  output_dir = "outputs"
)

message("[PHASE 2] Pathway enrichment (fgsea, Hallmark, ranked by t-statistic)...")
pathway_result <- run_pathway_enrichment(de_result, output_dir = "outputs")

if (is.null(pathway_result$sig) || nrow(pathway_result$sig) == 0)
  stop("[ERROR] No Hallmark pathways significant at FDR <= 0.05; cannot build signatures. Try a more lenient fdr_cutoff.")

message("[PHASE 2] Extracting pathway signatures...")
signatures <- extract_signatures(pathway_result, de_result, output_dir = "outputs")
gene_sets  <- signatures_to_genesets(signatures)

message("[PHASE 3] ssGSEA scoring (train/test) on pathway signatures...")
gse_test_h <- harmonized$test
ssgsea_scores <- run_universal_scoring(harmonized$train, gse_test_h, gene_sets, output_dir = "outputs")

message("[PHASE 3] Extracting gene-level leading-edge feature matrix...")
leading_edge <- extract_leading_edge_matrix(
  pathway_result, harmonized$train, gse_test_h,
  top_n_pathways = 5, output_dir = "outputs"
)

message("[PHASE 3] Scoring transcription factor activity (DoRothEA + VIPER)...")
tf_scores <- score_tf_activity_all_cohorts(harmonized$train, gse_test_h, output_dir = "outputs")

message("[PHASE 3] Building the AROHI Signature Gene Panel...")
gene_panel <- build_signature_panel(
  de_result           = de_result,
  enrich_result       = pathway_result,
  leading_edge_genes  = unique(unlist(gene_sets)),
  tf_names            = colnames(tf_scores$train),
  output_dir          = "outputs"
)

plot_panel_heatmap(harmonized$train, labels_train, gene_panel, top_n = 50, output_dir = "outputs")

message("[PIPELINE] Done through Phase 3.")