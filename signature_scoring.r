if (!exists(".ensure_pkg", mode = "function")) {
  .ensure_pkg <- function(pkgs, bioc = FALSE) {
    need <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
    if (!length(need)) return(invisible())
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(need, ask = FALSE, update = FALSE)
    } else {
      install.packages(need, repos = "https://cloud.r-project.org")
    }
  }
}

.ensure_pkg("GSVA", bioc = TRUE)
suppressPackageStartupMessages(library(GSVA))

run_ssgsea <- function(expr_mat, gene_sets, min_size = 5L, max_size = 500L) {
  expr_mat <- as.matrix(expr_mat)

  if (exists("ssgseaParam", where = asNamespace("GSVA"))) {
    param  <- GSVA::ssgseaParam(expr_mat, gene_sets, minSize = min_size, maxSize = max_size)
    scores <- GSVA::gsva(param, verbose = FALSE)
  } else {
    scores <- GSVA::gsva(
      expr = expr_mat, gset.idx.list = gene_sets, method = "ssgsea",
      min.sz = min_size, max.sz = max_size, verbose = FALSE
    )
  }

  as.matrix(scores)   # signatures x samples
}

run_universal_scoring <- function(gse_train, gse_test, gene_sets, output_dir = "outputs") {

  message(sprintf("[SSGSEA] Scoring train set (%d samples) against %d signature(s)...",
                   ncol(gse_train), length(gene_sets)))
  train_scores <- run_ssgsea(gse_train, gene_sets)

  message(sprintf("[SSGSEA] Scoring test set (%d samples) against %d signature(s)...",
                   ncol(gse_test), length(gene_sets)))
  test_scores <- run_ssgsea(gse_test, gene_sets)

  scores_dir <- file.path(output_dir, "scores")
  dir.create(scores_dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(train_scores, file.path(scores_dir, "gse_train_ssgsea.csv"))
  write.csv(test_scores,  file.path(scores_dir, "gse_test_ssgsea.csv"))

  message(sprintf(
    "[SSGSEA] Saved %s, %s",
    file.path(scores_dir, "gse_train_ssgsea.csv"),
    file.path(scores_dir, "gse_test_ssgsea.csv")
  ))

  list(train = train_scores, test = test_scores)
}
