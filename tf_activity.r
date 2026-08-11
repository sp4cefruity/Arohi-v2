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

.ensure_pkg(c("dorothea", "viper"), bioc = TRUE)

TARGET_TFS <- c("TFAP2C", "TFAP2A", "SP1", "MYC", "E2F1", "ESR1")

compute_tf_activity <- function(expr_mat, cohort_label, output_dir = "outputs") {
  suppressPackageStartupMessages({
    library(dorothea)
    library(viper)
  })

  regulon_df <- dorothea::dorothea_hs
  regulon_df <- regulon_df[regulon_df$confidence %in% c("A", "B", "C"), ]
  regulons   <- dorothea::df2regulon(regulon_df)

  message(sprintf(
    "[TF] [%s] Running VIPER (method='scale') on %d samples, %d regulons (confidence A/B/C)...",
    cohort_label, ncol(expr_mat), length(regulons)
  ))

  viper_res <- viper::viper(as.matrix(expr_mat), regulons, method = "scale", verbose = FALSE)

  missing_tfs <- setdiff(TARGET_TFS, rownames(viper_res))
  if (length(missing_tfs)) {
    warning(sprintf(
      "[TF] [%s] TF(s) missing from regulon output, filled with 0: %s",
      cohort_label, paste(missing_tfs, collapse = ", ")
    ))
  }

  tf_mat <- matrix(
    0, nrow = length(TARGET_TFS), ncol = ncol(viper_res),
    dimnames = list(TARGET_TFS, colnames(viper_res))
  )
  present_tfs <- intersect(TARGET_TFS, rownames(viper_res))
  tf_mat[present_tfs, ] <- viper_res[present_tfs, , drop = FALSE]

  tf_df <- as.data.frame(t(tf_mat))[, TARGET_TFS, drop = FALSE]

  scores_dir <- file.path(output_dir, "scores")
  dir.create(scores_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(scores_dir, sprintf("%s_tf_scores.csv", cohort_label))
  write.csv(tf_df, out_file)

  message(sprintf("[TF] [%s] Saved TF activity scores -> %s", cohort_label, out_file))

  tf_df
}

score_tf_activity_all_cohorts <- function(expr_train, expr_test, output_dir = "outputs") {
  train_scores <- compute_tf_activity(expr_train, cohort_label = "gse_train", output_dir = output_dir)
  test_scores  <- compute_tf_activity(expr_test,  cohort_label = "gse_test",  output_dir = output_dir)
  list(train = train_scores, test = test_scores)
}
