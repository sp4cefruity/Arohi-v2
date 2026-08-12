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

extract_leading_edge_matrix <- function(enrich_result, expr_train, expr_test,
                                         top_n_pathways = 5, output_dir = "outputs") {

  sig <- enrich_result$sig
  if (is.null(sig) || nrow(sig) == 0)
    stop("[ERROR] enrich_result$sig has no significant pathways to draw leading-edge genes from.")

  top_idx      <- order(-abs(sig$NES))[seq_len(min(top_n_pathways, nrow(sig)))]
  top_pathways <- sig[top_idx, ]

  message(sprintf(
    "[LEADING-EDGE] Top %d pathway(s) by |NES|: %s",
    nrow(top_pathways), paste(top_pathways$pathway, collapse = ", ")
  ))

  gene_lists <- lapply(top_pathways$leadingEdge, function(s) {
    g <- strsplit(s, ";")[[1]]
    trimws(g[nzchar(trimws(g))])
  })

  union_genes <- unique(unlist(gene_lists))
  message(sprintf(
    "[LEADING-EDGE] Union of leading-edge genes across selected pathways: %d genes.",
    length(union_genes)
  ))

  present_train <- intersect(union_genes, rownames(expr_train))
  present_test  <- intersect(union_genes, rownames(expr_test))
  common_genes  <- intersect(present_train, present_test)

  missing <- setdiff(union_genes, common_genes)
  if (length(missing)) {
    warning(sprintf(
      "[LEADING-EDGE] %d leading-edge gene(s) missing from train and/or test matrix and were dropped: %s",
      length(missing), paste(missing, collapse = ", ")
    ))
  }
  if (!length(common_genes))
    stop("[ERROR] No leading-edge genes present in both train and test matrices.")

  train_mat <- expr_train[common_genes, , drop = FALSE]
  test_mat  <- expr_test[common_genes, , drop = FALSE]

  enrich_dir <- file.path(output_dir, "enrichment")
  scores_dir <- file.path(output_dir, "scores")
  dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(scores_dir, recursive = TRUE, showWarnings = FALSE)

  .ensure_pkg("jsonlite", bioc = FALSE)
  jsonlite::write_json(
    common_genes, file.path(enrich_dir, "leading_edge_genes.json"),
    pretty = TRUE, auto_unbox = FALSE
  )

  write.csv(train_mat, file.path(scores_dir, "gse_train_leading_edge.csv"))
  write.csv(test_mat,  file.path(scores_dir, "gse_test_leading_edge.csv"))

  message(sprintf(
    "[LEADING-EDGE] Saved %d-gene matrix -> %s, %s",
    length(common_genes),
    file.path(scores_dir, "gse_train_leading_edge.csv"),
    file.path(scores_dir, "gse_test_leading_edge.csv")
  ))

  list(train = train_mat, test = test_mat, genes = common_genes)
}
