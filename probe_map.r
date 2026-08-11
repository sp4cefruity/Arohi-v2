
.ensure_pkg <- function(pkgs, bioc = FALSE) {
  need <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (!length(need)) {
    return(invisible())
  }
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(need, ask = FALSE, update = FALSE)
  } else {
    install.packages(need, repos = "https://cloud.r-project.org")
  }
}
 
collapse_probes_to_symbols <- function(expr_mat, chip = "hgu133plus2.db"
                                       , method = "maxIQR") {
  .ensure_pkg(c("AnnotationDbi", chip), bioc = TRUE)
  suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(package = chip, character.only = TRUE)
  })
 
  db <- get(chip)
 
  message(sprintf("[MAP] Looking up %d probe IDs in %s...",
                  nrow(expr_mat), chip))
 
  symbol_map <- tryCatch(
    AnnotationDbi::mapIds(db, keys = rownames(expr_mat),
                          column = "SYMBOL",
                          keytype = "PROBEID",
                          multiVals = "first"),
    error = function(e) {
      stop(sprintf(
        "[ERROR] Annotation lookup failed : %s",
        conditionMessage(e)
      ))
    }
  )
 
  # filter
  keep <- !is.na(symbol_map) & nchar(trimws(symbol_map)) > 0
  n_total <- length(symbol_map)
  n_annotated <- sum(keep)
  message(sprintf("[MAP] Annotated: %d / %d probes (%.1f%%).",
                  n_annotated, n_total, 100 * n_annotated / n_total))
 
  expr_sub <- expr_mat[keep, , drop = FALSE]
  symbols <- symbol_map[keep]
 
  # collapse
  if (method == "maxIQR") {
    probe_iqr <- apply(expr_sub, 1L, IQR)
    gene_groups <- split(seq_len(nrow(expr_sub)), symbols)
    gene_mat <- do.call(rbind, lapply(gene_groups, function(idx) {
      best <- idx[which.max(probe_iqr[idx])]
      expr_sub[best, , drop = FALSE]
    }))
    rownames(gene_mat) <- names(gene_groups)
  } else if (method == "mean") {
    gene_groups <- split(seq_len(nrow(expr_sub)), symbols)
    gene_mat <- do.call(rbind, lapply(gene_groups, function(idx) {
      colMeans(expr_sub[idx, , drop = FALSE])
    }))
  } else {
    stop(sprintf("[ERROR] Unknown collapse method : '%s'. 
                  Use 'maxIQR' or 'mean'.", method))
  }
  gene_mat <- gene_mat[order(rownames(gene_mat)), ]
  message(sprintf("[MAP] Collapsed to %d unique Gene Symbols.", nrow(gene_mat)))
  gene_mat
}