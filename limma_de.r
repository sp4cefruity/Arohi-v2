
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

.ensure_pkg("limma", bioc = TRUE)
suppressPackageStartupMessages(library(limma))

check_matrix_for_limma <- function(expr_mat) {
  vals <- as.numeric(expr_mat)
  vals <- vals[is.finite(vals)]
  rng  <- range(vals)
  prop_over_100 <- mean(vals > 100)

  message(sprintf("[LIMMA-QC] Expression range: [%.3f, %.3f]", rng[1], rng[2]))

  if (prop_over_100 > 0.10) {
    warning(sprintf(
      "[WARN] %.1f%% of values exceed 100 -- this looks like raw counts, not log2 data. Check upstream normalisation before running limma.",
      100 * prop_over_100
    ))
  }

  invisible(list(range = rng, prop_over_100 = prop_over_100))
}

run_limma_de <- function(expr_mat, labels, meta = NULL, adj_method = "BH",
                          fdr_cutoff = 0.05, lfc_cutoff = 0.3,
                          output_dir = "outputs", voom_weights = NULL) {

  check_matrix_for_limma(expr_mat)

  samp_ids <- colnames(expr_mat)
  labels   <- droplevels(labels[samp_ids])
  labels   <- relevel(labels, ref = "RD")   # positive logFC == higher in pCR

  if (!is.null(meta) && "gse_id" %in% names(meta)) {
    meta_sub <- meta[match(samp_ids, meta$sample_id), ]
    gse_id   <- factor(meta_sub$gse_id)
    message(sprintf(
      "[LIMMA] Batch-correcting for dataset-of-origin (%d datasets) via model.matrix(~ gse_id + labels).",
      nlevels(gse_id)
    ))
    design <- model.matrix(~ gse_id + labels)
  } else {
    message("[LIMMA] No gse_id column in meta; fitting design without a batch covariate.")
    design <- model.matrix(~ labels)
  }

  if (!is.null(voom_weights)) {
    message("[LIMMA] Using voom precision weights in lmFit().")
    fit <- lmFit(expr_mat[, samp_ids], design, weights = voom_weights[, samp_ids])
  } else {
    fit <- lmFit(expr_mat[, samp_ids], design)
  }

  fit <- eBayes(fit, trend = TRUE, robust = TRUE)

  if (!"labelspCR" %in% colnames(design))
    stop("[ERROR] Expected coefficient 'labelspCR' not found in design matrix columns: ",
         paste(colnames(design), collapse = ", "))

  tt_full <- topTable(fit, coef = "labelspCR", number = Inf,
                       adjust.method = adj_method, sort.by = "none")
  tt_full$gene <- rownames(tt_full)

  ranked_by_t <- tt_full[order(-tt_full$t), ]

  tt_sig <- tt_full[!is.na(tt_full$adj.P.Val) &
                       tt_full$adj.P.Val <= fdr_cutoff &
                       abs(tt_full$logFC) >= lfc_cutoff, ]

  message(sprintf(
    "[LIMMA] %d / %d genes significant at FDR <= %.2f & |logFC| >= %.2f.",
    nrow(tt_sig), nrow(tt_full), fdr_cutoff, lfc_cutoff
  ))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(tt_full, file.path(output_dir, "de_results_full.csv"), row.names = FALSE)
  write.csv(tt_sig,  file.path(output_dir, "de_results_sig.csv"),  row.names = FALSE)

  
  volcano_file <- file.path(output_dir, "de_volcano.png")
  tt_full$sig <- with(tt_full, !is.na(adj.P.Val) & adj.P.Val <= fdr_cutoff & abs(logFC) >= lfc_cutoff)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ggplot2))
    p <- ggplot(tt_full, aes(x = logFC, y = -log10(P.Value), color = sig)) +
      geom_point(alpha = 0.6, size = 1) +
      scale_color_manual(values = c(`TRUE` = "#de2d26", `FALSE` = "grey70")) +
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
      labs(
        title = "Differential expression: pCR vs RD", x = "logFC (pCR vs RD)",
        y = "-log10(P value)",
        color = sprintf("FDR<=%.2f & |logFC|>=%.2f", fdr_cutoff, lfc_cutoff)
      ) +
      theme_minimal()
    ggsave(volcano_file, plot = p, width = 7, height = 6, dpi = 150)
  } else {
    png(volcano_file, width = 700, height = 600)
    with(tt_full, plot(logFC, -log10(P.Value), col = ifelse(sig, "red", "grey70"),
                        pch = 16, main = "Differential expression: pCR vs RD"))
    abline(v = c(-lfc_cutoff, lfc_cutoff), lty = 2)
    dev.off()
  }
  message(sprintf("[LIMMA] Saved volcano plot -> %s", volcano_file))

  list(
    table_full  = tt_full,
    table_sig   = tt_sig,
    gene_list   = tt_sig$gene,
    ranked_by_t = ranked_by_t
  )
}