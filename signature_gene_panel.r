build_signature_panel <- function(de_result, enrich_result, leading_edge_genes,
                                   tf_names, output_dir = "outputs") {

  de_tab <- de_result$table_full
  de_gene_col <- if ("gene" %in% names(de_tab)) de_tab$gene else rownames(de_tab)
  rownames(de_tab) <- de_gene_col

  # (a) top 50 DE genes by |t|
  top50 <- head(de_gene_col[order(-abs(de_tab$t))], 50)

  sig <- enrich_result$sig
  pathway_membership <- list()
  if (!is.null(sig) && nrow(sig) > 0) {
    for (i in seq_len(nrow(sig))) {
      genes_i <- strsplit(sig$leadingEdge[i], ";")[[1]]
      genes_i <- trimws(genes_i)
      genes_i <- genes_i[nzchar(genes_i)]
      for (g in genes_i) pathway_membership[[g]] <- c(pathway_membership[[g]], sig$pathway[i])
    }
  }

  # (b) union of leading-edge genes from significant pathways
  le_all <- unique(c(leading_edge_genes, names(pathway_membership)))

  # (c) all TF names
  panel_genes <- unique(c(top50, le_all, tf_names))

  panel <- data.frame(gene_symbol = panel_genes, stringsAsFactors = FALSE)
  panel$logFC     <- de_tab[panel$gene_symbol, "logFC"]
  panel$t_stat    <- de_tab[panel$gene_symbol, "t"]
  panel$adj_p_val <- de_tab[panel$gene_symbol, "adj.P.Val"]

  panel$leading_edge_pathways <- vapply(panel$gene_symbol, function(g) {
    p <- pathway_membership[[g]]
    if (is.null(p)) "" else paste(unique(p), collapse = ";")
  }, character(1))

  panel$is_tf <- panel$gene_symbol %in% tf_names

  panel$n_pathway_lists <- vapply(panel$gene_symbol, function(g) {
    p <- pathway_membership[[g]]
    if (is.null(p)) 0L else length(unique(p))
  }, integer(1))

  panel <- panel[order(-abs(panel$t_stat)), ]
  rownames(panel) <- NULL

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(output_dir, "AROHI_signature_gene_panel.csv")
  write.csv(panel, out_file, row.names = FALSE)

  message(sprintf("[PANEL] AROHI Signature Gene Panel: %d genes -> %s", nrow(panel), out_file))

  cat("\n==================== Top 20 genes (by |t|) ====================\n")
  top20 <- head(panel, 20)
  for (i in seq_len(nrow(top20))) {
    r <- top20[i, ]
    cat(sprintf(
      "%2d. %-12s logFC=%6.2f  t=%7.2f  adj.P=%.2e  TF=%-5s  pathways=%d\n",
      i, r$gene_symbol, r$logFC, r$t_stat, r$adj_p_val, r$is_tf, r$n_pathway_lists
    ))
  }
  cat("\n")

  panel
}

plot_panel_heatmap <- function(expr_train, labels_train, panel_genes, top_n = 50, output_dir = "outputs") {
  .ensure_pkg("pheatmap", bioc = TRUE)
  suppressPackageStartupMessages(library(pheatmap))

  if (is.data.frame(panel_genes)) {
    genes_sorted <- panel_genes$gene_symbol
    logfc_lookup <- setNames(panel_genes$logFC, panel_genes$gene_symbol)
  } else {
    genes_sorted <- panel_genes
    logfc_lookup <- NULL
  }

  top_genes <- intersect(head(genes_sorted, top_n), rownames(expr_train))
  if (length(top_genes) < min(top_n, length(genes_sorted))) {
    warning(sprintf(
      "[HEATMAP] Only %d / %d requested panel genes found in expr_train.",
      length(top_genes), min(top_n, length(genes_sorted))
    ))
  }

  mat <- expr_train[top_genes, , drop = FALSE]

  samp_ids     <- colnames(mat)
  labels_train <- droplevels(labels_train[samp_ids])
  ord          <- order(labels_train)
  mat          <- mat[, ord, drop = FALSE]
  labels_ord   <- labels_train[ord]

  col_anno <- data.frame(Response = as.character(labels_ord), row.names = colnames(mat))

  if (!is.null(logfc_lookup)) {
    row_dir <- ifelse(logfc_lookup[top_genes] >= 0, "Up_in_pCR", "Up_in_RD")
  } else {
    grp_means <- t(apply(mat, 1, function(x) tapply(x, labels_ord, mean, na.rm = TRUE)))
    row_dir <- ifelse(grp_means[, "pCR"] - grp_means[, "RD"] >= 0, "Up_in_pCR", "Up_in_RD")
  }
  row_anno <- data.frame(logFC_direction = row_dir, row.names = top_genes)

  anno_colors <- list(
    Response         = c(pCR = "#2b8cbe", RD = "#de2d26"),
    logFC_direction  = c(Up_in_pCR = "#2b8cbe", Up_in_RD = "#de2d26")
  )

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(output_dir, "AROHI_panel_heatmap.png")
  bwr <- colorRampPalette(c("blue", "white", "red"))(100)

  pheatmap::pheatmap(
    mat,
    color             = bwr,
    scale             = "row",
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    annotation_col    = col_anno,
    annotation_row    = row_anno,
    annotation_colors = anno_colors,
    show_colnames     = FALSE,
    fontsize_row      = 6,
    main              = sprintf("AROHI Signature Panel (top %d genes)", length(top_genes)),
    filename          = out_file,
    width             = 10,
    height            = max(6, 0.15 * length(top_genes))
  )

  message(sprintf("[HEATMAP] Saved -> %s", out_file))
  invisible(out_file)
}