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

.build_ranked_list <- function(de_table_full) {
  if (!"t" %in% names(de_table_full))
    stop("[ERROR] de_table_full must contain a 't' (t-statistic) column.")

  genes <- if ("gene" %in% names(de_table_full)) de_table_full$gene else rownames(de_table_full)

  ranked <- setNames(de_table_full$t, genes)
  ranked <- ranked[!is.na(ranked)]
  sort(ranked, decreasing = TRUE)
}

.fetch_gene_sets <- function(collection = "H", subcollection = NULL) {
  .ensure_pkg("msigdbr", bioc = FALSE)
  suppressPackageStartupMessages(library(msigdbr))

  args <- list(species = "Homo sapiens", collection = collection)
  if (!is.null(subcollection)) args$subcollection <- subcollection

  sets_df    <- do.call(msigdbr::msigdbr, args)
  gene_sets  <- split(sets_df$gene_symbol, sets_df$gs_name)

  message(sprintf(
    "[GSEA] Fetched %d gene sets from MSigDB collection '%s'%s.",
    length(gene_sets), collection,
    if (!is.null(subcollection)) sprintf(" (%s)", subcollection) else ""
  ))
  gene_sets
}

.run_fgsea <- function(ranked_genes, gene_sets, min_size = 15L, max_size = 500L,
                        n_perm = 1000L, seed = 42L) {
  .ensure_pkg("fgsea", bioc = TRUE)
  suppressPackageStartupMessages(library(fgsea))
  set.seed(seed)

  res <- fgsea::fgsea(
    pathways    = gene_sets,
    stats       = ranked_genes,
    minSize     = min_size,
    maxSize     = max_size,
    nPermSimple = n_perm
  )

  res$padj        <- p.adjust(res$pval, method = "BH")
  res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse = ";"), character(1))

  res <- as.data.frame(res)
  res[order(-res$NES), ]
}


run_pathway_enrichment <- function(de_result, fdr_cutoff = 0.05, output_dir = "outputs") {
  ranked        <- .build_ranked_list(de_result$table_full)
  hallmark_sets <- .fetch_gene_sets(collection = "H")

  hallmark_res <- .run_fgsea(ranked, hallmark_sets)
  hallmark_res$direction <- ifelse(hallmark_res$NES > 0, "RESPONSE_pCR", "RESISTANCE_RD")

  sig_res <- hallmark_res[!is.na(hallmark_res$padj) & hallmark_res$padj <= fdr_cutoff, ]

  message(sprintf(
    "[GSEA] %d / %d Hallmark pathways significant at FDR <= %.2f.",
    nrow(sig_res), nrow(hallmark_res), fdr_cutoff
  ))

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(hallmark_res, file.path(output_dir, "fgsea_hallmark.csv"), row.names = FALSE)
  write.csv(sig_res,      file.path(output_dir, "fgsea_sig_pathways.csv"), row.names = FALSE)

  barplot_file <- file.path(output_dir, "enrichment_barplot.png")
  if (nrow(sig_res) > 0) {
    plot_df <- sig_res[order(sig_res$NES), ]
    plot_df$pathway <- factor(plot_df$pathway, levels = plot_df$pathway)

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      suppressPackageStartupMessages(library(ggplot2))
      p <- ggplot(plot_df, aes(x = pathway, y = NES, fill = direction)) +
        geom_col() + coord_flip() +
        labs(title = "Significant Hallmark pathways", x = NULL, y = "Normalized Enrichment Score") +
        theme_minimal()
      ggsave(barplot_file, plot = p, width = 9, height = max(4, 0.3 * nrow(plot_df)), dpi = 150)
    } else {
      png(barplot_file, width = 900, height = max(400, 30 * nrow(plot_df)))
      barplot(plot_df$NES, names.arg = plot_df$pathway, horiz = TRUE, las = 1,
              col = ifelse(plot_df$NES > 0, "#de2d26", "#2b8cbe"),
              main = "Significant Hallmark pathways")
      dev.off()
    }
    message(sprintf("[GSEA] Saved barplot -> %s", barplot_file))
  } else {
    message("[GSEA] No significant pathways at chosen FDR cutoff; skipping barplot.")
  }

  list(hallmark = hallmark_res, sig = sig_res)
}

extract_signatures <- function(enrich_result, de_result, use_leading_edge = TRUE, output_dir = "outputs") {
  sig_pathways <- enrich_result$sig
  if (is.null(sig_pathways) || nrow(sig_pathways) == 0) {
    warning("[SIGNATURES] No significant pathways in enrich_result$sig; returning empty list.")
    return(list())
  }

  sig_dir <- file.path(output_dir, "enrichment", "signatures")
  dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)

  de_tab   <- de_result$table_full
  de_genes <- if ("gene" %in% names(de_tab)) de_tab$gene else rownames(de_tab)
  rownames(de_tab) <- de_genes

  full_sets <- NULL
  if (!use_leading_edge) {
    full_sets <- tryCatch(.fetch_gene_sets(collection = "H"), error = function(e) {
      warning(sprintf(
        "[SIGNATURES] Could not fetch full gene sets (%s); falling back to leading-edge only.",
        conditionMessage(e)
      ))
      NULL
    })
  }

  signatures <- list()

  for (i in seq_len(nrow(sig_pathways))) {
    row       <- sig_pathways[i, ]
    pathway   <- row$pathway
    nes       <- row$NES
    direction <- row$direction

    le_genes <- strsplit(row$leadingEdge, ";")[[1]]
    le_genes <- trimws(le_genes)
    le_genes <- le_genes[nzchar(le_genes)]

    if (use_leading_edge || is.null(full_sets)) {
      candidate_genes <- le_genes
    } else {
      pathway_all <- full_sets[[pathway]]
      candidate_genes <- if (!is.null(pathway_all)) union(le_genes, pathway_all) else le_genes
    }

    present <- intersect(candidate_genes, de_genes)
    if (!length(present)) {
      warning(sprintf(
        "[SIGNATURES] Pathway '%s': no leading-edge genes found in de_result$table_full; skipping.", pathway
      ))
      next
    }

    df <- de_tab[present, c("logFC", "t", "adj.P.Val")]
    df$gene            <- present
    df$in_leading_edge <- present %in% le_genes
    df$pathway         <- pathway
    df$NES             <- nes
    df$direction       <- direction

    df <- df[order(-abs(df$t)), ]
    df <- df[, c("gene", "logFC", "t", "adj.P.Val", "in_leading_edge", "pathway", "NES", "direction")]
    rownames(df) <- NULL

    safe_name <- gsub("[^A-Za-z0-9_]+", "_", pathway)
    write.csv(df, file.path(sig_dir, sprintf("sig_%s.csv", safe_name)), row.names = FALSE)

    signatures[[pathway]] <- df
  }

  if (!length(signatures)) {
    warning("[SIGNATURES] No signatures extracted.")
    return(list())
  }

  master <- do.call(rbind, signatures)
  write.csv(master, file.path(output_dir, "enrichment", "signatures_summary.csv"), row.names = FALSE)

  response_pw   <- names(signatures)[sapply(signatures, function(d) d$direction[1] == "RESPONSE_pCR")]
  resistance_pw <- names(signatures)[sapply(signatures, function(d) d$direction[1] == "RESISTANCE_RD")]

  cat("\n==================== RESPONSE (pCR) pathways ====================\n")
  for (p in response_pw) {
    d <- signatures[[p]]
    cat(sprintf("- %s (NES=%.2f): %s\n", p, d$NES[1], paste(head(d$gene, 8), collapse = ", ")))
  }
  cat("\n==================== RESISTANCE (RD) pathways ====================\n")
  for (p in resistance_pw) {
    d <- signatures[[p]]
    cat(sprintf("- %s (NES=%.2f): %s\n", p, d$NES[1], paste(head(d$gene, 8), collapse = ", ")))
  }
  cat("\n")

  message(sprintf(
    "[SIGNATURES] Extracted %d pathway signature(s) -> %s (+ signatures_summary.csv)",
    length(signatures), sig_dir
  ))

  signatures
}

signatures_to_genesets <- function(signatures) {
  lapply(signatures, function(df) unique(df$gene))
}
