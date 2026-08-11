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

map_ensembl_to_symbol <- function(count_mat, max_retries = 3L, retry_wait_sec = 5) {
  .ensure_pkg("biomaRt", bioc = TRUE)
  suppressPackageStartupMessages(library(biomaRt))

  rownames(count_mat) <- sub("\\..*$", "", rownames(count_mat))  # strip version suffixes

  connect_mart <- function() {
    biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  }

  mart <- NULL
  attempt <- 0L
  while (attempt < max_retries && is.null(mart)) {
    attempt <- attempt + 1L
    mart <- tryCatch(connect_mart(), error = function(e) {
      warning(sprintf("[WARN] biomaRt connection attempt %d/%d failed: %s",
                       attempt, max_retries, conditionMessage(e)))
      NULL
    })
    if (is.null(mart) && attempt < max_retries) Sys.sleep(retry_wait_sec)
  }
  if (is.null(mart)) stop("[ERROR] Could not connect to Ensembl BioMart after retries.")

  message(sprintf("[ENSEMBL] Querying hgnc_symbol for %d Ensembl gene IDs...", nrow(count_mat)))

  query_biomart <- function() {
    biomaRt::getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters    = "ensembl_gene_id",
      values     = rownames(count_mat),
      mart       = mart
    )
  }

  map_df <- tryCatch(
    query_biomart(),
    error = function(e) {
      warning(sprintf("[WARN] getBM() failed: %s -- retrying once.", conditionMessage(e)))
      Sys.sleep(retry_wait_sec)
      query_biomart()
    }
  )

  map_df <- map_df[!is.na(map_df$hgnc_symbol) & nchar(trimws(map_df$hgnc_symbol)) > 0, ]
  map_df <- map_df[!duplicated(map_df$ensembl_gene_id), ]

  keep <- rownames(count_mat) %in% map_df$ensembl_gene_id
  message(sprintf("[ENSEMBL] Mapped %d / %d Ensembl IDs to HGNC symbols.", sum(keep), nrow(count_mat)))
  if (!any(keep)) stop("[ERROR] No Ensembl IDs mapped to an HGNC symbol.")

  mat_sub <- count_mat[keep, , drop = FALSE]
  symbols <- map_df$hgnc_symbol[match(rownames(mat_sub), map_df$ensembl_gene_id)]

  probe_iqr   <- apply(mat_sub, 1L, IQR)
  gene_groups <- split(seq_len(nrow(mat_sub)), symbols)
  gene_mat <- do.call(rbind, lapply(gene_groups, function(idx) {
    best <- idx[which.max(probe_iqr[idx])]
    mat_sub[best, , drop = FALSE]
  }))
  rownames(gene_mat) <- names(gene_groups)
  gene_mat <- gene_mat[order(rownames(gene_mat)), ]

  message(sprintf("[ENSEMBL] Collapsed to %d unique HGNC gene symbols.", nrow(gene_mat)))
  gene_mat
}