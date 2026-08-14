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

split_train_test <- function(gene_mat, labels, meta, train_prop = 0.70, seed = 42L) {
  set.seed(seed)

  samp_ids <- colnames(gene_mat)
  labels   <- labels[samp_ids]
  meta     <- meta[match(samp_ids, meta$sample_id), ]

  if (!is.null(meta$patient_id) && any(!is.na(meta$patient_id) & nzchar(meta$patient_id))) {
    dedup_keep <- rep(TRUE, length(samp_ids))
    for (g in unique(meta$gse_id)) {
      g_idx <- which(meta$gse_id == g)
      g_pid <- meta$patient_id[g_idx]
      dup <- duplicated(g_pid, incomparables = NA_character_)
      if (any(dup)) {
        message(sprintf("[DEDUP] %s: removing %d duplicate-patient sample(s).", g, sum(dup)))
        dedup_keep[g_idx[dup]] <- FALSE
      }
    }
    if (any(!dedup_keep)) {
      samp_ids <- samp_ids[dedup_keep]
      labels   <- labels[samp_ids]
      meta     <- meta[match(samp_ids, meta$sample_id), ]
    }
  }

  stratum   <- interaction(meta$gse_id, labels, drop = TRUE)
  train_idx <- logical(length(samp_ids))

  for (s in levels(stratum)) {
    idx <- which(stratum == s)
    n   <- length(idx)
    if (n == 0) next
    if (n == 1) {
      train_idx[idx] <- TRUE
      warning(sprintf("[SPLIT] Stratum '%s' has only 1 sample; assigned to train.", s))
      next
    }
    n_train <- max(1, min(n - 1, round(train_prop * n)))
    train_idx[sample(idx, n_train)] <- TRUE
  }

  
  for (g in unique(meta$gse_id)) {
    g_idx <- which(meta$gse_id == g)
    if (length(g_idx) >= 2 && (all(train_idx[g_idx]) || all(!train_idx[g_idx]))) {
      flip <- sample(g_idx, 1)
      train_idx[flip] <- !train_idx[flip]
      warning(sprintf(
        "[SPLIT] Dataset %s was entirely in one partition; reassigned 1 sample to guarantee both-side representation.",
        g
      ))
    }
  }

  train_ids <- samp_ids[train_idx]
  test_ids  <- samp_ids[!train_idx]

  train_mat <- gene_mat[, train_ids, drop = FALSE]
  test_mat  <- gene_mat[, test_ids,  drop = FALSE]

  labels_train <- droplevels(labels[train_ids])
  labels_test  <- droplevels(labels[test_ids])

  meta_train <- meta[match(train_ids, meta$sample_id), ]
  meta_test  <- meta[match(test_ids,  meta$sample_id), ]

  message(sprintf(
    "[SPLIT] Train: %d samples / %d datasets | Test: %d samples / %d datasets.",
    length(train_ids), length(unique(meta_train$gse_id)),
    length(test_ids),  length(unique(meta_test$gse_id))
  ))

  dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
  write.csv(test_mat, "outputs/gse_test_locked.csv")
  message("[SPLIT] Locked raw test partition -> outputs/gse_test_locked.csv")

  result <- list(
    train        = train_mat,
    test         = test_mat,
    labels_train = labels_train,
    labels_test  = labels_test,
    meta_train   = meta_train,
    meta_test    = meta_test
  )

  
  rm(train_mat, test_mat)
  gc(verbose = FALSE)

  result
}

harmonize_zscore <- function(train_mat, ...) {
  extra <- list(...)

  gene_mean <- rowMeans(train_mat, na.rm = TRUE)
  gene_sd   <- apply(train_mat, 1, sd, na.rm = TRUE)

  low_sd <- gene_sd < 1e-8
  if (any(low_sd)) {
    message(sprintf("[HARMONIZE] %d gene(s) had SD < 1e-8; replaced with SD = 1.", sum(low_sd)))
    gene_sd[low_sd] <- 1
  }

  apply_z <- function(mat) {
    mat <- mat[names(gene_mean), , drop = FALSE]
    sweep(sweep(mat, 1, gene_mean, "-"), 1, gene_sd, "/")
  }

  out <- list(train = apply_z(train_mat))
  for (nm in names(extra)) out[[nm]] <- apply_z(extra[[nm]])
  out
}

harmonize_tdm <- function(train_mat, query_mat, label = "query") {
  .ensure_pkg("preprocessCore", bioc = TRUE)
  suppressPackageStartupMessages(library(preprocessCore))

  common_genes <- intersect(rownames(train_mat), rownames(query_mat))
  train_sub <- as.matrix(train_mat[common_genes, , drop = FALSE])
  query_sub <- as.matrix(query_mat[common_genes, , drop = FALSE])

  
  target <- preprocessCore::normalize.quantiles.determine.target(
    train_sub, target.length = nrow(train_sub)
  )

  query_norm <- preprocessCore::normalize.quantiles.use.target(
    query_sub, target, copy = TRUE
  )
  dimnames(query_norm) <- dimnames(query_sub)

  message(sprintf(
    "[TDM] Mapped %s (%d samples, %d shared genes) onto training reference distribution.",
    label, ncol(query_sub), length(common_genes)
  ))

  query_norm
}

choose_harmonisation <- function(train_mat, test_mat, meta) {
  all_ids  <- c(colnames(train_mat), colnames(test_mat))
  meta_sub <- meta[match(all_ids, meta$sample_id), ]

  platforms <- unique(na.omit(meta_sub$platform))
  has_microarray <- any(grepl(
    "affy|microarray|agilent|illumina|processed|unknown",
    platforms, ignore.case = TRUE
  ))
  has_rnaseq     <- any(grepl("rnaseq|rna_seq|counts|voom|vst", platforms, ignore.case = TRUE))

  if (has_microarray && has_rnaseq) {
    message(sprintf(
      "[HARMONIZE] Both microarray/processed-expression and RNA-seq platforms present (%s) -> using harmonize_tdm() to correct cross-platform dynamic-range differences.",
      paste(platforms, collapse = ", ")
    ))
    train_h <- harmonize_tdm(train_mat, train_mat, label = "train")
    test_h  <- harmonize_tdm(train_mat, test_mat,  label = "test")
  } else {
    message(sprintf(
      "[HARMONIZE] Single platform type present (%s) -> using harmonize_zscore().",
      paste(platforms, collapse = ", ")
    ))
    z <- harmonize_zscore(train_mat, test = test_mat)
    train_h <- z$train
    test_h  <- z$test
  }

  list(train = train_h, test = test_h)
}