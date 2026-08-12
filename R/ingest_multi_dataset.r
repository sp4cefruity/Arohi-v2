suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

DATASET_LIST <- unique(c(
  "GSE25066", "GSE194040", "GSE16716", "GSE41998", "GSE180962",
  "GSE20271",  "GSE34138",  "GSE50948", "GSE149322", "GSE22226",
  "GSE32603",  "GSE22358",  "GSE32646", "GSE16446",  "GSE130788",
  "GSE231629", "GSE123845", "GSE4779",  "GSE173839", "GSE22093",
  "GSE21997",  "GSE42822",  "GSE66399", "GSE181574", "GSE23988",
  "GSE41656",  "GSE8465",   "GSE122630","GSE21974",  "GSE207248",
  "GSE191127"
))

download_and_normalise <- function(gse_id, cache_dir = "data/raw/") {

  cache_out_dir <- "outputs/cache"
  dir.create(cache_out_dir, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(cache_out_dir, sprintf("%s_normalised.csv", gse_id))

  message(sprintf("[INGEST] === %s ===", gse_id))

  gse_meta <- tryCatch(
    GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE),
    error = function(e) {
      warning(sprintf("[WARN] Failed to fetch metadata for %s: %s", gse_id, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(gse_meta)) return(NULL)
  eset_meta <- if (is.list(gse_meta)) gse_meta[[1]] else gse_meta
  pheno_df  <- Biobase::pData(eset_meta)

  platform <- detect_platform(gse_id)

  if (file.exists(cache_file)) {
    message(sprintf("[INGEST] Cache hit -> %s (skipping download).", cache_file))
    expr_mat <- as.matrix(read.csv(cache_file, row.names = 1, check.names = FALSE))
    return(list(expr_mat = expr_mat, pheno_df = pheno_df, platform = platform, gse_id = gse_id))
  }

  message(sprintf("[INGEST] No cache found; downloading raw data for %s (%s)...", gse_id, platform))
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  raw_dir <- file.path(cache_dir, gse_id)
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

  data_path <- NULL

  if (platform == "affy_microarray") {
    supp <- tryCatch(
      GEOquery::getGEOSuppFiles(gse_id, baseDir = cache_dir, makeDirectory = TRUE),
      error = function(e) {
        warning(sprintf("[WARN] getGEOSuppFiles() failed for %s: %s", gse_id, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(supp)) return(NULL)

    tar_files <- list.files(raw_dir, pattern = "\\.tar$", full.names = TRUE)
    for (tf in tar_files) untar(tf, exdir = raw_dir)

    gz_files <- list.files(raw_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE)
    if (length(gz_files)) {
      .ensure_pkg("R.utils", bioc = FALSE)
      for (gz in gz_files) R.utils::gunzip(gz, overwrite = TRUE, remove = FALSE)
    }
    data_path <- raw_dir

  } else {
    
    supp <- tryCatch(
      GEOquery::getGEOSuppFiles(gse_id, baseDir = cache_dir, makeDirectory = TRUE),
      error = function(e) {
        warning(sprintf("[WARN] getGEOSuppFiles() failed for %s: %s", gse_id, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(supp)) return(NULL)

    count_candidates <- list.files(raw_dir, pattern = "count|matrix", full.names = TRUE, ignore.case = TRUE)
    count_candidates <- count_candidates[grepl("\\.(csv|txt|tsv)(\\.gz)?$", count_candidates, ignore.case = TRUE)]

    if (!length(count_candidates)) {
      warning(sprintf(
        "[WARN] No count-matrix-like file found for %s among supplementary files; skipping.", gse_id
      ))
      return(NULL)
    }

    count_file <- count_candidates[1]
    if (grepl("\\.gz$", count_file, ignore.case = TRUE)) {
      .ensure_pkg("R.utils", bioc = FALSE)
      R.utils::gunzip(count_file, overwrite = TRUE, remove = FALSE)
      count_file <- sub("\\.gz$", "", count_file, ignore.case = TRUE)
    }
    data_path <- count_file
  }

  norm_result <- tryCatch(
    run_normalization(data_path = data_path, gse_id = gse_id, sample_meta = pheno_df),
    error = function(e) {
      warning(sprintf("[WARN] Normalisation failed for %s: %s", gse_id, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(norm_result)) return(NULL)

  expr_mat <- norm_result$expr_mat
  
  if (!is.matrix(expr_mat)) {
    warning(sprintf("[WARN] Unexpected non-matrix normalisation result for %s; skipping.", gse_id))
    return(NULL)
  }

  write.csv(expr_mat, cache_file)
  message(sprintf("[INGEST] Cached normalised matrix -> %s", cache_file))

  list(expr_mat = expr_mat, pheno_df = pheno_df, platform = norm_result$platform, gse_id = gse_id)
}


extract_labels <- function(pheno_df, gse_id) {
  keywords     <- c("pathologic", "pcr", "residual", "response")
  keyword_rgx  <- paste(keywords, collapse = "|")

  
  name_hits <- names(pheno_df)[
    sapply(names(pheno_df), function(cn) grepl(keyword_rgx, cn, ignore.case = TRUE))
  ]

  candidate_cols <- name_hits

  
  if (!length(candidate_cols)) {
    candidate_cols <- names(pheno_df)[
      sapply(pheno_df, function(col) any(grepl(keyword_rgx, as.character(col), ignore.case = TRUE)))
    ]
  }

  if (!length(candidate_cols)) {
    warning(sprintf("[WARN] No pCR/RD outcome column found for %s; returning NULL.", gse_id))
    return(NULL)
  }

  match_col <- candidate_cols[1]
  raw_vals  <- as.character(pheno_df[[match_col]])

  
  rd_pattern  <- "residual|non[ _-]?pcr|no[ _-]?pcr|(?<![a-z])rd(?![a-z])"
  pcr_pattern <- "complete response|(?<![a-z])pcr(?![a-z])|(?<![a-z])cr(?![a-z])"

  is_rd  <- grepl(rd_pattern,  raw_vals, ignore.case = TRUE, perl = TRUE)
  is_pcr <- grepl(pcr_pattern, raw_vals, ignore.case = TRUE, perl = TRUE) & !is_rd

  label <- rep(NA_character_, length(raw_vals))
  label[is_pcr] <- "pCR"
  label[is_rd]  <- "RD"

  result <- factor(label, levels = c("pCR", "RD"))
  names(result) <- rownames(pheno_df)

  message(sprintf(
    "[LABELS] %s: matched column '%s' -> %d pCR, %d RD, %d unmapped.",
    gse_id, match_col, sum(result == "pCR", na.rm = TRUE),
    sum(result == "RD", na.rm = TRUE), sum(is.na(result))
  ))

  result
}

merge_datasets <- function(dataset_list) {
  dir.create("outputs/qc", showWarnings = FALSE, recursive = TRUE)
  dir.create("data/raw",   showWarnings = FALSE, recursive = TRUE)

  loaded  <- list()
  skipped <- character(0)

  for (gse_id in dataset_list) {
    res <- tryCatch(
      download_and_normalise(gse_id),
      error = function(e) {
        warning(sprintf("[WARN] %s failed during download/normalise: %s", gse_id, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res)) { skipped <- c(skipped, gse_id); next }

    labels <- tryCatch(
      extract_labels(res$pheno_df, gse_id),
      error = function(e) {
        warning(sprintf("[WARN] Label extraction failed for %s: %s", gse_id, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(labels)) { skipped <- c(skipped, gse_id); next }

    res$labels <- labels
    loaded[[gse_id]] <- res
  }

  if (length(skipped)) {
    writeLines(skipped, "outputs/qc/skipped_datasets.txt")
    message(sprintf(
      "[MERGE] Skipped %d dataset(s) (see outputs/qc/skipped_datasets.txt): %s",
      length(skipped), paste(skipped, collapse = ", ")
    ))
  }

  if (!length(loaded)) stop("[ERROR] No datasets successfully loaded; cannot merge.")

  
  gene_sets     <- lapply(loaded, function(d) rownames(d$expr_mat))
  shared_genes  <- Reduce(intersect, gene_sets)
  message(sprintf(
    "[MERGE] Gene symbol intersection across %d loaded dataset(s): %d genes.",
    length(loaded), length(shared_genes)
  ))
  if (!length(shared_genes)) stop("[ERROR] Gene intersection across datasets is empty.")

  expr_list   <- lapply(loaded, function(d) d$expr_mat[shared_genes, , drop = FALSE])
  merged_expr <- do.call(cbind, expr_list)

  get_field <- function(pheno, samp_ids, patterns) {
    rgx <- paste(patterns, collapse = "|")
    hit <- names(pheno)[sapply(names(pheno), function(cn) grepl(rgx, cn, ignore.case = TRUE))]
    if (!length(hit)) return(rep(NA_character_, length(samp_ids)))
    as.character(pheno[[hit[1]]])[match(samp_ids, rownames(pheno))]
  }

  meta_rows <- list()
  label_vec <- character(0)

  for (gse_id in names(loaded)) {
    d <- loaded[[gse_id]]
    samp_ids <- colnames(d$expr_mat)
    lab      <- as.character(d$labels[samp_ids])

    meta_rows[[gse_id]] <- data.frame(
      sample_id   = samp_ids,
      gse_id      = gse_id,
      platform    = d$platform,
      response    = lab,
      age         = get_field(d$pheno_df, samp_ids, c("age")),
      stage       = get_field(d$pheno_df, samp_ids, c("stage")),
      er_status   = get_field(d$pheno_df, samp_ids, c("er[ _-]?status", "estrogen")),
      pr_status   = get_field(d$pheno_df, samp_ids, c("pr[ _-]?status", "progesterone")),
      her2_status = get_field(d$pheno_df, samp_ids, c("her2")),
      grade       = get_field(d$pheno_df, samp_ids, c("grade")),
      stringsAsFactors = FALSE
    )

    names(lab) <- samp_ids
    label_vec  <- c(label_vec, lab)
  }

  combined_meta <- do.call(rbind, meta_rows)
  rownames(combined_meta) <- NULL

  keep <- !is.na(label_vec)
  n_dropped <- sum(!keep)
  if (n_dropped) message(sprintf("[MERGE] Dropping %d sample(s) with NA/unmapped response.", n_dropped))

  merged_expr   <- merged_expr[, keep, drop = FALSE]
  combined_meta <- combined_meta[keep, ]
  final_labels  <- factor(label_vec[keep], levels = c("pCR", "RD"))
  names(final_labels) <- names(label_vec)[keep]

  n_pcr   <- sum(final_labels == "pCR")
  n_rd    <- sum(final_labels == "RD")
  pct_pcr <- 100 * n_pcr / (n_pcr + n_rd)
  message(sprintf(
    "[MERGE] Final class balance: n_pCR = %d, n_RD = %d (%.1f%% pCR), total n = %d",
    n_pcr, n_rd, pct_pcr, n_pcr + n_rd
  ))

  dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
  write.csv(merged_expr, "outputs/merged_expression.csv")
  write.csv(
    data.frame(sample_id = names(final_labels), response = as.character(final_labels)),
    "data/raw/combined_labels.csv", row.names = FALSE
  )
  write.csv(combined_meta, "data/raw/combined_metadata.csv", row.names = FALSE)

  message("[MERGE] Saved outputs/merged_expression.csv, data/raw/combined_labels.csv, data/raw/combined_metadata.csv")

  list(expr_mat = merged_expr, labels = final_labels, metadata = combined_meta)
}

plot_dataset_summary <- function(metadata) {
  dir.create("outputs/qc", showWarnings = FALSE, recursive = TRUE)
  out_file <- "outputs/qc/dataset_composition.png"

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ggplot2))
    p <- ggplot(metadata, aes(x = gse_id, fill = response)) +
      geom_bar(position = "stack") +
      labs(
        x = "Dataset (GSE)", y = "Number of samples", fill = "Response",
        title = "Sample composition per dataset"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    ggsave(out_file, plot = p, width = 12, height = 6, dpi = 150)
  } else {
    tab <- table(metadata$response, metadata$gse_id)
    png(out_file, width = 1400, height = 700)
    barplot(
      tab, beside = FALSE, legend.text = TRUE, las = 2,
      col = c("#2b8cbe", "#de2d26"),
      main = "Sample composition per dataset", ylab = "Number of samples"
    )
    dev.off()
  }

  message(sprintf("[QC] Saved dataset composition plot -> %s", out_file))
  invisible(out_file)
}