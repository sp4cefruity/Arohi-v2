suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

DATASET_LIST <- unique(c(
  "GSE25066", "GSE194040", "GSE41998", "GSE180962",
  "GSE20271",  "GSE34138",  "GSE50948", "GSE149322", "GSE22226",
  "GSE32603",  "GSE22358",  "GSE32646", "GSE16446",  "GSE130788",
  "GSE231629", "GSE123845", "GSE4779",  "GSE173839", "GSE22093",
  "GSE21997",  "GSE42822",  "GSE66399", "GSE181574", "GSE23988",
  "GSE41656",  "GSE8465",   "GSE122630","GSE21974",  "GSE207248",
  "GSE191127"
))

MIN_LABELS_PER_CLASS <- 5L    # datasets with fewer labeled samples per class are skipped
MIN_ALIGN_FRAC       <- 0.50  # expression-matrix columns must match pData rownames at >= this rate

.geo_fetch <- function(fun, max_attempts = 3L, base_delay = 1) {
  delay <- base_delay
  for (attempt in seq_len(max_attempts)) {
    res <- tryCatch(fun(), error = function(e) e)
    if (!inherits(res, "error")) return(res)
    warning(sprintf(
      "[GEO] Attempt %d/%d failed: %s; retrying in %ds...",
      attempt, max_attempts, conditionMessage(res), delay * 2
    ))
    Sys.sleep(delay * 2)
    delay <- delay * 2
  }
  stop(sprintf("[GEO] All %d attempts failed: %s", max_attempts, conditionMessage(res)))
}

download_supplementary_files <- function(gse_id, cache_dir = "data/raw/") {
  raw_dir <- file.path(cache_dir, gse_id)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  useful_regex <- "\\.(tar|cel)(\\.gz)?$|\\.(txt|csv|tsv)(\\.gz)?$"
  file_df <- tryCatch(
    .geo_fetch(function() GEOquery::getGEOSuppFiles(
      GEO = gse_id, baseDir = cache_dir, makeDirectory = TRUE,
      fetch_files = FALSE, filter_regex = useful_regex
    )),
    error = function(e) NULL
  )
  if (is.null(file_df) || !nrow(file_df)) return(FALSE)

  ok <- TRUE
  for (i in seq_len(nrow(file_df))) {
    destfile <- file.path(raw_dir, file_df$fname[i])
    if (file.exists(destfile)) {
      message(sprintf("[GEO] Using cached supplementary file: %s", destfile))
      next
    }
    message(sprintf("[GEO] Downloading %s ...", file_df$url[i]))
    res <- tryCatch(
      download.file(file_df$url[i], destfile, mode = "wb", method = "libcurl", quiet = TRUE),
      error = function(e) {
        warning(sprintf("[GEO] Download failed for %s: %s", file_df$url[i], conditionMessage(e)))
        NA_integer_
      }
    )
    if (is.na(res) || res != 0L) ok <- FALSE
  }
  ok
}

ingest_skip <- function(gse_id, reason) {
  structure(list(gse_id = gse_id, reason = reason), class = "ingest_skip")
}

find_count_file <- function(raw_dir) {
  candidates <- list.files(raw_dir, pattern = "count|matrix", full.names = TRUE,
                           ignore.case = TRUE, recursive = TRUE)
  candidates <- candidates[grepl("\\.(csv|txt|tsv)(\\.gz)?$", candidates, ignore.case = TRUE)]
  if (!length(candidates)) return(NULL)
  candidates[1]
}

find_processed_matrix_file <- function(raw_dir) {
  candidates <- list.files(
    raw_dir,
    pattern = "gene.*level|expression|expr|matrix|normalized|processed",
    full.names = TRUE, ignore.case = TRUE, recursive = TRUE
  )
  candidates <- candidates[grepl("\\.(csv|txt|tsv)(\\.gz)?$", candidates, ignore.case = TRUE)]
  if (!length(candidates)) return(NULL)
  pref <- candidates[grepl("gene.*level|mean|geneexp|expdat", basename(candidates), ignore.case = TRUE)]
  pref <- pref[!grepl("probe.?level|probe.?mat", basename(pref), ignore.case = TRUE)]
  if (length(pref)) candidates <- pref
  candidates[which.max(file.info(candidates)$size)]
}

matrix_from_series_eset <- function(eset_meta, gpl_id = NULL) {
  mat <- as.matrix(Biobase::exprs(eset_meta))
  if (!is.matrix(mat) || ncol(mat) == 0L || nrow(mat) == 0L)
    stop("[ERROR] Series matrix has no expression values.")

  dup <- duplicated(rownames(mat))
  if (sum(dup)) {
    message(sprintf("[INGEST] Series matrix: dropped %d duplicated feature row(s).", sum(dup)))
    mat <- mat[!dup, , drop = FALSE]
  }

  keep_row <- !apply(is.na(mat), 1L, all)
  mat <- mat[keep_row, , drop = FALSE]

  row_sd <- apply(mat, 1L, sd, na.rm = TRUE)
  mat <- mat[!is.na(row_sd) & row_sd > 1e-8, , drop = FALSE]
  if (!nrow(mat)) stop("[ERROR] No nonzero-variance rows in series matrix.")

  if (.looks_like_ensembl(rownames(mat)) && exists("map_ensembl_to_symbol", mode = "function")) {
    mat <- map_ensembl_to_symbol(mat)
  } else if (.looks_like_affy_probes(rownames(mat))) {
    chip <- .affy_chip_for_gpl(gpl_id)
    if (!is.null(chip) && exists("collapse_probes_to_symbols", mode = "function")) {
      message(sprintf("[INGEST] Series-matrix rows are Affy probes; collapsing via %s.", chip))
      mat <- collapse_probes_to_symbols(mat, chip = chip)
    } else if (exists("map_probes_to_symbols", mode = "function") && .try_gpl_map(mat, gpl_id)) {
      mapped <- map_probes_to_symbols(mat, gpl_id)
      if (!is.null(mapped)) mat <- mapped
    } else {
      warning("[WARN] Series-matrix rows look like Affy probes but no annotation package is mapped; keeping probe IDs.")
    }
  } else if (exists("map_probes_to_symbols", mode = "function") &&
             !is.null(gpl_id) && nzchar(gpl_id) &&
             .feature_space(rownames(mat)) %in% c("illumina_probe", "agilent_probe", "probe", "numeric")) {
    message(sprintf("[INGEST] Series-matrix rows are %s IDs; mapping via GPL %s table.",
                    .feature_space(rownames(mat)), gpl_id))
    mapped <- map_probes_to_symbols(mat, gpl_id)
    if (!is.null(mapped)) mat <- mapped
    else warning("[WARN] Could not map series-matrix probe IDs to gene symbols; keeping as-is.")
  }

  message(sprintf("[NORM] Series-matrix fallback: %d features x %d samples.", nrow(mat), ncol(mat)))
  mat
}

download_raw_and_normalise <- function(gse_id, platform, raw_dir, pheno_df, gpl_id) {
  if (platform == "affy_microarray") {
    tar_files <- list.files(raw_dir, pattern = "\\.(tar|tar\\.gz)$", full.names = TRUE,
                            ignore.case = TRUE, recursive = TRUE)
    for (tf in tar_files) {
      if (grepl("\\.tar\\.gz$", tf, ignore.case = TRUE)) {
        .ensure_pkg("R.utils", bioc = FALSE)
        R.utils::gunzip(tf, overwrite = TRUE, remove = FALSE)
        tf <- sub("\\.gz$", "", tf, ignore.case = TRUE)
      }
      untar(tf, exdir = raw_dir)
    }

    gz_files <- list.files(raw_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE,
                           ignore.case = TRUE, recursive = TRUE)
    if (length(gz_files)) {
      .ensure_pkg("R.utils", bioc = FALSE)
      for (gz in gz_files) R.utils::gunzip(gz, overwrite = TRUE, remove = FALSE)
    }

    data_path <- raw_dir
  } else {
    count_file <- find_count_file(raw_dir)
    if (is.null(count_file)) {
      warning(sprintf(
        "[WARN] No count-matrix-like file found for %s among supplementary files.", gse_id
      ))
      return(NULL)
    }
    if (grepl("\\.gz$", count_file, ignore.case = TRUE)) {
      .ensure_pkg("R.utils", bioc = FALSE)
      R.utils::gunzip(count_file, overwrite = TRUE, remove = FALSE)
      count_file <- sub("\\.gz$", "", count_file, ignore.case = TRUE)
    }
    data_path <- count_file

    if (grepl("tpm|fpkm|expr|normalized|processed|log", basename(data_path), ignore.case = TRUE)) {
      message(sprintf("[INGEST] '%s' looks like an already-processed expression matrix; loading as-is.",
                      basename(data_path)))
      platform <- "processed"
    }
  }

  norm_result <- run_normalization(
    data_path   = data_path,
    gse_id      = gse_id,
    sample_meta = pheno_df,
    platform    = platform,
    gpl_id      = gpl_id
  )
  norm_result$source <- if (platform == "affy_microarray") "cel" else "counts"
  norm_result
}

.align_frac <- function(mat, sample_meta) {
  if (is.null(mat) || !is.matrix(mat) || !ncol(mat)) return(0)
  if (is.null(sample_meta)) return(0)
  sum(colnames(mat) %in% rownames(sample_meta)) / ncol(mat)
}

.align_columns_to_pheno <- function(mat, sample_meta) {
  if (is.null(mat) || is.null(sample_meta)) return(mat)
  if (.align_frac(mat, sample_meta) >= MIN_ALIGN_FRAC) return(mat)
  for (cn in intersect(c("title", "source_name_ch1"), names(sample_meta))) {
    vals <- as.character(sample_meta[[cn]])
    idx  <- match(colnames(mat), vals)
    if (sum(!is.na(idx)) >= MIN_ALIGN_FRAC * ncol(mat)) {
      new_names <- rownames(sample_meta)[idx]
      new_names[is.na(idx)] <- colnames(mat)[is.na(idx)]
      colnames(mat) <- make.unique(new_names)
      message(sprintf("[INGEST] Realigned matrix columns to pData rownames via '%s'.", cn))
      return(mat)
    }
  }
  mat
}

.try_gpl_map <- function(mat, gpl_id) {
  !is.null(gpl_id) && nzchar(gpl_id) &&
    .feature_space(rownames(mat)) %in% c("affy_probe", "illumina_probe", "agilent_probe", "probe", "numeric")
}

.pick_main_eset <- function(esets) {
  if (methods::is(esets, "ExpressionSet")) return(esets)
  if (!is.list(esets) || length(esets) == 1L) return(esets[[1L]])
  rows <- vapply(esets, function(e) {
    tryCatch(nrow(Biobase::exprs(e)), error = function(x) 0L)
  }, integer(1))
  esets[[which.max(rows)]]
}

.geo_fetch_series <- function(gse_id) {
  res <- tryCatch(
    .geo_fetch(function() GEOquery::getGEO(GEO = gse_id, GSEMatrix = TRUE, getGPL = FALSE,
                                            destdir = .geo_cache_dir())),
    error = function(e) NULL
  )
  if (!is.null(res)) return(res)
  message(sprintf("[GEO] Network fetch failed for %s; loading cached series-matrix file(s).", gse_id))
  pat  <- sprintf("^%s(-GPL[0-9]+)?_series_matrix\\.txt\\.gz$", gse_id)
  files <- list.files(.geo_cache_dir(), pattern = pat, full.names = TRUE)
  if (!length(files)) {
    stop(sprintf("could not fetch GEO series matrix for %s (network failed and no local cache)", gse_id))
  }
  esets <- lapply(files, function(f) {
    message(sprintf("[GEO] Loading cached series matrix -> %s", f))
    GEOquery::getGEO(filename = f, getGPL = FALSE)
  })
  if (length(esets) == 1L) esets[[1L]] else esets
}

download_and_normalise <- function(gse_id, cache_dir = "data/raw/") {

  cache_out_dir <- "outputs/cache"
  dir.create(cache_out_dir, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(cache_out_dir, sprintf("%s_normalised.csv", gse_id))
  cache_gz   <- paste0(cache_file, ".gz")

  .cache_path <- function() if (file.exists(cache_gz)) cache_gz else cache_file
  .read_cache <- function(path) {
    con <- gzfile(path, "rt")
    on.exit(close(con), add = TRUE)
    as.matrix(read.csv(con, row.names = 1, check.names = FALSE))
  }
  .write_cache <- function(mat) {
    con <- gzfile(cache_gz, "wt")
    on.exit(close(con), add = TRUE)
    write.csv(mat, con)
  }
  .remove_cache <- function() {
    for (f in c(cache_file, cache_gz)) if (file.exists(f)) file.remove(f)
  }

  message(sprintf("[INGEST] === %s ===", gse_id))

  gse_meta <- tryCatch(
    .geo_fetch_series(gse_id),
    error = function(e) {
      warning(sprintf("[WARN] Failed to fetch metadata for %s: %s", gse_id, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(gse_meta)) {
    return(ingest_skip(gse_id, "could not fetch GEO series matrix"))
  }
  eset_meta <- .pick_main_eset(gse_meta)
  pheno_df  <- Biobase::pData(eset_meta)
  gpl_id    <- tryCatch(Biobase::annotation(eset_meta), error = function(e) "")

  platform <- detect_platform(gse_id, eset = eset_meta, gpl_id = gpl_id)

  if (platform == "cnv_array") {
    return(ingest_skip(gse_id, sprintf(
      "GPL %s is a SNP/copy-number (CNV) array, not gene expression; cannot merge with expression datasets",
      gpl_id
    )))
  }

  if (file.exists(cache_file) || file.exists(cache_gz)) {
    cache_use <- .cache_path()
    message(sprintf("[INGEST] Cache hit -> %s (skipping download).", cache_use))
    expr_mat <- .read_cache(cache_use)
    space <- .feature_space(rownames(expr_mat))
    if (space != "symbol" && .try_gpl_map(expr_mat, gpl_id) &&
        exists("map_probes_to_symbols", mode = "function")) {
      message(sprintf("[INGEST] Upgrading cached matrix from '%s' space to gene symbols...", space))
      mapped <- map_probes_to_symbols(expr_mat, gpl_id)
      if (!is.null(mapped)) {
        expr_mat <- mapped
        .write_cache(expr_mat)
      }
    } else if (platform == "affy_microarray" && .looks_like_affy_probes(rownames(expr_mat))) {
      chip <- .affy_chip_for_gpl(gpl_id)
      if (!is.null(chip) && exists("collapse_probes_to_symbols", mode = "function")) {
        message("[INGEST] Upgrading cached affy matrix from probe IDs to gene symbols...")
        expr_mat <- collapse_probes_to_symbols(expr_mat, chip = chip)
        .write_cache(expr_mat)
      }
    }
    expr_mat <- .align_columns_to_pheno(expr_mat, pheno_df)
    if (.align_frac(expr_mat, pheno_df) >= MIN_ALIGN_FRAC) {
      return(list(expr_mat = expr_mat, pheno_df = pheno_df, platform = platform,
                  source = "cache", gse_id = gse_id, gpl_id = gpl_id))
    }
    message(sprintf(
      "[INGEST] Cached matrix for %s does not align with pData (%.0f%%); removing cache and regenerating.",
      gse_id, 100 * .align_frac(expr_mat, pheno_df)
    ))
    .remove_cache()
  }

  message(sprintf("[INGEST] No cache found; downloading raw data for %s (%s)...", gse_id, platform))
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  raw_dir <- file.path(cache_dir, gse_id)
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

  supp <- NULL
  if (length(list.files(raw_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
    message(sprintf("[INGEST] Supplementary files already present in '%s'; skipping download.", raw_dir))
    supp <- TRUE
  } else {
    message(sprintf("[INGEST] Downloading supplementary files for %s...", gse_id))
    supp <- tryCatch(
      download_supplementary_files(gse_id, cache_dir),
      error = function(e) {
        warning(sprintf("[WARN] Supplementary download failed for %s: %s", gse_id, conditionMessage(e)))
        FALSE
      }
    )
    if (isFALSE(supp) && length(list.files(raw_dir, all.files = TRUE, no.. = TRUE)) == 0L) {
      message(sprintf(
        "[INGEST] No supplementary files obtained for %s; will try processed/series-matrix fallback.", gse_id
      ))
    }
  }

  norm_result <- NULL
  if (platform %in% c("affy_microarray", "rnaseq_counts")) {
    norm_result <- tryCatch(
      download_raw_and_normalise(gse_id, platform, raw_dir, pheno_df, gpl_id),
      error = function(e) {
        warning(sprintf("[WARN] Raw download/normalisation failed for %s: %s", gse_id, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(norm_result)) {
      norm_result$expr_mat <- .align_columns_to_pheno(norm_result$expr_mat, pheno_df)
      if (.align_frac(norm_result$expr_mat, pheno_df) < MIN_ALIGN_FRAC) {
        message(sprintf(
          "[INGEST] %s: raw-matrix sample IDs do not align with pData; trying series-matrix fallback.", gse_id
        ))
        norm_result <- NULL
      }
    }
  }

  if (is.null(norm_result)) {
    message(sprintf("[INGEST] No usable raw data for %s; trying processed-matrix fallback...", gse_id))
    proc_file <- find_processed_matrix_file(raw_dir)
    if (!is.null(proc_file)) {
      norm_result <- tryCatch(
        run_normalization(
          data_path   = proc_file,
          gse_id      = gse_id,
          sample_meta = pheno_df,
          platform    = "processed",
          gpl_id      = gpl_id
        ),
        error = function(e) {
          warning(sprintf("[WARN] Processed-matrix load failed for %s: %s", gse_id, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(norm_result)) {
        norm_result$source <- "processed"
        norm_result$expr_mat <- .align_columns_to_pheno(norm_result$expr_mat, pheno_df)
        if (.align_frac(norm_result$expr_mat, pheno_df) < MIN_ALIGN_FRAC) {
          message(sprintf(
            "[INGEST] %s: processed-matrix sample IDs do not align with pData; trying series-matrix fallback.", gse_id
          ))
          norm_result <- NULL
        }
      }
    }

    if (is.null(norm_result)) {
      message(sprintf("[INGEST] No supplementary matrix for %s; falling back to the GSE series matrix...", gse_id))
      series_mat <- tryCatch(
        matrix_from_series_eset(eset_meta, gpl_id = gpl_id),
        error = function(e) {
          warning(sprintf("[WARN] Series-matrix fallback failed for %s: %s", gse_id, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(series_mat)) {
        return(ingest_skip(gse_id, sprintf(
          "no CEL, count, or processed matrix file found under '%s' and series-matrix fallback failed", raw_dir
        )))
      }
      norm_result <- list(expr_mat = series_mat, platform = "processed", source = "series_matrix")
    }
  }

  expr_mat <- norm_result$expr_mat

  if (!is.matrix(expr_mat)) {
    return(ingest_skip(gse_id, "normalisation produced a non-matrix result"))
  }

  .write_cache(expr_mat)
  message(sprintf("[INGEST] Cached normalised matrix -> %s", cache_gz))

  list(expr_mat = expr_mat, pheno_df = pheno_df, platform = norm_result$platform,
       source = norm_result$source, gse_id = gse_id, gpl_id = gpl_id)
}


.classify_response_value <- function(s) {
  lc <- trimws(tolower(s))
  if (is.na(lc) || !nzchar(lc)) return(NA_character_)

  if (lc %in% c("na", "n/a", "nan", "-", "--", "---", "-1", "?", "none", "unavailable",
                "unknown", "not available", "unable to determine", "not evaluable",
                "ne", "not evaluated", "not assessed")) return(NA_character_)

  if (grepl("^[0-9.]+$", lc)) {
    if (lc == "1") return("pCR")
    if (lc == "0") return("RD")
    return(NA_character_)
  }

  if (lc %in% c("yes", "y")) return("pCR")
  if (lc %in% c("no", "n"))  return("RD")

  # near-CR / combined responder groups are intentionally left unlabeled
  if (grepl("near|npcr|\\bpcr\\b\\s*[+/]?\\s*npcr", lc)) return(NA_character_)
  # combined/ambiguous token such as "CR+kl"
  if (grepl("\\bcr\\b\\s*[+/]", lc)) return(NA_character_)

  rd_patterns <- c(
    "residual",
    "(no|non|not)[ _-]?(patholog(ic|ical)[ _-]?)?complete response",
    "(no|non|not)[ _-]?pcr",
    "no[ _-]?cr",
    "\\bnocr\\b",
    "\\bncr\\b",
    "(no|non|not)[ _-]?response",
    "(no|non|not)[ _-]?responder",
    "partial[ _-]?response",
    "progressive disease",
    "stable disease",
    "pr[+ ]nr",
    "primary refrac",
    "\\bpd\\b",
    "\\bnr\\b",
    "\\bpr\\b",
    "\\brd\\b"
  )
  if (any(grepl(paste(rd_patterns, collapse = "|"), lc, perl = TRUE))) return("RD")

  pcr_patterns <- c(
    "complete response",
    "complete[ _-]?pathologic|pathologic[ _-]?complete",
    "\\bpcr\\b",
    "\\bcr\\b"
  )
  if (any(grepl(paste(pcr_patterns, collapse = "|"), lc, perl = TRUE))) return("pCR")

  NA_character_
}

.parse_response_col <- function(x) {
  v <- as.character(x)
  lab <- vapply(v, .classify_response_value, character(1))
  factor(unname(lab), levels = c("pCR", "RD"))
}

.content_score <- function(x) {
  lab <- .parse_response_col(x)
  total <- sum(!is.na(as.character(x)) & nzchar(trimws(as.character(x))))
  if (!total) return(0)
  sum(!is.na(lab)) / total
}

.score_col_name <- function(cn) {
  lc <- tolower(cn)
  s <- 0
  if (grepl("pcr|patholog", lc))             s <- s + 100
  if (grepl("complete response", lc))        s <- s + 100
  if (grepl("residual", lc))                 s <- s + 60
  if (grepl("response", lc))                 s <- s + 20
  if (grepl("pred|predict|dlda|model|signature|score", lc)) s <- s - 200
  if (grepl(paste0(
    "rcb|recurrence|rfs|dfs|survival|surv|event|death|metasta|grade|stage|",
    "er[ _]?status|her2|pr[ _]?status|estrogen|progester|age|histolog|subtype|",
    "molecular|tumou?r|size|biopsy|tissue|lymph|node|cycle|visit|timepoint|",
    "time point|batch|treatment|chemo|drug|dose|side.?effect|toxicity|",
    "radiotherap|radiation|distant"
  ), lc)) s <- s - 100
  s
}

extract_labels <- function(pheno_df, gse_id) {
  cols <- names(pheno_df)

  name_score <- vapply(cols, .score_col_name, numeric(1))
  content    <- vapply(cols, function(cn) .content_score(pheno_df[[cn]]), numeric(1))

  primary  <- cols[name_score > 0]
  fallback <- cols[name_score <= 0 & content >= 0.5]

  candidate_cols <- c(primary, fallback)
  if (!length(candidate_cols)) {
    warning(sprintf("[WARN] No pCR/RD outcome column found for %s; returning NULL.", gse_id))
    return(NULL)
  }

  ord <- order(name_score[candidate_cols], content[candidate_cols], decreasing = TRUE, method = "radix")
  candidate_cols <- candidate_cols[ord]

  match_col <- NULL
  lab <- NULL
  for (cn in candidate_cols) {
    l <- .parse_response_col(pheno_df[[cn]])
    if (sum(!is.na(l)) > 0) { match_col <- cn; lab <- l; break }
  }
  if (is.null(match_col)) {
    warning(sprintf("[WARN] No interpretable pCR/RD values in any candidate column for %s; returning NULL.", gse_id))
    return(NULL)
  }

  names(lab) <- rownames(pheno_df)

  message(sprintf(
    "[LABELS] %s: matched column '%s' -> %d pCR, %d RD, %d unmapped.",
    gse_id, match_col, sum(lab == "pCR", na.rm = TRUE),
    sum(lab == "RD", na.rm = TRUE), sum(is.na(lab))
  ))

  lab
}

ensure_symbol_space <- function(res) {
  mat   <- res$expr_mat
  gpl   <- res$gpl_id
  space <- .feature_space(rownames(mat))
  message(sprintf("[MERGE] %s: feature space = '%s' (%d features).", res$gse_id, space, nrow(mat)))

  if (space == "ensembl" && exists("map_ensembl_to_symbol", mode = "function")) {
    mat <- tryCatch(map_ensembl_to_symbol(mat), error = function(e) {
      warning(sprintf("[WARN] Ensembl mapping failed for %s: %s", res$gse_id, conditionMessage(e)))
      mat
    })
  }

  if (.feature_space(rownames(mat)) != "symbol" && exists("map_probes_to_symbols", mode = "function")) {
    mapped <- tryCatch(map_probes_to_symbols(mat, gpl), error = function(e) {
      warning(sprintf("[WARN] GPL mapping failed for %s: %s", res$gse_id, conditionMessage(e)))
      NULL
    })
    if (!is.null(mapped)) mat <- mapped
  }

  space2 <- .feature_space(rownames(mat))
  if (space2 != "symbol") {
    return(list(skip = TRUE, reason = sprintf(
      "features are in '%s' space and could not be mapped to gene symbols", space2
    )))
  }
  res$expr_mat <- mat
  res
}

infer_patient_id <- function(pheno, samp_ids) {
  if (is.null(pheno) || is.null(samp_ids)) return(rep(NA_character_, length(samp_ids)))
  hit <- names(pheno)[sapply(names(pheno), function(cn)
    grepl("patient|subject|case|participant|sample.?id", cn, ignore.case = TRUE))]
  if (length(hit)) {
    v <- sub("^.*id[: ]?", "", as.character(pheno[[hit[1]]]))
    v <- sub("^ISPY2[_ ]?", "", v)
    v <- trimws(v)
    return(v[match(samp_ids, rownames(pheno))])
  }
  if ("title" %in% names(pheno)) {
    v <- sub("_[0-9]+$", "", as.character(pheno$title))
    return(v[match(samp_ids, rownames(pheno))])
  }
  rep(NA_character_, length(samp_ids))
}

merge_datasets <- function(dataset_list) {
  dir.create("outputs/qc", showWarnings = FALSE, recursive = TRUE)
  dir.create("data/raw",   showWarnings = FALSE, recursive = TRUE)

  loaded       <- list()
  skip_reasons <- list()

  for (gse_id in dataset_list) {
    res <- tryCatch(
      download_and_normalise(gse_id),
      error = function(e) {
        warning(sprintf("[WARN] %s failed during download/normalise: %s", gse_id, conditionMessage(e)))
        NULL
      }
    )
    if (inherits(res, "ingest_skip")) {
      skip_reasons[[res$gse_id]] <- res$reason
      message(sprintf("[SKIP] %s: %s", res$gse_id, res$reason))
      next
    }
    if (is.null(res)) {
      skip_reasons[[gse_id]] <- "unexpected failure (see warnings above)"
      message(sprintf("[SKIP] %s: unexpected failure (see warnings above)", gse_id))
      next
    }

    labels <- tryCatch(
      extract_labels(res$pheno_df, gse_id),
      error = function(e) {
        warning(sprintf("[WARN] Label extraction failed for %s: %s", gse_id, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(labels)) {
      skip_reasons[[gse_id]] <- "label extraction failed (no pCR/RD outcome column matched)"
      message(sprintf(
        "[SKIP] %s: label extraction failed (no pCR/RD outcome column matched)", gse_id
      ))
      next
    }

    n_pcr <- sum(labels == "pCR", na.rm = TRUE)
    n_rd  <- sum(labels == "RD", na.rm = TRUE)
    if (n_pcr < MIN_LABELS_PER_CLASS || n_rd < MIN_LABELS_PER_CLASS) {
      reason <- sprintf(
        "too few labeled samples per class (pCR=%d, RD=%d); need >= %d of each",
        n_pcr, n_rd, MIN_LABELS_PER_CLASS
      )
      skip_reasons[[gse_id]] <- reason
      message(sprintf("[SKIP] %s: %s", gse_id, reason))
      next
    }

    n_aligned <- sum(!is.na(labels[colnames(res$expr_mat)]))
    if (n_aligned < MIN_LABELS_PER_CLASS) {
      reason <- sprintf(
        "expression-matrix sample IDs do not align with pData labels (%d aligned labeled samples)",
        n_aligned
      )
      skip_reasons[[gse_id]] <- reason
      message(sprintf("[SKIP] %s: %s", gse_id, reason))
      next
    }

    res$labels <- labels
    res <- ensure_symbol_space(res)
    if (isTRUE(res$skip)) {
      skip_reasons[[gse_id]] <- res$reason
      message(sprintf("[SKIP] %s: %s", gse_id, res$reason))
      next
    }

    loaded[[gse_id]] <- res
  }

  if (length(skip_reasons)) {
    skip_df <- data.frame(
      gse_id = names(skip_reasons),
      reason = unlist(skip_reasons, use.names = FALSE),
      stringsAsFactors = FALSE
    )
    write.csv(skip_df, "outputs/qc/skipped_datasets.csv", row.names = FALSE)
    write.csv(skip_df, "outputs/qc/skipped_datasets.txt", row.names = FALSE)
  }

  if (!length(loaded)) stop("[ERROR] No datasets successfully loaded; cannot merge.")

  MIN_SHARED_GENES <- 1000L

  gene_sets    <- lapply(loaded, function(d) rownames(d$expr_mat))
  shared_genes <- Reduce(intersect, gene_sets)
  message(sprintf(
    "[MERGE] Gene symbol intersection across %d loaded dataset(s): %d genes.",
    length(loaded), length(shared_genes)
  ))

  if (length(shared_genes) < MIN_SHARED_GENES) {
    message(sprintf(
      "[MERGE] Full intersection too small (%d genes); identifying dataset(s) that share too few genes with the cohort...",
      length(shared_genes)
    ))

    gids  <- names(gene_sets)
    n_ids <- length(gids)

    overlap_ok <- matrix(FALSE, n_ids, n_ids, dimnames = list(gids, gids))
    for (i in seq_len(n_ids)) {
      for (j in seq_len(i - 1L)) {
        if (length(intersect(gene_sets[[gids[i]]], gene_sets[[gids[j]]])) >= MIN_SHARED_GENES) {
          overlap_ok[i, j] <- overlap_ok[j, i] <- TRUE
        }
      }
    }

    largest_component <- function(adj) {
      best_idx <- NULL
      best_n   <- 0L
      for (start in seq_len(nrow(adj))) {
        in_comp <- rep(FALSE, nrow(adj))
        in_comp[start] <- TRUE
        repeat {
          expand <- rowSums(adj[, in_comp, drop = FALSE]) > 0
          added  <- expand & !in_comp
          if (!any(added)) break
          in_comp <- in_comp | expand
        }
        if (sum(in_comp) > best_n) { best_n <- sum(in_comp); best_idx <- in_comp }
      }
      best_idx
    }

    keep <- largest_component(overlap_ok)
    for (g in gids[!keep]) {
      skip_reasons[[g]] <- sprintf(
        "shares < %d gene symbols with other loaded datasets (incompatible feature space); dropped at merge",
        MIN_SHARED_GENES
      )
      message(sprintf(
        "[MERGE] Dropping %s from merge: shares < %d gene symbols with the retained cohort.",
        g, MIN_SHARED_GENES
      ))
    }

    loaded       <- loaded[gids[keep]]
    gene_sets    <- gene_sets[gids[keep]]
    shared_genes <- Reduce(intersect, gene_sets)
    message(sprintf(
      "[MERGE] After dropping incompatible dataset(s): intersection across %d retained dataset(s) = %d genes.",
      length(loaded), length(shared_genes)
    ))

    if (length(skip_reasons)) {
      skip_df <- data.frame(
        gse_id = names(skip_reasons),
        reason = unlist(skip_reasons, use.names = FALSE),
        stringsAsFactors = FALSE
      )
      write.csv(skip_df, "outputs/qc/skipped_datasets.csv", row.names = FALSE)
      write.csv(skip_df, "outputs/qc/skipped_datasets.txt", row.names = FALSE)
    }
  }

  if (length(shared_genes) < MIN_SHARED_GENES || length(loaded) < 3L) {
    stop(sprintf(
      "[ERROR] Gene intersection across %d retained dataset(s) is too small (%d genes) for a robust merge (%s). See outputs/qc/skipped_datasets.csv.",
      length(loaded), length(shared_genes), paste(names(loaded), collapse = ", ")
    ))
  }

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
      patient_id  = infer_patient_id(d$pheno_df, samp_ids),
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

  dedup_keep <- rep(TRUE, nrow(combined_meta))
  pid <- combined_meta$patient_id
  if (any(!is.na(pid) & nzchar(pid))) {
    for (g in unique(combined_meta$gse_id)) {
      g_idx <- which(combined_meta$gse_id == g)
      g_pid <- pid[g_idx]
      dup <- duplicated(g_pid, incomparables = NA_character_)
      if (any(dup)) {
        message(sprintf(
          "[DEDUP] %s: removing %d duplicate-patient sample(s) (keeping first per patient).",
          g, sum(dup)
        ))
        dedup_keep[g_idx[dup]] <- FALSE
      }
    }
  }
  if (any(!dedup_keep)) {
    merged_expr   <- merged_expr[, dedup_keep, drop = FALSE]
    combined_meta <- combined_meta[dedup_keep, ]
    pre_names     <- names(final_labels)
    final_labels  <- factor(as.character(final_labels)[dedup_keep], levels = c("pCR", "RD"))
    names(final_labels) <- pre_names[dedup_keep]
    message(sprintf("[MERGE] Patient-level dedup removed %d sample(s); %d remain.",
                    sum(!dedup_keep), length(final_labels)))
  }

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