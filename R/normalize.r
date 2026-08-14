.ensure_pkg <- function(pkgs, bioc = FALSE) {
  need <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (!length(need)) return(invisible())

  if (!exists("check_project_dependencies", mode = "function")) {
    source(file.path("R", "dependencies.R"))
  }

  if (bioc) {
    message(sprintf("[DEPENDENCIES] Missing BioConductor package(s): %s", paste(need, collapse = ", ")))
    install_project_dependencies(need)
  } else {
    install.packages(need, repos = "https://cloud.r-project.org")
  }
  invisible(TRUE)
}

suppressPackageStartupMessages({
  library("Biobase")
  library("limma")
  library("edgeR")
  library("DESeq2")
  library("GEOquery")
})

# NCBI GEO is served over HTTPS; R's default 'wininet' download method stalls
# on long transfers. libcurl is reliable and honours options(timeout) below.
options(download.file.method = "libcurl")
options(timeout = 300)

DATASET_REGISTRY <- list(
  affy_microarray = list(
    label    = "Affymetrix microarray (RMA) -- limma downstream",
    packages = list(bioc = c("affy", "affyio", "limma")),
    fn       = "normalize_microarray_affy",
    id_type  = "probe"
  ),

  rnaseq_counts = list(
    label    = "RNA-seq raw counts (DESeq2 VST)",
    packages = list(bioc = c("DESeq2", "edgeR")),
    fn       = "normalize_rnaseq_vst",
    id_type  = "gene"
  ),

  rnaseq_voom = list(
    label    = "RNA-seq raw counts (limma-voom + TMM)",
    packages = list(bioc = c("limma", "edgeR")),
    fn       = "normalize_rnaseq_voom",
    id_type  = "gene"
  ),

  agilent = list(
    label    = "Agilent microarray -- processed matrix (as-is)",
    packages = list(),
    fn       = "load_processed_matrix",
    id_type  = "gene"
  ),

  illumina = list(
    label    = "Illumina microarray -- processed matrix (as-is)",
    packages = list(),
    fn       = "load_processed_matrix",
    id_type  = "gene"
  ),

  processed = list(
    label    = "Processed expression matrix (as-is)",
    packages = list(),
    fn       = "load_processed_matrix",
    id_type  = "gene"
  )
)

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

.robust_fetch <- function(fun) {
  if (exists(".geo_fetch", mode = "function", inherits = TRUE)) {
    .geo_fetch(fun)
  } else {
    fun()
  }
}

.geo_cache_dir <- function() {
  dir <- file.path("outputs", "cache", "geo")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

.gpl_desc <- local({
  cache_file <- file.path("outputs", "cache", "geo", "gpl_desc.csv")

  .read <- function() {
    if (!file.exists(cache_file)) return(NULL)
    df <- read.csv(cache_file, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c("gpl_id", "description") %in% names(df))) {
      setNames(df$description, df$gpl_id)
    } else NULL
  }

  .write <- function(map) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    write.csv(
      data.frame(gpl_id = names(map), description = unname(map),
                 stringsAsFactors = FALSE),
      cache_file, row.names = FALSE
    )
  }

  function(gpl_id) {
    cache <- .read()
    if (!is.null(cache) && gpl_id %in% names(cache)) {
      return(unname(cache[[gpl_id]]))
    }

    desc <- ""
    gpl_obj <- tryCatch(
      .robust_fetch(function() GEOquery::getGEO(gpl_id, destdir = .geo_cache_dir())),
      error = function(e) NULL
    )
    if (!is.null(gpl_obj)) {
      meta <- tryCatch(GEOquery::Meta(gpl_obj), error = function(e) list())
      first_string <- function(x) {
        if (is.null(x)) return("")
        v <- unlist(x)
        if (!length(v)) return("")
        paste(trimws(v), collapse = " ")
      }
      desc <- paste(
        first_string(meta$title),
        first_string(meta$description),
        first_string(meta$technology),
        collapse = " "
      )
    }

    cache <- c(cache, setNames(desc, gpl_id))
    .write(cache)
    desc
  }
})

detect_platform <- function(gse_id, eset = NULL, gpl_id = NULL) {
  suppressPackageStartupMessages(library(GEOquery))

  known_affy <- c("GPL96", "GPL570", "GPL571", "GPL201")
  known_cnv  <- c("GPL6801", "GPL3720", "GPL5175", "GPL19138", "GPL19730",
                  "GPL21337", "GPL2823", "GPL11664", "GPL3921", "GPL8775")

  if (is.null(eset)) {
    message(sprintf("[PLATFORM] Fetching GEO metadata for %s...", gse_id))
    gse_obj <- tryCatch(
      .robust_fetch(function() GEOquery::getGEO(GEO = gse_id, GSEMatrix = TRUE,
                                                 getGPL = FALSE, destdir = .geo_cache_dir())),
      error = function(e) {
        warning(sprintf(
          "[WARN] Could not fetch GEO metadata for %s (%s).",
          gse_id, conditionMessage(e)
        ))
        NULL
      }
    )
    if (is.null(gse_obj)) return("unknown")
    eset <- if (is.list(gse_obj)) gse_obj[[1]] else gse_obj
  }

  if (is.null(gpl_id) || !nzchar(gpl_id)) {
    gpl_id <- tryCatch(Biobase::annotation(eset), error = function(e) "")
  }

  if (!nzchar(gpl_id)) {
    warning(sprintf("[WARN] No GPL annotation found for %s.", gse_id))
    return("unknown")
  }

  message(sprintf("[PLATFORM] %s -> %s", gse_id, gpl_id))

  if (gpl_id %in% known_affy) {
    return("affy_microarray")
  }

  desc <- tolower(.gpl_desc(gpl_id))

  if (gpl_id %in% known_cnv ||
      grepl("genomewide?snp|snp[ _]?6|copy.?number|cnv|cytoscan|genotyping|snp array", desc)) {
    message(sprintf(
      "[PLATFORM] %s -> %s detected as a SNP/copy-number (CNV) platform; not gene-expression.",
      gse_id, gpl_id
    ))
    return("cnv_array")
  }

  if (grepl("agilent", desc))            return("agilent")
  if (grepl("illumina", desc))           return("illumina")
  if (grepl("affymetrix|genechip", desc)) return("affy_microarray")
  if (grepl("rna-seq|rna seq|counts", desc)) return("rnaseq_counts")

  warning(sprintf(
    "[WARN] Could not confidently classify platform for %s (GPL %s, description: '%s'); will try processed-matrix fallback.",
    gse_id, gpl_id, substr(desc, 1, 200)
  ))
  "unknown"
}

`%||%` <- function(a, b) if (is.null(a) || !nzchar(a)) b else a


run_normalization <- function(data_path, gse_id, sample_meta = NULL,
                              platform = NULL, gpl_id = NULL, ...) {

  if (is.null(platform)) platform <- detect_platform(gse_id)

  if (!platform %in% names(DATASET_REGISTRY))
    stop(sprintf(
      "[ERROR] detect_platform() returned unknown platform: '%s'\nRegistered: %s",
      platform, paste(names(DATASET_REGISTRY), collapse = ", ")
    ))

  cfg <- DATASET_REGISTRY[[platform]]
  message(sprintf("[NORM] Dataset : %s", gse_id))
  message(sprintf("[NORM] Type    : %s", cfg$label))
  message(sprintf("[NORM] Path    : %s", data_path))

  if (!is.null(cfg$packages$bioc)) .ensure_pkg(cfg$packages$bioc, bioc = TRUE)
  if (!is.null(cfg$packages$cran)) .ensure_pkg(cfg$packages$cran, bioc = FALSE)

  norm_fn  <- match.fun(cfg$fn)
  expr_mat <- norm_fn(data_path = data_path, sample_meta = sample_meta, gpl_id = gpl_id, ...)

  message(sprintf(
    "[NORM] Output  : %d features x %d samples.", nrow(expr_mat), ncol(expr_mat)
  ))

  list(
    expr_mat    = expr_mat,
    platform    = platform,
    voom_weights = NULL   # only normalize_rnaseq_voom() (called explicitly,
                           # not via this auto-dispatch path) produces weights
  )
}


.filter_low_counts <- function(count_mat, min_count = 10, min_prop = 0.20) {
  n_min <- max(2L, round(min_prop * ncol(count_mat)))
  keep  <- rowSums(count_mat >= min_count) >= n_min
  mat   <- count_mat[keep, , drop = FALSE]
  message(sprintf(
    "[NORM] Low-count filter: kept %d / %d genes (>= %d counts in >= %.0f%% of samples).",
    nrow(mat), nrow(count_mat), min_count, 100 * min_prop
  ))
  if (!nrow(mat)) stop("[ERROR] All genes removed -- relax min_count / min_prop.")
  mat
}

.load_count_csv <- function(path) {
  if (!file.exists(path))
    stop(sprintf("[ERROR] Count file not found: '%s'", path))
  mat <- as.matrix(read.csv(path, row.names = 1L, check.names = FALSE))
  storage.mode(mat) <- "integer"
  if (any(mat < 0L, na.rm = TRUE))
    stop("[ERROR] Negative values found. Counts must be non-negative integers.")
  mat
}

.looks_like_ensembl <- function(ids) {
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) return(FALSE)
  mean(grepl("^ENSG[0-9]+", ids)) > 0.5
}

.looks_like_affy_probes <- function(ids) {
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) return(FALSE)
  mean(grepl("_at$", ids)) > 0.5
}

.feature_space <- function(ids) {
  ids <- ids[!is.na(ids) & nzchar(trimws(ids))]
  n <- length(ids)
  if (!n) return("unknown")
  f <- function(x) mean(x)
  if (f(grepl("^ENSG[0-9]+", ids)) > 0.5)                         "ensembl"
  else if (f(grepl("_at$", ids)) > 0.5)                           "affy_probe"
  else if (f(grepl("^ILMN_", ids)) > 0.5)                         "illumina_probe"
  else if (f(grepl("^A_[0-9]+_P[0-9]+$", ids)) > 0.5)             "agilent_probe"
  else if (f(grepl("^[A-Za-z]{1,3}[0-9]{4,}$", ids)) > 0.5)       "probe"
  else if (f(grepl("^[0-9]+$", ids)) > 0.5)                       "numeric"
  else if (f(grepl("^[A-Za-z][A-Za-z0-9_.-]*$", ids)) > 0.8)      "symbol"
  else "unknown"
}

.affy_chip_for_gpl <- function(gpl_id) {
  map <- c(
    "GPL96"   = "hgu133a.db",
    "GPL571"  = "hgu133a.db",
    "GPL570"  = "hgu133plus2.db",
    "GPL201"  = "hgu133b.db",
    "GPL1352" = "hgu133a.db"
  )
  if (is.null(gpl_id) || !nzchar(gpl_id) || !gpl_id %in% names(map)) return(NULL)
  unname(map[gpl_id])
}

load_processed_matrix <- function(data_path, sample_meta = NULL, min_sd = 1e-8,
                                  gpl_id = NULL, ...) {
  if (!file.exists(data_path))
    stop(sprintf("[ERROR] Processed matrix file not found: '%s'", data_path))

  .ensure_pkg("data.table", bioc = FALSE)
  suppressPackageStartupMessages(library(data.table))

  raw <- tryCatch(
    data.table::fread(data_path, header = TRUE, check.names = FALSE,
                      na.strings = c("NA", "NaN", "na", "N/A", "")),
    error = function(e) {
      stop(sprintf("[ERROR] Could not read processed matrix '%s': %s", data_path, conditionMessage(e)))
    }
  )
  raw <- as.data.frame(raw)
  if (ncol(raw) < 2L)
    stop(sprintf("[ERROR] Processed matrix '%s' has fewer than 2 columns.", data_path))

  is_num_col <- vapply(raw, is.numeric, logical(1))
  id_col <- if (!any(is_num_col)) 1L else which(!is_num_col)[1]
  if (is.na(id_col)) id_col <- 1L

  ids <- as.character(raw[[id_col]])
  if (all(is.na(ids)) || any(!nzchar(trimws(ids[!is.na(ids)]))))
    stop(sprintf("[ERROR] Could not identify a feature-ID column in '%s'.", data_path))

  num_cols <- setdiff(seq_len(ncol(raw)), id_col)
  mat <- as.matrix(raw[, num_cols, drop = FALSE])
  for (j in seq_len(ncol(mat))) {
    mat[, j] <- suppressWarnings(as.numeric(as.character(mat[, j])))
  }
  keep_col <- !apply(is.na(mat), 2L, all)
  if (!any(keep_col)) stop(sprintf("[ERROR] No numeric sample columns in '%s'.", data_path))
  if (sum(!keep_col))
    warning(sprintf("[WARN] Dropped %d non-numeric column(s) from processed matrix.", sum(!keep_col)))
  mat <- mat[, keep_col, drop = FALSE]

  rownames(mat) <- ids
  mat <- mat[
    !is.na(rownames(mat)) &
      nzchar(trimws(rownames(mat))) &
      !duplicated(rownames(mat)),
    ,
    drop = FALSE
  ]

  row_sd <- apply(mat, 1L, sd, na.rm = TRUE)
  mat <- mat[!is.na(row_sd) & row_sd > min_sd, , drop = FALSE]
  if (!nrow(mat)) stop("[ERROR] No rows with nonzero variance in processed matrix.")

  if (.looks_like_ensembl(rownames(mat)) && exists("map_ensembl_to_symbol", mode = "function")) {
    mat <- map_ensembl_to_symbol(mat)
  } else if (.looks_like_affy_probes(rownames(mat))) {
    chip <- .affy_chip_for_gpl(gpl_id)
    if (!is.null(chip) && exists("collapse_probes_to_symbols", mode = "function")) {
      message(sprintf("[NORM] Processed-matrix rows are Affy probes; collapsing via %s.", chip))
      mat <- collapse_probes_to_symbols(mat, chip = chip)
    } else if (exists("map_probes_to_symbols", mode = "function") &&
               !is.null(gpl_id) && nzchar(gpl_id) &&
               .feature_space(rownames(mat)) %in% c("affy_probe", "illumina_probe", "agilent_probe", "probe")) {
      mapped <- map_probes_to_symbols(mat, gpl_id)
      if (!is.null(mapped)) mat <- mapped
    } else {
      warning("[WARN] Processed matrix rows look like Affymetrix probe IDs; they may not merge with gene-symbol matrices.")
    }
  } else if (exists("map_probes_to_symbols", mode = "function") &&
             !is.null(gpl_id) && nzchar(gpl_id) &&
             .feature_space(rownames(mat)) %in% c("illumina_probe", "agilent_probe", "probe", "numeric")) {
    message(sprintf("[NORM] Processed-matrix rows are %s IDs; mapping via GPL %s table.",
                    .feature_space(rownames(mat)), gpl_id))
    mapped <- map_probes_to_symbols(mat, gpl_id)
    if (!is.null(mapped)) mat <- mapped
    else warning("[WARN] Could not map processed-matrix probe IDs to gene symbols; keeping as-is.")
  }

  if (!is.null(sample_meta)) {
    pheno_ids <- rownames(sample_meta)
    col_ids   <- colnames(mat)
    mapped    <- intersect(col_ids, pheno_ids)
    if (length(mapped) < 0.5 * length(col_ids)) {
      warning(sprintf(
        "[WARN] Only %d / %d processed-matrix columns matched pData sample IDs; column order retained as-is.",
        length(mapped), length(col_ids)
      ))
    }
    if (length(mapped)) mat <- mat[, mapped, drop = FALSE]
  }

  message(sprintf("[NORM] Processed matrix: %d features x %d samples.", nrow(mat), ncol(mat)))
  mat
}


normalize_microarray_affy <- function(data_path, sample_meta = NULL, gpl_id = NULL, ...) {

  # limma is loaded here so it is available to every downstream step that
  # calls this function (as specified).
  suppressPackageStartupMessages({
    library(affy)
    library(affyio)
    library(limma)
  })

  cel_dir <- if (dir.exists(data_path)) data_path else getwd()
  old_wd  <- setwd(cel_dir)
  on.exit(setwd(old_wd), add = TRUE)

  cel_files <- list.files(".", pattern = "\\.CEL$", full.names = TRUE, ignore.case = TRUE,
                          recursive = TRUE)
  if (!length(cel_files))
    stop(sprintf("[ERROR] No .CEL files found in '%s'.", cel_dir))

  
  valid <- Filter(function(fp) {
    tryCatch({
      affyio::read.celfile.header(fp)
      TRUE
    }, error = function(e) {
      warning(sprintf("[WARN] Quarantining corrupt CEL file: %s (%s)",
                       basename(fp), conditionMessage(e)))
      FALSE
    })
  }, cel_files)

  message(sprintf("[NORM] Valid CEL files: %d / %d.", length(valid), length(cel_files)))
  if (!length(valid)) stop("[ERROR] No valid CEL files remain after validation.")

  raw  <- ReadAffy(filenames = valid)
  eset <- rma(raw)

  clean_names <- sub("_.*$", "", basename(sampleNames(eset)))
  sampleNames(eset) <- clean_names

  mat <- exprs(eset)   # probes x samples, log2 scale

  chip <- .affy_chip_for_gpl(gpl_id)
  if (!is.null(chip) && exists("collapse_probes_to_symbols", mode = "function")) {
    mat <- collapse_probes_to_symbols(mat, chip = chip)
  } else {
    warning(sprintf(
      "[WARN] No annotation package mapped for GPL '%s'; keeping probe IDs (may not merge with gene-symbol matrices).",
      gpl_id
    ))
  }

  mat
}


normalize_rnaseq_vst <- function(data_path,
                                  sample_meta = NULL,
                                  min_count   = 10,
                                  min_prop    = 0.20,
                                  blind       = TRUE,
                                  ...) {
  suppressPackageStartupMessages({
    library(DESeq2)
    library(edgeR)
  })

  count_mat <- .load_count_csv(data_path)
  count_mat <- .filter_low_counts(count_mat, min_count, min_prop)

  col_data <- data.frame(
    condition = factor(rep("unknown", ncol(count_mat))),
    row.names = colnames(count_mat)
  )

  
  dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = col_data, design = ~1)
  dds <- estimateSizeFactors(dds)

  vst_obj <- tryCatch(
    vst(dds, blind = blind),
    error = function(e) {
      message("[NORM] vst() failed; falling back to varianceStabilizingTransformation().")
      varianceStabilizingTransformation(dds, blind = blind)
    }
  )

  assay(vst_obj)   # approximately log2 scale
}


normalize_rnaseq_voom <- function(data_path,
                                   sample_meta = NULL,
                                   min_count   = 10,
                                   min_prop    = 0.20,
                                   ...) {
  suppressPackageStartupMessages({
    library(limma)
    library(edgeR)
  })

  count_mat <- .load_count_csv(data_path)
  count_mat <- .filter_low_counts(count_mat, min_count, min_prop)

  dge <- DGEList(counts = count_mat)
  dge <- calcNormFactors(dge, method = "TMM")

  message(sprintf(
    "[NORM] TMM factors: min=%.3f, max=%.3f.",
    min(dge$samples$norm.factors), max(dge$samples$norm.factors)
  ))

  dir.create("outputs/qc", showWarnings = FALSE, recursive = TRUE)
  png("outputs/qc/voom_trend.png", width = 800L, height = 600L)
  v <- voom(dge, plot = TRUE)
  dev.off()
  message("[NORM] voom mean-variance trend plot saved to outputs/qc/voom_trend.png.")

  list(
    log2_cpm = v$E,
    weights  = v$weights
  )
}