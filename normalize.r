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

suppressPackageStartupMessages({
  library("Biobase")
  library("limma")
  library("edgeR")
  library("DESeq2")
  library("GEOquery")
})

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
  )

detect_platform <- function(gse_id) {
  suppressPackageStartupMessages(library(GEOquery))

  known_affy <- c("GPL96", "GPL570", "GPL571", "GPL201")

  message(sprintf("[PLATFORM] Fetching GEO metadata for %s...", gse_id))
  gse_obj <- tryCatch(
    GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE),
    error = function(e) {
      warning(sprintf(
        "[WARN] Could not fetch GEO metadata for %s (%s); defaulting to affy_microarray.",
        gse_id, conditionMessage(e)
      ))
      NULL
    }
  )
  if (is.null(gse_obj)) return("affy_microarray")

  eset   <- if (is.list(gse_obj)) gse_obj[[1]] else gse_obj
  gpl_id <- tryCatch(Biobase::annotation(eset), error = function(e) "")

  if (!nzchar(gpl_id)) {
    warning(sprintf(
      "[WARN] No GPL annotation found for %s; defaulting to affy_microarray.", gse_id
    ))
    return("affy_microarray")
  }

  message(sprintf("[PLATFORM] %s -> %s", gse_id, gpl_id))

  if (gpl_id %in% known_affy) {
    return("affy_microarray")
  }

  
  gpl_obj <- tryCatch(GEOquery::getGEO(gpl_id), error = function(e) NULL)
  desc <- ""
  if (!is.null(gpl_obj)) {
    meta <- tryCatch(Biobase::Meta(gpl_obj), error = function(e) list())
    desc <- paste(
      meta$title %||% "",
      meta$description %||% "",
      meta$technology %||% "",
      collapse = " "
    )
  }
  desc_lower <- tolower(desc)

  if (grepl("rna-seq|rna seq|counts", desc_lower)) {
    return("rnaseq_counts")
  }

  warning(sprintf(
    "[WARN] Could not confidently classify platform for %s (GPL %s, description: '%s'); defaulting to affy_microarray.",
    gse_id, gpl_id, desc
  ))
  "affy_microarray"
}

`%||%` <- function(a, b) if (is.null(a) || !nzchar(a)) b else a


run_normalization <- function(data_path, gse_id, sample_meta = NULL, ...) {

  platform <- detect_platform(gse_id)

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
  expr_mat <- norm_fn(data_path = data_path, sample_meta = sample_meta, ...)

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


normalize_microarray_affy <- function(data_path, sample_meta = NULL, ...) {

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

  cel_files <- list.files(".", pattern = "\\.CEL$", full.names = TRUE, ignore.case = TRUE)
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

  exprs(eset)   # genes x samples, log2 scale
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