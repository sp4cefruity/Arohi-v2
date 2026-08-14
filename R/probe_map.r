
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
 
.collapse_gene_mat <- function(expr_sub, symbols, method = "maxIQR") {
  symbols <- trimws(symbols)
  keep <- !is.na(symbols) & nzchar(symbols)
  expr_sub <- expr_sub[keep, , drop = FALSE]
  symbols  <- symbols[keep]

  if (method == "maxIQR") {
    probe_iqr   <- apply(expr_sub, 1L, IQR, na.rm = TRUE)
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
    stop(sprintf("[ERROR] Unknown collapse method: '%s'. Use 'maxIQR' or 'mean'.", method))
  }
  gene_mat <- gene_mat[order(rownames(gene_mat)), , drop = FALSE]
  gene_mat
}

collapse_probes_to_symbols <- function(expr_mat, chip = "hgu133plus2.db",
                                       method = "maxIQR") {
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

  keep <- !is.na(symbol_map) & nchar(trimws(symbol_map)) > 0
  n_total <- length(symbol_map)
  n_annotated <- sum(keep)
  message(sprintf("[MAP] Annotated: %d / %d probes (%.1f%%).",
                  n_annotated, n_total, 100 * n_annotated / n_total))

  expr_sub <- expr_mat[keep, , drop = FALSE]
  symbols  <- symbol_map[keep]

  gene_mat <- .collapse_gene_mat(expr_sub, symbols, method)
  message(sprintf("[MAP] Collapsed to %d unique Gene Symbols.", nrow(gene_mat)))
  gene_mat
}

.geo_cache_dir_local <- function() {
  if (exists(".geo_cache_dir", mode = "function", inherits = TRUE)) {
    .geo_cache_dir()
  } else {
    dir <- file.path("outputs", "cache", "geo")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    dir
  }
}

.robust_fetch_local <- function(fun) {
  if (exists(".robust_fetch", mode = "function", inherits = TRUE)) .robust_fetch(fun) else fun()
}

.find_symbol_col <- function(tab) {
  nm <- names(tab)
  score <- vapply(nm, function(cn) {
    lc <- tolower(cn)
    if (grepl("gene.?symbol|hgnc", lc))      100
    else if (grepl("^symbol$", lc))          90
    else if (grepl("symbol", lc))            80
    else if (grepl("gene.?name", lc))        70
    else if (grepl("gene", lc) && grepl("name", lc)) 60
    else 0
  }, numeric(1))
  ord <- order(score, decreasing = TRUE, method = "radix")
  for (cn in nm[ord]) {
    if (score[[cn]] == 0) break
    vals <- as.character(tab[[cn]])
    vals <- vals[!is.na(vals) & nzchar(trimws(vals))]
    if (length(vals) < 0.1 * nrow(tab)) next
    if (mean(grepl("^[A-Z0-9][A-Za-z0-9_.-]*$", vals)) > 0.5) return(cn)
  }
  NULL
}

map_entrez_to_symbols <- function(expr_mat, method = "maxIQR") {
  .ensure_pkg(c("AnnotationDbi", "org.Hs.eg.db"), bioc = TRUE)
  suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Hs.eg.db)
  })

  rn <- rownames(expr_mat)
  ids <- sub("^0+", "", trimws(rn))
  ids <- ids[ids != ""]

  symbol_map <- tryCatch(
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL",
                          keytype = "ENTREZID", multiVals = "first"),
    error = function(e) {
      warning(sprintf("[MAP] Entrez->symbol lookup failed: %s", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(symbol_map)) return(NULL)

  valid <- !is.na(symbol_map) & nzchar(trimws(symbol_map))
  n_mapped <- sum(valid)
  message(sprintf("[MAP] org.Hs.eg.db: mapped %d / %d Entrez IDs to symbols.",
                  n_mapped, length(ids)))
  if (n_mapped < 0.5 * length(ids)) {
    warning(sprintf(
      "[MAP] Too few Entrez IDs mapped (%.0f%%); aborting mapping.",
      100 * n_mapped / length(ids)
    ))
    return(NULL)
  }

  rownames(expr_mat) <- ids
  expr_sub <- expr_mat[names(symbol_map[valid]), , drop = FALSE]
  symbols  <- unname(symbol_map[valid])
  gene_mat <- .collapse_gene_mat(expr_sub, symbols, method)
  message(sprintf("[MAP] Collapsed to %d unique Gene Symbols via org.Hs.eg.db.", nrow(gene_mat)))
  gene_mat
}


map_probes_to_symbols <- function(expr_mat, gpl_id, method = "maxIQR") {
  if (.feature_space(rownames(expr_mat)) == "numeric") {
    message(sprintf("[MAP] Rows are numeric Entrez IDs (%d); mapping via org.Hs.eg.db.",
                    nrow(expr_mat)))
    return(map_entrez_to_symbols(expr_mat, method))
  }
  if (is.null(gpl_id) || !nzchar(gpl_id)) {
    message("[MAP] No GPL id given; cannot map probe IDs to symbols.")
    return(NULL)
  }
  gpl_obj <- tryCatch(
    .robust_fetch_local(function() GEOquery::getGEO(gpl_id, getGPL = TRUE,
                                                    destdir = .geo_cache_dir_local())),
    error = function(e) {
      warning(sprintf("[MAP] GPL %s lookup failed: %s", gpl_id, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(gpl_obj)) return(NULL)
  tab <- tryCatch(GEOquery::Table(gpl_obj), error = function(e) NULL)
  if (is.null(tab) || !nrow(tab)) {
    warning(sprintf("[MAP] GPL %s has no annotation table.", gpl_id))
    return(NULL)
  }

  rn <- rownames(expr_mat)
  id_hits <- vapply(names(tab), function(cn) {
    vals <- as.character(tab[[cn]])
    sum(!is.na(vals) & vals %in% rn)
  }, numeric(1))
  id_col <- which.max(id_hits)
  if (id_hits[id_col] < 0.5 * length(rn)) {
    warning(sprintf(
      "[MAP] GPL %s: feature IDs overlap only %.0f%% of rows; not mapping.",
      gpl_id, 100 * id_hits[id_col] / length(rn)
    ))
    return(NULL)
  }

  ids  <- as.character(tab[[id_col]])
  sym_col <- .find_symbol_col(tab)
  if (is.null(sym_col)) {
    warning(sprintf("[MAP] GPL %s: no gene-symbol column found.", gpl_id))
    return(NULL)
  }
  syms <- as.character(tab[[sym_col]])

  valid <- !is.na(ids) & !is.na(syms) & nzchar(trimws(syms))
  id_map <- setNames(trimws(syms[valid]), ids[valid])
  id_map <- id_map[!duplicated(ids[valid])]

  probe_ids <- rn[rn %in% names(id_map)]
  message(sprintf("[MAP] GPL %s: mapped %d / %d probe IDs to gene symbols.",
                  gpl_id, length(probe_ids), length(rn)))
  if (length(probe_ids) < 0.5 * length(rn)) {
    warning(sprintf(
      "[MAP] Too few probes mapped (%.0f%%); aborting mapping.",
      100 * length(probe_ids) / length(rn)
    ))
    return(NULL)
  }

  expr_sub <- expr_mat[probe_ids, , drop = FALSE]
  symbols  <- unname(id_map[probe_ids])
  gene_mat <- .collapse_gene_mat(expr_sub, symbols, method)
  message(sprintf("[MAP] Collapsed to %d unique Gene Symbols via GPL %s.", nrow(gene_mat), gpl_id))
  gene_mat
}