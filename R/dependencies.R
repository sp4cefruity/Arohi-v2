REQUIRED_PACKAGES <- c(
  "here",
  "BiocManager",
  "Biobase",
  "limma",
  "edgeR",
  "DESeq2",
  "GEOquery",
  "affy",
  "affyio",
  "fgsea",
  "msigdbr",
  "GSVA",
  "pheatmap",
  "preprocessCore",
  "AnnotationDbi",
  "biomaRt",
  "dorothea",
  "viper",
  "ggplot2",
  "R.utils"
)

check_project_dependencies <- function(required = REQUIRED_PACKAGES) {
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0L) {
    stop(
      sprintf(
        "[DEPENDENCIES] Missing required R package(s): %s.\nInstall them in a terminal with:\n  source('R/dependencies.R'); install_project_dependencies()",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

install_project_dependencies <- function(required = REQUIRED_PACKAGES,
                                        repos = "https://cloud.r-project.org") {
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

  if (!length(missing)) {
    message("[DEPENDENCIES] All required R packages are already installed.")
    return(invisible(TRUE))
  }

  message(sprintf("[DEPENDENCIES] Installing missing packages: %s", paste(missing, collapse = ", ")))

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = repos)
  }

  options(repos = c(CRAN = repos))
  options(BioC_mirror = "https://bioconductor.org")

  BiocManager::install(missing, ask = FALSE, update = FALSE, INSTALL_opts = c("--no-test-load"))
  invisible(TRUE)
}
