# Helper utilities shared by all test files
# Loaded automatically by testthat before test files.

#' Build a tiny SingleCellExperiment with known dgCMatrix counts
#' @param nrow Number of genes (rows)
#' @param ncol Number of cells (columns)
#' @param density Fraction of non-zero entries
make_test_sce <- function(nrow = 30L, ncol = 20L, density = 0.3) {
  set.seed(42L)
  counts <- Matrix::rsparsematrix(nrow, ncol, density = density, repr = "C")
  # Ensure values are non-negative integers (typical for counts)
  counts <- abs(round(counts))
  counts <- methods::as(counts, "dgCMatrix")

  obs_df <- data.frame(
    cell_type = sample(c("A", "B"), ncol, replace = TRUE),
    row.names  = paste0("cell", seq_len(ncol))
  )
  var_df <- data.frame(
    gene_group = sample(c("protein_coding", "lncRNA"), nrow, replace = TRUE),
    row.names   = paste0("gene", seq_len(nrow))
  )

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays   = list(counts = counts),
    colData  = S4Vectors::DataFrame(obs_df),
    rowData  = S4Vectors::DataFrame(var_df)
  )
  sce
}

#' Skip if reticulate + scipy are not available without basilisk
skip_if_no_scipy <- function() {
  testthat::skip_if_not(
    reticulate::py_module_available("scipy") &&
      reticulate::py_module_available("numpy") &&
      reticulate::py_module_available("anndata"),
    message = "scipy/numpy/anndata Python modules not available; skipping."
  )
}

#' Skip if reticulate + mudata are not available
skip_if_no_mudata <- function() {
  testthat::skip_if_not(
    reticulate::py_module_available("mudata") &&
      reticulate::py_module_available("anndata") &&
      reticulate::py_module_available("scipy"),
    message = "mudata/anndata/scipy Python modules not available; skipping."
  )
}
