#' Convert a SingleCellExperiment to a Python AnnData object
#'
#' @description
#' Creates a Python `anndata.AnnData` object from a
#' [SingleCellExperiment::SingleCellExperiment] using a
#' `scipy.sparse.csc_matrix` built with [sce_assay_to_scipy_csc()].
#'
#' ## Copy semantics
#' * **X** (the main count matrix): the dgCMatrix is first converted to a
#'   `scipy.sparse.csc_matrix` (shape `genes × cells`) via
#'   [sce_assay_to_scipy_csc()], then transposed to
#'   `scipy.sparse.csr_matrix` (shape `cells × genes`) by `.T`.  The
#'   transpose shares all data buffers with the original CSC
#'   (`flags.owndata = FALSE`), so **no values are copied**.
#' * **obs / var metadata**: a pandas `DataFrame` is constructed from the R
#'   `data.frame`; reticulate performs one copy here, which is unavoidable for
#'   columnar R data frames → row-oriented Python objects.
#' * **layers / obsm / varm**: passed through as-is; the caller is responsible
#'   for any copy semantics of those objects.
#'
#' @param x A [SingleCellExperiment::SingleCellExperiment].
#' @param assay A single string: the assay to place in `X` (default
#'   `"counts"`).
#' @param obs `NULL` or a `data.frame` of cell-level metadata.  Defaults to
#'   `as.data.frame(colData(x))` with `rownames` set to `colnames(x)`.
#' @param var `NULL` or a `data.frame` of gene-level metadata.  Defaults to
#'   `as.data.frame(rowData(x))` with `rownames` set to `rownames(x)`.
#' @param layers A named list of additional assay matrices to place in
#'   `AnnData$layers`.
#' @param obsm A named list of matrices added to `AnnData$obsm`.
#' @param varm A named list of matrices added to `AnnData$varm`.
#'
#' @return A Python `anndata.AnnData` object (reticulate, `convert = FALSE`).
#'
#' @seealso [sce_assay_to_scipy_csc()], [sce_to_reticulate_anndata()]
#'
#' @examples
#' \dontrun{
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = Matrix::rsparsematrix(200, 50, 0.1, repr = "C"))
#' )
#' adata <- sce_to_anndata(sce)
#' cat("obs x var:", reticulate::py_str(adata$shape), "\n")
#' }
#'
#' @importFrom SummarizedExperiment colData rowData assayNames
#' @export
sce_to_anndata <- function(
    x,
    assay  = "counts",
    obs    = NULL,
    var    = NULL,
    layers = list(),
    obsm   = list(),
    varm   = list()) {

  # ---- validate x ----------------------------------------------------------
  if (!is(x, "SingleCellExperiment")) {
    stop(
      "'x' must be a SingleCellExperiment object, not '",
      paste(class(x), collapse = "', '"), "'."
    )
  }

  # ---- build main sparse matrix (zero-copy where possible) -----------------
  # sce_assay_to_scipy_csc returns csc_matrix with shape (genes, cells).
  # AnnData expects X with shape (obs=cells, vars=genes), so we transpose.
  # scipy CSC.T yields a CSR matrix sharing the same data buffers (no copy).
  csc <- sce_assay_to_scipy_csc(x, assay = assay)
  X   <- csc$T  # csr_matrix shape (n_cells, n_genes), zero-copy

  # ---- build obs (cell metadata) -------------------------------------------
  if (is.null(obs)) {
    obs <- as.data.frame(SummarizedExperiment::colData(x))
    rownames(obs) <- colnames(x)      # ensure cell barcodes are row names
  } else {
    stopifnot(is.data.frame(obs), nrow(obs) == ncol(x))
  }
  obs <- .coerce_df_to_character_rownames(obs, default_prefix = "cell")

  # ---- build var (gene metadata) -------------------------------------------
  if (is.null(var)) {
    var <- as.data.frame(SummarizedExperiment::rowData(x))
    rownames(var) <- rownames(x)      # ensure gene names are row names
  } else {
    stopifnot(is.data.frame(var), nrow(var) == nrow(x))
  }
  var <- .coerce_df_to_character_rownames(var, default_prefix = "var")

  # ---- Python imports (convert = FALSE) ------------------------------------
  pd <- reticulate::import("pandas", convert = FALSE)
  ad <- reticulate::import("anndata", convert = FALSE)

  # Build pandas DataFrames; index must be the obs/var names.
  py_obs <- .df_to_pandas(obs, pd)
  py_var <- .df_to_pandas(var, pd)

  # ---- construct AnnData ---------------------------------------------------
  py_adata <- ad$AnnData(
    X   = X,
    obs = py_obs,
    var = py_var
  )

  # ---- optional layers / obsm / varm ---------------------------------------
  if (length(layers) > 0L) {
    if (is.null(names(layers)) || any(names(layers) == "")) {
      stop("'layers' must be a *named* list.")
    }
    for (nm in names(layers)) {
      # Use __setitem__ to mutate the Python Layers mapping in-place.
      # R's `[[<-` on a nested Python attribute does not reliably propagate
      # back through reticulate's attribute proxy chain.
      py_adata$layers$`__setitem__`(nm, reticulate::r_to_py(layers[[nm]]))
    }
  }

  if (length(obsm) > 0L) {
    if (is.null(names(obsm)) || any(names(obsm) == "")) {
      stop("'obsm' must be a *named* list.")
    }
    for (nm in names(obsm)) {
      py_adata$obsm$`__setitem__`(nm, reticulate::r_to_py(obsm[[nm]]))
    }
  }

  if (length(varm) > 0L) {
    if (is.null(names(varm)) || any(names(varm) == "")) {
      stop("'varm' must be a *named* list.")
    }
    for (nm in names(varm)) {
      py_adata$varm$`__setitem__`(nm, reticulate::r_to_py(varm[[nm]]))
    }
  }

  py_adata
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Ensure a data.frame has character row names
#'
#' @param df A data.frame.
#' @param default_prefix Prefix used when row names are missing.
#' @return The same data.frame with valid character row names.
#' @keywords internal
.coerce_df_to_character_rownames <- function(df, default_prefix = "obs") {
  rn <- rownames(df)
  if (is.null(rn) || all(rn == "") || all(rn == as.character(seq_len(nrow(df))))) {
    rownames(df) <- paste0(default_prefix, seq_len(nrow(df)))
  }
  df
}

#' Convert an R data.frame to a pandas DataFrame with index = rownames
#'
#' @param df An R `data.frame` with row names set.
#' @param pd A reticulate handle to the `pandas` module.
#' @return A Python pandas `DataFrame`.
#' @keywords internal
.df_to_pandas <- function(df, pd) {
  idx  <- rownames(df)
  # r_to_py converts the data.frame; then we set the index explicitly.
  py_df <- pd$DataFrame(reticulate::r_to_py(df))
  py_df$index <- reticulate::r_to_py(idx)
  py_df
}
