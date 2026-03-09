#' Convert a SummarizedExperiment to a Python AnnData object
#'
#' @description
#' Creates a Python `anndata.AnnData` object from a
#' [SummarizedExperiment::SummarizedExperiment] (or any subclass, including
#' [SingleCellExperiment::SingleCellExperiment]).  The main matrix (`X`) is
#' built by [se_assay_to_python_matrix()], which dispatches to the
#' appropriate zero-copy NumPy / SciPy converter for the assay's matrix class.
#'
#' ## Copy semantics
#' * **X** (the main count matrix): converted via [se_assay_to_python_matrix()]
#'   and then transposed (`$T`) for AnnData orientation.  See that function
#'   for per-class zero-copy details.
#' * **layers**: each element is converted via
#'   [.matrix_to_python_array_or_sparse()] if it is an R object, or
#'   forwarded as-is if it is already a Python object.
#' * **obsm / varm**: converted via the same dispatcher with
#'   `prefer_sparse = FALSE` (dense NumPy preferred for embedding matrices);
#'   already-Python objects are forwarded unchanged.
#' * **obs / var metadata**: one unavoidable copy per `data.frame` column
#'   (R columnar → Python row-oriented).
#'
#' @param x A [SummarizedExperiment::SummarizedExperiment] or subclass.
#' @param assay A single string: the assay to place in `X` (default
#'   `"counts"`).
#' @param obs `NULL` or a `data.frame` of cell-level metadata.  Defaults to
#'   `as.data.frame(colData(x))` with `rownames` set to `colnames(x)`.
#' @param var `NULL` or a `data.frame` of feature-level metadata.  Defaults
#'   to `as.data.frame(rowData(x))` with `rownames` set to `rownames(x)`.
#' @param layers A named list of additional assay matrices for
#'   `AnnData$layers`.  Each element may be an R matrix-like object or an
#'   already-constructed Python object.
#' @param obsm A named list of matrices for `AnnData$obsm`.
#' @param varm A named list of matrices for `AnnData$varm`.
#'
#' @return A Python `anndata.AnnData` object (reticulate, `convert = FALSE`).
#'
#' @seealso [se_assay_to_python_matrix()], [sce_assay_to_scipy_csc()],
#'   [sce_to_reticulate_anndata()]
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
  .validate_summarized_experiment(x)

  # ---- build main matrix (zero-copy where possible) ------------------------
  # se_assay_to_python_matrix returns shape (n_features, n_samples).
  # AnnData expects X with shape (obs=cells, vars=genes), so we transpose.
  # For sparse matrices, .T yields a view sharing all data buffers (no copy).
  X <- se_assay_to_python_matrix(x, assay = assay, prefer_sparse = TRUE)
  X <- X$T  # (n_cells, n_genes), zero-copy for sparse; view for dense

  # ---- build obs (cell metadata) -------------------------------------------
  if (is.null(obs)) {
    obs <- as.data.frame(SummarizedExperiment::colData(x))
    rownames(obs) <- colnames(x)      # ensure cell barcodes are row names
  } else {
    stopifnot(is.data.frame(obs), nrow(obs) == ncol(x))
  }
  obs <- .coerce_df_to_character_rownames(obs, default_prefix = "cell")

  # ---- build var (feature metadata) ----------------------------------------
  if (is.null(var)) {
    var <- as.data.frame(SummarizedExperiment::rowData(x))
    rownames(var) <- rownames(x)      # ensure feature names are row names
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

  # ---- optional layers -----------------------------------------------------
  # Each element is dispatched through the matrix converter if it is an R
  # object, or forwarded as-is if it is already a Python object.
  # R assay matrices have shape (features, cells); AnnData layers must match
  # X which is (obs=cells, vars=features), so R objects are transposed.
  if (length(layers) > 0L) {
    if (is.null(names(layers)) || any(names(layers) == "")) {
      stop("'layers' must be a *named* list.")
    }
    for (nm in names(layers)) {
      mat_py <- if (.is_python_object(layers[[nm]])) {
        layers[[nm]]   # already Python: caller controls orientation
      } else {
        # Use __setitem__ to mutate the Python Layers mapping in-place.
        # R's `[[<-` on a nested Python attribute does not reliably propagate
        # back through reticulate's attribute proxy chain.
        .matrix_to_python_array_or_sparse(layers[[nm]], prefer_sparse = TRUE)$T
      }
      py_adata$layers$`__setitem__`(nm, mat_py)
    }
  }

  # ---- optional obsm -------------------------------------------------------
  # Embeddings / coordinates are typically dense; prefer_sparse = FALSE.
  if (length(obsm) > 0L) {
    if (is.null(names(obsm)) || any(names(obsm) == "")) {
      stop("'obsm' must be a *named* list.")
    }
    for (nm in names(obsm)) {
      mat_py <- if (.is_python_object(obsm[[nm]])) {
        obsm[[nm]]
      } else {
        .matrix_to_python_array_or_sparse(obsm[[nm]], prefer_sparse = FALSE)
      }
      py_adata$obsm$`__setitem__`(nm, mat_py)
    }
  }

  # ---- optional varm -------------------------------------------------------
  if (length(varm) > 0L) {
    if (is.null(names(varm)) || any(names(varm) == "")) {
      stop("'varm' must be a *named* list.")
    }
    for (nm in names(varm)) {
      mat_py <- if (.is_python_object(varm[[nm]])) {
        varm[[nm]]
      } else {
        .matrix_to_python_array_or_sparse(varm[[nm]], prefer_sparse = FALSE)
      }
      py_adata$varm$`__setitem__`(nm, mat_py)
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
