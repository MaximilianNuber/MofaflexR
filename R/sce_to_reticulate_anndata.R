#' Convert a SingleCellExperiment to an anndataR ReticulateAnnData
#'
#' @description
#' Converts a [SingleCellExperiment::SingleCellExperiment] to an
#' `anndataR::ReticulateAnnData` object backed by a Python
#' `anndata.AnnData`.  The main matrix (`X`) is a
#' `scipy.sparse.csc_matrix` built without densifying the data; see
#' [sce_assay_to_scipy_csc()] for zero-copy semantics and lifetime caveats.
#'
#' The returned object behaves as an `anndataR` AnnData: it can be passed to
#' scanpy functions, written to HDF5/`.h5ad`, or inspected with the standard
#' `$` accessor.
#'
#' @param x A [SingleCellExperiment::SingleCellExperiment].
#' @param assay A single string: the assay to place in `X` (default
#'   `"counts"`).
#' @param ... Additional arguments forwarded to [sce_to_anndata()].
#'
#' @return An `anndataR::ReticulateAnnData` R6 object backed by a Python
#'   `anndata.AnnData`.
#'
#' @seealso [sce_assay_to_scipy_csc()], [sce_to_anndata()]
#'
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' library(Matrix)
#'
#' counts <- rsparsematrix(200, 50, density = 0.1, repr = "C")
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#'
#' rada <- sce_to_reticulate_anndata(sce)
#' print(rada)      # shows n_obs x n_vars
#'
#' # Access the underlying Python object
#' py_adata <- rada$py_anndata()
#' cat("Shape:", reticulate::py_str(py_adata$shape), "\n")
#' cat("X is sparse:", reticulate::py_bool(
#'   reticulate::import("scipy.sparse")$issparse(py_adata$X)
#' ), "\n")
#' }
#'
#' @export
sce_to_reticulate_anndata <- function(x, assay = "counts", ...) {
  # ---- validate x ----------------------------------------------------------
  if (!is(x, "SingleCellExperiment")) {
    stop(
      "'x' must be a SingleCellExperiment object, not '",
      paste(class(x), collapse = "', '"), "'."
    )
  }

  # ---- build Python AnnData ------------------------------------------------
  py_adata <- sce_to_anndata(x, assay = assay, ...)

  # ---- wrap in ReticulateAnnData -------------------------------------------
  # anndataR:::ReticulateAnnData is an R6 class; its $initialize() accepts
  # a 'py_anndata' argument that is a Python anndata.AnnData object.
  rada <- get_ReticulateAnnData_constructor()
  rada$new(py_anndata = py_adata)
}

# ---------------------------------------------------------------------------
# Helper: access ReticulateAnnData constructor regardless of export status
# ---------------------------------------------------------------------------

#' @keywords internal
get_ReticulateAnnData_constructor <- function() {
  # ReticulateAnnData is an R6 generator defined in anndataR.
  # It is accessible via the namespace regardless of export status.
  get("ReticulateAnnData", envir = asNamespace("anndataR"), inherits = FALSE)
}
