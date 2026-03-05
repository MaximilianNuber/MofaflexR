#' Convert a SingleCellExperiment assay to a scipy CSC sparse matrix
#'
#' @description
#' Extracts an assay from a [SingleCellExperiment::SingleCellExperiment] object,
#' coerces it to a [Matrix::dgCMatrix-class], and constructs a Python
#' `scipy.sparse.csc_matrix` by passing [reticulate::np_array()] views of the
#' underlying R memory (true zero-copy: `flags.owndata = FALSE`).
#'
#' ## Slot mapping
#' | R (`dgCMatrix`) | Python / scipy CSC | dtype |
#' |---|---|---|
#' | `x@x` | `data` | `float64` |
#' | `x@i` | `indices` | `int32` |
#' | `x@p` | `indptr` | `int32` |
#' | `x@Dim` | `shape` | – |
#'
#' ## Zero-copy semantics
#' [reticulate::np_array()] uses the C buffer protocol to expose R's
#' contiguous memory to NumPy **without copying** (`flags.owndata = FALSE`).
#' This works because plain R vectors are already C-contiguous:
#' * R numeric (`double`) maps directly to NumPy `float64`.
#' * R integer (`int32_t`) maps directly to NumPy `int32`.
#'
#' **Note on `np.array(..., copy=False)`**: NumPy >= 2.0 changed this flag to
#' raise `ValueError` when a copy would be needed.  Using
#' [reticulate::np_array()] avoids this entirely and is the canonical
#' reticulate API for zero-copy array creation from R.
#'
#' **Lifetime caveat**: the returned Python object holds a borrowed reference
#' to the underlying R vectors.  Do *not* modify, reallocate, or allow the
#' `mat@x`, `mat@i`, `mat@p` vectors to be garbage-collected while the Python
#' object is alive.
#'
#' ## Why CSC?
#' `dgCMatrix` stores columns contiguously (Compressed Sparse Column).
#' We map directly to `csc_matrix` without reordering values.  Building a CSR
#' matrix would require transposing the data (= always a copy).
#'
#' ## int32 for indptr
#' R integer vectors are 32-bit (`int32_t`).  scipy CSC accepts `int32` for
#' both `indices` and `indptr`, supporting up to \eqn{2^{31}-1} non-zero
#' entries per matrix.
#'
#' @param x A [SingleCellExperiment::SingleCellExperiment].
#' @param assay A single string: the name of the assay to extract
#'   (default `"counts"`).
#'
#' @return A Python `scipy.sparse.csc_matrix` (a reticulate Python object;
#'   not auto-converted to R).
#'
#' @seealso [sce_to_anndata()], [sce_to_reticulate_anndata()]
#'
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' library(Matrix)
#' counts <- rsparsematrix(200, 50, density = 0.1, repr = "C")
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#' csc <- sce_assay_to_scipy_csc(sce, assay = "counts")
#' sp <- reticulate::import("scipy.sparse", convert = FALSE)
#' shape_r <- reticulate::py_to_r(csc$shape)
#' cat("Shape:", shape_r[[1]], "x", shape_r[[2]], "\n")
#' cat("Sparse:", reticulate::py_to_r(sp$issparse(csc)), "\n")
#' cat("owndata:", reticulate::py_to_r(csc$data$flags$owndata), "\n")
#' }
#'
#' @importFrom SummarizedExperiment assay assayNames
#' @export
sce_assay_to_scipy_csc <- function(x, assay = "counts") {
  # ---- input validation ----------------------------------------------------
  if (!is(x, "SingleCellExperiment")) {
    stop(
      "'x' must be a SingleCellExperiment object, not '",
      paste(class(x), collapse = "', '"), "'."
    )
  }
  if (!is.character(assay) || length(assay) != 1L) {
    stop("'assay' must be a single string.")
  }
  avail <- SummarizedExperiment::assayNames(x)
  if (!assay %in% avail) {
    stop(
      "Assay '", assay, "' not found. Available assays: ",
      paste0("'", avail, "'", collapse = ", "), "."
    )
  }

  mat <- SummarizedExperiment::assay(x, assay)

  # Coerce to dgCMatrix if needed (e.g. lgCMatrix, ngCMatrix, dgTMatrix …).
  # Route through "dMatrix" for robustness with recent Matrix deprecations.
  if (!is(mat, "dgCMatrix")) {
    message(
      "Assay '", assay, "' is a '", class(mat),
      "'; coercing to dgCMatrix."
    )
    mat <- tryCatch(
      methods::as(mat, "dgCMatrix"),
      error = function(e) methods::as(methods::as(mat, "dMatrix"), "dgCMatrix")
    )
  }

  # ---- import scipy (convert = FALSE: keep all Python objects as-is) -------
  sp <- reticulate::import("scipy.sparse", convert = FALSE)

  # ---- extract dgCMatrix slots as plain R vectors --------------------------
  # @x  = stored values     (R numeric = C double   → float64)
  # @i  = row indices       (R integer = C int32_t  → int32, already 0-based)
  # @p  = column pointers   (R integer = C int32_t  → int32)
  # @Dim = c(nrow, ncol)
  data_r    <- mat@x
  indices_r <- mat@i
  indptr_r  <- mat@p
  nr        <- mat@Dim[[1L]]
  nc        <- mat@Dim[[2L]]

  # ---- build zero-copy NumPy views into R memory ---------------------------
  # reticulate::np_array() uses the C buffer protocol: R's contiguous memory
  # is exposed to NumPy directly; flags.owndata = FALSE, no heap allocation.
  # Contrast with np$array(x, copy=FALSE) which in NumPy ≥ 2.0 raises
  # ValueError whenever a copy would be needed.
  py_data    <- reticulate::np_array(data_r,    dtype = "float64")
  py_indices <- reticulate::np_array(indices_r, dtype = "int32")
  py_indptr  <- reticulate::np_array(indptr_r,  dtype = "int32")

  shape <- reticulate::tuple(
    reticulate::r_to_py(nr),
    reticulate::r_to_py(nc)
  )

  # ---- construct scipy CSC matrix ------------------------------------------
  # scipy.sparse.csc_matrix((data, indices, indptr), shape=shape)
  csc <- sp$csc_matrix(
    reticulate::tuple(py_data, py_indices, py_indptr),
    shape = shape
  )

  csc
}
