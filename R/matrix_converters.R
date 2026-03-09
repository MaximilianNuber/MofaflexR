# ===========================================================================
# R → Python NumPy / SciPy matrix conversion utilities
# ===========================================================================
#
# Internal helpers to convert common R matrix representations to Python
# objects via reticulate (convert = FALSE throughout).
#
# Supported input types:
#   * base R matrix          → numpy.ndarray (zero-copy via r_to_py)
#   * Matrix::dgeMatrix      → numpy.ndarray (zero-copy via @x + reshape)
#   * Matrix::dgCMatrix      → scipy.sparse.csc_matrix (zero-copy slot views)
#   * Matrix::dgRMatrix      → scipy.sparse.csr_matrix (zero-copy slot views)
#   * SparseArray::COO_SparseArray (2-D) → scipy.sparse.coo_matrix
#
# Zero-copy summary:
#   - sparse slot vectors (@x, @i/@j, @p) are exposed to NumPy as zero-copy
#     views via reticulate::np_array() (C buffer protocol, owndata = FALSE).
#   - dgeMatrix: @x slot → 1-D zero-copy np_array, then reshaped Fortran-
#     order → 2-D zero-copy view (reshape returns a view, not a copy).
#   - base matrix: reticulate::r_to_py() uses NumPy's buffer protocol and
#     produces a Fortran-order float64 ndarray with owndata = FALSE (observed
#     behaviour; not formally guaranteed by the reticulate API).
#   - COO coordinates: @nzcoo is 1-based and must be decremented to 0-based
#     before passing to SciPy.  This requires a copy of the index arrays.
#     @nzdata values are still a zero-copy view.
#
# Lifetime caveat: reticulate::np_array() holds a *borrowed* C reference
# into the backing R vector.  Do NOT allow the source R objects to be
# modified or garbage-collected while the Python result is alive.
# ===========================================================================


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

#' Check that x is a SummarizedExperiment (or subclass)
#' @keywords internal
.validate_summarized_experiment <- function(x) {
  if (!methods::is(x, "SummarizedExperiment")) {
    stop(
      "'x' must be a SummarizedExperiment (or subclass), not '",
      paste(class(x), collapse = "', '"), "'."
    )
  }
  invisible(x)
}

#' Validate that assay is a single string present in x
#' @keywords internal
.validate_assay_name <- function(x, assay) {
  if (!is.character(assay) || length(assay) != 1L) {
    stop("'assay' must be a single string.")
  }
  avail <- SummarizedExperiment::assayNames(x)
  if (!assay %in% avail) {
    stop(
      "Assay '", assay, "' not found in object. ",
      "Available: ", paste0("'", avail, "'", collapse = ", "), "."
    )
  }
  invisible(assay)
}

#' Build a Python tuple suitable as a scipy `shape` argument
#' @keywords internal
.python_shape_tuple <- function(nr, nc) {
  reticulate::tuple(
    reticulate::r_to_py(as.integer(nr)),
    reticulate::r_to_py(as.integer(nc))
  )
}


# ---------------------------------------------------------------------------
# NumPy view helpers
# ---------------------------------------------------------------------------
#
# reticulate::np_array() uses NumPy's C buffer protocol to expose an R
# vector's memory to Python WITHOUT copying (flags.owndata = FALSE), as long
# as the R vector is already stored as the requested dtype:
#   * R numeric (double, 64-bit) → dtype = "float64"  — always zero-copy
#   * R integer (int32_t)        → dtype = "int32"    — always zero-copy
# If the R vector has the wrong type a copy is unavoidable (as.double /
# as.integer below), but the subsequent np_array call is still zero-copy.

#' Zero-copy NumPy view of an R numeric vector
#' @param x A numeric R vector.
#' @param dtype NumPy dtype string (default `"float64"`).
#' @return A Python `numpy.ndarray` (owndata = FALSE when x is already double).
#' @keywords internal
.as_numpy_view_numeric <- function(x, dtype = "float64") {
  if (!is.double(x)) x <- as.double(x)   # copy if coercion needed
  reticulate::np_array(x, dtype = dtype)  # buffer-protocol view, no copy
}

#' Zero-copy NumPy view of an R integer vector
#' @param x An integer R vector.
#' @param dtype NumPy dtype string (default `"int32"`).
#' @return A Python `numpy.ndarray` (owndata = FALSE when x is already integer).
#' @keywords internal
.as_numpy_view_integer <- function(x, dtype = "int32") {
  if (!is.integer(x)) x <- as.integer(x) # copy if coercion needed
  reticulate::np_array(x, dtype = dtype)  # buffer-protocol view, no copy
}


# ---------------------------------------------------------------------------
# Python object detection
# ---------------------------------------------------------------------------

#' Return TRUE if x is already a reticulate Python proxy object
#'
#' Useful for deciding whether layers / obsm / varm entries need conversion,
#' or whether they are already Python objects that should be forwarded as-is.
#'
#' @param x Any R object.
#' @return Logical scalar.
#' @keywords internal
.is_python_object <- function(x) {
  inherits(x, "python.builtin.object")
}


# ---------------------------------------------------------------------------
# Dense matrix converters
# ---------------------------------------------------------------------------

#' Convert a dense R matrix or dgeMatrix to a NumPy ndarray (zero-copy route)
#'
#' @description
#' ## Memory layout and zero-copy guarantees
#'
#' R stores matrices **column-major** (Fortran order); NumPy accepts both
#' layouts.
#'
#' * **`dgeMatrix`**: the `@x` slot is a plain R `double` vector laid out
#'   column-major.  `reticulate::np_array(@x)` creates a 1-D zero-copy view
#'   (`owndata = FALSE`).  Calling `.reshape((nr, nc), order = "F")` on that
#'   1-D array produces a Fortran-contiguous 2-D *view* (still `owndata =
#'   FALSE`) because the memory layout already matches — no copy occurs.
#'
#' * **base `matrix`**: `reticulate::r_to_py(mat)` exposes R's column-major
#'   memory to NumPy via the C buffer protocol, producing a Fortran-order
#'   `float64` ndarray with `owndata = FALSE`.  This is observed behaviour;
#'   the reticulate API does not formally guarantee it, but it is stable in
#'   practice.
#'
#' @param x A base R `matrix` or a `Matrix::dgeMatrix`.
#' @return A Python `numpy.ndarray` with dtype `float64` and shape
#'   `(nrow(x), ncol(x))`.
#' @keywords internal
.dense_matrix_to_numpy <- function(x) {
  if (methods::is(x, "dgeMatrix")) {
    # Zero-copy path via @x slot + Fortran-order reshape.
    # @x is the column-major data vector; @Dim is c(nrow, ncol).
    data_r <- if (is.double(x@x)) x@x else as.double(x@x)
    arr1d  <- reticulate::np_array(data_r, dtype = "float64") # owndata=FALSE
    nr     <- x@Dim[[1L]]
    nc     <- x@Dim[[2L]]
    # reshape(order='F') returns a Fortran-contiguous 2-D *view* because the
    # 1-D buffer is already laid out in column-major order. No copy occurs.
    return(arr1d$reshape(
      reticulate::tuple(as.integer(nr), as.integer(nc)),
      order = "F"
    ))
  }

  if (!is.matrix(x)) {
    # Coerce unknown dense types to base matrix with one R-level copy.
    message("Coercing '", class(x)[[1L]], "' to base matrix for NumPy conversion.")
    x <- as.matrix(x)
  }

  # Ensure double storage (no-op when already double).
  if (!is.double(x)) storage.mode(x) <- "double"

  # r_to_py on a base R matrix:
  #   Uses NumPy's buffer protocol; produces a Fortran-order float64 array
  #   that views R's column-major memory (owndata = FALSE, no copy).
  reticulate::r_to_py(x)
}


# ---------------------------------------------------------------------------
# Sparse matrix converters
# ---------------------------------------------------------------------------

#' Convert a dgCMatrix to scipy.sparse.csc_matrix (zero-copy slot views)
#'
#' @description
#' Slot mapping (all indices in R and scipy CSC are 0-based):
#'
#' | `dgCMatrix` slot | SciPy CSC field | dtype   |
#' |---|---|---|
#' | `@x`             | `data`          | float64 |
#' | `@i`             | `indices`       | int32   |
#' | `@p`             | `indptr`        | int32   |
#' | `@Dim`           | `shape`         | –       |
#'
#' `@i` (row indices) and `@p` (column pointers) are R `integer` (int32_t)
#' and are already 0-based → zero-copy `np_array` views.
#'
#' SciPy does not copy the backing NumPy arrays during `csc_matrix`
#' construction from the `(data, indices, indptr)` tuple — it stores
#' references to them.
#'
#' @param x A `Matrix::dgCMatrix`.
#' @return A Python `scipy.sparse.csc_matrix`.
#' @keywords internal
.dgCMatrix_to_scipy_csc <- function(x) {
  sp <- reticulate::import("scipy.sparse", convert = FALSE)

  py_data    <- .as_numpy_view_numeric(x@x, dtype = "float64") # owndata=FALSE
  py_indices <- .as_numpy_view_integer(x@i, dtype = "int32")   # 0-based, owndata=FALSE
  py_indptr  <- .as_numpy_view_integer(x@p, dtype = "int32")   # 0-based, owndata=FALSE
  shape      <- .python_shape_tuple(x@Dim[[1L]], x@Dim[[2L]])

  sp$csc_matrix(
    reticulate::tuple(py_data, py_indices, py_indptr),
    shape = shape
  )
}

#' Convert a dgRMatrix to scipy.sparse.csr_matrix (zero-copy slot views)
#'
#' @description
#' Slot mapping for `dgRMatrix` (CSR / Compressed Sparse Row format):
#'
#' | `dgRMatrix` slot | SciPy CSR field | dtype   |
#' |---|---|---|
#' | `@x`             | `data`          | float64 |
#' | `@j`             | `indices`       | int32   |
#' | `@p`             | `indptr`        | int32   |
#' | `@Dim`           | `shape`         | –       |
#'
#' `@j` stores **column** indices (0-based) and `@p` stores row pointers
#' (0-based), mirroring the CSR layout exactly.  Both are R `integer`
#' (int32_t) → zero-copy `np_array` views.
#'
#' @param x A `Matrix::dgRMatrix`.
#' @return A Python `scipy.sparse.csr_matrix`.
#' @keywords internal
.dgRMatrix_to_scipy_csr <- function(x) {
  sp <- reticulate::import("scipy.sparse", convert = FALSE)

  # @j = column indices (0-based); @p = row pointers (0-based).
  py_data    <- .as_numpy_view_numeric(x@x, dtype = "float64")
  py_indices <- .as_numpy_view_integer(x@j, dtype = "int32")   # col indices
  py_indptr  <- .as_numpy_view_integer(x@p, dtype = "int32")   # row pointers
  shape      <- .python_shape_tuple(x@Dim[[1L]], x@Dim[[2L]])

  sp$csr_matrix(
    reticulate::tuple(py_data, py_indices, py_indptr),
    shape = shape
  )
}

#' Convert a COO_SparseArray (2-D) to scipy.sparse.coo_matrix
#'
#' @description
#' `COO_SparseArray` (from the `SparseArray` package) stores coordinates in
#' `@nzcoo`, a 1-based integer matrix of shape `(nnz, ndim)`.  For 2-D
#' arrays:
#'
#' * `@nzcoo[, 1]` — row indices (1-based → subtract 1 for 0-based)
#' * `@nzcoo[, 2]` — col indices (1-based → subtract 1 for 0-based)
#' * `@nzdata`     — non-zero values
#' * `@dim`        — `c(nrow, ncol)`
#'
#' **Zero-copy is not possible for the coordinate arrays** because the
#' 1 → 0 base conversion requires allocating new integer vectors.  The value
#' array `@nzdata` is still a zero-copy NumPy view.
#'
#' @param x A `SparseArray::COO_SparseArray` (or `COO_SparseMatrix`) with
#'   exactly 2 dimensions.
#' @return A Python `scipy.sparse.coo_matrix`.
#' @keywords internal
.coo_matrix_to_scipy_coo <- function(x) {
  if (length(dim(x)) != 2L) {
    stop(
      ".coo_matrix_to_scipy_coo requires a 2-D COO_SparseArray; ",
      "got ", length(dim(x)), " dimensions."
    )
  }

  sp <- reticulate::import("scipy.sparse", convert = FALSE)

  # @nzdata: non-zero values → zero-copy numeric view.
  py_data <- .as_numpy_view_numeric(x@nzdata, dtype = "float64")

  # @nzcoo: 1-based coordinates (nnz x 2 integer matrix).
  # Subtracting 1 produces new R vectors (unavoidable copy for the 0-base
  # conversion).  We then view them zero-copy in NumPy.
  row_0 <- x@nzcoo[, 1L] - 1L  # copy here (1→0 base conversion)
  col_0 <- x@nzcoo[, 2L] - 1L  # copy here
  py_row <- .as_numpy_view_integer(row_0, dtype = "int32")
  py_col <- .as_numpy_view_integer(col_0, dtype = "int32")

  shape <- .python_shape_tuple(x@dim[[1L]], x@dim[[2L]])

  sp$coo_matrix(
    reticulate::tuple(py_data, reticulate::tuple(py_row, py_col)),
    shape = shape
  )
}


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

#' Dispatch an R matrix-like object to the appropriate Python converter
#'
#' @description
#' Dispatch order (first match wins):
#'
#' 1. `base matrix`     → [.dense_matrix_to_numpy()] → `numpy.ndarray`
#' 2. `dgeMatrix`       → [.dense_matrix_to_numpy()] → `numpy.ndarray`
#' 3. `dgCMatrix`       → [.dgCMatrix_to_scipy_csc()] → `scipy.sparse.csc_matrix`
#' 4. `dgRMatrix`       → [.dgRMatrix_to_scipy_csr()] → `scipy.sparse.csr_matrix`
#' 5. `COO_SparseArray` → [.coo_matrix_to_scipy_coo()] → `scipy.sparse.coo_matrix`
#' 6. Fallback: if `prefer_sparse = TRUE`, coerce to `dgCMatrix` with a message
#'    then apply path 3; otherwise densify to base matrix with a message then
#'    apply path 1.
#'
#' @param x An R matrix-like object.
#' @param prefer_sparse Logical.  Controls the fallback coercion direction.
#'   Default `TRUE`.
#' @return A Python `numpy.ndarray` or `scipy.sparse` matrix
#'   (`convert = FALSE`).
#' @keywords internal
.matrix_to_python_array_or_sparse <- function(x, prefer_sparse = TRUE) {
  if (is.matrix(x)) {
    return(.dense_matrix_to_numpy(x))
  }
  if (methods::is(x, "dgeMatrix")) {
    return(.dense_matrix_to_numpy(x))
  }
  if (methods::is(x, "dgCMatrix")) {
    return(.dgCMatrix_to_scipy_csc(x))
  }
  if (methods::is(x, "dgRMatrix")) {
    return(.dgRMatrix_to_scipy_csr(x))
  }
  # COO_SparseArray is the abstract parent; COO_SparseMatrix is the 2-D
  # concrete subclass.  Both are handled by .coo_matrix_to_scipy_coo().
  if (methods::is(x, "COO_SparseArray")) {
    return(.coo_matrix_to_scipy_coo(x))
  }

  # Fallback: coerce with an informative message.
  if (prefer_sparse) {
    message(
      "No direct converter for '", class(x)[[1L]], "'; ",
      "coercing to dgCMatrix for SciPy CSC conversion."
    )
    mat <- tryCatch(
      methods::as(x, "dgCMatrix"),
      error = function(e) methods::as(methods::as(x, "dMatrix"), "dgCMatrix")
    )
    return(.dgCMatrix_to_scipy_csc(mat))
  } else {
    message(
      "No direct converter for '", class(x)[[1L]], "'; ",
      "densifying to base matrix for NumPy conversion."
    )
    return(.dense_matrix_to_numpy(as.matrix(x)))
  }
}


# ---------------------------------------------------------------------------
# Container-level helper (exported)
# ---------------------------------------------------------------------------

#' Extract a SummarizedExperiment assay and convert to a Python matrix
#'
#' @description
#' Extracts the named assay from a
#' [SummarizedExperiment::SummarizedExperiment] (or any subclass, including
#' [SingleCellExperiment::SingleCellExperiment]), validates it, and dispatches
#' to the appropriate Python converter via
#' [.matrix_to_python_array_or_sparse()].
#'
#' The returned Python object has shape `(n_features, n_samples)` — the same
#' orientation as the R assay (rows = features, columns = samples).  For
#' AnnData — which expects `(n_obs, n_vars) = (cells, features)` — the caller
#' is responsible for transposing: `X <- se_assay_to_python_matrix(...); X$T`.
#'
#' ## Supported assay matrix types
#' | R class           | Python result              | Zero-copy? |
#' |---|---|---|
#' | `base matrix`     | `numpy.ndarray`            | yes (r_to_py) |
#' | `dgeMatrix`       | `numpy.ndarray`            | yes (@x + reshape) |
#' | `dgCMatrix`       | `scipy.sparse.csc_matrix`  | yes (slot views) |
#' | `dgRMatrix`       | `scipy.sparse.csr_matrix`  | yes (slot views) |
#' | `COO_SparseArray` | `scipy.sparse.coo_matrix`  | values yes; coords no |
#' | other             | `csc_matrix` or `ndarray`  | no (coercion) |
#'
#' @param x A [SummarizedExperiment::SummarizedExperiment] or subclass.
#' @param assay A single string: the assay name to extract (default
#'   `"counts"`).
#' @param prefer_sparse Logical.  If the assay matrix has no direct converter,
#'   coerce to `dgCMatrix` when `TRUE` (default), otherwise densify.
#'
#' @return A Python `numpy.ndarray` or `scipy.sparse` matrix
#'   (`convert = FALSE`) with shape `(n_features, n_samples)`.
#'
#' @seealso [sce_to_anndata()], [sce_assay_to_scipy_csc()]
#'
#' @importFrom SummarizedExperiment assayNames assay
#' @export
se_assay_to_python_matrix <- function(x, assay = "counts", prefer_sparse = TRUE) {
  .validate_summarized_experiment(x)
  .validate_assay_name(x, assay)
  mat <- SummarizedExperiment::assay(x, assay)
  .matrix_to_python_array_or_sparse(mat, prefer_sparse = prefer_sparse)
}
