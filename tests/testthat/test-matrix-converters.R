# tests/testthat/test-matrix-converters.R
#
# Tests for the matrix conversion utilities in R/matrix_converters.R.
#
# Covers:
#   * base dense matrix         → numpy.ndarray
#   * Matrix::dgeMatrix         → numpy.ndarray
#   * Matrix::dgCMatrix         → scipy.sparse.csc_matrix
#   * Matrix::dgRMatrix         → scipy.sparse.csr_matrix
#   * SparseArray::COO_SparseMatrix → scipy.sparse.coo_matrix
#   * se_assay_to_python_matrix() dispatcher (SummarizedExperiment)
#   * sce_to_anndata() with dense assay, sparse assay, extra layer, obsm
#
# Numerical round-trip checks use small matrices with known values.

library(Matrix)
library(methods)


# ---------------------------------------------------------------------------
# Small helper fixtures
# ---------------------------------------------------------------------------

# 4 × 3 numeric matrix with known values (column-major in R)
make_small_dense <- function() {
  matrix(c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0),
         nrow = 4L, ncol = 3L)
}

make_small_dge <- function() {
  methods::as(make_small_dense(), "dgeMatrix")
}

# 5 × 4 sparse matrix with a few known non-zeros
make_small_dgC <- function() {
  # row 0, col 0: 1.0;  row 2, col 1: 2.5;  row 4, col 3: 3.0
  i <- c(0L, 2L, 4L)
  j <- c(0L, 1L, 3L)
  x <- c(1.0, 2.5, 3.0)
  methods::as(
    Matrix::sparseMatrix(i = i + 1L, j = j + 1L, x = x, dims = c(5L, 4L)),
    "dgCMatrix"
  )
}

make_small_dgR <- function() {
  # as(..., "RsparseMatrix") is the supported coercion to dgRMatrix.
  methods::as(make_small_dgC(), "RsparseMatrix")
}

make_small_coo <- function() {
  skip_if_not_installed("SparseArray")
  methods::as(make_small_dgC(), "COO_SparseMatrix")
}


# ---------------------------------------------------------------------------
# 1. Validation helpers
# ---------------------------------------------------------------------------

test_that(".validate_summarized_experiment rejects non-SE objects", {
  expect_error(
    .validate_summarized_experiment(list()),
    regexp = "SummarizedExperiment"
  )
  expect_error(
    .validate_summarized_experiment(data.frame(x = 1)),
    regexp = "SummarizedExperiment"
  )
})

test_that(".validate_summarized_experiment accepts SE and subclasses", {
  se  <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = make_small_dgC())
  )
  expect_silent(.validate_summarized_experiment(se))

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = make_small_dgC())
  )
  expect_silent(.validate_summarized_experiment(sce))
})

test_that(".validate_assay_name rejects missing assay", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = make_small_dgC())
  )
  expect_error(.validate_assay_name(se, "logcounts"), regexp = "logcounts")
})

test_that(".validate_assay_name rejects non-string", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = make_small_dgC())
  )
  expect_error(.validate_assay_name(se, 1L), regexp = "single string")
})


# ---------------------------------------------------------------------------
# 2. .is_python_object
# ---------------------------------------------------------------------------

test_that(".is_python_object returns FALSE for plain R objects", {
  expect_false(.is_python_object(matrix(1)))
  expect_false(.is_python_object(list()))
  expect_false(.is_python_object(NULL))
})

test_that(".is_python_object returns TRUE for Python objects", {
  skip_if_no_scipy()
  np  <- reticulate::import("numpy", convert = FALSE)
  arr <- reticulate::np_array(c(1.0, 2.0))
  expect_true(.is_python_object(arr))
})


# ---------------------------------------------------------------------------
# 3. Dense matrix → NumPy ndarray
# ---------------------------------------------------------------------------

test_that(".dense_matrix_to_numpy: base matrix — correct shape and values", {
  skip_if_no_scipy()
  m    <- make_small_dense()
  arr  <- .dense_matrix_to_numpy(m)

  expect_true(reticulate::is_py_object(arr))

  shape <- reticulate::py_to_r(arr$shape)
  expect_equal(as.integer(shape[[1L]]), nrow(m))
  expect_equal(as.integer(shape[[2L]]), ncol(m))

  # Round-trip: convert back to R and compare
  m_rt <- reticulate::py_to_r(arr)
  expect_equal(m_rt, m)
})

test_that(".dense_matrix_to_numpy: base matrix — owndata = FALSE (zero-copy observed)", {
  skip_if_no_scipy()
  m   <- make_small_dense()
  arr <- .dense_matrix_to_numpy(m)
  # r_to_py on a base R matrix is observed to produce owndata=FALSE.
  # This is not guaranteed by the reticulate API, so we report rather than
  # hard-fail if it changes.
  owndata <- reticulate::py_to_r(arr$flags$owndata)
  expect_false(owndata,
    label = "base matrix: r_to_py should produce owndata=FALSE (zero-copy)"
  )
})

test_that(".dense_matrix_to_numpy: dgeMatrix — correct shape and values", {
  skip_if_no_scipy()
  m   <- make_small_dge()
  arr <- .dense_matrix_to_numpy(m)

  expect_true(reticulate::is_py_object(arr))

  shape <- reticulate::py_to_r(arr$shape)
  expect_equal(as.integer(shape[[1L]]), nrow(m))
  expect_equal(as.integer(shape[[2L]]), ncol(m))

  m_rt <- reticulate::py_to_r(arr)
  expect_equal(m_rt, as.matrix(m))
})

test_that(".dense_matrix_to_numpy: dgeMatrix — 1D slot view is owndata=FALSE", {
  skip_if_no_scipy()
  m   <- make_small_dge()
  # The 1D np_array of @x should have owndata=FALSE (buffer-protocol view).
  arr1d <- reticulate::np_array(m@x, dtype = "float64")
  owndata <- reticulate::py_to_r(arr1d$flags$owndata)
  expect_false(owndata,
    label = "np_array of @x vector: owndata should be FALSE"
  )
})


# ---------------------------------------------------------------------------
# 4. dgCMatrix → scipy.sparse.csc_matrix
# ---------------------------------------------------------------------------

test_that(".dgCMatrix_to_scipy_csc: correct Python class", {
  skip_if_no_scipy()
  m   <- make_small_dgC()
  csc <- .dgCMatrix_to_scipy_csc(m)

  sp       <- reticulate::import("scipy.sparse", convert = FALSE)
  is_csc   <- reticulate::py_to_r(sp$isspmatrix_csc(csc))
  expect_true(is_csc, label = "result must be scipy csc_matrix")
})

test_that(".dgCMatrix_to_scipy_csc: shape matches input", {
  skip_if_no_scipy()
  m     <- make_small_dgC()
  csc   <- .dgCMatrix_to_scipy_csc(m)
  shape <- reticulate::py_to_r(csc$shape)
  expect_equal(as.integer(shape[[1L]]), nrow(m))
  expect_equal(as.integer(shape[[2L]]), ncol(m))
})

test_that(".dgCMatrix_to_scipy_csc: nnz and data match @x slot", {
  skip_if_no_scipy()
  m     <- make_small_dgC()
  csc   <- .dgCMatrix_to_scipy_csc(m)
  nnz   <- as.integer(reticulate::py_to_r(csc$nnz))
  expect_equal(nnz, length(m@x))

  # py_to_r on a 1-D NumPy array may return with a dim attribute; as.vector
  # strips it so the comparison is against a plain R vector.
  data_py <- as.vector(reticulate::py_to_r(csc$data))
  expect_equal(sort(data_py), sort(as.vector(m@x)))
})

test_that(".dgCMatrix_to_scipy_csc: data array is owndata=FALSE (zero-copy)", {
  skip_if_no_scipy()
  m   <- make_small_dgC()
  csc <- .dgCMatrix_to_scipy_csc(m)
  owndata <- reticulate::py_to_r(csc$data$flags$owndata)
  expect_false(owndata, label = "CSC data slot: owndata should be FALSE")
})

test_that(".dgCMatrix_to_scipy_csc: roundtrip values correct", {
  skip_if_no_scipy()
  m   <- make_small_dgC()
  csc <- .dgCMatrix_to_scipy_csc(m)
  # Convert to dense for value comparison
  dense_py <- reticulate::py_to_r(csc$toarray())
  dense_r  <- as.matrix(m)
  expect_equal(dense_py, dense_r)
})


# ---------------------------------------------------------------------------
# 5. dgRMatrix → scipy.sparse.csr_matrix
# ---------------------------------------------------------------------------

test_that(".dgRMatrix_to_scipy_csr: correct Python class", {
  skip_if_no_scipy()
  m   <- make_small_dgR()
  csr <- .dgRMatrix_to_scipy_csr(m)

  sp     <- reticulate::import("scipy.sparse", convert = FALSE)
  is_csr <- reticulate::py_to_r(sp$isspmatrix_csr(csr))
  expect_true(is_csr, label = "result must be scipy csr_matrix")
})

test_that(".dgRMatrix_to_scipy_csr: shape matches input", {
  skip_if_no_scipy()
  m     <- make_small_dgR()
  csr   <- .dgRMatrix_to_scipy_csr(m)
  shape <- reticulate::py_to_r(csr$shape)
  expect_equal(as.integer(shape[[1L]]), nrow(m))
  expect_equal(as.integer(shape[[2L]]), ncol(m))
})

test_that(".dgRMatrix_to_scipy_csr: roundtrip values correct", {
  skip_if_no_scipy()
  m       <- make_small_dgR()
  csr     <- .dgRMatrix_to_scipy_csr(m)
  dense_py <- reticulate::py_to_r(csr$toarray())
  dense_r  <- as.matrix(m)
  expect_equal(dense_py, dense_r)
})

test_that(".dgRMatrix_to_scipy_csr: data array is owndata=FALSE", {
  skip_if_no_scipy()
  m   <- make_small_dgR()
  csr <- .dgRMatrix_to_scipy_csr(m)
  owndata <- reticulate::py_to_r(csr$data$flags$owndata)
  expect_false(owndata, label = "CSR data slot: owndata should be FALSE")
})


# ---------------------------------------------------------------------------
# 6. COO_SparseMatrix → scipy.sparse.coo_matrix
# ---------------------------------------------------------------------------

test_that(".coo_matrix_to_scipy_coo: correct Python class", {
  skip_if_no_scipy()
  skip_if_not_installed("SparseArray")
  m   <- make_small_coo()
  coo <- .coo_matrix_to_scipy_coo(m)

  sp     <- reticulate::import("scipy.sparse", convert = FALSE)
  is_coo <- reticulate::py_to_r(sp$isspmatrix_coo(coo))
  expect_true(is_coo, label = "result must be scipy coo_matrix")
})

test_that(".coo_matrix_to_scipy_coo: shape matches input", {
  skip_if_no_scipy()
  skip_if_not_installed("SparseArray")
  m     <- make_small_coo()
  coo   <- .coo_matrix_to_scipy_coo(m)
  shape <- reticulate::py_to_r(coo$shape)
  expect_equal(as.integer(shape[[1L]]), m@dim[[1L]])
  expect_equal(as.integer(shape[[2L]]), m@dim[[2L]])
})

test_that(".coo_matrix_to_scipy_coo: roundtrip values and nnz correct", {
  skip_if_no_scipy()
  skip_if_not_installed("SparseArray")
  m       <- make_small_coo()
  coo     <- .coo_matrix_to_scipy_coo(m)
  nnz     <- as.integer(reticulate::py_to_r(coo$nnz))
  expect_equal(nnz, length(m@nzdata))

  dense_py <- reticulate::py_to_r(coo$toarray())
  # Reference: convert the original dgCMatrix to dense
  dense_r  <- as.matrix(make_small_dgC())
  expect_equal(dense_py, dense_r)
})

test_that(".coo_matrix_to_scipy_coo: 0-based coordinates are correct", {
  skip_if_no_scipy()
  skip_if_not_installed("SparseArray")
  m   <- make_small_coo()
  coo <- .coo_matrix_to_scipy_coo(m)

  # @nzcoo is 1-based; subtract 1 should match coo.row and coo.col
  expected_row <- m@nzcoo[, 1L] - 1L
  expected_col <- m@nzcoo[, 2L] - 1L

  py_row <- as.vector(reticulate::py_to_r(coo$row))
  py_col <- as.vector(reticulate::py_to_r(coo$col))

  # Sort both by (row, col) for comparison
  ord_r  <- order(expected_row, expected_col)
  ord_py <- order(py_row, py_col)
  expect_equal(expected_row[ord_r], as.integer(py_row[ord_py]))
  expect_equal(expected_col[ord_r], as.integer(py_col[ord_py]))
})

test_that(".coo_matrix_to_scipy_coo: data value array is owndata=FALSE", {
  skip_if_no_scipy()
  skip_if_not_installed("SparseArray")
  m   <- make_small_coo()
  coo <- .coo_matrix_to_scipy_coo(m)
  # @nzdata is exposed via np_array (zero-copy); coordinates are copied.
  owndata <- reticulate::py_to_r(coo$data$flags$owndata)
  expect_false(owndata, label = "COO data: owndata should be FALSE")
})


# ---------------------------------------------------------------------------
# 7. Dispatcher: .matrix_to_python_array_or_sparse
# ---------------------------------------------------------------------------

test_that("dispatcher routes base matrix to numpy", {
  skip_if_no_scipy()
  m    <- make_small_dense()
  py   <- .matrix_to_python_array_or_sparse(m, prefer_sparse = FALSE)
  np   <- reticulate::import("numpy", convert = FALSE)
  expect_true(reticulate::py_to_r(np$ndim(py)) == 2L)
})

test_that("dispatcher routes dgeMatrix to numpy", {
  skip_if_no_scipy()
  py   <- .matrix_to_python_array_or_sparse(make_small_dge())
  np   <- reticulate::import("numpy", convert = FALSE)
  expect_true(reticulate::py_to_r(np$ndim(py)) == 2L)
})

test_that("dispatcher routes dgCMatrix to scipy csc", {
  skip_if_no_scipy()
  sp  <- reticulate::import("scipy.sparse", convert = FALSE)
  py  <- .matrix_to_python_array_or_sparse(make_small_dgC())
  expect_true(reticulate::py_to_r(sp$isspmatrix_csc(py)))
})

test_that("dispatcher routes dgRMatrix to scipy csr", {
  skip_if_no_scipy()
  sp  <- reticulate::import("scipy.sparse", convert = FALSE)
  py  <- .matrix_to_python_array_or_sparse(make_small_dgR())
  expect_true(reticulate::py_to_r(sp$isspmatrix_csr(py)))
})

test_that("dispatcher routes COO_SparseMatrix to scipy coo", {
  skip_if_no_scipy()
  skip_if_not_installed("SparseArray")
  sp  <- reticulate::import("scipy.sparse", convert = FALSE)
  py  <- .matrix_to_python_array_or_sparse(make_small_coo())
  expect_true(reticulate::py_to_r(sp$isspmatrix_coo(py)))
})

test_that("dispatcher prefer_sparse=FALSE densifies unknown sparse type", {
  skip_if_no_scipy()
  # lgCMatrix is not directly handled: should trigger fallback message + densify
  lmat <- methods::as(
    Matrix::rsparsematrix(6L, 5L, 0.4, repr = "C") > 0,
    "lgCMatrix"
  )
  np <- reticulate::import("numpy", convert = FALSE)
  expect_message(
    py <- .matrix_to_python_array_or_sparse(lmat, prefer_sparse = FALSE),
    regexp = "densifying|Coercing"
  )
  expect_true(reticulate::py_to_r(np$ndim(py)) == 2L)
})


# ---------------------------------------------------------------------------
# 8. se_assay_to_python_matrix container helper
# ---------------------------------------------------------------------------

test_that("se_assay_to_python_matrix rejects non-SE input", {
  expect_error(
    se_assay_to_python_matrix(list()),
    regexp = "SummarizedExperiment"
  )
})

test_that("se_assay_to_python_matrix rejects bad assay name", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = make_small_dgC())
  )
  expect_error(
    se_assay_to_python_matrix(se, assay = "xyz"),
    regexp = "xyz"
  )
})

test_that("se_assay_to_python_matrix: plain SE with dgCMatrix assay", {
  skip_if_no_scipy()
  m  <- make_small_dgC()
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))
  py <- se_assay_to_python_matrix(se, assay = "counts")

  # Shape must be (n_features, n_samples) = (nrow, ncol) of assay
  sp    <- reticulate::import("scipy.sparse", convert = FALSE)
  shape <- reticulate::py_to_r(py$shape)
  expect_equal(as.integer(shape[[1L]]), nrow(m))
  expect_equal(as.integer(shape[[2L]]), ncol(m))
  expect_true(reticulate::py_to_r(sp$isspmatrix_csc(py)))
})

test_that("se_assay_to_python_matrix: plain SE with dense assay", {
  skip_if_no_scipy()
  m  <- make_small_dense()
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))
  np <- reticulate::import("numpy", convert = FALSE)
  py <- se_assay_to_python_matrix(se, assay = "counts", prefer_sparse = FALSE)

  shape <- reticulate::py_to_r(py$shape)
  expect_equal(as.integer(shape[[1L]]), nrow(m))
  expect_equal(as.integer(shape[[2L]]), ncol(m))
  expect_true(reticulate::py_to_r(np$ndim(py)) == 2L)
})

test_that("se_assay_to_python_matrix: SCE subclass also accepted", {
  skip_if_no_scipy()
  m   <- make_small_dgC()
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
  py  <- se_assay_to_python_matrix(sce, assay = "counts")
  expect_true(reticulate::is_py_object(py))
})

test_that("se_assay_to_python_matrix: orientation is (n_features, n_samples)", {
  skip_if_no_scipy()
  m  <- make_small_dgC()  # 5 rows (features), 4 cols (samples)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))
  py <- se_assay_to_python_matrix(se)

  shape <- reticulate::py_to_r(py$shape)
  expect_equal(as.integer(shape[[1L]]), 5L, label = "dim 0 = n_features")
  expect_equal(as.integer(shape[[2L]]), 4L, label = "dim 1 = n_samples")
})


# ---------------------------------------------------------------------------
# 9. sce_to_anndata integration tests with new helpers
# ---------------------------------------------------------------------------

test_that("sce_to_anndata accepts plain SummarizedExperiment (not only SCE)", {
  skip_if_no_scipy()
  m     <- make_small_dgC()
  se    <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))
  adata <- sce_to_anndata(se)
  expect_true(reticulate::is_py_object(adata))
  shape <- reticulate::py_to_r(adata$shape)
  expect_equal(as.integer(shape[[1L]]), ncol(m), label = "n_obs = n_cells")
  expect_equal(as.integer(shape[[2L]]), nrow(m), label = "n_vars = n_features")
})

test_that("sce_to_anndata: sparse assay produces sparse X", {
  skip_if_no_scipy()
  m     <- make_small_dgC()
  se    <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))
  adata <- sce_to_anndata(se)
  sp    <- reticulate::import("scipy.sparse", convert = FALSE)
  expect_true(reticulate::py_to_r(sp$issparse(adata$X)))
})

test_that("sce_to_anndata: dense assay produces ndarray X", {
  skip_if_no_scipy()
  m     <- make_small_dense()
  se    <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))
  adata <- sce_to_anndata(se)
  np    <- reticulate::import("numpy", convert = FALSE)
  expect_true(reticulate::py_to_r(np$ndim(adata$X)) == 2L)
  shape <- reticulate::py_to_r(adata$X$shape)
  expect_equal(as.integer(shape[[1L]]), ncol(m))
  expect_equal(as.integer(shape[[2L]]), nrow(m))
})

test_that("sce_to_anndata: extra sparse layer round-trips correctly", {
  skip_if_no_scipy()
  m      <- make_small_dgC()
  layer  <- make_small_dgC()   # same shape, re-use fixture
  se     <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))
  adata  <- sce_to_anndata(se, layers = list(raw = layer))
  sp     <- reticulate::import("scipy.sparse", convert = FALSE)

  # The layer should exist and be sparse
  py_layer <- adata$layers$`__getitem__`("raw")
  expect_true(reticulate::py_to_r(sp$issparse(py_layer)))

  # Values should round-trip (compare dense representations; note AnnData
  # stores layers in obs × vars orientation, i.e. transposed vs assay)
  dense_py <- reticulate::py_to_r(py_layer$toarray())
  dense_r  <- t(as.matrix(layer))   # transpose to (cells, features)
  expect_equal(dense_py, dense_r)
})

test_that("sce_to_anndata: dense obsm entry is stored correctly", {
  skip_if_no_scipy()
  set.seed(42L)
  m_assay <- make_small_dgC()   # 5 features × 4 cells
  # obsm entry must have n_obs rows, i.e. ncol(assay) rows
  emb     <- matrix(rnorm(4L * 2L), nrow = 4L, ncol = 2L)

  se    <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = m_assay)
  )
  adata <- sce_to_anndata(se, obsm = list(X_pca = emb))
  np    <- reticulate::import("numpy", convert = FALSE)

  py_emb <- adata$obsm$`__getitem__`("X_pca")
  expect_true(reticulate::py_to_r(np$ndim(py_emb)) == 2L)

  shape <- reticulate::py_to_r(py_emb$shape)
  expect_equal(as.integer(shape[[1L]]), 4L, label = "obsm rows = n_obs")
  expect_equal(as.integer(shape[[2L]]), 2L, label = "obsm cols = embedding dim")

  # Round-trip values
  emb_rt <- reticulate::py_to_r(py_emb)
  expect_equal(emb_rt, emb)
})

test_that("sce_to_anndata: already-Python layer forwarded without re-conversion", {
  skip_if_no_scipy()
  m     <- make_small_dgC()
  se    <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))

  # Build a Python CSR layer in (n_obs × n_vars) = (cells × features) shape,
  # matching AnnData's expected layer orientation. Passing this as a Python
  # object means the function forwards it without further conversion.
  py_layer <- .dgCMatrix_to_scipy_csc(make_small_dgC())$T   # (4 cells × 5 features)

  adata  <- sce_to_anndata(se, layers = list(prebuilt = py_layer))
  sp     <- reticulate::import("scipy.sparse", convert = FALSE)
  result <- adata$layers$`__getitem__`("prebuilt")
  expect_true(reticulate::py_to_r(sp$issparse(result)))

  # Shape must be (n_obs, n_vars)
  shape <- reticulate::py_to_r(result$shape)
  expect_equal(as.integer(shape[[1L]]), ncol(m), label = "layer rows = n_obs")
  expect_equal(as.integer(shape[[2L]]), nrow(m), label = "layer cols = n_vars")
})
