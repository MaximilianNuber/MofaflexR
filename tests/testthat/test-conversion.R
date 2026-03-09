library(Matrix)
library(SingleCellExperiment)
library(methods)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

skip_if_no_python <- function() {
  testthat::skip_if_not(
    tryCatch({
      reticulate::py_available(initialize = TRUE) &&
        reticulate::py_module_available("numpy") &&
        reticulate::py_module_available("scipy") &&
        reticulate::py_module_available("anndata")
    }, error = function(e) FALSE),
    message = paste(
      "Python with numpy/scipy/anndata not available;",
      "set RETICULATE_PYTHON to a suitable environment to run these tests."
    )
  )
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

make_sce <- function(nr = 30L, nc = 20L, density = 0.3) {
  set.seed(42L)
  counts <- Matrix::rsparsematrix(nr, nc, density = density, repr = "C")
  counts <- abs(round(counts))
  counts <- methods::as(counts, "dgCMatrix")

  obs_df <- data.frame(
    cell_type = sample(c("A", "B"), nc, replace = TRUE),
    row.names  = paste0("cell", seq_len(nc))
  )
  var_df <- data.frame(
    gene_group = c("coding", "ncRNA")[sample(1:2, nr, replace = TRUE)],
    row.names   = paste0("gene", seq_len(nr))
  )
  SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = counts),
    colData = S4Vectors::DataFrame(obs_df),
    rowData = S4Vectors::DataFrame(var_df)
  )
}

# ---------------------------------------------------------------------------
# 1. Input validation
# ---------------------------------------------------------------------------

test_that("sce_assay_to_scipy_csc rejects non-SCE input", {
  expect_error(
    sce_assay_to_scipy_csc(list(), assay = "counts"),
    regexp = "SummarizedExperiment",
    info   = "Should mention expected class in error message"
  )
})

test_that("sce_assay_to_scipy_csc rejects missing assay name", {
  sce <- make_sce()
  expect_error(
    sce_assay_to_scipy_csc(sce, assay = "logcounts"),
    regexp = "logcounts",
    info   = "Error should mention the missing assay name"
  )
})

test_that("sce_assay_to_scipy_csc rejects non-character assay argument", {
  sce <- make_sce()
  expect_error(
    sce_assay_to_scipy_csc(sce, assay = 1L),
    regexp = "single string",
    info   = "Error should mention 'single string'"
  )
})

test_that("sce_to_anndata rejects non-SCE input", {
  expect_error(
    sce_to_anndata(data.frame(x = 1:3)),
    regexp = "SummarizedExperiment"
  )
})

test_that("sce_to_reticulate_anndata rejects non-SCE input", {
  expect_error(
    sce_to_reticulate_anndata("not a sce"),
    regexp = "SummarizedExperiment"
  )
})

# ---------------------------------------------------------------------------
# 2. Non-dgCMatrix assay is coerced (not rejected)
# ---------------------------------------------------------------------------

test_that("sce_assay_to_scipy_csc coerces lgCMatrix to dgCMatrix with message", {
  skip_if_no_python()

  set.seed(42L)
  lmat <- methods::as(
    Matrix::rsparsematrix(10L, 8L, 0.3, repr = "C") > 0,
    "lgCMatrix"
  )
  sce_lg <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = lmat)
  )
  # Should emit a message but succeed
  expect_message(
    csc <- sce_assay_to_scipy_csc(sce_lg),
    regexp = "coercing to dgCMatrix"
  )
  # Result is still a Python object
  expect_true(reticulate::is_py_object(csc))
})

# ---------------------------------------------------------------------------
# 3. scipy.sparse.csc_matrix shape and nnz
# ---------------------------------------------------------------------------

test_that("sce_assay_to_scipy_csc returns correct shape and nnz", {
  skip_if_no_python()

  nr  <- 30L
  nc  <- 20L
  sce <- make_sce(nr = nr, nc = nc)
  mat <- methods::as(SummarizedExperiment::assay(sce, "counts"), "dgCMatrix")
  expected_nnz <- length(mat@x)

  csc <- sce_assay_to_scipy_csc(sce, assay = "counts")

  # Shape
  shape_r <- reticulate::py_to_r(csc$shape)
  expect_equal(as.integer(shape_r[[1L]]), nr)
  expect_equal(as.integer(shape_r[[2L]]), nc)

  # nnz
  actual_nnz <- reticulate::py_to_r(csc$nnz)
  expect_equal(as.integer(actual_nnz), expected_nnz)
})

# ---------------------------------------------------------------------------
# 3.5 Zero-copy verification: owndata = FALSE
# ---------------------------------------------------------------------------

test_that("sce_assay_to_scipy_csc: data array is a zero-copy view (owndata=FALSE)", {
  skip_if_no_python()

  sce <- make_sce()
  csc <- sce_assay_to_scipy_csc(sce, assay = "counts")

  # reticulate::np_array() exposes R memory via the C buffer protocol;
  # the resulting numpy array does NOT own its data (owndata=False).
  owndata <- reticulate::py_to_r(csc$data$flags$owndata)
  expect_false(owndata, label = "csc.data.flags.owndata should be FALSE (zero-copy)")
})

# ---------------------------------------------------------------------------
# 4. Slot lengths: data, indices, indptr
# ---------------------------------------------------------------------------

test_that("sce_assay_to_scipy_csc: data/indices/indptr lengths match dgCMatrix slots", {
  skip_if_no_python()

  sce <- make_sce()
  mat <- methods::as(SummarizedExperiment::assay(sce, "counts"), "dgCMatrix")
  csc <- sce_assay_to_scipy_csc(sce, assay = "counts")

  data_len    <- reticulate::py_to_r(reticulate::py_len(csc$data))
  indices_len <- reticulate::py_to_r(reticulate::py_len(csc$indices))
  indptr_len  <- reticulate::py_to_r(reticulate::py_len(csc$indptr))

  expect_equal(as.integer(data_len),    length(mat@x))
  expect_equal(as.integer(indices_len), length(mat@i))
  expect_equal(as.integer(indptr_len),  length(mat@p))
})

# ---------------------------------------------------------------------------
# 5. No densification: result must be sparse
# ---------------------------------------------------------------------------

test_that("sce_assay_to_scipy_csc: result is scipy sparse (not dense)", {
  skip_if_no_python()

  sce <- make_sce()
  csc <- sce_assay_to_scipy_csc(sce, assay = "counts")

  sp      <- reticulate::import("scipy.sparse", convert = FALSE)
  is_spar <- reticulate::py_to_r(sp$issparse(csc))
  expect_true(is_spar, label = "scipy.sparse.issparse(X) should be TRUE")
})

test_that("sce_to_anndata: adata$X is scipy sparse (not dense)", {
  skip_if_no_python()

  sce    <- make_sce()
  adata  <- sce_to_anndata(sce)
  sp     <- reticulate::import("scipy.sparse", convert = FALSE)
  is_spar <- reticulate::py_to_r(sp$issparse(adata$X))
  expect_true(is_spar, label = "adata.X should remain sparse")
})

# ---------------------------------------------------------------------------
# 6. AnnData dimensions match SCE
# ---------------------------------------------------------------------------

test_that("sce_to_anndata: AnnData obs/var dimensions match SCE", {
  skip_if_no_python()

  nr  <- 30L
  nc  <- 20L
  sce <- make_sce(nr = nr, nc = nc)
  adata <- sce_to_anndata(sce)

  shape_r <- reticulate::py_to_r(adata$shape)
  expect_equal(as.integer(shape_r[[1L]]), nc) # AnnData: rows = cells = n_obs
  expect_equal(as.integer(shape_r[[2L]]), nr) # AnnData: cols = genes = n_vars
})

# ---------------------------------------------------------------------------
# 7. obs and var names are preserved
# ---------------------------------------------------------------------------

test_that("sce_to_anndata: obs_names and var_names match SCE names", {
  skip_if_no_python()

  sce     <- make_sce(nr = 10L, nc = 8L)
  adata   <- sce_to_anndata(sce)

  obs_names_py <- reticulate::py_to_r(adata$obs_names$tolist())
  var_names_py <- reticulate::py_to_r(adata$var_names$tolist())

  expect_equal(obs_names_py, colnames(sce))
  expect_equal(var_names_py, rownames(sce))
})

# ---------------------------------------------------------------------------
# 8. ReticulateAnnData smoke test
# ---------------------------------------------------------------------------

test_that("sce_to_reticulate_anndata returns a ReticulateAnnData-like object", {
  skip_if_no_python()

  sce  <- make_sce()
  rada <- sce_to_reticulate_anndata(sce)

  # Must be an R6 object with $py_anndata method
  expect_true(
    R6::is.R6(rada),
    label = "Returned object should be R6"
  )
  expect_true(
    is.function(rada$py_anndata),
    label = "ReticulateAnnData should expose $py_anndata()"
  )

  py_obj <- rada$py_anndata()
  expect_true(
    inherits(py_obj, "anndata._core.anndata.AnnData"),
    label = "py_anndata() should return a Python AnnData object"
  )
})

test_that("sce_to_reticulate_anndata: underlying AnnData$X is still sparse", {
  skip_if_no_python()

  sce    <- make_sce()
  rada   <- sce_to_reticulate_anndata(sce)
  py_obj <- rada$py_anndata()
  sp     <- reticulate::import("scipy.sparse", convert = FALSE)
  is_spar <- reticulate::py_to_r(sp$issparse(py_obj$X))
  expect_true(is_spar, label = "X should still be sparse after ReticulateAnnData wrapping")
})

# ---------------------------------------------------------------------------
# 9. mofaflex_basilisk_env returns a BasiliskEnvironment
# ---------------------------------------------------------------------------

test_that("mofaflex_basilisk_env returns a BasiliskEnvironment", {
  env <- mofaflex_basilisk_env()
  expect_s4_class(env, "BasiliskEnvironment")
})
