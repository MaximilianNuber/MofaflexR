# tests/testthat/test-ReticulateMuData.R
# Unit tests for the ReticulateMuData R6 wrapper.

# ---------------------------------------------------------------------------
# Helper: build a minimal MuData in Python and return the Python object.
# Called inside individual tests so that reticulate is already initialised.
# ---------------------------------------------------------------------------

.make_py_mudata <- function() {
  mu <- reticulate::import("mudata",        convert = FALSE)
  ad <- reticulate::import("anndata",       convert = FALSE)
  sp <- reticulate::import("scipy.sparse",  convert = FALSE)
  np <- reticulate::import("numpy",         convert = FALSE)

  # 3 obs × 4 vars (rna), 3 obs × 2 vars (prot)
  X1 <- sp$csr_matrix(reticulate::r_to_py(matrix(1:12, 3L, 4L)))
  X2 <- sp$csr_matrix(reticulate::r_to_py(matrix(1:6,  3L, 2L)))

  obs_names <- reticulate::r_to_py(c("obs1", "obs2", "obs3"))

  a1 <- ad$AnnData(X = X1)
  a2 <- ad$AnnData(X = X2)
  a1$obs_names <- obs_names
  a2$obs_names <- obs_names

  py_dict <- reticulate::py_dict(
    keys   = list("rna", "prot"),
    values = list(a1, a2),
    convert = FALSE
  )
  mu$MuData(py_dict)
}


# ===========================================================================
# Construction tests
# ===========================================================================

test_that("ReticulateMuData construction from Python MuData works", {
  skip_if_no_mudata()

  py_md <- .make_py_mudata()
  mdata <- ReticulateMuData$new(py_md)

  expect_s3_class(mdata, "ReticulateMuData")
})

test_that("ReticulateMuData construction is idempotent (accepts existing wrapper)", {
  skip_if_no_mudata()

  py_md  <- .make_py_mudata()
  mdata1 <- ReticulateMuData$new(py_md)
  mdata2 <- ReticulateMuData$new(mdata1)  # wrap of wrap

  expect_s3_class(mdata2, "ReticulateMuData")
  # Both should reference the same underlying Python object
  expect_identical(mdata1$py_mudata(), mdata2$py_mudata())
})

test_that("ReticulateMuData rejects invalid input", {
  expect_error(ReticulateMuData$new(42L),         "mudata.MuData")
  expect_error(ReticulateMuData$new("not-mudata"), "mudata.MuData")
  expect_error(ReticulateMuData$new(list()),       "mudata.MuData")
})


# ===========================================================================
# Helper-function tests
# ===========================================================================

test_that("is_reticulate_mudata returns correct logical", {
  skip_if_no_mudata()

  py_md <- .make_py_mudata()
  mdata <- ReticulateMuData$new(py_md)

  expect_true(is_reticulate_mudata(mdata))
  expect_false(is_reticulate_mudata(42L))
  expect_false(is_reticulate_mudata(list()))
  expect_false(is_reticulate_mudata(NULL))
})

test_that("as_reticulate_mudata coerces Python MuData", {
  skip_if_no_mudata()

  py_md <- .make_py_mudata()
  mdata <- as_reticulate_mudata(py_md)
  expect_s3_class(mdata, "ReticulateMuData")
})

test_that("as_reticulate_mudata is idempotent on ReticulateMuData", {
  skip_if_no_mudata()

  py_md  <- .make_py_mudata()
  mdata  <- ReticulateMuData$new(py_md)
  mdata2 <- as_reticulate_mudata(mdata)

  expect_identical(mdata, mdata2)
})

test_that("reticulate_mudata() constructor alias works", {
  skip_if_no_mudata()

  py_md <- .make_py_mudata()
  mdata <- reticulate_mudata(py_md)
  expect_s3_class(mdata, "ReticulateMuData")
})

test_that("py_mudata() module import returns Python mudata module", {
  skip_if_no_mudata()

  mu <- py_mudata()
  # The returned object should be a Python module, not an R object
  expect_true(inherits(mu, "python.builtin.module"))
})


# ===========================================================================
# Accessor tests
# ===========================================================================

test_that("modality_names() returns character vector of mod keys", {
  skip_if_no_mudata()

  mdata <- ReticulateMuData$new(.make_py_mudata())
  nms   <- mdata$modality_names()

  expect_type(nms, "character")
  expect_setequal(nms, c("rna", "prot"))
})

test_that("n_obs() returns correct count", {
  skip_if_no_mudata()

  mdata <- ReticulateMuData$new(.make_py_mudata())
  expect_equal(mdata$n_obs(), 3L)
})

test_that("obs_names active binding returns character vector of length n_obs", {
  skip_if_no_mudata()

  mdata <- ReticulateMuData$new(.make_py_mudata())
  on    <- mdata$obs_names

  expect_type(on, "character")
  expect_length(on, 3L)
  expect_equal(on, c("obs1", "obs2", "obs3"))
})

test_that("py_mudata() returns Python object with correct class", {
  skip_if_no_mudata()

  mdata  <- ReticulateMuData$new(.make_py_mudata())
  py_obj <- mdata$py_mudata()

  expect_true(inherits(py_obj, "mudata._core.mudata.MuData"))
})

test_that("print() runs without error and returns mdata invisibly", {
  skip_if_no_mudata()

  mdata <- ReticulateMuData$new(.make_py_mudata())
  out   <- capture.output(result <- mdata$print())

  expect_s3_class(result, "ReticulateMuData")
  expect_true(any(grepl("ReticulateMuData", out)))
  expect_true(any(grepl("rna", out)))
  expect_true(any(grepl("prot", out)))
})


# ===========================================================================
# [[ accessor – critical invariant: always returns ReticulateAnnData
# ===========================================================================

test_that("[[ returns ReticulateAnnData for each modality", {
  skip_if_no_mudata()

  mdata <- ReticulateMuData$new(.make_py_mudata())

  rna_r  <- mdata[["rna"]]
  prot_r <- mdata[["prot"]]

  # Critical invariant: must be ReticulateAnnData (anndataR R6 class)
  expect_s3_class(rna_r,  "ReticulateAnnData")
  expect_s3_class(prot_r, "ReticulateAnnData")
})

test_that("[[ returns ReticulateAnnData with correct dimensions", {
  skip_if_no_mudata()

  mdata <- ReticulateMuData$new(.make_py_mudata())

  rna  <- mdata[["rna"]]
  prot <- mdata[["prot"]]

  # rna: 3 obs × 4 vars; prot: 3 obs × 2 vars
  expect_equal(rna$n_obs(),   3L)
  expect_equal(rna$n_vars(),  4L)
  expect_equal(prot$n_obs(),  3L)
  expect_equal(prot$n_vars(), 2L)
})

test_that("[[ result is NOT auto-converted to R list or data.frame", {
  skip_if_no_mudata()

  mdata   <- ReticulateMuData$new(.make_py_mudata())
  rna_r   <- mdata[["rna"]]
  py_rna  <- rna_r$py_anndata()

  # The underlying Python object must still be a Python object, not
  # auto-converted to an R list.
  expect_true(inherits(py_rna, "python.builtin.object"))
  expect_false(is.list(py_rna))
})


# ===========================================================================
# [[<- setter tests
# ===========================================================================

test_that("[[<- works with ReticulateAnnData input", {
  skip_if_no_mudata()

  mdata <- ReticulateMuData$new(.make_py_mudata())

  # Build a new small AnnData (3 obs × 5 vars) and wrap it
  ad   <- reticulate::import("anndata",      convert = FALSE)
  sp   <- reticulate::import("scipy.sparse", convert = FALSE)
  X_new <- sp$csr_matrix(reticulate::r_to_py(matrix(1:15, 3L, 5L)))
  new_py_ad <- ad$AnnData(
    X = X_new,
    obs = reticulate::r_to_py(data.frame(row.names = c("obs1","obs2","obs3")))
  )
  RAD <- get("ReticulateAnnData", envir = asNamespace("anndataR"), inherits = FALSE)
  new_rada <- RAD$new(py_anndata = new_py_ad)

  # Assign new modality via wrapper
  mdata[["rna"]] <- new_rada

  # Verify update through [[
  updated <- mdata[["rna"]]
  expect_s3_class(updated, "ReticulateAnnData")
  expect_equal(updated$n_vars(), 5L)
})

test_that("[[<- works with raw Python AnnData input", {
  skip_if_no_mudata()

  mdata <- ReticulateMuData$new(.make_py_mudata())

  ad   <- reticulate::import("anndata",      convert = FALSE)
  sp   <- reticulate::import("scipy.sparse", convert = FALSE)
  X_new <- sp$csr_matrix(reticulate::r_to_py(matrix(1:9, 3L, 3L)))
  new_py_ad <- ad$AnnData(
    X = X_new,
    obs = reticulate::r_to_py(data.frame(row.names = c("obs1","obs2","obs3")))
  )

  # Assign raw Python AnnData
  mdata[["prot"]] <- new_py_ad

  updated <- mdata[["prot"]]
  expect_s3_class(updated, "ReticulateAnnData")
  expect_equal(updated$n_vars(), 3L)
})


# ===========================================================================
# py_call() – centralised Python-return interception
# ===========================================================================

test_that("py_call() wraps MuData return as ReticulateMuData", {
  skip_if_no_mudata()

  mdata  <- ReticulateMuData$new(.make_py_mudata())
  result <- mdata$py_call("copy")

  # copy() returns a new MuData -> should be wrapped
  expect_s3_class(result, "ReticulateMuData")
  expect_equal(result$n_obs(), 3L)
  expect_setequal(result$modality_names(), c("rna", "prot"))
})

test_that("py_call() wraps AnnData return as ReticulateAnnData", {
  skip_if_no_mudata()

  # MuData.__getitem__("rna") returns an AnnData
  mdata  <- ReticulateMuData$new(.make_py_mudata())
  result <- mdata$py_call("__getitem__", reticulate::r_to_py("rna"))

  expect_s3_class(result, "ReticulateAnnData")
  expect_equal(result$n_vars(), 4L)
})

test_that("py_call() returns non-AnnData/MuData Python objects unchanged", {
  skip_if_no_mudata()

  mdata  <- ReticulateMuData$new(.make_py_mudata())
  # n_obs is an attribute; __len__ is a method call -> returns integer
  result <- mdata$py_call("__len__")

  expect_false(inherits(result, "ReticulateMuData"))
  expect_false(inherits(result, "ReticulateAnnData"))
})
