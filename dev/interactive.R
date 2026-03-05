# dev/interactive.R
# Interactive exploration script for MofaflexR.
# Source this file or run sections manually in an R session.
#
# Prerequisites:
#   - MofaflexR installed (R CMD INSTALL .)
#   - Python env with numpy / scipy / anndata available
#     (set RETICULATE_PYTHON or run via:
#      RETICULATE_PYTHON=/home/maximiliann/miniforge3/envs/mofaflex_test/bin/python Rscript dev/interactive.R)

# ---------------------------------------------------------------------------
# 0. Load packages
# ---------------------------------------------------------------------------
library(MofaflexR)
library(SingleCellExperiment)
library(Matrix)
library(reticulate)

# ---------------------------------------------------------------------------
# 1. Build test data
# ---------------------------------------------------------------------------

set.seed(42L)

N_GENES  <- 200L   # rows  (features)
N_CELLS  <- 80L    # cols  (observations)
DENSITY  <- 0.08   # ~8 % non-zero

counts <- rsparsematrix(N_GENES, N_CELLS, density = DENSITY, repr = "C")
counts <- abs(round(counts))                       # non-negative integers
counts <- methods::as(counts, "dgCMatrix")

# Cell metadata
cell_meta <- data.frame(
  cell_type  = sample(c("T_cell", "B_cell", "NK"), N_CELLS, replace = TRUE),
  condition  = sample(c("ctrl", "treat"),           N_CELLS, replace = TRUE),
  total_umi  = Matrix::colSums(counts),
  row.names  = paste0("cell", seq_len(N_CELLS))
)

# Gene metadata
gene_meta <- data.frame(
  biotype        = sample(c("protein_coding", "lncRNA", "pseudogene"),
                          N_GENES, replace = TRUE),
  mean_expression = Matrix::rowMeans(counts),
  row.names       = paste0("gene", seq_len(N_GENES))
)

sce <- SingleCellExperiment(
  assays  = list(counts = counts),
  colData = S4Vectors::DataFrame(cell_meta),
  rowData = S4Vectors::DataFrame(gene_meta)
)

cat("--- SCE summary ---\n")
print(sce)
cat("Assay dims (genes x cells):", dim(counts), "\n")
cat("nnz:", length(counts@x), "\n\n")

# ---------------------------------------------------------------------------
# 2. sce_assay_to_scipy_csc  —  raw CSC (genes x cells)
# ---------------------------------------------------------------------------

cat("--- sce_assay_to_scipy_csc ---\n")
csc <- sce_assay_to_scipy_csc(sce, assay = "counts")

sp       <- import("scipy.sparse", convert = FALSE)
shape_r  <- py_to_r(csc$shape)

cat("Python type  :", class(csc)[[1L]], "\n")
cat("Shape        :", shape_r[[1L]], "x", shape_r[[2L]], "\n")
cat("nnz          :", as.integer(py_to_r(csc$nnz)), "\n")
cat("issparse     :", py_to_r(sp$issparse(csc)), "\n")
cat("data.owndata :", py_to_r(csc$data$flags$owndata),  "  <- FALSE = zero-copy\n")
cat("data.dtype   :", py_to_r(csc$data$dtype$name), "\n")
cat("indices dtype:", py_to_r(csc$indices$dtype$name), "\n")
cat("indptr dtype :", py_to_r(csc$indptr$dtype$name), "\n\n")

# Verify slot lengths match R
stopifnot(as.integer(py_to_r(py_len(csc$data)))    == length(counts@x))
stopifnot(as.integer(py_to_r(py_len(csc$indices))) == length(counts@i))
stopifnot(as.integer(py_to_r(py_len(csc$indptr)))  == length(counts@p))
cat("Slot length checks: OK\n\n")

# ---------------------------------------------------------------------------
# 3. sce_to_anndata  —  Python AnnData (cells x genes)
# ---------------------------------------------------------------------------

cat("--- sce_to_anndata ---\n")
adata    <- sce_to_anndata(sce, assay = "counts")
shape_ad <- py_to_r(adata$shape)

cat("Python type  :", class(adata)[[1L]], "\n")
cat("Shape (obs x vars):", shape_ad[[1L]], "x", shape_ad[[2L]], "\n")
cat("X type       :", class(adata$X)[[1L]], "\n")
cat("X issparse   :", py_to_r(sp$issparse(adata$X)), "\n")
cat("X owndata    :", py_to_r(adata$X$data$flags$owndata), "\n")
cat("obs_names[1:3]:", paste(py_to_r(adata$obs_names$tolist())[1:3], collapse=", "), "\n")
cat("var_names[1:3]:", paste(py_to_r(adata$var_names$tolist())[1:3], collapse=", "), "\n\n")

# ---------------------------------------------------------------------------
# 4. sce_to_reticulate_anndata  —  ReticulateAnnData
# ---------------------------------------------------------------------------

cat("--- sce_to_reticulate_anndata ---\n")
rada <- sce_to_reticulate_anndata(sce, assay = "counts")
print(rada)

py_obj <- rada$py_anndata()
shape_py <- py_to_r(py_obj$shape)
cat("Underlying Python AnnData shape:", shape_py[[1L]], "x", shape_py[[2L]], "\n")
cat("X issparse:", py_to_r(sp$issparse(py_obj$X)), "\n")
cat("X owndata :", py_to_r(py_obj$X$data$flags$owndata), "\n\n")

# Access obs/var metadata via ReticulateAnnData
cat("First 3 rows of rada$obs:\n")
print(head(rada$obs, 3))

cat("\nFirst 3 rows of rada$var:\n")
print(head(rada$var, 3))

# ---------------------------------------------------------------------------
# 5. Optional: extra assay as a layer
# ---------------------------------------------------------------------------

cat("\n--- sce_to_anndata with extra layer ---\n")
# For layers, pass a plain base R matrix (cells x genes) so reticulate
# converts it to a numpy array without issue.
logcounts_mat <- base::t(base::log1p(base::as.matrix(counts)))  # (cells x genes)

adata2 <- sce_to_anndata(
  sce,
  assay  = "counts",
  layers = list(logcounts = logcounts_mat)
)

# Check the layer was stored using Python's __contains__
layer_in    <- reticulate::py_to_r(adata2$layers$`__contains__`("logcounts"))
cat("'logcounts' layer present:", layer_in, "\n")
layer_shape <- reticulate::py_to_r(adata2$layers[["logcounts"]]$shape)
cat("logcounts layer shape    :", layer_shape[[1L]], "x", layer_shape[[2L]], "\n")

cat("\nAll done. No errors.\n")
