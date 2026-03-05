# MofaflexR

**Fit MOFAFlex models from R via basilisk + reticulate.**

MofaflexR provides an R-native bridge between
[SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment/)
objects and the Python [MOFAFlex](https://github.com/bioFAM/mofaflex)
package. The key contribution is a **zero-copy sparse bridge**: count
matrices stored as `Matrix::dgCMatrix` are exposed to Python
`scipy.sparse` without densifying or allocating a second copy of the
data.

The Python environment (numpy, scipy, anndata, mofaflex) is managed
automatically by
[basilisk](https://bioconductor.org/packages/basilisk/); users do
**not** need to install Python manually.

------------------------------------------------------------------------

## Installation

``` r
# Bioconductor + CRAN dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "basilisk", "basilisk.utils",
                        "anndataR"))

# Install MofaflexR from source (once on CRAN/Bioconductor, use that)
# remotes::install_github("yourorg/MofaflexR")
```

------------------------------------------------------------------------

## Quick start

``` r
library(MofaflexR)
library(SingleCellExperiment)
library(Matrix)

# ----- 1. Build a toy SCE -----
set.seed(42)
counts <- rsparsematrix(200, 50, density = 0.1, repr = "C")
counts <- abs(round(counts))   # non-negative integer counts
sce <- SingleCellExperiment(assays = list(counts = counts))
colnames(sce) <- paste0("cell", seq_len(ncol(sce)))
rownames(sce) <- paste0("gene", seq_len(nrow(sce)))

# ----- 2. Convert to ReticulateAnnData (zero-copy sparse bridge) -----
rada <- sce_to_reticulate_anndata(sce)
print(rada)
#> ReticulateAnnData object with n_obs × n_vars = 50 × 200

# ----- 3. Access the Python AnnData directly -----
py_adata <- rada$py_anndata()

sp <- reticulate::import("scipy.sparse", convert = FALSE)
cat("Shape:   ", reticulate::py_to_r(py_adata$shape)[[1]], "x",
    reticulate::py_to_r(py_adata$shape)[[2]], "\n")
cat("Sparse:  ", reticulate::py_to_r(sp$issparse(py_adata$X)), "\n")
cat("owndata: ", reticulate::py_to_r(py_adata$X$data$flags$owndata), "\n")
# Shape:    50 x 200
# Sparse:   TRUE
# owndata:  FALSE   ← no copy of the data array
```

### Low-level usage: CSC matrix only

``` r
# Returns a Python scipy.sparse.csc_matrix (shape = genes × cells)
csc <- sce_assay_to_scipy_csc(sce, assay = "counts")

shape_r <- reticulate::py_to_r(csc$shape)
cat("Shape:", shape_r[[1]], "x", shape_r[[2]], "\n")  # 200 x 50
cat("nnz:  ", reticulate::py_to_r(csc$nnz), "\n")
```

### Using the basilisk-managed environment

``` r
result <- with_mofaflex_env({
  mf <- reticulate::import("mofaflex", convert = FALSE)
  # ... fit model, return R-friendly result
})
```

------------------------------------------------------------------------

## Fit a MOFA-FLEX model

[`fit_mofaflex()`](reference/fit_mofaflex.md) converts an SCE (or an
existing `ReticulateAnnData` / `ReticulateMuData`) to the required
Python data structure, builds and trains a MOFA-FLEX model, and
optionally saves the result in MOFA2-compatible HDF5 format.

``` r
library(MofaflexR)
library(muscData)   # BiocManager::install("muscData")
library(scran);  library(scuttle)

# Build SCE and normalise (see vignette for full QC)
sce <- Kang18_8vs8()
sce <- sce[, sce$multiplets == "singlet" & !is.na(sce$cell)]
sce <- scuttle::logNormCounts(sce)
hvgs <- scran::getTopHVGs(scran::modelGeneVar(sce), n = 3000)
sce_hvg <- sce[hvgs, ]

# Zero-copy bridge to Python AnnData
adata <- sce_to_reticulate_anndata(sce_hvg, assay = "logcounts")

# Train MOFA-FLEX with two groups (ctrl / stim)
# Passing a single AnnData with group_by automatically wraps it in MuData.
model_path <- tempfile(fileext = ".hdf5")
model <- fit_mofaflex(
  data             = adata,
  data_options     = list(group_by           = "stim",
                          scale_per_group    = TRUE,
                          subset_var         = NULL,
                          plot_data_overview = FALSE),
  model_options    = list(n_factors    = 10L,
                          weight_prior = "Horseshoe"),
  training_options = list(max_epochs             = 5000L,
                          early_stopper_patience = 50L,
                          seed                   = 42L,
                          device                 = "cpu",
                          save_path              = model_path),
  mofa_compat      = "full"   # MOFA2-readable HDF5
)

# Load into MOFA2 for downstream analysis
mofa <- load_mofaflex_model(model_path)
MOFA2::plot_variance_explained(mofa)
```

> **Tip.** If `group_by` is specified and `data` is a single `AnnData`,
> [`fit_mofaflex()`](reference/fit_mofaflex.md) automatically wraps it
> in a `mudata.MuData` and promotes the grouping column to the top-level
> `obs`. The wrapped modality is named **`rna`** automatically.

------------------------------------------------------------------------

## ReticulateMuData: wrapping Python MuData objects

[mudata](https://mudata.readthedocs.io/) is the Python multi-modal data
container used by the scverse ecosystem. `ReticulateMuData` is an R6
class that wraps a Python `mudata.MuData` object, mirroring the design
of
[`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html).

### Critical invariant

> **Any `AnnData` extracted from a `ReticulateMuData` is returned as an
> [`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html).**

This means scanpy/muon functions receiving `rna$py_anndata()` work
without any extra conversion step.

### Quick example

``` r
library(MofaflexR)
library(reticulate)

# Import Python modules (convert = FALSE keeps objects Python-native)
mu  <- py_mudata()                             # mudata module
ad  <- import("anndata",       convert = FALSE)
sp  <- import("scipy.sparse",  convert = FALSE)

# Build two AnnData objects
X1  <- sp$csr_matrix(r_to_py(matrix(1:12, 3L, 4L)))
X2  <- sp$csr_matrix(r_to_py(matrix(1:6,  3L, 2L)))
a1  <- ad$AnnData(X = X1)
a2  <- ad$AnnData(X = X2)
a1$obs_names <- r_to_py(c("obs1", "obs2", "obs3"))
a2$obs_names <- r_to_py(c("obs1", "obs2", "obs3"))

# Combine into MuData
py_md <- mu$MuData(py_dict(list("rna" = a1, "prot" = a2), convert = FALSE))

# Wrap in ReticulateMuData
mdata <- ReticulateMuData$new(py_md)
print(mdata)
# ReticulateMuData
#   n_obs     : 3
#   modalities: rna, prot
#     rna:           3 obs x 4 vars
#     prot:          3 obs x 2 vars

# Access a modality – returns ReticulateAnnData automatically
rna <- mdata[["rna"]]
inherits(rna, "ReticulateAnnData")   # TRUE
rna$n_vars()                         # 4L

# Assign a new modality (no data copy)
mdata[["rna"]] <- new_rna_anndata   # accepts ReticulateAnnData or Python AnnData

# Call a Python method with automatic return-type wrapping
mdata_copy <- mdata$py_call("copy")
inherits(mdata_copy, "ReticulateMuData")  # TRUE

# Access the raw Python object for muon/scanpy functions
py_obj <- mdata$py_mudata()
# muon.tl.mofa(py_obj, ...)
```

### ReticulateMuData API reference

| Method / function                       | Description                                                     |
|-----------------------------------------|-----------------------------------------------------------------|
| `ReticulateMuData$new(x)`               | Wrap a Python `MuData` or existing wrapper                      |
| `mdata[["mod"]]`                        | Return modality as `ReticulateAnnData`                          |
| `mdata[["mod"]] <- adata`               | Assign modality (accepts `ReticulateAnnData` or Python AnnData) |
| `mdata$py_mudata()`                     | Raw Python `mudata.MuData` object                               |
| `mdata$modality_names()`                | Character vector of modality names                              |
| `mdata$n_obs()`                         | Number of observations                                          |
| `mdata$obs_names()`                     | Observation names                                               |
| `mdata$py_call(method, ...)`            | Call Python method with wrapped return                          |
| `is_reticulate_mudata(x)`               | Logical test                                                    |
| `as_reticulate_mudata(x)`               | Coerce Python MuData or existing wrapper                        |
| `reticulate_mudata(py_obj)`             | Constructor alias                                               |
| [`py_mudata()`](reference/py_mudata.md) | Import the Python `mudata` module (`convert = FALSE`)           |

------------------------------------------------------------------------

## Zero-copy sparse bridge: how it works and what can force copies

### How it works

`dgCMatrix` stores its non-zero data in three contiguous C arrays:

| Slot | Content               | R type    | NumPy dtype |
|------|-----------------------|-----------|-------------|
| `@x` | non-zero values       | `double`  | `float64`   |
| `@i` | row indices (0-based) | `integer` | `int32`     |
| `@p` | column pointers       | `integer` | `int32`     |

`reticulate::np_array(<R vector>, dtype = ...)` exposes the underlying R
memory to NumPy via the **C buffer protocol** — no data allocation
occurs (`numpy.ndarray.flags.owndata == False`).

These three arrays are then passed directly to
`scipy.sparse.csc_matrix((data, indices, indptr), shape=...)`.

`scipy` itself does not copy the buffers at construction time; it wraps
them by reference.

The resulting `csc_matrix` (shape `genes × cells`) is then transposed
via `.T` to obtain a `csr_matrix` (shape `cells × genes`) suitable for
AnnData’s `obs × vars` layout. The `scipy` transpose is also
**zero-copy**: it reuses the same C arrays, only changing strides and
shape metadata.

### What CAN force a copy

| Trigger                                                                                                            | Why                                                    |
|--------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------|
| Coercion from a non-`dgCMatrix` assay                                                                              | `as(mat, "dgCMatrix")` may allocate a new object       |
| `adata$X = some_dense_array`                                                                                       | AnnData will densify if you assign a dense numpy array |
| Reading `adata$X` back to R via [`py_to_r()`](https://rstudio.github.io/reticulate/reference/r-py-conversion.html) | reticulate converts to a dense R matrix                |
| `adata.X.toarray()` or `adata.X.todense()` in Python                                                               | Explicit densification                                 |
| `.copy()` on any sparse matrix                                                                                     | Explicit copy                                          |
| Operations that require sorted indices (e.g. `adata.X.sort_indices()`)                                             | Creates new arrays                                     |

### Lifetime caveat

The Python `csc_matrix` holds a **borrowed reference** to the R vectors
`mat@x`, `mat@i`, and `mat@p`. You must ensure that:

1.  The `dgCMatrix` object (and therefore its slots) remains alive while
    Python is using the matrix.
2.  No R code modifies the vectors in-place (copy-on-modify semantics
    usually prevent this for normal R code).
3.  The R garbage collector has not freed the underlying memory.

In practice, as long as the `MofaflexR` call site holds a reference to
`sce` (or the extracted `dgCMatrix`) within the same R session, the
memory is safe.

------------------------------------------------------------------------

## API reference

| Function                                                                              | Description                                                                                          |
|---------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|
| `sce_assay_to_scipy_csc(x, assay)`                                                    | Extract assay → `scipy.sparse.csc_matrix` (genes × cells, zero-copy)                                 |
| `sce_to_anndata(x, assay, obs, var, ...)`                                             | Full SCE → Python `anndata.AnnData` (cells × genes)                                                  |
| `sce_to_reticulate_anndata(x, assay, ...)`                                            | SCE → [`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html) |
| `fit_mofaflex(data, data_options, model_options, training_options, mofa_compat, ...)` | Fit a MOFA-FLEX model; plain `AnnData` is auto-wrapped in `MuData`                                   |
| `load_mofaflex_model(path, ...)`                                                      | Load a mofaflex HDF5 into MOFA2, patching HDF5 format differences                                    |
| `with_mofaflex_env(expr)`                                                             | Evaluate `expr` inside the managed Python env                                                        |
| `ReticulateMuData$new(x)`                                                             | Wrap Python `mudata.MuData` → `ReticulateMuData`                                                     |
| `is_reticulate_mudata(x)`                                                             | Logical: is `x` a `ReticulateMuData`?                                                                |
| `as_reticulate_mudata(x)`                                                             | Coerce Python MuData or existing wrapper                                                             |
| `reticulate_mudata(py_obj)`                                                           | Constructor alias for `ReticulateMuData$new`                                                         |
| [`py_mudata()`](reference/py_mudata.md)                                               | Import Python `mudata` module (`convert = FALSE`)                                                    |

------------------------------------------------------------------------

## Session info / reproducibility

``` r
sessionInfo()
reticulate::py_config()
```
