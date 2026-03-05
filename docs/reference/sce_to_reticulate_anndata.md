# Convert a SingleCellExperiment to an anndataR ReticulateAnnData

Converts a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
to an
[`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html)
object backed by a Python `anndata.AnnData`. The main matrix (`X`) is a
`scipy.sparse.csc_matrix` built without densifying the data; see
[`sce_assay_to_scipy_csc()`](sce_assay_to_scipy_csc.md) for zero-copy
semantics and lifetime caveats.

The returned object behaves as an `anndataR` AnnData: it can be passed
to scanpy functions, written to HDF5/`.h5ad`, or inspected with the
standard `$` accessor.

## Usage

``` r
sce_to_reticulate_anndata(x, assay = "counts", ...)
```

## Arguments

- x:

  A
  [SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).

- assay:

  A single string: the assay to place in `X` (default `"counts"`).

- ...:

  Additional arguments forwarded to
  [`sce_to_anndata()`](sce_to_anndata.md).

## Value

An
[`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html)
R6 object backed by a Python `anndata.AnnData`.

## See also

[`sce_assay_to_scipy_csc()`](sce_assay_to_scipy_csc.md),
[`sce_to_anndata()`](sce_to_anndata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(SingleCellExperiment)
library(Matrix)

counts <- rsparsematrix(200, 50, density = 0.1, repr = "C")
sce <- SingleCellExperiment(assays = list(counts = counts))

rada <- sce_to_reticulate_anndata(sce)
print(rada)      # shows n_obs x n_vars

# Access the underlying Python object
py_adata <- rada$py_anndata()
cat("Shape:", reticulate::py_str(py_adata$shape), "\n")
cat("X is sparse:", reticulate::py_bool(
  reticulate::import("scipy.sparse")$issparse(py_adata$X)
), "\n")
} # }
```
