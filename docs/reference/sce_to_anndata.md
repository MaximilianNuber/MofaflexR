# Convert a SingleCellExperiment to a Python AnnData object

Creates a Python `anndata.AnnData` object from a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
using a `scipy.sparse.csc_matrix` built with
[`sce_assay_to_scipy_csc()`](sce_assay_to_scipy_csc.md).

### Copy semantics

- **X** (the main count matrix): the dgCMatrix is first converted to a
  `scipy.sparse.csc_matrix` (shape `genes × cells`) via
  [`sce_assay_to_scipy_csc()`](sce_assay_to_scipy_csc.md), then
  transposed to `scipy.sparse.csr_matrix` (shape `cells × genes`) by
  `.T`. The transpose shares all data buffers with the original CSC
  (`flags.owndata = FALSE`), so **no values are copied**.

- **obs / var metadata**: a pandas `DataFrame` is constructed from the R
  `data.frame`; reticulate performs one copy here, which is unavoidable
  for columnar R data frames → row-oriented Python objects.

- **layers / obsm / varm**: passed through as-is; the caller is
  responsible for any copy semantics of those objects.

## Usage

``` r
sce_to_anndata(
  x,
  assay = "counts",
  obs = NULL,
  var = NULL,
  layers = list(),
  obsm = list(),
  varm = list()
)
```

## Arguments

- x:

  A
  [SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).

- assay:

  A single string: the assay to place in `X` (default `"counts"`).

- obs:

  `NULL` or a `data.frame` of cell-level metadata. Defaults to
  `as.data.frame(colData(x))` with `rownames` set to `colnames(x)`.

- var:

  `NULL` or a `data.frame` of gene-level metadata. Defaults to
  `as.data.frame(rowData(x))` with `rownames` set to `rownames(x)`.

- layers:

  A named list of additional assay matrices to place in
  `AnnData$layers`.

- obsm:

  A named list of matrices added to `AnnData$obsm`.

- varm:

  A named list of matrices added to `AnnData$varm`.

## Value

A Python `anndata.AnnData` object (reticulate, `convert = FALSE`).

## See also

[`sce_assay_to_scipy_csc()`](sce_assay_to_scipy_csc.md),
[`sce_to_reticulate_anndata()`](sce_to_reticulate_anndata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = Matrix::rsparsematrix(200, 50, 0.1, repr = "C"))
)
adata <- sce_to_anndata(sce)
cat("obs x var:", reticulate::py_str(adata$shape), "\n")
} # }
```
