# Convert a SingleCellExperiment assay to a scipy CSC sparse matrix

Extracts an assay from a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object, coerces it to a
[Matrix::dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html),
and constructs a Python `scipy.sparse.csc_matrix` by passing
[`reticulate::np_array()`](https://rstudio.github.io/reticulate/reference/np_array.html)
views of the underlying R memory (true zero-copy:
`flags.owndata = FALSE`).

### Slot mapping

|                 |                    |           |
|-----------------|--------------------|-----------|
| R (`dgCMatrix`) | Python / scipy CSC | dtype     |
| `x@x`           | `data`             | `float64` |
| `x@i`           | `indices`          | `int32`   |
| `x@p`           | `indptr`           | `int32`   |
| `x@Dim`         | `shape`            | –         |

### Zero-copy semantics

[`reticulate::np_array()`](https://rstudio.github.io/reticulate/reference/np_array.html)
uses the C buffer protocol to expose R's contiguous memory to NumPy
**without copying** (`flags.owndata = FALSE`). This works because plain
R vectors are already C-contiguous:

- R numeric (`double`) maps directly to NumPy `float64`.

- R integer (`int32_t`) maps directly to NumPy `int32`.

**Note on `np.array(..., copy=False)`**: NumPy \>= 2.0 changed this flag
to raise `ValueError` when a copy would be needed. Using
[`reticulate::np_array()`](https://rstudio.github.io/reticulate/reference/np_array.html)
avoids this entirely and is the canonical reticulate API for zero-copy
array creation from R.

**Lifetime caveat**: the returned Python object holds a borrowed
reference to the underlying R vectors. Do *not* modify, reallocate, or
allow the `mat@x`, `mat@i`, `mat@p` vectors to be garbage-collected
while the Python object is alive.

### Why CSC?

`dgCMatrix` stores columns contiguously (Compressed Sparse Column). We
map directly to `csc_matrix` without reordering values. Building a CSR
matrix would require transposing the data (= always a copy).

### int32 for indptr

R integer vectors are 32-bit (`int32_t`). scipy CSC accepts `int32` for
both `indices` and `indptr`, supporting up to \\2^{31}-1\\ non-zero
entries per matrix.

## Usage

``` r
sce_assay_to_scipy_csc(x, assay = "counts")
```

## Arguments

- x:

  A
  [SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).

- assay:

  A single string: the name of the assay to extract (default
  `"counts"`).

## Value

A Python `scipy.sparse.csc_matrix` (a reticulate Python object; not
auto-converted to R).

## See also

[`sce_to_anndata()`](sce_to_anndata.md),
[`sce_to_reticulate_anndata()`](sce_to_reticulate_anndata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(SingleCellExperiment)
library(Matrix)
counts <- rsparsematrix(200, 50, density = 0.1, repr = "C")
sce <- SingleCellExperiment(assays = list(counts = counts))
csc <- sce_assay_to_scipy_csc(sce, assay = "counts")
sp <- reticulate::import("scipy.sparse", convert = FALSE)
shape_r <- reticulate::py_to_r(csc$shape)
cat("Shape:", shape_r[[1]], "x", shape_r[[2]], "\n")
cat("Sparse:", reticulate::py_to_r(sp$issparse(csc)), "\n")
cat("owndata:", reticulate::py_to_r(csc$data$flags$owndata), "\n")
} # }
```
