# ReticulateMuData

An R6 wrapper around a Python **mudata** `MuData` object accessed via
reticulate. All data remain in Python; the R object is a thin proxy that
forwards attribute reads/writes to the underlying Python object.

Modalities (sub-`AnnData` objects) are always returned as
[`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html)
wrappers, so the critical invariant holds: **any `AnnData` extracted
from a `ReticulateMuData` is returned as an
[`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html).**

## Constructor

    mdata <- ReticulateMuData$new(x)

where `x` is either

- a Python `mudata.MuData` object (from `reticulate::import("mudata")`),
  or

- an existing `ReticulateMuData` (wrapping is idempotent).

## Modality access

    rna <- mdata[["rna"]]          # returns ReticulateAnnData
    mdata[["rna"]] <- new_rna_adata  # accepts ReticulateAnnData or raw Python AnnData

## Python object access

    py_obj <- mdata$py_mudata()    # the raw Python MuData

Use
[`py_call()`](https://rstudio.github.io/reticulate/reference/py_call.html)
to invoke Python methods with automatic wrapping of the result:

    result <- mdata$py_call("copy")  # returns ReticulateMuData if result is MuData

## See also

[`is_reticulate_mudata()`](is_reticulate_mudata.md),
[`as_reticulate_mudata()`](as_reticulate_mudata.md),
[`reticulate_mudata()`](reticulate_mudata.md)

## Active bindings

- `obs`:

  Global observation data frame (pandas DataFrame ↔ R data.frame).

- `var`:

  Global variable data frame (pandas DataFrame ↔ R data.frame).

- `obs_names`:

  Global observation names (character vector).

- `var_names`:

  Global variable names (character vector).

- `shape`:

  Integer vector `c(n_obs, n_vars)` (read-only).

- `axis`:

  MuData axis: 0 = shared obs, 1 = shared var (read-only).

- `obsm`:

  Multi-dimensional observation annotation (named R list).

- `varm`:

  Multi-dimensional variable annotation (named R list).

- `obsp`:

  Pairwise observation annotation (named R list of square matrices).

- `varp`:

  Pairwise variable annotation (named R list of square matrices).

- `uns`:

  Unstructured annotation (named R list).

- `mod`:

  Named R list of modalities; each value is a
  [`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html).
  Use `mdata[["name"]] <- adata` to update individual modalities.

## Methods

### Public methods

- [`ReticulateMuData$new()`](#method-ReticulateMuData-new)

- [`ReticulateMuData$py_mudata()`](#method-ReticulateMuData-py_mudata)

- [`ReticulateMuData$modality_names()`](#method-ReticulateMuData-modality_names)

- [`ReticulateMuData$n_obs()`](#method-ReticulateMuData-n_obs)

- [`ReticulateMuData$n_vars()`](#method-ReticulateMuData-n_vars)

- [`ReticulateMuData$obs_keys()`](#method-ReticulateMuData-obs_keys)

- [`ReticulateMuData$var_keys()`](#method-ReticulateMuData-var_keys)

- [`ReticulateMuData$obsm_keys()`](#method-ReticulateMuData-obsm_keys)

- [`ReticulateMuData$varm_keys()`](#method-ReticulateMuData-varm_keys)

- [`ReticulateMuData$uns_keys()`](#method-ReticulateMuData-uns_keys)

- [`ReticulateMuData$update()`](#method-ReticulateMuData-update)

- [`ReticulateMuData$update_obs()`](#method-ReticulateMuData-update_obs)

- [`ReticulateMuData$update_var()`](#method-ReticulateMuData-update_var)

- [`ReticulateMuData$pull_obs()`](#method-ReticulateMuData-pull_obs)

- [`ReticulateMuData$pull_var()`](#method-ReticulateMuData-pull_var)

- [`ReticulateMuData$push_obs()`](#method-ReticulateMuData-push_obs)

- [`ReticulateMuData$push_var()`](#method-ReticulateMuData-push_var)

- [`ReticulateMuData$copy()`](#method-ReticulateMuData-copy)

- [`ReticulateMuData$py_call()`](#method-ReticulateMuData-py_call)

- [`ReticulateMuData$print()`](#method-ReticulateMuData-print)

------------------------------------------------------------------------

### Method [`new()`](https://rdrr.io/r/methods/new.html)

Create a new `ReticulateMuData` wrapper.

#### Usage

    ReticulateMuData$new(x)

#### Arguments

- `x`:

  A Python `mudata.MuData` object (class `"mudata._core.mudata.MuData"`)
  **or** an existing `ReticulateMuData`. Any other type raises an error.

------------------------------------------------------------------------

### Method [`py_mudata()`](py_mudata.md)

Return the underlying Python `mudata.MuData` object.

#### Usage

    ReticulateMuData$py_mudata()

#### Returns

A Python `mudata.MuData` object (`convert = FALSE`).

------------------------------------------------------------------------

### Method `modality_names()`

Return the names of available modalities.

#### Usage

    ReticulateMuData$modality_names()

#### Returns

A character vector of modality names.

------------------------------------------------------------------------

### Method `n_obs()`

Number of observations (shared across modalities).

#### Usage

    ReticulateMuData$n_obs()

#### Returns

A single integer.

------------------------------------------------------------------------

### Method `n_vars()`

Total number of variables (features across all modalities).

#### Usage

    ReticulateMuData$n_vars()

#### Returns

A single integer.

------------------------------------------------------------------------

### Method `obs_keys()`

Keys of global observation annotation (`obs` columns).

#### Usage

    ReticulateMuData$obs_keys()

#### Returns

A character vector.

------------------------------------------------------------------------

### Method `var_keys()`

Keys of global variable annotation (`var` columns).

#### Usage

    ReticulateMuData$var_keys()

#### Returns

A character vector.

------------------------------------------------------------------------

### Method `obsm_keys()`

Keys in `obsm` (multi-dimensional observation annotation).

#### Usage

    ReticulateMuData$obsm_keys()

#### Returns

A character vector.

------------------------------------------------------------------------

### Method `varm_keys()`

Keys in `varm` (multi-dimensional variable annotation).

#### Usage

    ReticulateMuData$varm_keys()

#### Returns

A character vector.

------------------------------------------------------------------------

### Method `uns_keys()`

Keys in `uns` (unstructured annotation).

#### Usage

    ReticulateMuData$uns_keys()

#### Returns

A character vector.

------------------------------------------------------------------------

### Method [`update()`](https://rdrr.io/r/stats/update.html)

Update global `obs` / `var` indices from individual modalities.
Equivalent to Python `mdata.update()`.

#### Usage

    ReticulateMuData$update()

#### Returns

`self` invisibly.

------------------------------------------------------------------------

### Method `update_obs()`

Update global `obs_names` from individual modalities. Equivalent to
Python `mdata.update_obs()`.

#### Usage

    ReticulateMuData$update_obs()

#### Returns

`self` invisibly.

------------------------------------------------------------------------

### Method `update_var()`

Update global `var_names` from individual modalities. Equivalent to
Python `mdata.update_var()`.

#### Usage

    ReticulateMuData$update_var()

#### Returns

`self` invisibly.

------------------------------------------------------------------------

### Method `pull_obs()`

Copy columns from modality-level `.obs` to the global `.obs`. Wraps
Python `mdata.pull_obs(...)`.

#### Usage

    ReticulateMuData$pull_obs(...)

#### Arguments

- `...`:

  Arguments forwarded to the Python method.

#### Returns

`self` invisibly.

------------------------------------------------------------------------

### Method `pull_var()`

Copy columns from modality-level `.var` to the global `.var`. Wraps
Python `mdata.pull_var(...)`.

#### Usage

    ReticulateMuData$pull_var(...)

#### Arguments

- `...`:

  Arguments forwarded to the Python method.

#### Returns

`self` invisibly.

------------------------------------------------------------------------

### Method `push_obs()`

Push columns from global `.obs` back to individual modalities. Wraps
Python `mdata.push_obs(...)`.

#### Usage

    ReticulateMuData$push_obs(...)

#### Arguments

- `...`:

  Arguments forwarded to the Python method.

#### Returns

`self` invisibly.

------------------------------------------------------------------------

### Method `push_var()`

Push columns from global `.var` back to individual modalities. Wraps
Python `mdata.push_var(...)`.

#### Usage

    ReticulateMuData$push_var(...)

#### Arguments

- `...`:

  Arguments forwarded to the Python method.

#### Returns

`self` invisibly.

------------------------------------------------------------------------

### Method `copy()`

Create a deep copy of this MuData.

#### Usage

    ReticulateMuData$copy()

#### Returns

A new ReticulateMuData wrapping the Python copy.

------------------------------------------------------------------------

### Method [`py_call()`](https://rstudio.github.io/reticulate/reference/py_call.html)

Invoke any Python method on the underlying MuData and return the result
with automatic wrapping:

- `anndata.AnnData` →
  [`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html)

- `mudata.MuData` → `ReticulateMuData`

- anything else → returned as-is (Python object, `convert = FALSE`)

#### Usage

    ReticulateMuData$py_call(method, ...)

#### Arguments

- `method`:

  A string: the Python method name.

- `...`:

  Arguments forwarded to the Python method.

#### Returns

The (possibly wrapped) result.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print a human-readable summary.

#### Usage

    ReticulateMuData$print(...)

#### Arguments

- `...`:

  Ignored.

## Examples

``` r
if (FALSE) { # \dontrun{
reticulate::use_condaenv("mofaflex_test", required = TRUE)
mu   <- reticulate::import("mudata",  convert = FALSE)
ad   <- reticulate::import("anndata", convert = FALSE)
sp   <- reticulate::import("scipy.sparse", convert = FALSE)
np   <- reticulate::import("numpy",   convert = FALSE)

X1 <- sp$csr_matrix(reticulate::r_to_py(matrix(1:12, 3L, 4L)))
X2 <- sp$csr_matrix(reticulate::r_to_py(matrix(1:6,  3L, 2L)))
a1  <- ad$AnnData(X = X1)
a2  <- ad$AnnData(X = X2)
md  <- mu$MuData(reticulate::py_dict(list("rna" = a1, "prot" = a2)))

mdata <- ReticulateMuData$new(md)
print(mdata)

rna   <- mdata[["rna"]]   # ReticulateAnnData
print(rna)
} # }
```
