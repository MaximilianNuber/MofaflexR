# Assign a modality into a ReticulateMuData

Assigns an AnnData modality into the MuData. `value` may be either an
[`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html)
(its underlying Python object is extracted) or a raw Python
`anndata.AnnData` object. No data is copied.

## Usage

``` r
# S3 method for class 'ReticulateMuData'
x[[i]] <- value
```

## Arguments

- x:

  A [ReticulateMuData](ReticulateMuData.md).

- i:

  A single string: the modality name.

- value:

  An
  [`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html)
  **or** a Python `anndata.AnnData` object.

## Value

`x` (invisibly, for chaining; the Python-side mutation happens
in-place).

## Examples

``` r
if (FALSE) { # \dontrun{
mdata[["rna"]] <- new_rna_reticulateanndata
mdata[["rna"]] <- raw_python_anndata_object
} # }
```
