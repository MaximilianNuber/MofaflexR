# Constructor alias for ReticulateMuData

Thin alias for `ReticulateMuData$new(py_obj)`. Equivalent to calling the
constructor directly but reads more naturally in a pipeline.

## Usage

``` r
reticulate_mudata(py_obj)
```

## Arguments

- py_obj:

  A Python `mudata.MuData` object.

## Value

A [ReticulateMuData](ReticulateMuData.md).

## Examples

``` r
if (FALSE) { # \dontrun{
mdata <- reticulate_mudata(py_mdata_from_python)
} # }
```
