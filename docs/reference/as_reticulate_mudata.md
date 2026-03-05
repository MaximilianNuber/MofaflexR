# Coerce to ReticulateMuData

If `x` is already a [ReticulateMuData](ReticulateMuData.md) it is
returned unchanged. If it is a Python `mudata.MuData` object it is
wrapped. Otherwise an error is raised.

## Usage

``` r
as_reticulate_mudata(x)
```

## Arguments

- x:

  A [ReticulateMuData](ReticulateMuData.md) or a Python `mudata.MuData`
  object.

## Value

A [ReticulateMuData](ReticulateMuData.md).

## Examples

``` r
if (FALSE) { # \dontrun{
mdata2 <- as_reticulate_mudata(py_mdata_obj)
} # }
```
