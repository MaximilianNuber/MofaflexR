# Wrap a Python object in the appropriate R wrapper

Central dispatch: if `obj` is a Python `anndata.AnnData` it is wrapped
in
[`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html);
if it is a `mudata.MuData` it is wrapped in
[ReticulateMuData](ReticulateMuData.md). Any other Python object (or
plain R object) is returned unchanged.

## Usage

``` r
.wrap_py_return(obj)
```

## Arguments

- obj:

  Any R or Python object.

## Value

Wrapped R object or `obj` unchanged.
