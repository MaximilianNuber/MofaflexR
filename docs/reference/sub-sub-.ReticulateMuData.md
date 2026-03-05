# Extract a modality from a ReticulateMuData

Returns the named modality as an
[`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html)
wrapper. The invariant **"any AnnData extracted from MuData is a
ReticulateAnnData"** is enforced here.

## Usage

``` r
# S3 method for class 'ReticulateMuData'
x[[i]]
```

## Arguments

- x:

  A [ReticulateMuData](ReticulateMuData.md).

- i:

  A single string: the modality name (e.g. `"rna"`).

## Value

An
[`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html)
object.

## Examples

``` r
if (FALSE) { # \dontrun{
rna <- mdata[["rna"]]
inherits(rna, "ReticulateAnnData")  # TRUE
} # }
```
