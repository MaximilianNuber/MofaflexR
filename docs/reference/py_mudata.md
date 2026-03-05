# Import the Python `mudata` module (convert = FALSE)

Returns the Python `mudata` module imported with
[`reticulate::import()`](https://rstudio.github.io/reticulate/reference/import.html)
and `convert = FALSE`. The result is cached so the import only happens
once per R session.

The module is useful for creating `MuData` objects inside tests or
interactive sessions before wrapping them with
[ReticulateMuData](ReticulateMuData.md):

    mu    <- py_mudata()
    mdata <- mu$MuData(dict_of_adatas)

## Usage

``` r
py_mudata()
```

## Value

The Python `mudata` module (`convert = FALSE`).

## See also

[ReticulateMuData](ReticulateMuData.md),
[`reticulate_mudata()`](reticulate_mudata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
mu <- py_mudata()
mu$`__version__`
} # }
```
