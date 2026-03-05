# Evaluate an expression inside the MOFAFlex Python environment

A thin wrapper around
[`basilisk::basiliskRun()`](https://rdrr.io/pkg/basilisk/man/basiliskStart.html)
that evaluates `expr` inside the managed Python environment returned by
`mofaflex_basilisk_env()`. Use this wrapper to ensure that all Python
imports resolve to the correct environment.

## Usage

``` r
with_mofaflex_env(expr, ..., .packages = "MofaflexR")
```

## Arguments

- expr:

  An R expression (unquoted) to evaluate.

## Value

The value of `expr`.

## Examples

``` r
if (FALSE) { # \dontrun{
with_mofaflex_env({
  np <- reticulate::import("numpy", convert = FALSE)
  np$array(1:5)
})
} # }
```
