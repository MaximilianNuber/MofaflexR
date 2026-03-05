# Load a MOFA-FLEX model into MOFA2

A thin wrapper around
[`MOFA2::load_model()`](https://rdrr.io/pkg/MOFA2/man/load_model.html)
that transparently handles an incompatibility between the HDF5 files
written by **mofaflex** and the HDF5-group detection logic used by MOFA2
1.x.

**Background.** MOFA2's `load_model()` calls
[`rhdf5::h5ls()`](https://rdrr.io/pkg/rhdf5/man/h5ls.html) recursively
and checks whether `"covariates"` appears *anywhere* in the resulting
name column — including in deeply-nested paths such as
`/mofaflex/state/data_opts/covariates`. When the string is found it
unconditionally calls `h5read(file, "covariates")`, which fails because
no *top-level* `/covariates` group exists in mofaflex output files.

This function detects that situation and writes an empty top-level
`/covariates` entry into a temporary copy of the file before passing it
to MOFA2, so that `load_model()` parses it without error.

## Usage

``` r
load_mofaflex_model(path, ...)
```

## Arguments

- path:

  Character scalar. Path to the HDF5 file written by
  [`fit_mofaflex()`](fit_mofaflex.md) (with `mofa_compat = "full"`).

- ...:

  Additional arguments forwarded verbatim to
  [`MOFA2::load_model()`](https://rdrr.io/pkg/MOFA2/man/load_model.html).

## Value

A trained `MOFA` object as returned by
[`MOFA2::load_model()`](https://rdrr.io/pkg/MOFA2/man/load_model.html).

## See also

[`fit_mofaflex()`](fit_mofaflex.md),
[`MOFA2::load_model()`](https://rdrr.io/pkg/MOFA2/man/load_model.html)
