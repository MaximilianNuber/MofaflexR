# MOFAFlex basilisk Python environment

Defines the managed Python environment used by **MofaflexR** for all
Python operations. The environment is created on first use by
[`basilisk::BasiliskEnvironment()`](https://rdrr.io/pkg/basilisk/man/BasiliskEnvironment-class.html)
and reused afterward; the user does **not** need to install Python
manually.

The environment includes:

- **numpy** – used for zero-copy array views of R vectors

- **scipy** – provides `scipy.sparse.csc_matrix`

- **anndata** – AnnData container

- **mofaflex** – the MOFAFlex model-fitting library (installed via pip)

## Usage

``` r
.mofaflex_env
```

## Format

An object of class `BasiliskEnvironment` of length 1.

## Value

A
[basilisk::BasiliskEnvironment](https://rdrr.io/pkg/basilisk/man/BasiliskEnvironment-class.html)
object.

## Examples

``` r
env <- mofaflex_basilisk_env()
#> Error in mofaflex_basilisk_env(): could not find function "mofaflex_basilisk_env"
class(env)
#> Error: object 'env' not found
```
