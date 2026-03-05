# Convert an R data.frame to a pandas DataFrame with index = rownames

Convert an R data.frame to a pandas DataFrame with index = rownames

## Usage

``` r
.df_to_pandas(df, pd)
```

## Arguments

- df:

  An R `data.frame` with row names set.

- pd:

  A reticulate handle to the `pandas` module.

## Value

A Python pandas `DataFrame`.
