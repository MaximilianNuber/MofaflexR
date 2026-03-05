# Wrap a bare Python AnnData in a single-modality MuData

MOFAFLEX only accepts `MuData` or a nested dict. When the user passes a
single `AnnData` (e.g.\\ the output of `sce_to_reticulate_anndata`),
this helper transparently wraps it in `mudata.MuData({"rna": adata})`.

## Usage

``` r
.maybe_wrap_anndata_in_mudata(py_obj, group_by = NULL)
```

## Arguments

- py_obj:

  A Python object extracted from the data argument.

- group_by:

  Character vector of `obs` column names used for grouping, or `NULL`.

## Value

Either `py_obj` unchanged (if it is already a MuData / dict), or a new
`mudata.MuData` wrapping it.

## Details

When `group_by` is supplied, the corresponding `obs` column(s) are
promoted to the top-level `mdata.obs` (without a modality prefix)—which
is required for `MOFAFLEX(group_by=...)` to find the column.
