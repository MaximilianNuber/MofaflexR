# Convert a MultiAssayExperiment to a ReticulateMuData

Converts each experiment in a
[`MultiAssayExperiment::MultiAssayExperiment()`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html)
to a Python `anndata.AnnData` modality, then assembles them into a
[ReticulateMuData](ReticulateMuData.md).

### Required package

This function requires the **MultiAssayExperiment** Bioconductor
package, which is listed under `Suggests`. Install it with:

    BiocManager::install("MultiAssayExperiment")

### Experiment handling

- [SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  and
  [SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  experiments are converted via [`sce_to_anndata()`](sce_to_anndata.md).

- Experiments of any other type emit a
  [`warning()`](https://rdrr.io/r/base/warning.html) and are
  **skipped**.

### Sample (column) metadata

The MAE's `colData` is joined to each modality's `obs` before building
the MuData, so that all per-sample annotations are available globally
via `mdata$obs` once `mdata$update()` is called.

## Usage

``` r
mae_to_reticulate_mudata(mae, assay = NULL, ...)
```

## Arguments

- mae:

  A
  [MultiAssayExperiment::MultiAssayExperiment](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html).

- assay:

  Character string or `NULL`. Name of the assay to use as `X` for
  **every** experiment. If `NULL` (default) the first available assay is
  used for each experiment individually.

- ...:

  Additional arguments forwarded to
  [`sce_to_anndata()`](sce_to_anndata.md) for each experiment modality.

## Value

A [ReticulateMuData](ReticulateMuData.md) wrapping a Python
`mudata.MuData`.

## See also

[`sce_to_reticulate_mudata()`](sce_to_reticulate_mudata.md),
[`sce_to_anndata()`](sce_to_anndata.md),
[ReticulateMuData](ReticulateMuData.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# ---- Tiny MAE with two experiments ----
library(MultiAssayExperiment)
library(SingleCellExperiment)
library(Matrix)
library(S4Vectors)

set.seed(1)
n_samples <- 50
sample_ids <- paste0("sample_", seq_len(n_samples))

rna_mat  <- Matrix::rsparsematrix(200, n_samples, 0.1, repr = "C")
rownames(rna_mat) <- paste0("gene_", seq_len(200))
colnames(rna_mat) <- sample_ids

adt_mat  <- Matrix::rsparsematrix(30, n_samples, 0.4, repr = "C")
rownames(adt_mat) <- paste0("protein_", seq_len(30))
colnames(adt_mat) <- sample_ids

sce_rna  <- SingleCellExperiment(assays = list(counts = rna_mat))
sce_adt  <- SingleCellExperiment(assays = list(counts = adt_mat))

sample_df <- S4Vectors::DataFrame(
  condition = sample(c("ctrl", "treat"), n_samples, TRUE),
  row.names = sample_ids
)

el <- MultiAssayExperiment::ExperimentList(rna = sce_rna, adt = sce_adt)
cl <- MultiAssayExperiment::listToMap(
        setNames(lapply(c("rna","adt"), function(nm) {
          data.frame(primary = sample_ids, colname = sample_ids)
        }), c("rna","adt"))
      )
mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = el,
  colData     = sample_df,
  sampleMap   = cl
)

mdata <- mae_to_reticulate_mudata(mae)
print(mdata)
mdata$n_obs()
mdata[["rna"]]
mdata[["adt"]]
} # }
```
