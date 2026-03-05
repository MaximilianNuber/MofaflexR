# Convert a SingleCellExperiment to a ReticulateMuData

Converts a
[SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
(possibly with alternative experiments attached via
[`SingleCellExperiment::altExps()`](https://rdrr.io/pkg/SingleCellExperiment/man/altExps.html))
to a [ReticulateMuData](ReticulateMuData.md).

The main SCE is mapped to the modality named `main_modality`. Each
alternative experiment (accessed with
[`SingleCellExperiment::altExpNames()`](https://rdrr.io/pkg/SingleCellExperiment/man/altExps.html))
is mapped to its own modality using the altExp name as the key.

All matrix transfers are **zero-copy** (via
[`sce_assay_to_scipy_csc()`](sce_assay_to_scipy_csc.md)).

## Usage

``` r
sce_to_reticulate_mudata(
  x,
  main_modality = "rna",
  assay = "counts",
  altexp_assay = NULL,
  obs = list(),
  obsm = list(),
  ...
)
```

## Arguments

- x:

  A
  [SingleCellExperiment::SingleCellExperiment](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html).

- main_modality:

  Character string; name for the main SCE modality. Default: `"rna"`.

- assay:

  Character string; assay name to use as `X` in the main modality
  (passed to [`sce_to_anndata()`](sce_to_anndata.md)). Default:
  `"counts"`.

- altexp_assay:

  Character string or `NULL`. Assay to use as `X` in each alternative
  experiment modality. If `NULL` (default) the first available assay is
  used for each altExp.

- obs:

  Named list of additional columns to add to the modalities' `obs` data
  frames. Each element must be a vector/factor of length `ncol(x)`.
  Forwarded to [`sce_to_anndata()`](sce_to_anndata.md) for each
  modality.

- obsm:

  Named list of embedding matrices (cells × dims) to attach to each
  modality's `obsm`. Forwarded to
  [`sce_to_anndata()`](sce_to_anndata.md) for the main SCE.

- ...:

  Additional arguments forwarded to
  [`sce_to_anndata()`](sce_to_anndata.md) for **every** modality
  (including altExps). Arguments that differ between modalities can be
  set via `obs`, `obsm`, etc., which only apply to the main modality.

## Value

A [ReticulateMuData](ReticulateMuData.md) wrapping a Python
`mudata.MuData`.

## See also

[`sce_to_anndata()`](sce_to_anndata.md),
[`mae_to_reticulate_mudata()`](mae_to_reticulate_mudata.md),
[ReticulateMuData](ReticulateMuData.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# ---- Minimal synthetic CITE-seq example ----
library(SingleCellExperiment)
library(Matrix)

set.seed(42)
n_cells <- 100; n_genes <- 500; n_proteins <- 30

rna_counts <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.1,
                                    repr = "C")
adt_counts <- Matrix::rsparsematrix(n_proteins, n_cells, density = 0.5,
                                    repr = "C")

sce <- SingleCellExperiment(
  assays  = list(counts = rna_counts),
  colData = S4Vectors::DataFrame(cell_type = sample(c("T","B"), n_cells, TRUE))
)
rownames(sce) <- paste0("gene_", seq_len(n_genes))
colnames(sce) <- paste0("cell_", seq_len(n_cells))

adt_sce <- SingleCellExperiment(assays = list(counts = adt_counts))
rownames(adt_sce) <- paste0("protein_", seq_len(n_proteins))
colnames(adt_sce) <- colnames(sce)

SingleCellExperiment::altExp(sce, "adt") <- adt_sce

mdata <- sce_to_reticulate_mudata(sce, main_modality = "rna", assay = "counts")
print(mdata)
mdata$obs_names   # shared cell barcodes
mdata[["rna"]]    # ReticulateAnnData for RNA
mdata[["adt"]]    # ReticulateAnnData for ADT
} # }
```
