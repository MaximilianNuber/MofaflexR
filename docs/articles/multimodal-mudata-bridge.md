# Bridging Single-Cell Multi-Omics: SCE and MAE to Python MuData

## Introduction

**MofaflexR** provides matrix-aware bridges between Bioconductor
single-cell objects and the Python `mudata` / `anndata` ecosystem.
Common dense and sparse assay matrix types are converted to NumPy arrays
or SciPy sparse matrices, as appropriate, and zero-copy is used where
possible for supported representations. This vignette demonstrates:

1.  **[`sce_to_reticulate_mudata()`](../reference/sce_to_reticulate_mudata.md)**:
    Convert a `SingleCellExperiment` – with alternative experiments
    (altExps) encoding additional modalities such as surface-protein
    (ADT/CITE-seq) or chromatin-accessibility (ATAC) – into a
    `ReticulateMuData`.

2.  **[`mae_to_reticulate_mudata()`](../reference/mae_to_reticulate_mudata.md)**:
    Convert a `MultiAssayExperiment` with multiple experiments into a
    `ReticulateMuData`.

3.  **`ReticulateMuData` active bindings and methods**: Access Python
    MuData fields directly from R without any data copies.

Matrix transfers use the same dispatch machinery as
[`sce_to_anndata()`](../reference/sce_to_anndata.md): each modality’s
assay is routed to the appropriate converter based on its class
(`dgCMatrix` → `csc_matrix`, `dgRMatrix` → `csr_matrix`, `dgeMatrix` /
base `matrix` → `ndarray`, `COO_SparseMatrix` → `coo_matrix`). For the
sparse formats used in typical single-cell count matrices, zero-copy is
used for the backing slot arrays. No Python ↔︎ R serialisation takes
place unless you explicitly call
[`reticulate::py_to_r()`](https://rstudio.github.io/reticulate/reference/r-py-conversion.html).

### Prerequisites

``` r
# Ensure the mofaflex_test conda environment is available (or your own)
reticulate::use_condaenv("mofaflex_test", required = TRUE)

# Or configure via ~/.Renviron:
#   RETICULATE_PYTHON=/path/to/python
```

``` r
suppressPackageStartupMessages({
  library(MofaflexR)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(S4Vectors)
})
```

------------------------------------------------------------------------

## Synthetic CITE-seq Dataset

We construct a minimal **CITE-seq** dataset that faithfully mimics real
single-cell data: 300 cells, 1 000 genes (RNA), and 40 surface proteins
(ADT). Both RNA and ADT are represented as sparse `dgCMatrix` objects,
which is typical for single-cell count data. The same conversion
machinery also handles dense assay matrices (e.g. `dgeMatrix`, base
`matrix`).

``` r
set.seed(42)

n_cells    <- 300L
n_genes    <- 1000L
n_proteins <- 40L

cell_ids    <- paste0("cell_",    seq_len(n_cells))
gene_ids    <- paste0("gene_",    seq_len(n_genes))
protein_ids <- paste0("protein_", seq_len(n_proteins))

# ---- RNA counts (sparse, ~5% density) -----------------------------------
rna_mat <- Matrix::rsparsematrix(
  nrow    = n_genes,
  ncol    = n_cells,
  density = 0.05,
  repr    = "C"
)
rna_mat <- round(abs(rna_mat) * 100)   # non-negative integer counts
rownames(rna_mat) <- gene_ids
colnames(rna_mat) <- cell_ids

# ---- ADT / protein counts (~50% density, lower diversity) ---------------
adt_mat <- Matrix::rsparsematrix(
  nrow    = n_proteins,
  ncol    = n_cells,
  density = 0.50,
  repr    = "C"
)
adt_mat <- round(abs(adt_mat) * 500)
rownames(adt_mat) <- protein_ids
colnames(adt_mat) <- cell_ids

# ---- Main SCE (RNA) -----------------------------------------------------
cell_types <- sample(c("CD4_T", "CD8_T", "B_cell", "Monocyte", "NK"),
                     n_cells, replace = TRUE)
batches    <- sample(c("batch1", "batch2"), n_cells, replace = TRUE)

sce <- SingleCellExperiment(
  assays  = list(counts = rna_mat),
  colData = DataFrame(
    cell_type = cell_types,
    batch     = batches,
    n_counts  = Matrix::colSums(rna_mat)
  ),
  rowData = DataFrame(
    gene_length = sample(500:5000, n_genes, replace = TRUE),
    highly_variable = sample(c(TRUE, FALSE), n_genes,
                             replace = TRUE, prob = c(0.1, 0.9))
  )
)

# ---- ADT altExp --------------------------------------------------------
adt_sce <- SingleCellExperiment(
  assays  = list(counts = adt_mat),
  rowData = DataFrame(
    isotype = sample(c("IgG", "IgM", "IgA"), n_proteins, replace = TRUE)
  )
)

SingleCellExperiment::altExp(sce, "adt") <- adt_sce

sce
#> class: SingleCellExperiment 
#> dim: 1000 300 
#> metadata(0):
#> assays(1): counts
#> rownames(1000): gene_1 gene_2 ... gene_999 gene_1000
#> rowData names(2): gene_length highly_variable
#> colnames(300): cell_1 cell_2 ... cell_299 cell_300
#> colData names(3): cell_type batch n_counts
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(1): adt
```

``` r
cat("Main SCE dimensions  : ", dim(sce)[1], "genes x", dim(sce)[2], "cells\n")
#> Main SCE dimensions  :  1000 genes x 300 cells
cat("altExp names         : ", altExpNames(sce), "\n")
#> altExp names         :  adt
cat("ADT dimensions       : ", dim(altExp(sce, "adt")), "\n")
#> ADT dimensions       :  40 300
```

------------------------------------------------------------------------

## SCE → ReticulateMuData

### Convert with `sce_to_reticulate_mudata()`

[`sce_to_reticulate_mudata()`](../reference/sce_to_reticulate_mudata.md)
maps the main SCE to one modality and each `altExp` to its own modality.
Each modality’s assay is converted independently using the same
matrix-aware dispatch as
[`sce_to_anndata()`](../reference/sce_to_anndata.md): for sparse count
matrices stored in a supported format (e.g. `dgCMatrix`), zero-copy is
used for the underlying data arrays; dense matrices are converted to
NumPy arrays via the buffer protocol. The resulting Python objects share
memory with R for the duration of the session.

``` r
mdata <- sce_to_reticulate_mudata(
  x             = sce,
  main_modality = "rna",
  assay         = "counts",
  altexp_assay  = "counts"
)
mdata
#> ReticulateMuData
#>   n_obs x n_vars : 300 x 1040
#>   modalities (2) : rna, adt
#>     rna:           300 obs x 1000 vars
#>     adt:           300 obs x 40 vars
#>   obs keys       : rna:cell_type, rna:batch, rna:n_counts
#>   obsm keys      : rna, adt
```

`mdata` is a thin R6 proxy around a Python `mudata.MuData` object. All
data live in Python; R just holds a pointer.

### Active bindings – access Python fields directly

#### `$shape` – dimensions

``` r
mdata$shape          # c(n_obs, n_vars_total)
#> [1]  300 1040
```

#### `$obs_names` – shared observation names

``` r
head(mdata$obs_names, 10)
#>  [1] "cell_1"  "cell_2"  "cell_3"  "cell_4"  "cell_5"  "cell_6"  "cell_7" 
#>  [8] "cell_8"  "cell_9"  "cell_10"
length(mdata$obs_names)
#> [1] 300
```

#### `$var_names` – all variable names

``` r
head(mdata$var_names, 6)
#> [1] "gene_1" "gene_2" "gene_3" "gene_4" "gene_5" "gene_6"
tail(mdata$var_names, 6)   # protein names appear at the end
#> [1] "protein_35" "protein_36" "protein_37" "protein_38" "protein_39"
#> [6] "protein_40"
```

#### `$obs` – global observation metadata

After construction, `mudata.MuData` automatically pulls
observation-level metadata from the individual modalities into the
global `obs` table.

``` r
str(mdata$obs)
#> 'data.frame':    300 obs. of  3 variables:
#>  $ rna:cell_type: chr  "B_cell" "B_cell" "Monocyte" "B_cell" ...
#>  $ rna:batch    : chr  "batch2" "batch2" "batch2" "batch1" ...
#>  $ rna:n_counts : num  4973 3290 3945 3131 4498 ...
#>  - attr(*, "pandas.index")=Index(['cell_1', 'cell_2', 'cell_3', 'cell_4', 'cell_5', 'cell_6', 'cell_7',
#>        'cell_8', 'cell_9', 'cell_10',
#>        ...
#>        'cell_291', 'cell_292', 'cell_293', 'cell_294', 'cell_295', 'cell_296',
#>        'cell_297', 'cell_298', 'cell_299', 'cell_300'],
#>       dtype='object', length=300)
head(mdata$obs)
#>        rna:cell_type rna:batch rna:n_counts
#> cell_1        B_cell    batch2         4973
#> cell_2        B_cell    batch2         3290
#> cell_3      Monocyte    batch2         3945
#> cell_4        B_cell    batch1         3131
#> cell_5            NK    batch2         4498
#> cell_6         CD8_T    batch2         4645
```

#### `$mod` – extracted modalities

The `$mod` active binding returns a named R list of
[`anndataR::ReticulateAnnData`](https://anndataR.scverse.org/reference/ReticulateAnnData.html)
wrappers.

``` r
names(mdata$mod)
#> [1] "rna" "adt"
rna_adata <- mdata$mod$rna
rna_adata
#> ReticulateAnnData object with n_obs × n_vars = 300 × 1000
#>     obs: 'cell_type', 'batch', 'n_counts'
#>     var: 'gene_length', 'highly_variable'
```

Or equivalently using `[[`:

``` r
adt_adata <- mdata[["adt"]]
adt_adata
#> ReticulateAnnData object with n_obs × n_vars = 300 × 40
#>     var: 'isotype'
```

Both wrappers remain zero-copy: they point to the same Python `AnnData`
objects that are held inside the `MuData`.

#### `$obsm` – multi-dimensional observation annotations

The `obsm` slot is empty initially; we add a synthetic PCA embedding to
illustrate the round-trip.

``` r
# Create a synthetic PCA result (cells × 10 PCs)
pca_result <- matrix(rnorm(n_cells * 10), nrow = n_cells, ncol = 10,
                     dimnames = list(cell_ids, paste0("PC", 1:10)))

# Store directly via the Python AnnData binding
rna_adata$obsm[["X_pca"]] <- pca_result
rna_adata$obsm_keys()
#> [1] "X_pca"
```

After adding embeddings to a modality, call `mdata$update_obs()` to
synchronise the global obs index:

``` r
mdata$update_obs()
cat("Global n_obs after update:", mdata$n_obs(), "\n")
#> Global n_obs after update: 300
```

#### obs / var key introspection

``` r
mdata$obs_keys()     # global obs column names
#> [1] "rna:cell_type" "rna:batch"     "rna:n_counts"
mdata$var_keys()     # global var column names (if any)
#> [1] "rna:gene_length"     "rna:highly_variable" "adt:isotype"
mdata$uns_keys()     # unstructured annotation (empty here)
#> list()
```

------------------------------------------------------------------------

### Python method forwarding

Any Python `MuData` method can be called directly from R via
`$<method>()`. The return value is automatically wrapped when
appropriate.

#### `$copy()` – deep copy

``` r
mdata_copy <- mdata$copy()
inherits(mdata_copy, "ReticulateMuData")   # TRUE
#> [1] TRUE
identical(mdata_copy$py_mudata(), mdata$py_mudata())   # FALSE – new object
#> [1] FALSE
```

#### `$update()` – sync global obs / var from modalities

``` r
mdata$update()
mdata$n_obs()
#> [1] 300
mdata$n_vars()
#> [1] 1040
```

#### `$py_call()` – generic Python method with auto-wrapping

``` r
# Any Python method not directly surfaced as an R method can be called via
# py_call().  The return is auto-wrapped into ReticulateMuData / ReticulateAnnData
# when applicable.
result <- mdata$py_call("copy")
class(result)
#> [1] "ReticulateMuData" "R6"
```

------------------------------------------------------------------------

## MultiAssayExperiment → ReticulateMuData

### Construct a tiny MAE

``` r
library(MultiAssayExperiment)

# Reuse the same RNA and protein matrices from above
sce_rna <- SingleCellExperiment(
  assays  = list(counts = rna_mat),
  colData = DataFrame(cell_type = cell_types, batch = batches)
)

sce_adt <- SingleCellExperiment(
  assays  = list(counts = adt_mat)
)

# Sample map: each experiment has identical primary ↔ colname
make_sm <- function(primary, colname, assay) {
  DataFrame(assay = assay, primary = primary, colname = colname)
}

smap <- rbind(
  make_sm(cell_ids, cell_ids, "rna"),
  make_sm(cell_ids, cell_ids, "adt")
)

col_df <- DataFrame(
  donor     = sample(paste0("donor_", 1:3), n_cells, replace = TRUE),
  row.names = cell_ids
)

mae <- MultiAssayExperiment(
  experiments = ExperimentList(rna = sce_rna, adt = sce_adt),
  colData     = col_df,
  sampleMap   = smap
)

mae
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] rna: SingleCellExperiment with 1000 rows and 300 columns
#>  [2] adt: SingleCellExperiment with 40 rows and 300 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```

### Convert with `mae_to_reticulate_mudata()`

``` r
mdata_mae <- mae_to_reticulate_mudata(mae, assay = "counts")
mdata_mae
#> ReticulateMuData
#>   n_obs x n_vars : 300 x 1040
#>   modalities (2) : rna, adt
#>     rna:           300 obs x 1000 vars
#>     adt:           300 obs x 40 vars
#>   obs keys       : rna:cell_type, rna:batch
#>   obsm keys      : rna, adt
```

``` r
mdata_mae$modality_names()
#> [1] "rna" "adt"
mdata_mae$n_obs()
#> [1] 300
mdata_mae$n_vars()
#> [1] 1040

# Check that MAE colData was propagated to global obs
str(mdata_mae$obs)
#> 'data.frame':    300 obs. of  2 variables:
#>  $ rna:cell_type: chr  "B_cell" "B_cell" "Monocyte" "B_cell" ...
#>  $ rna:batch    : chr  "batch2" "batch2" "batch2" "batch1" ...
#>  - attr(*, "pandas.index")=Index(['cell_1', 'cell_2', 'cell_3', 'cell_4', 'cell_5', 'cell_6', 'cell_7',
#>        'cell_8', 'cell_9', 'cell_10',
#>        ...
#>        'cell_291', 'cell_292', 'cell_293', 'cell_294', 'cell_295', 'cell_296',
#>        'cell_297', 'cell_298', 'cell_299', 'cell_300'],
#>       dtype='object', length=300)
```

``` r
rna_from_mae <- mdata_mae[["rna"]]
adt_from_mae <- mdata_mae[["adt"]]

cat("RNA modality: ", rna_from_mae$n_obs(), "obs x", rna_from_mae$n_vars(), "vars\n")
#> RNA modality:  300 obs x 1000 vars
cat("ADT modality: ", adt_from_mae$n_obs(), "obs x", adt_from_mae$n_vars(), "vars\n")
#> ADT modality:  300 obs x 40 vars
```

------------------------------------------------------------------------

## Real Dataset Examples

The sections above use synthetic data so the vignette can be rendered
offline. The following examples run the same bridge on real datasets.

### SCE with altExps: Buettner et al. ESC data

[`BuettnerESCData()`](https://rdrr.io/pkg/scRNAseq/man/BuettnerESCData.html)
from the `scRNAseq` package provides 288 mouse embryonic stem cells with
ERCC spike-in counts stored in an `altExp`. This is a compact real
dataset that exercises the multi-modality path.

``` r
library(scRNAseq)

sce_buettner <- BuettnerESCData()
altExpNames(sce_buettner)   # "ERCC"
#> [1] "ERCC"

mdata_buettner <- sce_to_reticulate_mudata(
  x             = sce_buettner,
  main_modality = "rna",
  assay         = "counts"
)
mdata_buettner$modality_names()   # "rna" and "ERCC"
#> [1] "rna"  "ERCC"
mdata_buettner$shape
#> [1]   288 38385
```

### MultiAssayExperiment: ACC data

`miniACC` ships with the `MultiAssayExperiment` package and contains
five assays from an adrenocortical carcinoma study. `Mutations` (plain
`matrix`) is skipped because it is not a `SummarizedExperiment`;
`gistict` is skipped because its assay is unnamed. The remaining three
SE assays are converted normally.

``` r
data(miniACC)
names(experiments(miniACC))
#> [1] "RNASeq2GeneNorm" "gistict"         "RPPAArray"       "Mutations"      
#> [5] "miRNASeqGene"

mdata_acc <- mae_to_reticulate_mudata(miniACC)
mdata_acc$modality_names()
#> [1] "RNASeq2GeneNorm" "RPPAArray"       "miRNASeqGene"
mdata_acc$n_obs()
#> [1] 205
```

------------------------------------------------------------------------

## Interoperability Round-trip

A `ReticulateMuData` wraps a live Python `MuData` object. The raw Python
object is accessible via `$py_mudata()` and can be passed to any Python
code in the mudata ecosystem. Here we write `mdata` (the synthetic
CITE-seq object from earlier) to HDF5 and read it back.

``` r
mu <- reticulate::import("mudata", convert = FALSE)

h5mu_path <- tempfile(fileext = ".h5mu")
invisible(mu$write_h5mu(h5mu_path, mdata$py_mudata()))

py_loaded  <- mu$read_h5mu(h5mu_path)
mdata_back <- ReticulateMuData$new(py_loaded)

head(mdata_back$obs_names, 4)
#> [1] "cell_1" "cell_2" "cell_3" "cell_4"
mdata_back$shape
#> [1]  300 1040
```

------------------------------------------------------------------------

## Session Info

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Berlin
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] ensembldb_2.34.0            AnnotationFilter_1.34.0    
#>  [3] GenomicFeatures_1.62.0      AnnotationDbi_1.72.0       
#>  [5] scRNAseq_2.24.0             MultiAssayExperiment_1.36.1
#>  [7] Matrix_1.7-4                SingleCellExperiment_1.32.0
#>  [9] SummarizedExperiment_1.40.0 Biobase_2.70.0             
#> [11] GenomicRanges_1.62.1        Seqinfo_1.0.0              
#> [13] IRanges_2.44.0              S4Vectors_0.48.0           
#> [15] BiocGenerics_0.56.0         generics_0.1.4             
#> [17] MatrixGenerics_1.22.0       matrixStats_1.5.0          
#> [19] MofaflexR_0.0.1             BiocStyle_2.38.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] DBI_1.3.0                bitops_1.0-9             httr2_1.2.2             
#>  [4] anndataR_1.0.2           rlang_1.1.7              magrittr_2.0.4          
#>  [7] otel_0.2.0               gypsum_1.6.0             compiler_4.5.2          
#> [10] RSQLite_2.4.6            png_0.1-8                systemfonts_1.3.1       
#> [13] vctrs_0.7.1              ProtGenerics_1.42.0      pkgconfig_2.0.3         
#> [16] crayon_1.5.3             fastmap_1.2.0            dbplyr_2.5.2            
#> [19] XVector_0.50.0           Rsamtools_2.26.0         rmarkdown_2.30          
#> [22] UCSC.utils_1.6.1         ragg_1.5.0               purrr_1.2.1             
#> [25] bit_4.6.0                xfun_0.56                cachem_1.1.0            
#> [28] cigarillo_1.0.0          GenomeInfoDb_1.46.2      jsonlite_2.0.0          
#> [31] blob_1.3.0               rhdf5filters_1.22.0      DelayedArray_0.36.0     
#> [34] Rhdf5lib_1.32.0          BiocParallel_1.44.0      parallel_4.5.2          
#> [37] R6_2.6.1                 bslib_0.10.0             reticulate_1.45.0       
#> [40] rtracklayer_1.70.1       jquerylib_0.1.4          Rcpp_1.1.1              
#> [43] bookdown_0.46            knitr_1.51               BiocBaseUtils_1.12.0    
#> [46] tidyselect_1.2.1         abind_1.4-8              yaml_2.3.12             
#> [49] codetools_0.2-20         curl_7.0.0               alabaster.sce_1.10.0    
#> [52] lattice_0.22-7           tibble_3.3.1             withr_3.0.2             
#> [55] KEGGREST_1.50.0          evaluate_1.0.5           desc_1.4.3              
#> [58] BiocFileCache_3.0.0      alabaster.schemas_1.10.0 ExperimentHub_3.0.0     
#> [61] Biostrings_2.78.0        pillar_1.11.1            BiocManager_1.30.27     
#> [64] filelock_1.0.3           RCurl_1.98-1.17          BiocVersion_3.22.0      
#> [67] alabaster.base_1.10.0    alabaster.ranges_1.10.0  glue_1.8.0              
#> [70] lazyeval_0.2.2           alabaster.matrix_1.10.0  tools_4.5.2             
#> [73] AnnotationHub_4.0.0      BiocIO_1.20.0            GenomicAlignments_1.46.0
#> [76] fs_1.6.6                 XML_3.99-0.22            rhdf5_2.54.1            
#> [79] grid_4.5.2               HDF5Array_1.38.0         restfulr_0.0.16         
#> [82] cli_3.6.5                rappdirs_0.3.4           textshaping_1.0.4       
#> [85] S4Arrays_1.10.1          dplyr_1.2.0              alabaster.se_1.10.0     
#> [88] sass_0.4.10              digest_0.6.39            SparseArray_1.10.8      
#> [91] rjson_0.2.23             htmlwidgets_1.6.4        memoise_2.0.1           
#> [94] htmltools_0.5.9          pkgdown_2.2.0            lifecycle_1.0.5         
#> [97] h5mread_1.2.1            httr_1.4.8               bit64_4.6.0-1
```
