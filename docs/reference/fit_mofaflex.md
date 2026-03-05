# Train a MOFAFLEX model

Creates and trains a `MOFAFLEX` model using the Python `mofaflex`
package via `reticulate`. The user provides data and three named R lists
describing the data, model, and training options—corresponding directly
to the Python `mofaflex.DataOptions`, `mofaflex.ModelOptions`, and
`mofaflex.TrainingOptions` dataclasses.

The function always forces `mofa_compat` in `training_options` to the
value of the `mofa_compat` argument (default `"modelonly"`), so the
saved model file is readable by the MOFA2 R package without any extra
steps.

## Usage

``` r
fit_mofaflex(
  data,
  model_options = list(),
  data_options = list(),
  training_options = list(),
  mofa_compat = "modelonly",
  smooth_options = NULL
)
```

## Arguments

- data:

  One of:

  - A [ReticulateMuData](ReticulateMuData.md) object (created by
    [`sce_to_reticulate_mudata()`](sce_to_reticulate_mudata.md),
    [`mae_to_reticulate_mudata()`](mae_to_reticulate_mudata.md), or
    [`reticulate_mudata()`](reticulate_mudata.md)).

  - A `ReticulateAnnData` object from the **anndataR** package (created
    by [`sce_to_reticulate_anndata()`](sce_to_reticulate_anndata.md)).

  - A Python `mudata.MuData` object (with `convert = FALSE`).

  - A Python `anndata.AnnData` object (with `convert = FALSE`).

  - A named R list of named lists, where outer names are group names,
    inner names are view names, and each leaf is a Python `AnnData`
    object. Example:
    `list(ctrl = list(rna = adata_ctrl), stim = list(rna = adata_stim))`.

- model_options:

  Named R list of model options passed to `mofaflex.ModelOptions(...)`.
  See the *Model options* section for all available fields and their
  defaults.

- data_options:

  Named R list of data options passed to `mofaflex.DataOptions(...)`.
  See the *Data options* section.

- training_options:

  Named R list of training options passed to
  `mofaflex.TrainingOptions(...)`. The `mofa_compat` field is **always**
  overridden by the top-level `mofa_compat` argument. See the *Training
  options* section.

- mofa_compat:

  Controls MOFA2 compatibility of the saved model file. One of:

  - `"modelonly"` *(default)* — save only model weights, readable by
    MOFA2.

  - `"full"` — save the full model state.

  - `FALSE` — disable MOFA2-compatible saving entirely.

  This value is injected into `training_options$mofa_compat` before the
  Python options object is constructed, overriding any user-supplied
  value.

- smooth_options:

  Optional named R list of smooth options passed to
  `mofaflex.SmoothOptions(...)`, required for the MEFISTO / GP-factor
  model. `NULL` (default) omits the smooth term entirely. See the
  *Smooth options* section.

## Value

A Python `mofaflex.MOFAFLEX` object (with `convert = FALSE`). Use `$` to
call methods on the trained model. Key methods include:

- `model$get_factors()` — latent factor scores as a dict of pandas
  DataFrames or AnnData objects (controlled by `return_type`).

- `model$get_weights()` — feature weights as a dict of pandas
  DataFrames.

- `model$get_r2()` — variance explained as a pandas DataFrame.

- `model$get_annotations()` — annotation enrichment results.

- `model$get_significant_factor_annotations()` — significant factor-
  programme associations.

To convert individual results to R, use
[`reticulate::py_to_r()`](https://rstudio.github.io/reticulate/reference/r-py-conversion.html)
or the `convert` helpers in the **reticulate** package. To save or
reload the model use `model$...` methods or the
`mofaflex.MOFAFLEX.load()` class method.

## Model options

All fields are keyword arguments passed verbatim to
`mofaflex.ModelOptions`. Integer fields should be supplied as R integers
(e.g., `n_factors = 5L`).

|                            |            |                                                                                                                                                                                     |
|----------------------------|------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Field                      | Default    | Description                                                                                                                                                                         |
| `n_factors`                | `0`        | Number of factors. `0` means determine automatically.                                                                                                                               |
| `weight_prior`             | `"Normal"` | Prior distribution on feature weights. One of `"Normal"`, `"Horseshoe"`, `"Laplace"`, `"SnS"`. Can be a named list mapping view names to priors.                                    |
| `factor_prior`             | `"Normal"` | Prior distribution on latent factors. One of `"Normal"`, `"GP"`, `"Horseshoe"`, `"Laplace"`, `"SnS"`.                                                                               |
| `likelihoods`              | `NULL`     | Named character vector of likelihood per view (e.g., `list(rna = "NegativeBinomial")`). Inferred automatically when `NULL`. One of `"Normal"`, `"NegativeBinomial"`, `"Bernoulli"`. |
| `nonnegative_weights`      | `FALSE`    | Enforce non-negative constraint on weights.                                                                                                                                         |
| `nonnegative_factors`      | `FALSE`    | Enforce non-negative constraint on factors.                                                                                                                                         |
| `annotation_confidence`    | `0.99`     | Confidence threshold for annotation-guided factors (used when `data_options$annotations_varm_key` is set).                                                                          |
| `guiding_vars_likelihoods` | `"Normal"` | Likelihood for guiding variable coefficients. One of `"Normal"`, `"Categorical"`, `"Bernoulli"`.                                                                                    |
| `guiding_vars_scales`      | `1.0`      | Prior scale for guiding variable coefficients.                                                                                                                                      |
| `init_factors`             | `"random"` | Factor initialisation strategy. One of `"random"`, `"orthogonal"`, `"pca"`, or a numeric value.                                                                                     |
| `init_scale`               | `0.1`      | Scale for random factor initialisation.                                                                                                                                             |

## Data options

All fields are keyword arguments passed verbatim to
`mofaflex.DataOptions`.

|                            |                     |                                                                                                                         |
|----------------------------|---------------------|-------------------------------------------------------------------------------------------------------------------------|
| Field                      | Default             | Description                                                                                                             |
| `group_by`                 | `NULL`              | Column(s) of `.obs` in `AnnData`/`MuData` used to define groups. Ignored for nested dict input.                         |
| `layer`                    | `NULL`              | Which layer to use. `NULL` → `.X`. A string applies to all views; a named list maps view names to layer names.          |
| `scale_per_group`          | `TRUE`              | Scale features within each group independently.                                                                         |
| `annotations_varm_key`     | `NULL`              | Key in `.varm` containing a binary feature-by-programme annotation matrix. When set, enables annotation-guided factors. |
| `covariates_obs_key`       | `NULL`              | `.obs` column(s) to include as fixed covariates.                                                                        |
| `covariates_obsm_key`      | `NULL`              | `.obsm` key to use as fixed covariates.                                                                                 |
| `guiding_vars_obs_keys`    | `NULL`              | `.obs` column(s) to use as guiding variables.                                                                           |
| `use_obs`                  | `"union"`           | How to align samples across views: `"union"` or `"intersection"`.                                                       |
| `use_var`                  | `"union"`           | How to align features across groups: `"union"` or `"intersection"`.                                                     |
| `subset_var`               | `"highly_variable"` | `.var` boolean column used to subset to highly variable features. `NULL` uses all features.                             |
| `plot_data_overview`       | `TRUE`              | Show a data-overview plot before training starts.                                                                       |
| `remove_constant_features` | `TRUE`              | Remove constant (zero-variance) features before fitting.                                                                |

## Training options

All fields are keyword arguments passed verbatim to
`mofaflex.TrainingOptions`. The `mofa_compat` field is **always**
overridden by the top-level `mofa_compat` argument.

|                          |          |                                                                                                                                                                          |
|--------------------------|----------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Field                    | Default  | Description                                                                                                                                                              |
| `device`                 | `"cuda"` | PyTorch device string: `"cuda"`, `"cpu"`, or `"mps"`. Falls back to CPU automatically when no GPU is available.                                                          |
| `batch_size`             | `0`      | Mini-batch size. `0` means full-batch (recommended default).                                                                                                             |
| `max_epochs`             | `10000`  | Maximum number of training epochs.                                                                                                                                       |
| `n_particles`            | `1`      | Number of particles for importance-weighted ELBO.                                                                                                                        |
| `lr`                     | `0.001`  | Base learning rate for the Adam optimizer.                                                                                                                               |
| `early_stopper_patience` | `100`    | Number of epochs without improvement before early stopping.                                                                                                              |
| `save_path`              | `NULL`   | File path to save the trained model. `NULL` does not save. When `mofa_compat` is set, the file is saved as an HDF5 readable by [MOFA2](https://biofam.github.io/MOFA2/). |
| `seed`                   | `NULL`   | Integer random seed. `NULL` → time-based seed.                                                                                                                           |
| `num_workers`            | `0`      | Number of `torch.utils.data.DataLoader` workers.                                                                                                                         |
| `pin_memory`             | `FALSE`  | Use pinned (page-locked) memory in the data loader.                                                                                                                      |

## Smooth options

Passed to `mofaflex.SmoothOptions(...)` only when `smooth_options` is
non-`NULL` (activates the MEFISTO / GP-factor extension).

|                            |                                                        |                                                            |
|----------------------------|--------------------------------------------------------|------------------------------------------------------------|
| Field                      | Default                                                | Description                                                |
| `n_inducing`               | `100`                                                  | Number of inducing points for the sparse GP approximation. |
| `kernel`                   | `"RBF"`                                                | GP kernel type. One of `"RBF"`, `"Matern"`.                |
| `mefisto_kernel`           | `TRUE`                                                 | Use the MEFISTO-style structured kernel.                   |
| `independent_lengthscales` | `FALSE`                                                | Fit an independent lengthscale per group.                  |
| `group_covar_rank`         | `1`                                                    | Rank of the inter-group covariance matrix.                 |
| `warp_groups`              | [`character()`](https://rdrr.io/r/base/character.html) | Group names to warp along the covariate axis.              |
| `warp_interval`            | `20`                                                   | Number of epochs between warping steps.                    |
| `warp_open_begin`          | `TRUE`                                                 | Allow temporal warping at the beginning of the time axis.  |
| `warp_open_end`            | `TRUE`                                                 | Allow temporal warping at the end of the time axis.        |
| `warp_reference_group`     | `NULL`                                                 | Name of the reference (unwarped) group.                    |

## See also

[`sce_to_reticulate_anndata()`](sce_to_reticulate_anndata.md),
[`sce_to_reticulate_mudata()`](sce_to_reticulate_mudata.md),
[`mae_to_reticulate_mudata()`](mae_to_reticulate_mudata.md),
[`sce_to_anndata()`](sce_to_anndata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## Minimal example with a single AnnData from a SingleCellExperiment
library(scran)
library(scuttle)

sce <- scRNAseq::ZeiselBrainData()
sce <- scuttle::logNormCounts(sce)
hvgs <- scran::getTopHVGs(scran::modelGeneVar(sce), n = 2000)
sce  <- sce[hvgs, ]

adata <- sce_to_anndata(sce, assay = "logcounts")

model <- fit_mofaflex(
  data             = list(group_1 = list(rna = adata)),
  model_options    = list(n_factors = 5L, weight_prior = "Horseshoe"),
  data_options     = list(scale_per_group = FALSE, subset_var = NULL,
                          plot_data_overview = FALSE),
  training_options = list(max_epochs = 2000L, seed = 42L, device = "cpu")
)

## Inspect variance explained
r2 <- reticulate::py_to_r(model$get_r2())
print(r2)

## Extract factor scores as AnnData
factors_ad <- model$get_factors(return_type = "anndata")
} # }
```
