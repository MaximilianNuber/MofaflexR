# ===========================================================================
# fit_mofaflex: train a MOFAFLEX model from R
# ===========================================================================


#' Train a MOFAFLEX model
#'
#' @description
#' Creates and trains a `MOFAFLEX` model using the Python `mofaflex` package
#' via `reticulate`. The user provides data and three named R lists describing
#' the data, model, and training options—corresponding directly to the Python
#' `mofaflex.DataOptions`, `mofaflex.ModelOptions`, and
#' `mofaflex.TrainingOptions` dataclasses.
#'
#' The function always forces `mofa_compat` in `training_options` to the value
#' of the `mofa_compat` argument (default `"modelonly"`), so the saved model
#' file is readable by the MOFA2 R package without any extra steps.
#'
#' @param data One of:
#'   * A [ReticulateMuData] object (created by [sce_to_reticulate_mudata()],
#'     [mae_to_reticulate_mudata()], or [reticulate_mudata()]).
#'   * A `ReticulateAnnData` object from the **anndataR** package (created by
#'     [sce_to_reticulate_anndata()]).
#'   * A Python `mudata.MuData` object (with `convert = FALSE`).
#'   * A Python `anndata.AnnData` object (with `convert = FALSE`).
#'   * A named R list of named lists, where outer names are group names,
#'     inner names are view names, and each leaf is a Python `AnnData` object.
#'     Example: `list(ctrl = list(rna = adata_ctrl), stim = list(rna = adata_stim))`.
#'
#' @param model_options Named R list of model options passed to
#'   `mofaflex.ModelOptions(...)`. See the *Model options* section for all
#'   available fields and their defaults.
#'
#' @param data_options Named R list of data options passed to
#'   `mofaflex.DataOptions(...)`. See the *Data options* section.
#'
#' @param training_options Named R list of training options passed to
#'   `mofaflex.TrainingOptions(...)`. The `mofa_compat` field is **always**
#'   overridden by the top-level `mofa_compat` argument. See the
#'   *Training options* section.
#'
#' @param mofa_compat Controls MOFA2 compatibility of the saved model file.
#'   One of:
#'   * `"modelonly"` *(default)* — save only model weights, readable by MOFA2.
#'   * `"full"` — save the full model state.
#'   * `FALSE` — disable MOFA2-compatible saving entirely.
#'
#'   This value is injected into `training_options$mofa_compat` before the
#'   Python options object is constructed, overriding any user-supplied value.
#'
#' @param smooth_options Optional named R list of smooth options passed to
#'   `mofaflex.SmoothOptions(...)`, required for the MEFISTO / GP-factor
#'   model. `NULL` (default) omits the smooth term entirely. See the
#'   *Smooth options* section.
#'
#' @section Model options:
#' All fields are keyword arguments passed verbatim to `mofaflex.ModelOptions`.
#' Integer fields should be supplied as R integers (e.g., `n_factors = 5L`).
#'
#' | Field | Default | Description |
#' |---|---|---|
#' | `n_factors` | `0` | Number of factors. `0` means determine automatically. |
#' | `weight_prior` | `"Normal"` | Prior distribution on feature weights. One of `"Normal"`, `"Horseshoe"`, `"Laplace"`, `"SnS"`. Can be a named list mapping view names to priors. |
#' | `factor_prior` | `"Normal"` | Prior distribution on latent factors. One of `"Normal"`, `"GP"`, `"Horseshoe"`, `"Laplace"`, `"SnS"`. |
#' | `likelihoods` | `NULL` | Named character vector of likelihood per view (e.g., `list(rna = "NegativeBinomial")`). Inferred automatically when `NULL`. One of `"Normal"`, `"NegativeBinomial"`, `"Bernoulli"`. |
#' | `nonnegative_weights` | `FALSE` | Enforce non-negative constraint on weights. |
#' | `nonnegative_factors` | `FALSE` | Enforce non-negative constraint on factors. |
#' | `annotation_confidence` | `0.99` | Confidence threshold for annotation-guided factors (used when `data_options$annotations_varm_key` is set). |
#' | `guiding_vars_likelihoods` | `"Normal"` | Likelihood for guiding variable coefficients. One of `"Normal"`, `"Categorical"`, `"Bernoulli"`. |
#' | `guiding_vars_scales` | `1.0` | Prior scale for guiding variable coefficients. |
#' | `init_factors` | `"random"` | Factor initialisation strategy. One of `"random"`, `"orthogonal"`, `"pca"`, or a numeric value. |
#' | `init_scale` | `0.1` | Scale for random factor initialisation. |
#'
#' @section Data options:
#' All fields are keyword arguments passed verbatim to `mofaflex.DataOptions`.
#'
#' | Field | Default | Description |
#' |---|---|---|
#' | `group_by` | `NULL` | Column(s) of `.obs` in `AnnData`/`MuData` used to define groups. Ignored for nested dict input. |
#' | `layer` | `NULL` | Which layer to use. `NULL` → `.X`. A string applies to all views; a named list maps view names to layer names. |
#' | `scale_per_group` | `TRUE` | Scale features within each group independently. |
#' | `annotations_varm_key` | `NULL` | Key in `.varm` containing a binary feature-by-programme annotation matrix. When set, enables annotation-guided factors. |
#' | `covariates_obs_key` | `NULL` | `.obs` column(s) to include as fixed covariates. |
#' | `covariates_obsm_key` | `NULL` | `.obsm` key to use as fixed covariates. |
#' | `guiding_vars_obs_keys` | `NULL` | `.obs` column(s) to use as guiding variables. |
#' | `use_obs` | `"union"` | How to align samples across views: `"union"` or `"intersection"`. |
#' | `use_var` | `"union"` | How to align features across groups: `"union"` or `"intersection"`. |
#' | `subset_var` | `"highly_variable"` | `.var` boolean column used to subset to highly variable features. `NULL` uses all features. |
#' | `plot_data_overview` | `TRUE` | Show a data-overview plot before training starts. |
#' | `remove_constant_features` | `TRUE` | Remove constant (zero-variance) features before fitting. |
#'
#' @section Training options:
#' All fields are keyword arguments passed verbatim to
#' `mofaflex.TrainingOptions`.  The `mofa_compat` field is **always**
#' overridden by the top-level `mofa_compat` argument.
#'
#' | Field | Default | Description |
#' |---|---|---|
#' | `device` | `"cuda"` | PyTorch device string: `"cuda"`, `"cpu"`, or `"mps"`. Falls back to CPU automatically when no GPU is available. |
#' | `batch_size` | `0` | Mini-batch size. `0` means full-batch (recommended default). |
#' | `max_epochs` | `10000` | Maximum number of training epochs. |
#' | `n_particles` | `1` | Number of particles for importance-weighted ELBO. |
#' | `lr` | `0.001` | Base learning rate for the Adam optimizer. |
#' | `early_stopper_patience` | `100` | Number of epochs without improvement before early stopping. |
#' | `save_path` | `NULL` | File path to save the trained model. `NULL` does not save. When `mofa_compat` is set, the file is saved as an HDF5 readable by [MOFA2](https://biofam.github.io/MOFA2/). |
#' | `seed` | `NULL` | Integer random seed. `NULL` → time-based seed. |
#' | `num_workers` | `0` | Number of `torch.utils.data.DataLoader` workers. |
#' | `pin_memory` | `FALSE` | Use pinned (page-locked) memory in the data loader. |
#'
#' @section Smooth options:
#' Passed to `mofaflex.SmoothOptions(...)` only when `smooth_options` is
#' non-`NULL` (activates the MEFISTO / GP-factor extension).
#'
#' | Field | Default | Description |
#' |---|---|---|
#' | `n_inducing` | `100` | Number of inducing points for the sparse GP approximation. |
#' | `kernel` | `"RBF"` | GP kernel type. One of `"RBF"`, `"Matern"`. |
#' | `mefisto_kernel` | `TRUE` | Use the MEFISTO-style structured kernel. |
#' | `independent_lengthscales` | `FALSE` | Fit an independent lengthscale per group. |
#' | `group_covar_rank` | `1` | Rank of the inter-group covariance matrix. |
#' | `warp_groups` | `character()` | Group names to warp along the covariate axis. |
#' | `warp_interval` | `20` | Number of epochs between warping steps. |
#' | `warp_open_begin` | `TRUE` | Allow temporal warping at the beginning of the time axis. |
#' | `warp_open_end` | `TRUE` | Allow temporal warping at the end of the time axis. |
#' | `warp_reference_group` | `NULL` | Name of the reference (unwarped) group. |
#'
#' @return A Python `mofaflex.MOFAFLEX` object (with `convert = FALSE`).
#'   Use `$` to call methods on the trained model. Key methods include:
#'
#'   * `model$get_factors()` — latent factor scores as a dict of pandas
#'     DataFrames or AnnData objects (controlled by `return_type`).
#'   * `model$get_weights()` — feature weights as a dict of pandas DataFrames.
#'   * `model$get_r2()` — variance explained as a pandas DataFrame.
#'   * `model$get_annotations()` — annotation enrichment results.
#'   * `model$get_significant_factor_annotations()` — significant factor-
#'     programme associations.
#'
#'   To convert individual results to R, use [reticulate::py_to_r()] or the
#'   `convert` helpers in the **reticulate** package.
#'   To save or reload the model use `model$...` methods or the
#'   `mofaflex.MOFAFLEX.load()` class method.
#'
#' @seealso
#'   [sce_to_reticulate_anndata()], [sce_to_reticulate_mudata()],
#'   [mae_to_reticulate_mudata()], [sce_to_anndata()]
#'
#' @examples
#' \dontrun{
#' ## Minimal example with a single AnnData from a SingleCellExperiment
#' library(scran)
#' library(scuttle)
#'
#' sce <- scRNAseq::ZeiselBrainData()
#' sce <- scuttle::logNormCounts(sce)
#' hvgs <- scran::getTopHVGs(scran::modelGeneVar(sce), n = 2000)
#' sce  <- sce[hvgs, ]
#'
#' adata <- sce_to_anndata(sce, assay = "logcounts")
#'
#' model <- fit_mofaflex(
#'   data             = list(group_1 = list(rna = adata)),
#'   model_options    = list(n_factors = 5L, weight_prior = "Horseshoe"),
#'   data_options     = list(scale_per_group = FALSE, subset_var = NULL,
#'                           plot_data_overview = FALSE),
#'   training_options = list(max_epochs = 2000L, seed = 42L, device = "cpu")
#' )
#'
#' ## Inspect variance explained
#' r2 <- reticulate::py_to_r(model$get_r2())
#' print(r2)
#'
#' ## Extract factor scores as AnnData
#' factors_ad <- model$get_factors(return_type = "anndata")
#' }
#'
#' @export
fit_mofaflex <- function(
    data,
    model_options    = list(),
    data_options     = list(),
    training_options = list(),
    mofa_compat      = "modelonly",
    smooth_options   = NULL
) {
  # Validate mofa_compat
  valid_compat <- c("modelonly", "full", FALSE)
  if (!identical(mofa_compat, FALSE) &&
      !(is.character(mofa_compat) && mofa_compat %in% c("modelonly", "full"))) {
    stop(
      "'mofa_compat' must be \"modelonly\", \"full\", or FALSE."
    )
  }

  mfl <- .mofaflex_module()

  # --------------------------------------------------------------------------
  # 1. Unwrap data from R wrappers to the underlying Python object
  # --------------------------------------------------------------------------
  py_data <- .extract_py_data_for_mofaflex(data)

  # MOFAFLEX.__init__ only accepts MuData | Mapping[str, Mapping[str, AnnData]].
  # A bare AnnData (common when the user passes sce_to_reticulate_anndata output)
  # must be wrapped in a single-modality MuData so group_by works correctly.
  group_by <- data_options[["group_by"]]
  py_data  <- .maybe_wrap_anndata_in_mudata(py_data, group_by)

  # --------------------------------------------------------------------------
  # 2. Force mofa_compat into training_options (always overrides user value)
  # --------------------------------------------------------------------------
  training_options[["mofa_compat"]] <- mofa_compat

  # --------------------------------------------------------------------------
  # 3. Build Python option objects
  # --------------------------------------------------------------------------
  data_opts  <- do.call(mfl$DataOptions,     data_options)
  model_opts <- do.call(mfl$ModelOptions,    model_options)
  train_opts <- do.call(mfl$TrainingOptions, training_options)

  # --------------------------------------------------------------------------
  # 4. Construct and train the model
  # --------------------------------------------------------------------------
  if (!is.null(smooth_options)) {
    smooth_opts <- do.call(mfl$SmoothOptions, smooth_options)
    model <- mfl$MOFAFLEX(py_data, data_opts, model_opts, train_opts, smooth_opts)
  } else {
    model <- mfl$MOFAFLEX(py_data, data_opts, model_opts, train_opts)
  }

  model
}


# ===========================================================================
# Internal helpers
# ===========================================================================

#' @keywords internal
.mofaflex_module <- function() {
  if (!exists("mofaflex", envir = .module_cache, inherits = FALSE)) {
    if (!reticulate::py_module_available("mofaflex")) {
      stop(
        "The Python package 'mofaflex' is not available in the current ",
        "Python environment.\n",
        "Install it with:\n",
        "  reticulate::py_install('mofaflex')\n",
        "or activate an environment that contains mofaflex."
      )
    }
    assign(
      "mofaflex",
      reticulate::import("mofaflex", convert = FALSE),
      envir = .module_cache
    )
  }
  get("mofaflex", envir = .module_cache, inherits = FALSE)
}


# #' @keywords internal
# .extract_py_data_for_mofaflex <- function(data) {
#   # ---- ReticulateMuData wrapper ------------------------------------------
#   if (inherits(data, "ReticulateMuData")) {
#     return(data$py_mudata())
#   }

#   # ---- anndataR ReticulateAnnData wrapper --------------------------------
#   # anndataR stores the Python object in private$py
#   if (inherits(data, "AnnDataR6")) {
#     py_obj <- tryCatch(
#       data$.__enclos_env__$private$py,
#       error = function(e) NULL
#     )
#     if (!is.null(py_obj)) {
#       return(py_obj)
#     }
#     stop(
#       "Could not extract the underlying Python AnnData from the supplied ",
#       "AnnDataR6 object. Please pass the Python AnnData directly."
#     )
#   }

#   # ---- Named R list (nested: group -> view -> Python AnnData) -------------
#   if (is.list(data) && !is.null(names(data))) {
#     return(.r_nested_list_to_py_dict(data))
#   }

#   # ---- Bare Python object (AnnData, MuData, etc.) — pass through ---------
#   data
# }

#' @keywords internal
.extract_py_data_for_mofaflex <- function(data) {

  # ---- ReticulateMuData wrapper ------------------------------------------
  if (inherits(data, "ReticulateMuData")) {
    # Pick ONE canonical accessor in your class design and stick to it.
    # Common: data$py (public field) or data$py_mudata()
    if (!is.null(data$py)) return(data$py)
    if (is.function(data$py_mudata)) return(data$py_mudata())
    stop("ReticulateMuData does not expose the underlying python object via `$py` or `$py_mudata()`.")
  }

  # ---- anndataR ReticulateAnnData wrapper --------------------------------
  # Accept both common class tags: "ReticulateAnnData" and/or "AnnDataR6"
  if (inherits(data, c("ReticulateAnnData"))) {

    # Try a few likely locations in decreasing brittleness.
    # 1) Public field (some wrappers use $py)
    py_obj <- tryCatch(data$py, error = function(e) NULL)
    if (!is.null(py_obj)) return(py_obj)

    # 2) Private field used by many R6 wrappers (what you already do)
    py_obj <- tryCatch(data$.__enclos_env__$private$py, error = function(e) NULL)
    if (!is.null(py_obj)) return(py_obj)

    # 3) Another common private name (seen sometimes)
    py_obj <- tryCatch(data$.__enclos_env__$private$adata, error = function(e) NULL)
    if (!is.null(py_obj)) return(py_obj)

    py_obj <- tryCatch(data$py_anndata(), error = function(e) NULL)
    if (!is.null(py_obj)) return(py_obj)

    stop(
      "Could not extract the underlying Python AnnData from the supplied anndataR wrapper. ",
      "Tried `$py`, `private$py`, `private$adata`. ",
      "As a workaround, pass the underlying python object directly."
    )
  }

  # ---- Named R list (nested: group -> view -> Python AnnData) -------------
  if (is.list(data) && !is.null(names(data))) {
    return(.r_nested_list_to_py_dict(data))
  }

  # ---- Bare Python object (AnnData, MuData, etc.) — pass through ----------
  # (Optional but strongly recommended) validate that it's AnnData/MuData if it's a python object
  if (inherits(data, "python.builtin.object")) {
    # best-effort type check without forcing heavy imports
    is_ok <- TRUE
    # If you want strictness, uncomment and implement `.py_is_anndata_or_mudata()`
    # is_ok <- .py_is_anndata_or_mudata(data)
    if (is_ok) return(data)
  }

  data
}


#' @keywords internal
.r_nested_list_to_py_dict <- function(x) {
  # Convert a named R list of named lists of Python objects to a Python dict
  # of dicts.  Each leaf must already be a Python object (e.g., AnnData).
  inner_dicts <- lapply(x, function(group) {
    if (!is.list(group) || is.null(names(group))) {
      stop(
        "When 'data' is a named list, each element must itself be a named ",
        "list mapping view names to Python AnnData objects.\n",
        "Example: list(group_1 = list(rna = adata1), group_2 = list(rna = adata2))"
      )
    }
    reticulate::py_dict(
      keys    = as.list(names(group)),
      values  = unname(group),
      convert = FALSE
    )
  })
  reticulate::py_dict(
    keys    = as.list(names(x)),
    values  = inner_dicts,
    convert = FALSE
  )
}

#' Wrap a bare Python AnnData in a single-modality MuData
#'
#' MOFAFLEX only accepts \code{MuData} or a nested dict.  When the user passes
#' a single \code{AnnData} (e.g.\ the output of \code{sce_to_reticulate_anndata}),
#' this helper transparently wraps it in \code{mudata.MuData(\{"rna": adata\})}.
#'
#' When \code{group_by} is supplied, the corresponding \code{obs} column(s) are
#' promoted to the top-level \code{mdata.obs} (without a modality prefix)—which
#' is required for \code{MOFAFLEX(group_by=...)} to find the column.
#'
#' @param py_obj A Python object extracted from the data argument.
#' @param group_by Character vector of \code{obs} column names used for
#'   grouping, or \code{NULL}.
#' @return Either \code{py_obj} unchanged (if it is already a MuData / dict),
#'   or a new \code{mudata.MuData} wrapping it.
#' @keywords internal
.maybe_wrap_anndata_in_mudata <- function(py_obj, group_by = NULL) {
  if (!inherits(py_obj, "python.builtin.object")) return(py_obj)

  # Check whether py_obj is an AnnData instance (via __class__.__name__).
  is_anndata <- tryCatch(
    identical(reticulate::py_to_r(py_obj$`__class__`$`__name__`), "AnnData"),
    error = function(e) FALSE
  )

  if (!isTRUE(is_anndata)) return(py_obj)  # already MuData or dict

  mu_mod <- reticulate::import("mudata", convert = FALSE)

  # Disable auto-pull so obs columns are NOT prefixed as "rna:<col>"
  mu_mod$set_options(pull_on_update = FALSE)

  mdata <- mu_mod$MuData(
    reticulate::py_dict(
      keys    = list("rna"),
      values  = list(py_obj),
      convert = FALSE
    )
  )

  # Promote group_by column(s) to top-level mdata.obs (no "rna:" prefix).
  # mofaflex DataOptions(group_by=...) looks up the column in mdata.obs.
  if (!is.null(group_by)) {
    mod_obs <- py_obj$obs          # adata.obs (individual modality)
    for (col in as.character(group_by)) {
      # Use pandas __setitem__ because R [[<- converts to attribute access
      mdata$obs$`__setitem__`(col, mod_obs$`__getitem__`(col))
    }
  }

  mdata
}

#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b