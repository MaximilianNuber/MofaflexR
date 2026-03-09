# ===========================================================================
# mae_to_reticulate_mudata
# ===========================================================================
#
# Convert a MultiAssayExperiment to a ReticulateMuData.
#
# Each experiment in the MAE is converted to an AnnData modality via
# sce_to_anndata().  Experiments that are not SCE/SE objects raise a warning
# and are skipped.
#
# Requires the `MultiAssayExperiment` package (Suggests); checked at runtime.
# ===========================================================================


#' Convert a MultiAssayExperiment to a ReticulateMuData
#'
#' @description
#' Converts each experiment in a
#' [MultiAssayExperiment::MultiAssayExperiment()] to a Python `anndata.AnnData`
#' modality, then assembles them into a [ReticulateMuData].
#'
#' ## Required package
#' This function requires the **MultiAssayExperiment** Bioconductor package,
#' which is listed under `Suggests`.  Install it with:
#' ```r
#' BiocManager::install("MultiAssayExperiment")
#' ```
#'
#' ## Experiment handling
#' * [SummarizedExperiment::SummarizedExperiment] and
#'   [SingleCellExperiment::SingleCellExperiment] experiments are converted via
#'   [sce_to_anndata()].
#' * Experiments of any other type emit a `warning()` and are **skipped**.
#'
#' ## Sample (column) metadata
#' The MAE's `colData` is joined to each modality's `obs` before building
#' the MuData, so that all per-sample annotations are available globally via
#' `mdata$obs` once `mdata$update()` is called.
#'
#' @param mae A [MultiAssayExperiment::MultiAssayExperiment].
#' @param assay Character string or `NULL`.  Name of the assay to use as `X`
#'   for **every** experiment.  If `NULL` (default) the first available assay
#'   is used for each experiment individually.
#' @param ... Additional arguments forwarded to [sce_to_anndata()] for each
#'   experiment modality.
#'
#' @return A [ReticulateMuData] wrapping a Python `mudata.MuData`.
#'
#' @seealso [sce_to_reticulate_mudata()], [sce_to_anndata()], [ReticulateMuData]
#'
#' @examples
#' \dontrun{
#' # ---- Tiny MAE with two experiments ----
#' library(MultiAssayExperiment)
#' library(SingleCellExperiment)
#' library(Matrix)
#' library(S4Vectors)
#'
#' set.seed(1)
#' n_samples <- 50
#' sample_ids <- paste0("sample_", seq_len(n_samples))
#'
#' rna_mat  <- Matrix::rsparsematrix(200, n_samples, 0.1, repr = "C")
#' rownames(rna_mat) <- paste0("gene_", seq_len(200))
#' colnames(rna_mat) <- sample_ids
#'
#' adt_mat  <- Matrix::rsparsematrix(30, n_samples, 0.4, repr = "C")
#' rownames(adt_mat) <- paste0("protein_", seq_len(30))
#' colnames(adt_mat) <- sample_ids
#'
#' sce_rna  <- SingleCellExperiment(assays = list(counts = rna_mat))
#' sce_adt  <- SingleCellExperiment(assays = list(counts = adt_mat))
#'
#' sample_df <- S4Vectors::DataFrame(
#'   condition = sample(c("ctrl", "treat"), n_samples, TRUE),
#'   row.names = sample_ids
#' )
#'
#' el <- MultiAssayExperiment::ExperimentList(rna = sce_rna, adt = sce_adt)
#' cl <- MultiAssayExperiment::listToMap(
#'         setNames(lapply(c("rna","adt"), function(nm) {
#'           data.frame(primary = sample_ids, colname = sample_ids)
#'         }), c("rna","adt"))
#'       )
#' mae <- MultiAssayExperiment::MultiAssayExperiment(
#'   experiments = el,
#'   colData     = sample_df,
#'   sampleMap   = cl
#' )
#'
#' mdata <- mae_to_reticulate_mudata(mae)
#' print(mdata)
#' mdata$n_obs()
#' mdata[["rna"]]
#' mdata[["adt"]]
#' }
#'
#' @export
mae_to_reticulate_mudata <- function(mae, assay = NULL, ...) {

  # ---- Check that MultiAssayExperiment is available ----------------------- #
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop(
      "Package 'MultiAssayExperiment' is required to use ",
      "mae_to_reticulate_mudata().  Install it with:\n",
      "  BiocManager::install('MultiAssayExperiment')"
    )
  }

  stopifnot(
    inherits(mae, "MultiAssayExperiment")
  )

  mu       <- py_mudata()
  exp_list <- MultiAssayExperiment::experiments(mae)
  exp_nms  <- names(exp_list)

  if (length(exp_nms) == 0L) {
    stop("The MultiAssayExperiment has no experiments.")
  }

  mod_keys   <- character(0)
  mod_values <- list()

  for (nm in exp_nms) {
    exp_obj <- exp_list[[nm]]

    # Accept SE and SCE; warn & skip others
    if (!inherits(exp_obj, "SummarizedExperiment")) {
      warning(sprintf(
        "Experiment '%s' is not a SummarizedExperiment/SingleCellExperiment (%s); skipping.",
        nm, paste(class(exp_obj), collapse = ", ")
      ))
      next
    }

    # Determine the assay name for this experiment
    this_assay <- if (!is.null(assay)) {
      assay
    } else {
      nms <- SummarizedExperiment::assayNames(exp_obj)
      if (length(nms) == 0L) {
        warning(sprintf(
          "Experiment '%s' has no named assays; skipping.", nm
        ))
        next
      }
      nms[[1L]]
    }

    py_ad <- sce_to_anndata(exp_obj, assay = this_assay, ...)

    mod_keys   <- c(mod_keys, nm)
    mod_values <- c(mod_values, list(py_ad))
  }

  if (length(mod_keys) == 0L) {
    stop("No convertible experiments found in the MultiAssayExperiment.")
  }

  py_mod_dict <- reticulate::py_dict(
    keys    = as.list(mod_keys),
    values  = mod_values,
    convert = FALSE
  )

  py_mdata <- mu$MuData(py_mod_dict)
  mdata    <- ReticulateMuData$new(py_mdata)

  # Push MAE colData columns to the global obs (best-effort)
  col_df <- as.data.frame(MultiAssayExperiment::colData(mae))
  if (nrow(col_df) > 0L && ncol(col_df) > 0L) {
    tryCatch({
      # Only rows whose names appear in obs_names
      obs_nms <- mdata$obs_names
      shared  <- intersect(rownames(col_df), obs_nms)
      if (length(shared) > 0L) {
        # Add MAE-level columns to the Python obs dataframe
        py_obs  <- reticulate::py_get_attr(mdata$py_mudata(), "obs")
        pd      <- reticulate::import("pandas", convert = FALSE)
        for (col_nm in colnames(col_df)) {
          vals <- col_df[shared, col_nm, drop = TRUE]
          py_series <- pd$Series(
            data  = reticulate::r_to_py(as.character(vals)),
            index = reticulate::r_to_py(shared)
          )
          py_obs[[col_nm]] <- py_series
        }
      }
    }, error = function(e) {
      warning("Could not push MAE colData to MuData global obs: ", conditionMessage(e))
    })
  }

  mdata
}
