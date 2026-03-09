# ===========================================================================
# sce_to_reticulate_mudata
# ===========================================================================
#
# Convert a SingleCellExperiment (with optional altExps) to a ReticulateMuData.
#
# Each modality is built via sce_to_anndata():
#   * The main SCE becomes the "main_modality" AnnData.
#   * Every altExp(x, nm) becomes a separate AnnData modality.
#
# Matrix transfers use se_assay_to_python_matrix(), which dispatches on the
# assay class (dgCMatrix → CSC, dgRMatrix → CSR, COO → COO, dense → NumPy)
# and applies zero-copy conversion where the class supports it.
# ===========================================================================


#' Convert a SingleCellExperiment to a ReticulateMuData
#'
#' @description
#' Converts a [SingleCellExperiment::SingleCellExperiment] (possibly with
#' alternative experiments attached via [SingleCellExperiment::altExps()]) to
#' a [ReticulateMuData].
#'
#' The main SCE is mapped to the modality named `main_modality`.  Each
#' alternative experiment (accessed with
#' [SingleCellExperiment::altExpNames()]) is mapped to its own modality using
#' the altExp name as the key.
#'
#' Matrix transfers use [se_assay_to_python_matrix()], dispatching by assay
#' matrix class (dgCMatrix → SciPy CSC, dgRMatrix → CSR, COO → COO, dense
#' → NumPy array) and applying zero-copy conversion where possible.
#'
#' @param x A [SingleCellExperiment::SingleCellExperiment].
#' @param main_modality Character string; name for the main SCE modality.
#'   Default: `"rna"`.
#' @param assay Character string; assay name to use as `X` in the main
#'   modality (passed to [sce_to_anndata()]).  Default: `"counts"`.
#' @param altexp_assay Character string or `NULL`.  Assay to use as `X` in
#'   each alternative experiment modality.  If `NULL` (default) the first
#'   available assay is used for each altExp.
#' @param obs Named list of additional columns to add to the modalities' `obs`
#'   data frames.  Each element must be a vector/factor of length
#'   `ncol(x)`. Forwarded to [sce_to_anndata()] for each modality.
#' @param obsm Named list of embedding matrices (cells × dims) to attach to
#'   each modality's `obsm`. Forwarded to [sce_to_anndata()] for the main SCE.
#' @param ... Additional arguments forwarded to [sce_to_anndata()] for
#'   **every** modality (including altExps).  Arguments that differ between
#'   modalities can be set via `obs`, `obsm`, etc., which only apply to the
#'   main modality.
#'
#' @return A [ReticulateMuData] wrapping a Python `mudata.MuData`.
#'
#' @seealso [sce_to_anndata()], [mae_to_reticulate_mudata()],
#'   [ReticulateMuData]
#'
#' @examples
#' \dontrun{
#' # ---- Minimal synthetic CITE-seq example ----
#' library(SingleCellExperiment)
#' library(Matrix)
#'
#' set.seed(42)
#' n_cells <- 100; n_genes <- 500; n_proteins <- 30
#'
#' rna_counts <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.1,
#'                                     repr = "C")
#' adt_counts <- Matrix::rsparsematrix(n_proteins, n_cells, density = 0.5,
#'                                     repr = "C")
#'
#' sce <- SingleCellExperiment(
#'   assays  = list(counts = rna_counts),
#'   colData = S4Vectors::DataFrame(cell_type = sample(c("T","B"), n_cells, TRUE))
#' )
#' rownames(sce) <- paste0("gene_", seq_len(n_genes))
#' colnames(sce) <- paste0("cell_", seq_len(n_cells))
#'
#' adt_sce <- SingleCellExperiment(assays = list(counts = adt_counts))
#' rownames(adt_sce) <- paste0("protein_", seq_len(n_proteins))
#' colnames(adt_sce) <- colnames(sce)
#'
#' SingleCellExperiment::altExp(sce, "adt") <- adt_sce
#'
#' mdata <- sce_to_reticulate_mudata(sce, main_modality = "rna", assay = "counts")
#' print(mdata)
#' mdata$obs_names   # shared cell barcodes
#' mdata[["rna"]]    # ReticulateAnnData for RNA
#' mdata[["adt"]]    # ReticulateAnnData for ADT
#' }
#'
#' @export
sce_to_reticulate_mudata <- function(
    x,
    main_modality = "rna",
    assay         = "counts",
    altexp_assay  = NULL,
    obs           = list(),
    obsm          = list(),
    ...) {

  stopifnot(
    inherits(x, "SingleCellExperiment"),
    is.character(main_modality), length(main_modality) == 1L,
    is.character(assay),         length(assay) == 1L
  )

  mu <- py_mudata()

  # Normalize: an empty list for obs/obsm means "derive from colData/nothing"
  # sce_to_anndata treats any non-NULL obs as a user-supplied data.frame.
  if (is.list(obs)  && !is.data.frame(obs)  && length(obs)  == 0L) obs  <- NULL
  if (is.list(obsm) && !is.data.frame(obsm) && length(obsm) == 0L) obsm <- list()

  # ------------------------------------------------------------------
  # Build Python AnnData for the main modality
  # ------------------------------------------------------------------
  py_main <- sce_to_anndata(x,
    assay = assay,
    obs   = obs,
    obsm  = obsm,
    ...
  )

  # ------------------------------------------------------------------
  # Build Python AnnData for each altExp modality
  # ------------------------------------------------------------------
  ae_names <- SingleCellExperiment::altExpNames(x)

  mod_keys   <- c(main_modality, ae_names)
  mod_values <- vector("list", length(mod_keys))
  mod_values[[1L]] <- py_main

  for (i in seq_along(ae_names)) {
    ae_nm  <- ae_names[[i]]
    ae_sce <- SingleCellExperiment::altExp(x, ae_nm)

    # Determine which assay to use for this altExp
    ae_assay <- if (!is.null(altexp_assay)) {
      altexp_assay
    } else {
      nms <- SummarizedExperiment::assayNames(ae_sce)
      if (length(nms) == 0L) {
        stop(sprintf(
          "altExp '%s' has no named assays.  Provide `altexp_assay`.",
          ae_nm
        ))
      }
      nms[[1L]]
    }

    mod_values[[i + 1L]] <- sce_to_anndata(ae_sce, assay = ae_assay, ...)
  }

  # ------------------------------------------------------------------
  # Build MuData from the modality dict
  # ------------------------------------------------------------------
  py_mod_dict <- reticulate::py_dict(
    keys   = as.list(mod_keys),
    values = mod_values,
    convert = FALSE
  )

  py_mdata <- mu$MuData(py_mod_dict)

  ReticulateMuData$new(py_mdata)
}
