# ===========================================================================
# load_mofaflex_model: load a MOFA-FLEX HDF5 model file into MOFA2
# ===========================================================================


#' Load a MOFA-FLEX model into MOFA2
#'
#' @description
#' A thin wrapper around [MOFA2::load_model()] that transparently handles an
#' incompatibility between the HDF5 files written by **mofaflex** and the
#' HDF5-group detection logic used by MOFA2 1.x.
#'
#' **Background.** MOFA2's `load_model()` calls `rhdf5::h5ls()` recursively
#' and checks whether `"covariates"` appears *anywhere* in the resulting name
#' column — including in deeply-nested paths such as
#' `/mofaflex/state/data_opts/covariates`. When the string is found it
#' unconditionally calls `h5read(file, "covariates")`, which fails because no
#' *top-level* `/covariates` group exists in mofaflex output files.
#'
#' This function detects that situation and writes an empty top-level
#' `/covariates` entry into a temporary copy of the file before passing it to
#' MOFA2, so that `load_model()` parses it without error.
#'
#' @param path Character scalar. Path to the HDF5 file written by
#'   [fit_mofaflex()] (with `mofa_compat = "full"`).
#' @param ... Additional arguments forwarded verbatim to
#'   [MOFA2::load_model()].
#'
#' @return A trained \code{MOFA} object as returned by [MOFA2::load_model()].
#'
#' @seealso [fit_mofaflex()], [MOFA2::load_model()]
#' @export
load_mofaflex_model <- function(path, ...) {
  if (!requireNamespace("MOFA2",  quietly = TRUE))
    stop("Package 'MOFA2' is required. Install with BiocManager::install('MOFA2').")
  if (!requireNamespace("rhdf5", quietly = TRUE))
    stop("Package 'rhdf5' is required. Install with BiocManager::install('rhdf5').")

  # ------------------------------------------------------------------
  # Check whether the false-positive covariates trigger is present.
  # MOFA2 uses h5ls(recursive = TRUE) and looks for ANY dataset named
  # "covariates". mofaflex stores one at /mofaflex/state/data_opts/covariates,
  # so MOFA2 falsely enters the covariates-reading branch and calls
  # h5read(file, "covariates"), which errors because no top-level
  # /covariates group exists.
  # ------------------------------------------------------------------
  all_names   <- rhdf5::h5ls(path, recursive = TRUE,  datasetinfo = FALSE)$name
  top_names   <- rhdf5::h5ls(path, recursive = FALSE, datasetinfo = FALSE)$name
  false_pos   <- ("covariates" %in% all_names) && (!"covariates" %in% top_names)

  if (false_pos) {
    # Write a minimal top-level /covariates group so MOFA2 can parse
    # it correctly (empty covariate names → no covariates in the model).
    tmp <- tempfile(fileext = ".hdf5")
    file.copy(path, tmp, overwrite = TRUE)
    on.exit(unlink(tmp), add = TRUE)

    fid <- rhdf5::H5Fopen(tmp)
    tryCatch({
      rhdf5::h5createGroup(fid, "covariates")
      # MOFA2 does: as.character(h5read(file, "covariates")[[1]])
      # An empty character dataset gives character(0), which is safe.
      rhdf5::h5write(character(0), fid, "covariates/covariate_names")
    }, finally = rhdf5::H5Fclose(fid))

    path <- tmp
  }

  MOFA2::load_model(path, ...)
}
