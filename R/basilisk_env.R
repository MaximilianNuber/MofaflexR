#' MOFAFlex basilisk Python environment
#'
#' @description
#' Defines the managed Python environment used by **MofaflexR** for all
#' Python operations.  The environment is created on first use by
#' [basilisk::BasiliskEnvironment()] and reused afterward; the user does
#' **not** need to install Python manually.
#'
#' The environment includes:
#' * **numpy** – used for zero-copy array views of R vectors
#' * **scipy** – provides `scipy.sparse.csc_matrix`
#' * **anndata** – AnnData container
#' * **mofaflex** – the MOFAFlex model-fitting library (installed via pip)
#'
#' @return A [basilisk::BasiliskEnvironment-class] object.
#'
#' @examples
#' env <- mofaflex_basilisk_env()
#' class(env)
#'
#' @export
# mofaflex_basilisk_env <- function() {
#   basilisk::BasiliskEnvironment(
#     envname  = "mofaflex_env_v1",
#     pkgname  = "MofaflexR",
#     packages = c(
#       "numpy>=1.24",
#       "scipy>=1.10",
#       "pandas>=1.5"
#     ),
#     pip = c(
#       "anndata>=0.10",
#       "mudata>=0.3",
#       "mofaflex>=0.1.0"
#     )
#   )
# }
.mofaflex_env <- basilisk::BasiliskEnvironment(
  envname = "mofaflexr_env",
  pkgname = "MofaflexR",
  packages = c(
    "python==3.11.9",
    "numpy==2.4.2",
    "scipy==1.17.1",
    "pandas==2.3.3",
    "h5py==3.15.1",
    "anndata==0.12.10",
    "mudata==0.3.3",
    "mofaflex==0.1.0.post1"
  )
)

#' Evaluate an expression inside the MOFAFlex Python environment
#'
#' @description
#' A thin wrapper around [basilisk::basiliskRun()] that evaluates `expr`
#' inside the managed Python environment returned by [mofaflex_basilisk_env()].
#' Use this wrapper to ensure that all Python imports resolve to the correct
#' environment.
#'
#' @param expr An R expression (unquoted) to evaluate.
#'
#' @return The value of `expr`.
#'
#' @examples
#' \dontrun{
#' with_mofaflex_env({
#'   np <- reticulate::import("numpy", convert = FALSE)
#'   np$array(1:5)
#' })
#' }
#'
#' @export
#' @export
with_mofaflex_env <- function(expr, ..., .packages = "MofaflexR") {
  env <- .mofaflex_env

  # Capture objects to send to the worker
  dots <- list(...)
  caller <- parent.frame()

  fun <- eval(bquote(function() {
    # Attach required R packages inside the worker
    for (pkg in .(as.list(.packages))) {
      suppressPackageStartupMessages(require(pkg, character.only = TRUE))
    }

    # Recreate exported objects in the worker's global environment
    if (length(.(as.list(names(dots)))) > 0) {
      nms <- .(as.list(names(dots)))
      vals <- .(as.list(dots))
      for (i in seq_along(nms)) assign(nms[[i]], vals[[i]], envir = .GlobalEnv)
    }

    .(substitute(expr))
  }), envir = caller)

  environment(fun) <- caller
  basilisk::basiliskRun(env = env, fun = fun)
}