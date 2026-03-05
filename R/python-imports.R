# ===========================================================================
# Cached Python module imports
# ===========================================================================
#
# These helpers import Python modules with `convert = FALSE` and cache the
# result in a package-level (function-local) environment so that the import
# call only happens once per R session.
# ===========================================================================


# Module cache environment (not exported, not documented externally)
.module_cache <- new.env(parent = emptyenv())


#' Import the Python \code{mudata} module (convert = FALSE)
#'
#' @description
#' Returns the Python `mudata` module imported with `reticulate::import()` and
#' `convert = FALSE`.  The result is cached so the import only happens once per
#' R session.
#'
#' The module is useful for creating `MuData` objects inside tests or
#' interactive sessions before wrapping them with [ReticulateMuData]:
#'
#' ```r
#' mu    <- py_mudata()
#' mdata <- mu$MuData(dict_of_adatas)
#' ```
#'
#' @return The Python `mudata` module (`convert = FALSE`).
#'
#' @seealso [ReticulateMuData], [reticulate_mudata()]
#'
#' @examples
#' \dontrun{
#' mu <- py_mudata()
#' mu$`__version__`
#' }
#'
#' @export
py_mudata <- function() {
  if (!exists("mudata", envir = .module_cache, inherits = FALSE)) {
    assign("mudata",
      reticulate::import("mudata", convert = FALSE),
      envir = .module_cache
    )
  }
  get("mudata", envir = .module_cache, inherits = FALSE)
}
