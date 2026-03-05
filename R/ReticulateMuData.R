# ===========================================================================
# ReticulateMuData – R6 wrapper for Python mudata.MuData objects
# ===========================================================================
#
# Design principles (mirrors anndataR::ReticulateAnnData):
#   * The underlying Python object is always stored with `convert = FALSE`.
#   * All Python access goes through reticulate with `convert = FALSE` to
#     avoid inadvertent densification or data copies.
#   * AnnData objects extracted from the MuData are always wrapped in
#     `anndataR::ReticulateAnnData`; MuData objects are wrapped in
#     `ReticulateMuData`.  This interception is centralised in the internal
#     helper `.wrap_py_return()`.
#   * `[[` / `[[<-` are provided as S3 methods so that users can index
#     modalities with familiar list syntax.
#
# Python class hierarchy used for type dispatch (reticulate S3 class names):
#   MuData  : "mudata._core.mudata.MuData"
#   AnnData : "anndata._core.anndata.AnnData"
#
# Active bindings (property-style, no `()`):
#   obs, var, obs_names, var_names, shape, axis
#   obsm, varm, obsp, varp, uns, mod
#
# Public methods (must be called with `()`):
#   py_mudata(), py_call(), print()
#   n_obs(), n_vars(), modality_names()
#   update(), update_obs(), update_var(), copy()
#   pull_obs(), pull_var(), push_obs(), push_var()
#   obs_keys(), var_keys(), obsm_keys(), varm_keys(), uns_keys()
# ===========================================================================


# ---------------------------------------------------------------------------
# Internal helper: centralised Python-return interception
# ---------------------------------------------------------------------------

#' Wrap a Python object in the appropriate R wrapper
#'
#' @description
#' Central dispatch: if `obj` is a Python `anndata.AnnData` it is wrapped in
#' `anndataR::ReticulateAnnData`; if it is a `mudata.MuData` it is wrapped in
#' [ReticulateMuData].  Any other Python object (or plain R object) is
#' returned unchanged.
#'
#' @param obj Any R or Python object.
#' @return Wrapped R object or `obj` unchanged.
#'
#' @keywords internal
.wrap_py_return <- function(obj) {
  if (!inherits(obj, "python.builtin.object")) {
    return(obj)
  }
  if (inherits(obj, "anndata._core.anndata.AnnData")) {
    RAD <- get("ReticulateAnnData",
      envir     = asNamespace("anndataR"),
      inherits  = FALSE
    )
    return(RAD$new(py_anndata = obj))
  }
  if (inherits(obj, "mudata._core.mudata.MuData")) {
    return(ReticulateMuData$new(obj))
  }
  obj
}


# ---------------------------------------------------------------------------
# Internal helpers for dict-like Python mapping → R named list
# ---------------------------------------------------------------------------

# Convert a Python mapping (obsm, varm, obsp, varp, uns, …) to a named R list.
# Each value is converted individually via py_to_r().
.py_mapping_to_list <- function(py_map) {
  bi   <- reticulate::import_builtins(convert = FALSE)
  keys <- reticulate::py_to_r(bi$list(py_map$keys()))
  out  <- vector("list", length(keys))
  names(out) <- keys
  for (k in keys) {
    out[[k]] <- tryCatch(
      reticulate::py_to_r(py_map[[k]]),
      error = function(e) py_map[[k]]   # fall back to raw Python object
    )
  }
  out
}


# ---------------------------------------------------------------------------
# R6 class definition
# ---------------------------------------------------------------------------

#' @title ReticulateMuData
#'
#' @description
#' An R6 wrapper around a Python **mudata** `MuData` object accessed via
#' \pkg{reticulate}.  All data remain in Python; the R object is a thin proxy
#' that forwards attribute reads/writes to the underlying Python object.
#'
#' Modalities (sub-`AnnData` objects) are always returned as
#' `anndataR::ReticulateAnnData` wrappers, so the critical invariant holds:
#' **any `AnnData` extracted from a `ReticulateMuData` is returned as an
#' `anndataR::ReticulateAnnData`.**
#'
#' @section Constructor:
#' ```r
#' mdata <- ReticulateMuData$new(x)
#' ```
#' where `x` is either
#' * a Python `mudata.MuData` object (from `reticulate::import("mudata")`), or
#' * an existing `ReticulateMuData` (wrapping is idempotent).
#'
#' @section Modality access:
#' ```r
#' rna <- mdata[["rna"]]          # returns ReticulateAnnData
#' mdata[["rna"]] <- new_rna_adata  # accepts ReticulateAnnData or raw Python AnnData
#' ```
#'
#' @section Python object access:
#' ```r
#' py_obj <- mdata$py_mudata()    # the raw Python MuData
#' ```
#'
#' Use `py_call()` to invoke Python methods with automatic wrapping of the
#' result:
#' ```r
#' result <- mdata$py_call("copy")  # returns ReticulateMuData if result is MuData
#' ```
#'
#' @seealso [is_reticulate_mudata()], [as_reticulate_mudata()],
#'   [reticulate_mudata()]
#'
#' @examples
#' \dontrun{
#' reticulate::use_condaenv("mofaflex_test", required = TRUE)
#' mu   <- reticulate::import("mudata",  convert = FALSE)
#' ad   <- reticulate::import("anndata", convert = FALSE)
#' sp   <- reticulate::import("scipy.sparse", convert = FALSE)
#' np   <- reticulate::import("numpy",   convert = FALSE)
#'
#' X1 <- sp$csr_matrix(reticulate::r_to_py(matrix(1:12, 3L, 4L)))
#' X2 <- sp$csr_matrix(reticulate::r_to_py(matrix(1:6,  3L, 2L)))
#' a1  <- ad$AnnData(X = X1)
#' a2  <- ad$AnnData(X = X2)
#' md  <- mu$MuData(reticulate::py_dict(list("rna" = a1, "prot" = a2)))
#'
#' mdata <- ReticulateMuData$new(md)
#' print(mdata)
#'
#' rna   <- mdata[["rna"]]   # ReticulateAnnData
#' print(rna)
#' }
#'
#' @export
ReticulateMuData <- R6::R6Class(          # nolint
  "ReticulateMuData",
  cloneable = FALSE,
  private   = list(
    .py_mudata = NULL,

    # Abort if the object has not been properly initialised.
    .check_valid = function() {
      if (is.null(private$.py_mudata)) {
        stop("ReticulateMuData: underlying Python MuData has not been initialised.")
      }
    }
  ),

  # -------------------------------------------------------------------------
  # Active bindings – Python properties exposed as R fields (no `()` needed)
  # -------------------------------------------------------------------------
  active = list(

    #' @field obs Global observation data frame (pandas DataFrame ↔ R data.frame).
    obs = function(value) {
      private$.check_valid()
      if (missing(value)) {
        reticulate::py_to_r(reticulate::py_get_attr(private$.py_mudata, "obs"))
      } else {
        reticulate::py_set_attr(private$.py_mudata, "obs", reticulate::r_to_py(value))
        invisible(self)
      }
    },

    #' @field var Global variable data frame (pandas DataFrame ↔ R data.frame).
    var = function(value) {
      private$.check_valid()
      if (missing(value)) {
        reticulate::py_to_r(reticulate::py_get_attr(private$.py_mudata, "var"))
      } else {
        reticulate::py_set_attr(private$.py_mudata, "var", reticulate::r_to_py(value))
        invisible(self)
      }
    },

    #' @field obs_names Global observation names (character vector).
    obs_names = function(value) {
      private$.check_valid()
      if (missing(value)) {
        bi <- reticulate::import_builtins(convert = FALSE)
        reticulate::py_to_r(bi$list(
          reticulate::py_get_attr(private$.py_mudata, "obs_names")
        ))
      } else {
        reticulate::py_set_attr(
          private$.py_mudata, "obs_names", reticulate::r_to_py(value)
        )
        invisible(self)
      }
    },

    #' @field var_names Global variable names (character vector).
    var_names = function(value) {
      private$.check_valid()
      if (missing(value)) {
        bi <- reticulate::import_builtins(convert = FALSE)
        reticulate::py_to_r(bi$list(
          reticulate::py_get_attr(private$.py_mudata, "var_names")
        ))
      } else {
        reticulate::py_set_attr(
          private$.py_mudata, "var_names", reticulate::r_to_py(value)
        )
        invisible(self)
      }
    },

    #' @field shape Integer vector `c(n_obs, n_vars)` (read-only).
    shape = function(value) {
      if (!missing(value)) stop("shape is read-only; adjust individual modalities instead.")
      private$.check_valid()
      as.integer(reticulate::py_to_r(
        reticulate::py_get_attr(private$.py_mudata, "shape")
      ))
    },

    #' @field axis MuData axis: 0 = shared obs, 1 = shared var (read-only).
    axis = function(value) {
      if (!missing(value)) stop("axis is read-only; set at construction time.")
      private$.check_valid()
      as.integer(reticulate::py_to_r(
        reticulate::py_get_attr(private$.py_mudata, "axis")
      ))
    },

    #' @field obsm Multi-dimensional observation annotation (named R list).
    obsm = function(value) {
      private$.check_valid()
      if (missing(value)) {
        .py_mapping_to_list(reticulate::py_get_attr(private$.py_mudata, "obsm"))
      } else {
        reticulate::py_set_attr(private$.py_mudata, "obsm", reticulate::r_to_py(value))
        invisible(self)
      }
    },

    #' @field varm Multi-dimensional variable annotation (named R list).
    varm = function(value) {
      private$.check_valid()
      if (missing(value)) {
        .py_mapping_to_list(reticulate::py_get_attr(private$.py_mudata, "varm"))
      } else {
        reticulate::py_set_attr(private$.py_mudata, "varm", reticulate::r_to_py(value))
        invisible(self)
      }
    },

    #' @field obsp Pairwise observation annotation (named R list of square matrices).
    obsp = function(value) {
      private$.check_valid()
      if (missing(value)) {
        .py_mapping_to_list(reticulate::py_get_attr(private$.py_mudata, "obsp"))
      } else {
        reticulate::py_set_attr(private$.py_mudata, "obsp", reticulate::r_to_py(value))
        invisible(self)
      }
    },

    #' @field varp Pairwise variable annotation (named R list of square matrices).
    varp = function(value) {
      private$.check_valid()
      if (missing(value)) {
        .py_mapping_to_list(reticulate::py_get_attr(private$.py_mudata, "varp"))
      } else {
        reticulate::py_set_attr(private$.py_mudata, "varp", reticulate::r_to_py(value))
        invisible(self)
      }
    },

    #' @field uns Unstructured annotation (named R list).
    uns = function(value) {
      private$.check_valid()
      if (missing(value)) {
        reticulate::py_to_r(reticulate::py_get_attr(private$.py_mudata, "uns"))
      } else {
        reticulate::py_set_attr(private$.py_mudata, "uns", reticulate::r_to_py(value))
        invisible(self)
      }
    },

    #' @field mod Named R list of modalities; each value is a
    #'   `anndataR::ReticulateAnnData`.  Use `mdata[["name"]] <- adata` to
    #'   update individual modalities.
    mod = function(value) {
      if (!missing(value)) {
        stop("Assign with `mdata[[\"name\"]] <- adata` to update a modality.")
      }
      private$.check_valid()
      bi     <- reticulate::import_builtins(convert = FALSE)
      py_mod <- reticulate::py_get_attr(private$.py_mudata, "mod")
      keys   <- reticulate::py_to_r(bi$list(py_mod$keys()))
      out    <- vector("list", length(keys))
      names(out) <- keys
      for (k in keys) {
        out[[k]] <- .wrap_py_return(py_mod[[k]])
      }
      out
    }
  ),

  # -------------------------------------------------------------------------
  # Public methods – must be called with `()`
  # -------------------------------------------------------------------------
  public = list(

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    #' @description
    #' Create a new `ReticulateMuData` wrapper.
    #'
    #' @param x A Python `mudata.MuData` object (class
    #'   `"mudata._core.mudata.MuData"`) **or** an existing
    #'   `ReticulateMuData`.  Any other type raises an error.
    initialize = function(x) {
      # Unwrap ReticulateMuData to get the underlying Python object.
      if (inherits(x, "ReticulateMuData")) {
        x <- x$py_mudata()
      }
      if (!inherits(x, "mudata._core.mudata.MuData")) {
        stop(
          "'x' must be a Python mudata.MuData object (got class: ",
          paste(class(x), collapse = ", "), ")."
        )
      }
      private$.py_mudata <- x
      invisible(self)
    },

    # ------------------------------------------------------------------
    # Core accessors
    # ------------------------------------------------------------------

    #' @description Return the underlying Python `mudata.MuData` object.
    #' @return A Python `mudata.MuData` object (`convert = FALSE`).
    py_mudata = function() {
      private$.check_valid()
      private$.py_mudata
    },

    #' @description Return the names of available modalities.
    #' @return A character vector of modality names.
    modality_names = function() {
      private$.check_valid()
      bi <- reticulate::import_builtins(convert = FALSE)
      reticulate::py_to_r(bi$list(private$.py_mudata$mod$keys()))
    },

    #' @description Number of observations (shared across modalities).
    #' @return A single integer.
    n_obs = function() {
      private$.check_valid()
      as.integer(reticulate::py_to_r(
        reticulate::py_get_attr(private$.py_mudata, "n_obs")
      ))
    },

    #' @description Total number of variables (features across all modalities).
    #' @return A single integer.
    n_vars = function() {
      private$.check_valid()
      as.integer(reticulate::py_to_r(
        reticulate::py_get_attr(private$.py_mudata, "n_vars")
      ))
    },

    # ------------------------------------------------------------------
    # Key-listing methods
    # ------------------------------------------------------------------

    #' @description Keys of global observation annotation (`obs` columns).
    #' @return A character vector.
    obs_keys = function() {
      private$.check_valid()
      reticulate::py_to_r(private$.py_mudata$obs_keys())
    },

    #' @description Keys of global variable annotation (`var` columns).
    #' @return A character vector.
    var_keys = function() {
      private$.check_valid()
      reticulate::py_to_r(private$.py_mudata$var_keys())
    },

    #' @description Keys in `obsm` (multi-dimensional observation annotation).
    #' @return A character vector.
    obsm_keys = function() {
      private$.check_valid()
      reticulate::py_to_r(private$.py_mudata$obsm_keys())
    },

    #' @description Keys in `varm` (multi-dimensional variable annotation).
    #' @return A character vector.
    varm_keys = function() {
      private$.check_valid()
      reticulate::py_to_r(private$.py_mudata$varm_keys())
    },

    #' @description Keys in `uns` (unstructured annotation).
    #' @return A character vector.
    uns_keys = function() {
      private$.check_valid()
      reticulate::py_to_r(private$.py_mudata$uns_keys())
    },

    # ------------------------------------------------------------------
    # Update / sync methods
    # ------------------------------------------------------------------

    #' @description
    #' Update global `obs` / `var` indices from individual modalities.
    #' Equivalent to Python `mdata.update()`.
    #' @return `self` invisibly.
    update = function() {
      private$.check_valid()
      private$.py_mudata$update()
      invisible(self)
    },

    #' @description
    #' Update global `obs_names` from individual modalities.
    #' Equivalent to Python `mdata.update_obs()`.
    #' @return `self` invisibly.
    update_obs = function() {
      private$.check_valid()
      private$.py_mudata$update_obs()
      invisible(self)
    },

    #' @description
    #' Update global `var_names` from individual modalities.
    #' Equivalent to Python `mdata.update_var()`.
    #' @return `self` invisibly.
    update_var = function() {
      private$.check_valid()
      private$.py_mudata$update_var()
      invisible(self)
    },

    # ------------------------------------------------------------------
    # Pull / push obs & var columns
    # ------------------------------------------------------------------

    #' @description
    #' Copy columns from modality-level `.obs` to the global `.obs`.
    #' Wraps Python `mdata.pull_obs(...)`.
    #' @param ... Arguments forwarded to the Python method.
    #' @return `self` invisibly.
    pull_obs = function(...) {
      private$.check_valid()
      private$.py_mudata$pull_obs(...)
      invisible(self)
    },

    #' @description
    #' Copy columns from modality-level `.var` to the global `.var`.
    #' Wraps Python `mdata.pull_var(...)`.
    #' @param ... Arguments forwarded to the Python method.
    #' @return `self` invisibly.
    pull_var = function(...) {
      private$.check_valid()
      private$.py_mudata$pull_var(...)
      invisible(self)
    },

    #' @description
    #' Push columns from global `.obs` back to individual modalities.
    #' Wraps Python `mdata.push_obs(...)`.
    #' @param ... Arguments forwarded to the Python method.
    #' @return `self` invisibly.
    push_obs = function(...) {
      private$.check_valid()
      private$.py_mudata$push_obs(...)
      invisible(self)
    },

    #' @description
    #' Push columns from global `.var` back to individual modalities.
    #' Wraps Python `mdata.push_var(...)`.
    #' @param ... Arguments forwarded to the Python method.
    #' @return `self` invisibly.
    push_var = function(...) {
      private$.check_valid()
      private$.py_mudata$push_var(...)
      invisible(self)
    },

    # ------------------------------------------------------------------
    # copy
    # ------------------------------------------------------------------

    #' @description
    #' Create a deep copy of this MuData.
    #' @return A new [ReticulateMuData] wrapping the Python copy.
    copy = function() {
      private$.check_valid()
      .wrap_py_return(private$.py_mudata$copy())
    },

    # ------------------------------------------------------------------
    # py_call – generic Python method forwarding
    # ------------------------------------------------------------------

    #' @description
    #' Invoke any Python method on the underlying MuData and return the result
    #' with automatic wrapping:
    #' * `anndata.AnnData` → `anndataR::ReticulateAnnData`
    #' * `mudata.MuData`   → `ReticulateMuData`
    #' * anything else     → returned as-is (Python object, `convert = FALSE`)
    #'
    #' @param method A string: the Python method name.
    #' @param ... Arguments forwarded to the Python method.
    #' @return The (possibly wrapped) result.
    py_call = function(method, ...) {
      private$.check_valid()
      fn     <- reticulate::py_get_attr(private$.py_mudata, method)
      result <- fn(...)
      .wrap_py_return(result)
    },

    # ------------------------------------------------------------------
    # print
    # ------------------------------------------------------------------

    #' @description Print a human-readable summary.
    #' @param ... Ignored.
    print = function(...) {
      private$.check_valid()
      bi   <- reticulate::import_builtins(convert = FALSE)
      mods <- reticulate::py_to_r(bi$list(private$.py_mudata$mod$keys()))
      shp  <- self$shape

      cat("ReticulateMuData\n")
      cat(sprintf("  n_obs x n_vars : %d x %d\n", shp[1L], shp[2L]))
      cat(sprintf("  modalities (%d) : %s\n", length(mods), paste(mods, collapse = ", ")))
      for (nm in mods) {
        py_ad <- private$.py_mudata$mod[[nm]]
        n_o   <- as.integer(reticulate::py_to_r(
          reticulate::py_get_attr(py_ad, "n_obs")
        ))
        n_v   <- as.integer(reticulate::py_to_r(
          reticulate::py_get_attr(py_ad, "n_vars")
        ))
        cat(sprintf("    %-14s %d obs x %d vars\n", paste0(nm, ":"), n_o, n_v))
      }
      obs_c <- tryCatch(self$obs_keys(), error = function(e) character(0))
      if (length(obs_c) > 0L) {
        cat(sprintf("  obs keys       : %s\n", paste(obs_c, collapse = ", ")))
      }
      obsm_c <- tryCatch(self$obsm_keys(), error = function(e) character(0))
      if (length(obsm_c) > 0L) {
        cat(sprintf("  obsm keys      : %s\n", paste(obsm_c, collapse = ", ")))
      }
      uns_c <- tryCatch(self$uns_keys(), error = function(e) character(0))
      if (length(uns_c) > 0L) {
        cat(sprintf("  uns keys       : %s\n", paste(uns_c, collapse = ", ")))
      }
      invisible(self)
    }
  )
)


# ---------------------------------------------------------------------------
# S3 methods for [[ and [[<-
# ---------------------------------------------------------------------------

#' Extract a modality from a ReticulateMuData
#'
#' @description
#' Returns the named modality as an `anndataR::ReticulateAnnData` wrapper.
#' The invariant **"any AnnData extracted from MuData is a ReticulateAnnData"**
#' is enforced here.
#'
#' @param x A [ReticulateMuData].
#' @param i A single string: the modality name (e.g. `"rna"`).
#'
#' @return An `anndataR::ReticulateAnnData` object.
#'
#' @examples
#' \dontrun{
#' rna <- mdata[["rna"]]
#' inherits(rna, "ReticulateAnnData")  # TRUE
#' }
#'
#' @export
`[[.ReticulateMuData` <- function(x, i) {
  py_adata <- x$py_mudata()$mod[[i]]
  .wrap_py_return(py_adata)
}

#' Assign a modality into a ReticulateMuData
#'
#' @description
#' Assigns an AnnData modality into the MuData.  `value` may be either an
#' `anndataR::ReticulateAnnData` (its underlying Python object is extracted)
#' or a raw Python `anndata.AnnData` object.  No data is copied.
#'
#' @param x A [ReticulateMuData].
#' @param i A single string: the modality name.
#' @param value An `anndataR::ReticulateAnnData` **or** a Python
#'   `anndata.AnnData` object.
#'
#' @return `x` (invisibly, for chaining; the Python-side mutation happens
#'   in-place).
#'
#' @examples
#' \dontrun{
#' mdata[["rna"]] <- new_rna_reticulateanndata
#' mdata[["rna"]] <- raw_python_anndata_object
#' }
#'
#' @export
`[[<-.ReticulateMuData` <- function(x, i, value) {
  # Unwrap ReticulateAnnData to access the raw Python AnnData.
  py_val <- if (inherits(value, "ReticulateAnnData")) {
    value$py_anndata()
  } else {
    value
  }
  # __setitem__ mutates the Python mapping in-place.
  x$py_mudata()$mod$`__setitem__`(i, py_val)
  x
}


# ---------------------------------------------------------------------------
# Exported helper functions
# ---------------------------------------------------------------------------

#' Test whether an object is a ReticulateMuData
#'
#' @param x Any R object.
#' @return `TRUE` if `x` inherits from `"ReticulateMuData"`, `FALSE` otherwise.
#'
#' @examples
#' is_reticulate_mudata(42)        # FALSE
#'
#' @export
is_reticulate_mudata <- function(x) {
  inherits(x, "ReticulateMuData")
}

#' Coerce to ReticulateMuData
#'
#' @description
#' If `x` is already a [ReticulateMuData] it is returned unchanged.  If it is
#' a Python `mudata.MuData` object it is wrapped.  Otherwise an error is raised.
#'
#' @param x A [ReticulateMuData] or a Python `mudata.MuData` object.
#' @return A [ReticulateMuData].
#'
#' @examples
#' \dontrun{
#' mdata2 <- as_reticulate_mudata(py_mdata_obj)
#' }
#'
#' @export
as_reticulate_mudata <- function(x) {
  if (inherits(x, "ReticulateMuData")) return(x)
  ReticulateMuData$new(x)
}

#' Constructor alias for ReticulateMuData
#'
#' @description
#' Thin alias for `ReticulateMuData$new(py_obj)`.  Equivalent to calling the
#' constructor directly but reads more naturally in a pipeline.
#'
#' @param py_obj A Python `mudata.MuData` object.
#' @return A [ReticulateMuData].
#'
#' @examples
#' \dontrun{
#' mdata <- reticulate_mudata(py_mdata_from_python)
#' }
#'
#' @export
reticulate_mudata <- function(py_obj) {
  ReticulateMuData$new(py_obj)
}
