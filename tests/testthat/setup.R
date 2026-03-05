# tests/testthat/setup.R
# Configure Python environment for tests.
# Uses the mofaflex_test conda environment if available, falling back to
# r-reticulate or whatever reticulate auto-discovers.

local({
  conda_bin <- Sys.which("conda")
  if (nchar(conda_bin) == 0L) {
    # Try miniforge3 path common on this system
    conda_bin <- "/home/maximiliann/miniforge3/bin/conda"
  }

  if (file.exists(conda_bin)) {
    py_path <- "/home/maximiliann/miniforge3/envs/mofaflex_test/bin/python"
    if (file.exists(py_path)) {
      Sys.setenv(RETICULATE_PYTHON = py_path)
    }
  }
})
