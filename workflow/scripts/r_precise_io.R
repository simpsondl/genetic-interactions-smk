# High-precision R workspace IO utilities
# Use version=3 RDS without compression and set numeric options for precision.

save_workspace <- function(obj, path, digits = 17, scipen = 999) {
  old_digits <- getOption("digits")
  old_scipen <- getOption("scipen")
  options(digits = digits, scipen = scipen)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, path, compress = FALSE, version = 3)
  options(digits = old_digits, scipen = old_scipen)
  invisible(path)
}

load_workspace <- function(path, digits = 17, scipen = 999) {
  options(digits = digits, scipen = scipen)
  readRDS(path)
}

validate_workspace <- function(ws, required_fields = character()) {
  missing <- setdiff(required_fields, names(ws))
  if (length(missing) > 0) {
    stop(sprintf("Workspace is missing required fields: %s", paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}
