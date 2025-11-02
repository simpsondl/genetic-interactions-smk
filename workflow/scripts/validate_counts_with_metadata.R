library(readr)
library(dplyr)

source("scripts/r_precise_io.R")

# Dual logging when running under Snakemake
if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
  source("scripts/dual_logging.R")
  .dual_cleanup <- setup_dual_logging(snakemake@log[[1]])
  on.exit({ .dual_cleanup() }, add = TRUE)
}

provided_path <- snakemake@input[["input_counts"]]

message(sprintf("[%s] validate_counts_with_metadata.R starting", Sys.time()))
message(sprintf("[%s] Provided counts: %s", Sys.time(), provided_path))

read_header <- function(path) {
  # if zip, find matching entry
  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    zlist <- utils::unzip(path, list = TRUE)
    target <- zlist$Name[grepl("_raw_counts\\.tsv$", zlist$Name, ignore.case = TRUE)]
    if (length(target) == 0) stop(sprintf("No raw counts TSV found inside zip: %s", path))
    tbl <- read_tsv(unz(path, target[1]), n_max = 0, show_col_types = FALSE)
  } else {
    tbl <- read_tsv(path, n_max = 0, show_col_types = FALSE)
  }
  colnames(tbl)
}

provided_cols <- read_header(provided_path)

expected_cols <- snakemake@params[["expected_count_columns"]]
message(sprintf("[%s] Using expected columns from snakemake params (EXPECTED_COUNT_COLUMNS)",
                Sys.time()))


# Normalize expected_cols: accept lists or comma-separated strings
if (!is.null(expected_cols) && !is.character(expected_cols)) {
  expected_cols <- unlist(expected_cols)
}
if (!is.null(expected_cols) && is.character(expected_cols) && length(expected_cols) == 1 && grepl(",", expected_cols)) {
  expected_cols <- strsplit(expected_cols, ",")[[1]]
  expected_cols <- trimws(expected_cols)
}

if (is.null(expected_cols) || length(expected_cols) == 0) {
  stop("Expected columns must be provided by config EXPECTED_COUNT_COLUMNS")
}

# Validate: provided must contain expected columns as its leading prefix (same names and order)
n_expected <- length(expected_cols)
if (length(provided_cols) < n_expected) {
  stop(sprintf("Provided counts file has fewer columns (%d) than expected (%d)", 
                length(provided_cols), n_expected))
}

prefix <- provided_cols[1:n_expected]
if (!all(prefix == expected_cols)) {
  diffs <- which(prefix != expected_cols)
  msg_lines <- sapply(diffs, 
                      function(i) {
                        sprintf("[%s] col %d: expected '%s' but got '%s'", 
                                Sys.time(), i, expected_cols[i], prefix[i])
                      })
  stop(paste(c("Provided counts columns do not match expected columns (prefix mismatch):", msg_lines), collapse = "\n"))
}

message(sprintf("[%s] Validation passed.", Sys.time()))

# Touch the marker output
out_path <- snakemake@output[["output_validation_marker"]]
file_connection <- file(out_path, "w")
writeLines(sprintf("[%s] validated: %s", Sys.time(), provided_path), con = file_connection)
close(file_connection)
message(sprintf("[%s] Wrote validation marker to %s", Sys.time(), out_path))
