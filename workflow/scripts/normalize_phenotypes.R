library(readr)

source("scripts/helper_functions.R")

# Dual logging when run under Snakemake
if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
  source("scripts/dual_logging.R")
  .dual_cleanup <- setup_dual_logging(snakemake@log[[1]])
  on.exit({ .dual_cleanup() }, add = TRUE)
}

input_path <- snakemake@input[["input_phenotypes"]]
message(sprintf("[%s] normalize_phenotypes.R reading %s", Sys.time(), input_path))
phenos <- read_tsv(input_path, show_col_types = FALSE)

# params
normalize_flag <- as.logical(snakemake@params[["normalize"]])
doublings <- snakemake@params[["doublings"]]
phenotype_prefix_map <- snakemake@params[["phenotype_prefix_map"]]

if (isTRUE(normalize_flag)) {
  message(sprintf("[%s] Normalization requested; running normalization", Sys.time()))
  phenos_out <- normalize_phenotypes(phenos, phenotype_prefix_map, doublings)
  message(sprintf("[%s] Normalization complete", Sys.time()))
} else {
  message(sprintf("[%s] Normalization disabled in config; copying through", Sys.time()))
  phenos_out <- phenos
}

message(sprintf("[%s] Writing normalized phenotypes to %s", Sys.time(), snakemake@output[["output_normalized_phenotypes"]]))
write_tsv(phenos_out, snakemake@output[["output_normalized_phenotypes"]])
