library(readr)

source("scripts/helper_functions.R")

# Dual logging when run under Snakemake
if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
  source("scripts/dual_logging.R")
  .dual_cleanup <- setup_dual_logging(snakemake@log[[1]])
  on.exit({ .dual_cleanup() }, add = TRUE)
}

input_path <- snakemake@input[["input_phenotypes_normalized"]]
message(sprintf("[%s] replicate_averaging.R reading %s", Sys.time(), input_path))
phenos <- read_tsv(input_path, show_col_types = FALSE)

phenotype_prefix_map <- snakemake@params[["phenotype_prefix_map"]]

# compute replicate averages (adds *.Avg columns)
message(sprintf("[%s] Computing replicate averages", Sys.time()))
phenos_avg <- compute_phenotype_averages(phenos, phenotype_prefix_map)

# compute orientation-independent phenotypes
message(sprintf("[%s] Computing orientation-independent phenotypes", Sys.time()))
orind <- calculate_averaged_phenotypes(phenos_avg)

# compute single-sgRNA phenotypes
message(sprintf("[%s] Computing single-sgRNA phenotypes", Sys.time()))
single_pheno <- calculate_single_sgrna_phenotypes(phenos_avg)

# write outputs
message(sprintf("[%s] Writing final phenotypes to %s", Sys.time(), snakemake@output[["output_phenotypes"]]))
write_tsv(phenos_avg, snakemake@output[["output_phenotypes"]])

message(sprintf("[%s] Writing orientation-independent phenotypes to %s", Sys.time(), snakemake@output[["output_orientation_indep_phenotypes"]]))
write_tsv(orind, snakemake@output[["output_orientation_indep_phenotypes"]])

message(sprintf("[%s] Writing single-sgRNA phenotypes to %s", Sys.time(), snakemake@output[["output_single_sgRNA_phenotypes"]]))
write_tsv(single_pheno, snakemake@output[["output_single_sgRNA_phenotypes"]])
