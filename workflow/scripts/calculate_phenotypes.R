library(readr)
library(dplyr)

source("scripts/helper_functions.R")

# Dual logging to both console and log file when running under Snakemake
if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
  source("scripts/dual_logging.R")
  .dual_cleanup <- setup_dual_logging(snakemake@log[[1]])
  on.exit({ 
    .dual_cleanup() 
  }, add = TRUE)
}

input_counts_path <- snakemake@input[["input_counts"]]

# Support either a plain TSV or a ZIP containing the TSV (expected to contain *_raw_counts.tsv)
if (grepl("\\.zip$", input_counts_path, ignore.case = TRUE)) {
    zlist <- utils::unzip(input_counts_path, list = TRUE)
    target_name <- zlist$Name[grepl("_raw_counts\\.tsv$", zlist$Name, ignore.case = TRUE)]
    if (length(target_name) == 0) {
        stop(sprintf("No file matching '*_raw_counts.tsv' found inside zip: %s", input_counts_path))
    }
    counts <- read_tsv(unz(input_counts_path, target_name[1]), show_col_types = FALSE)
} else {
    counts <- read_tsv(input_counts_path, show_col_types = FALSE)
}
to_filt <- read_tsv(snakemake@input[["input_filter_flags"]])

pc <- as.numeric(snakemake@params[["pseudocount"]])
norm <- as.logical(snakemake@params[["normalize"]])
doubs <- snakemake@params[["doublings"]]
cols_for_cond <- snakemake@params[["counts_cols"]]

message(sprintf("[%s] calculate_phenotypes.R starting; pseudocount=%s, normalize=%s", Sys.time(), pc, norm))
message(sprintf("[%s] doublings provided: %s", Sys.time(), paste(doubs, collapse = ", ")))
message(sprintf("[%s] columns used for conditions: %s", Sys.time(), 
                paste(colnames(counts)[cols_for_cond], collapse = ", ")))

# Define conditions in each arm
cond_df <- data.frame(Colname = colnames(counts[, cols_for_cond]),
                      Samplename = gsub("\\.R.*", "", colnames(counts[, cols_for_cond])),
                      Replicate = gsub(".*\\.", "", colnames(counts[, cols_for_cond])))

# Calculate phenotypes for non-flagged constructs
keep_constructs <- to_filt$ConstructID[is.na(to_filt$Flag)]
message(sprintf("[%s] Number of constructs to process: %d", Sys.time(), length(keep_constructs)))
message(sprintf("[%s] Starting calculate_phenotypes on %d constructs", Sys.time(), length(keep_constructs)))
raw_phenotypes <- calculate_phenotypes(counts = counts[counts$ConstructID %in% keep_constructs, ],
                                   conds = cond_df,
                                   pseudocount = pc,
                                   normalize = norm,
                                   doublings = doubs)
message(sprintf("[%s] calculate_phenotypes returned %d rows", Sys.time(), nrow(raw_phenotypes)))

# Calculate orientation independent, averaged phenotypes
message(sprintf("[%s] Starting calculate_averaged_phenotypes (orientation-independent)", Sys.time()))
orind_phenotypes <- calculate_averaged_phenotypes(phenos = raw_phenotypes)
message(sprintf("[%s] calculate_averaged_phenotypes returned %d rows", Sys.time(), nrow(orind_phenotypes)))

# Calculate single sgRNA phenotypes
message(sprintf("[%s] Starting calculate_single_sgRNA_phenotypes", Sys.time()))
single_pheno <- calculate_single_sgrna_phenotypes(raw_phenotypes)
message(sprintf("[%s] calculate_single_sgRNA_phenotypes returned %d rows", Sys.time(), nrow(single_pheno)))

# Save phenotype tables
message(sprintf("[%s] Writing raw phenotypes to %s", Sys.time(), snakemake@output[["output_phenotypes"]]))
write_tsv(raw_phenotypes, snakemake@output[["output_phenotypes"]])
message(sprintf("[%s] Writing orientation-independent phenotypes to %s", Sys.time(), 
                snakemake@output[["output_orientation_indep_phenotypes"]]))
write_tsv(orind_phenotypes, snakemake@output[["output_orientation_indep_phenotypes"]])
message(sprintf("[%s] Writing single-sgRNA phenotypes to %s", Sys.time(), 
                snakemake@output[["output_single_sgRNA_phenotypes"]]))
write_tsv(single_pheno, snakemake@output[["output_single_sgRNA_phenotypes"]])