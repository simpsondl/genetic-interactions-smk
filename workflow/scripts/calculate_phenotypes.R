library(readr)
library(dplyr)

source("scripts/helper_functions.R")

counts <- read_tsv(snakemake@input[["input_counts"]])
to_filt <- read_tsv(snakemake@input[["input_filter_flags"]])

pc <- as.numeric(snakemake@params[["pseudocount"]])
norm <- as.logical(snakemake@params[["normalize"]])
doubs <- snakemake@params[["doublings"]]
cols_for_cond <- snakemake@params[["counts_cols"]]

message(sprintf("[%s] calculate_phenotypes.R starting; pseudocount=%s, normalize=%s", Sys.time(), pc, norm))
message(sprintf("[%s] doublings provided: %s", Sys.time(), paste(doubs, collapse=", ")))
message(sprintf("[%s] columns used for conditions: %s", Sys.time(), 
                paste(colnames(counts)[cols_for_cond], collapse=", ")))

# Define conditions in each arm
cond_df <- data.frame(Colname = colnames(counts[, cols_for_cond]),
                      Samplename = gsub("\\.R.*", "", colnames(counts[, cols_for_cond])),
                      Replicate = gsub(".*\\.","", colnames(counts[, cols_for_cond])))

# Calculate phenotypes for non-flagged constructs
keep_constructs <- to_filt$ConstructID[is.na(to_filt$Flag)]
message(sprintf("[%s] Number of constructs to process: %d", Sys.time(), length(keep_constructs)))
message(sprintf("[%s] Starting calculate_phenotypes on %d constructs", Sys.time(), length(keep_constructs)))
raw_phenotypes <- calculate_phenotypes(counts = counts[counts$ConstructID %in% keep_constructs,],
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
single_pheno <- calculate_single_sgRNA_phenotypes(raw_phenotypes)
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