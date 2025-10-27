library(readr)
library(dplyr)

source("scripts/helper_functions.R")

counts <- read_tsv(snakemake@input[["input_counts"]])
to_filt <- read_tsv(snakemake@input[["input_filter_flags"]])

pc <- as.numeric(snakemake@params[["pseudocount"]])
norm <- as.logical(snakemake@params[["normalize"]])
doubs <- snakemake@params[["doublings"]]
cols_for_cond <- snakemake@params[["counts_cols"]]

# Define conditions in each arm
cond_df <- data.frame(Colname = colnames(counts[,cols_for_cond]),
                      Samplename = gsub("\\.R.*", "", colnames(counts[,cols_for_cond])),
                      Replicate = gsub(".*\\.","", colnames(counts[,cols_for_cond])))

# Calculate phenotypes for non-flagged constructs
raw_phenotypes <- calculate_phenotypes(counts = counts[counts$ConstructID %in% 
                                                        to_filt$ConstructID[is.na(to_filt$Flag)],],
                                   conds = cond_df,
                                   pseudocount = pc,
                                   normalize = norm,
                                   doublings = doubs)

# Calculate orientation independent, averaged phenotypes
orind_phenotypes <- calculate_averaged_phenotypes(phenos = raw_phenotypes)

# Calculate single sgRNA phenotypes
single_pheno <- calculate_single_sgRNA_phenotypes(raw_phenotypes)

# Save phenotype tables
write_tsv(raw_phenotypes, snakemake@output[["output_phenotypes"]])
write_tsv(orind_phenotypes, snakemake@output[["output_orientation_indep_phenotypes"]])
write_tsv(single_pheno, snakemake@output[["output_single_sgRNA_phenotypes"]])