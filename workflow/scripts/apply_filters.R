library(readr)
library(dplyr)

source("scripts/helper_functions.R")

counts <- read_tsv(snakemake@input[["input_counts"]], show_col_types = FALSE)

cols_for_cond <- snakemake@params[["counts_cols"]]

# Define conditions in each arm
cond.df <- data.frame(Colname = colnames(counts[,cols_for_cond]),
                      Samplename = gsub("\\.R.*", "", colnames(counts[,cols_for_cond])),
                      Replicate = gsub(".*\\.","", colnames(counts[,cols_for_cond])))
                      
# Identify individual sgRNAs which have low representation at T0 (med <= 35)
sgrna.filt <- filt_low_representation(counts = counts, 
                                      conds = cond.df, 
                                      filtersamp = "T0", 
                                      filterthresh = snakemake@params[["individual_sgRNA_median_threshold"]])

# Identify sgRNA combinations with low representation at T0 (count <= 50)
combination.filt <- filt_combinations(counts = counts, 
                                      conds = cond.df, 
                                      filtersamp = "T0", 
                                      filterthresh = snakemake@params[["combination_sgRNA_count_threshold"]])

# Keep track of filtered combinations and reason why
# These are the metadata columns, adjust if your data is different
to.filt <- counts[,1:11]
to.filt$Flag <- NA

# Flag combinations involving identified sgRNAs that have low median representation in either position
to.filt$Flag[to.filt$FirstPosition %in% sgrna.filt | to.filt$SecondPosition %in% sgrna.filt] <- 
  "Low median representation at T0"
# Flag individual combinations which have low representation at T0, 
# do not flag constructs which were previously removed
to.filt$Flag[to.filt$ConstructID %in% combination.filt & is.na(to.filt$Flag)] <- "Low representation at T0"

# Save
write_tsv(to.filt, snakemake@output[["output_filter_flags"]])