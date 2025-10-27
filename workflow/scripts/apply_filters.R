library(readr)
library(dplyr)

source("scripts/helper_functions.R")

counts <- read_tsv(snakemake@input[["input_counts"]], show_col_types = FALSE)

cols_for_cond <- snakemake@params[["counts_cols"]]

# Define conditions in each arm
cond_df <- data.frame(Colname = colnames(counts[,cols_for_cond]),
                      Samplename = gsub("\\.R.*", "", colnames(counts[,cols_for_cond])),
                      Replicate = gsub(".*\\.","", colnames(counts[,cols_for_cond])))
                      
# Log: columns used for condition parsing
message(sprintf("[%s] Columns used for conditions: %s", Sys.time(), 
        paste(colnames(counts)[cols_for_cond], collapse=", ")))

# Identify individual sgRNAs which have low representation at T0 (med <= threshold)
message(sprintf("[%s] Starting filt_low_representation (filtersamp=T0, threshold=%s)",
                Sys.time(), snakemake@params[["individual_sgRNA_median_threshold"]]))
sgrna_filt <- filt_low_representation(counts = counts, 
                                      conds = cond_df, 
                                      filtersamp = "T0", 
                                      filterthresh = snakemake@params[["individual_sgRNA_median_threshold"]])
message(sprintf("[%s] filt_low_representation returned %d sgRNA ids", Sys.time(), length(sgrna_filt)))

# Identify sgRNA combinations with low representation at T0 (count <= threshold)
message(sprintf("[%s] Starting filt_combinations (filtersamp=T0, threshold=%s)",
                Sys.time(), snakemake@params[["combination_sgRNA_count_threshold"]]))
combination_filt <- filt_combinations(counts = counts, 
                                      conds = cond_df, 
                                      filtersamp = "T0", 
                                      filterthresh = snakemake@params[["combination_sgRNA_count_threshold"]])
message(sprintf("[%s] filt_combinations returned %d construct ids", Sys.time(), length(combination_filt)))

# Keep track of filtered combinations and reason why
# These are the metadata columns, adjust if your data is different
to_filt <- counts[, 1:13]
to_filt$Flag <- NA

# Log the columns included in to_filt
message(sprintf("[%s] to_filt columns: %s", Sys.time(), paste(colnames(to_filt), collapse = ", ")))

# Flag combinations involving identified sgRNAs that have low median representation in either position
to_filt$Flag[to_filt$FirstPosition %in% sgrna_filt | to_filt$SecondPosition %in% sgrna_filt] <- 
  "Low median representation at T0"
# Flag individual combinations which have low representation at T0, 
# do not flag constructs which were previously removed
to_filt$Flag[to_filt$ConstructID %in% combination_filt & is.na(to_filt$Flag)] <- "Low representation at T0"

## Summarize flag counts (including non-flagged)
summary_counts <- to_filt %>%
  mutate(Flag = ifelse(is.na(Flag), "NOT_FLAGGED", Flag)) %>%
  group_by(Flag) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

message(sprintf("[%s] Flag summary:\n%s", Sys.time(), 
        paste(apply(summary_counts, 1, function(r) paste0(r[1], ": ", r[2])), collapse = "\n")))

# Save
write_tsv(to_filt, snakemake@output[["output_filter_flags"]])