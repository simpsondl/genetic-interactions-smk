library(readr)
library(dplyr)

source("scripts/helper_functions.R")

raw_phenotypes <- read_tsv(snakemake@input[["input_phenotypes"]])
orind_phenotypes <- read_tsv(snakemake@input[["input_orientation_indep_phenotypes"]])
single_phenotypes <- read_tsv(snakemake@input[["input_single_sgRNA_phenotypes"]])
to_filt <- read_tsv(snakemake@input[["input_filter_flags"]])

raw_phenotypes <- left_join(raw_phenotypes,
                            orind_phenotypes[, c("GuideCombinationID", "Gamma.OI.Avg", "Tau.OI.Avg")],
                            by = "GuideCombinationID")

# Identify single sgRNAs whose combinatorial and single phenotypes do not correlate
# Returns a list with four elements, first is a dataframe of all calculated correlations,
# second element is gamma sgRNAs which failed filter, and third is tau sgRNAs which failed filter
nocorr_filt <- filt_nocorrelation(combphenos = raw_phenotypes, 
                                  singlephenos = single_phenotypes,
                                  filterthresh = snakemake@params[["no_correlation_threshold"]])

# Pull out filteredphenotypes for interaction scores
phenos_filt <- raw_phenotypes[!(raw_phenotypes$FirstPosition %in% nocorr_filt[[2]]) &
                                  !(raw_phenotypes$SecondPosition %in% nocorr_filt[[2]]) &
                                    !(raw_phenotypes$FirstPosition %in% nocorr_filt[[3]]) &
                                      !(raw_phenotypes$SecondPosition %in% nocorr_filt[[3]]),]
single_pheno_filt <- single_phenotypes[!(single_phenotypes$sgRNA.ID %in% nocorr_filt[[2]]) &
                                    !(single_phenotypes$sgRNA.ID %in% nocorr_filt[[3]]),]

# Flag filtered constructs that were not previously filtered
to_filt$Flag[(to_filt$FirstPosition %in% nocorr_filt[[2]] | 
                to_filt$SecondPosition %in% nocorr_filt[[2]]) &
             is.na(to_filt$Flag)] <- "No correlation - Gamma phenotype"
# If an sgRNA was filtered due to having no gamma correlation, that sgRNA will not be flagged if
# it also had no tau correlation.
to_filt$Flag[(to_filt$FirstPosition %in% nocorr_filt[[3]] | 
                to_filt$SecondPosition %in% nocorr_filt[[3]]) &
             is.na(to_filt$Flag)] <- "No correlation - Tau phenotype"

# Save flags to file
write_tsv(to_filt, snakemake@output[["output_full_filter_flags"]])

# Save filtered phenotypes to file
write_tsv(phenos_filt, snakemake@output[["output_filtered_phenotypes"]])
write_tsv(single_pheno_filt, snakemake@output[["output_filtered_single_sgRNA_phenotypes"]])

# Save correlation results to file
write_tsv(nocorr_filt[[1]], snakemake@output[["output_correlation_results"]])