library(readr)
library(dplyr)

source("scripts/helper_functions.R")

# Dual logging to both console and log file when running under Snakemake
if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
  source("scripts/dual_logging.R")
  .dual_cleanup <- setup_dual_logging(snakemake@log[[1]])
  on.exit({ .dual_cleanup() }, add = TRUE)
}

raw_phenotypes <- read_tsv(snakemake@input[["input_phenotypes"]])
orind_phenotypes <- read_tsv(snakemake@input[["input_orientation_indep_phenotypes"]])
single_phenotypes <- read_tsv(snakemake@input[["input_single_sgRNA_phenotypes"]])
to_filt <- read_tsv(snakemake@input[["input_filter_flags"]])

raw_phenotypes <- left_join(raw_phenotypes,
                            orind_phenotypes[, c("GuideCombinationID", 
                                                 "Gamma.OI.R1", "Gamma.OI.R2", "Gamma.OI.Avg", 
                                                 "Tau.OI.R1", "Tau.OI.R2", "Tau.OI.Avg", 
                                                 "Rho.OI.R1", "Rho.OI.R2", "Rho.OI.Avg")],
                            by = "GuideCombinationID")

message(sprintf("[%s] apply_correlation_filter.R starting", Sys.time()))
message(sprintf("[%s] Inputs: phenotypes=%s, orientation_indep=%s, single=%s, filter_flags=%s",
                Sys.time(), snakemake@input[["input_phenotypes"]], 
                snakemake@input[["input_orientation_indep_phenotypes"]],
                snakemake@input[["input_single_sgRNA_phenotypes"]], 
                snakemake@input[["input_filter_flags"]]))
message(sprintf("[%s] Number of raw phenotype rows: %d", Sys.time(), nrow(raw_phenotypes)))
message(sprintf("[%s] Number of single-phenotype rows: %d", Sys.time(), nrow(single_phenotypes)))
message(sprintf("[%s] Number of filter-flag rows: %d", Sys.time(), nrow(to_filt)))

# Identify single sgRNAs whose combinatorial and single phenotypes do not correlate
# Returns a list with four elements, first is a dataframe of all calculated correlations,
# second element is gamma sgRNAs which failed filter, and third is tau sgRNAs which failed filter
message(sprintf("[%s] Starting filt_nocorrelation (threshold=%s)", Sys.time(), 
                snakemake@params[["no_correlation_threshold"]]))
nocorr_filt <- filt_nocorrelation(combphenos = raw_phenotypes,
                                  singlephenos = single_phenotypes,
                                  filterthresh = snakemake@params[["no_correlation_threshold"]])
message(sprintf("[%s] filt_nocorrelation completed: correlation table rows=%d", Sys.time(), nrow(nocorr_filt[[1]])))
message(sprintf("[%s] Number of Gamma sgRNAs failing correlation: %d", Sys.time(), length(nocorr_filt[[2]])))
message(sprintf("[%s] Number of Tau sgRNAs failing correlation: %d", Sys.time(), length(nocorr_filt[[3]])))

# show first few failing sgRNAs for quick inspection
if(length(nocorr_filt[[2]]) > 0) message(sprintf("[%s] Example Gamma failures: %s", 
                                               Sys.time(), paste(head(nocorr_filt[[2]], 5), collapse = ", ")))
if(length(nocorr_filt[[3]]) > 0) message(sprintf("[%s] Example Tau failures: %s", 
                                               Sys.time(), paste(head(nocorr_filt[[3]], 5), collapse = ", ")))

# Pull out filteredphenotypes for interaction scores
phenos_gamma_filt <- raw_phenotypes[!(raw_phenotypes$FirstPosition %in% nocorr_filt[[2]]) &
                                  !(raw_phenotypes$SecondPosition %in% nocorr_filt[[2]]), ]
single_gamma_pheno_filt <- single_phenotypes[!(single_phenotypes$sgRNA.ID %in% nocorr_filt[[2]]), ]


phenos_tau_filt <- raw_phenotypes[!(raw_phenotypes$FirstPosition %in% nocorr_filt[[2]]) &
                                  !(raw_phenotypes$SecondPosition %in% nocorr_filt[[2]]) &
                                    !(raw_phenotypes$FirstPosition %in% nocorr_filt[[3]]) &
                                      !(raw_phenotypes$SecondPosition %in% nocorr_filt[[3]]), ]
single_tau_pheno_filt <- single_phenotypes[!(single_phenotypes$sgRNA.ID %in% nocorr_filt[[2]]) &
                                    !(single_phenotypes$sgRNA.ID %in% nocorr_filt[[3]]), ]

# Flag filtered constructs that were not previously filtered
to_filt$Flag[(to_filt$FirstPosition %in% nocorr_filt[[2]] |
                to_filt$SecondPosition %in% nocorr_filt[[2]]) &
             is.na(to_filt$Flag)] <- "No correlation - Gamma phenotype - no gGI or tGI scores"
# If an sgRNA was filtered due to having no gamma correlation, that sgRNA will not be flagged if
# it also had no tau correlation.
to_filt$Flag[(to_filt$FirstPosition %in% nocorr_filt[[3]] |
                to_filt$SecondPosition %in% nocorr_filt[[3]]) &
             is.na(to_filt$Flag)] <- "No correlation - Tau phenotype - no tGI scores"

# Summarize flags after correlation filter
flag_summary <- to_filt %>% mutate(Flag = ifelse(is.na(Flag), "NOT_FLAGGED", Flag)) %>%
  group_by(Flag) %>% summarise(n = n()) %>% arrange(desc(n))
message(sprintf("[%s] Flag summary after correlation filter:\n%s", 
                Sys.time(), paste(apply(flag_summary, 1, function(r) paste0(r[1], ": ", r[2])), collapse = "\n")))

# Save flags to file
write_tsv(to_filt, snakemake@output[["output_full_filter_flags"]])

# Save filtered phenotypes to file
write_tsv(phenos_gamma_filt, snakemake@output[["output_filtered_gamma_phenotypes"]])
write_tsv(single_gamma_pheno_filt, snakemake@output[["output_filtered_gamma_single_sgRNA_phenotypes"]])
write_tsv(phenos_tau_filt, snakemake@output[["output_filtered_tau_phenotypes"]])
write_tsv(single_tau_pheno_filt, snakemake@output[["output_filtered_tau_single_sgRNA_phenotypes"]])


# Save correlation results to file
write_tsv(nocorr_filt[[1]], snakemake@output[["output_correlation_results"]])
message(sprintf("[%s] Wrote correlation results to %s", Sys.time(), snakemake@output[["output_correlation_results"]]))