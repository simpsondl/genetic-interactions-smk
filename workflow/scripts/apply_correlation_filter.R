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

raw_phenotypes <- read_tsv(snakemake@input[["input_phenotypes"]])
orind_phenotypes <- read_tsv(snakemake@input[["input_orientation_indep_phenotypes"]])
single_phenotypes <- read_tsv(snakemake@input[["input_single_sgRNA_phenotypes"]])
to_filt <- read_tsv(snakemake@input[["input_filter_flags"]])

# Determine phenotypes to process (passed from Snakemake params)
phenos_all <- snakemake@params[["phenotypes"]]
if (is.null(phenos_all) || length(phenos_all) == 0) {
  # default order if not provided
  phenos_all <- c("Gamma", "Tau", "Rho")
}

# Optional mapping from config telling which phenotypes should have the
# no-correlation filter applied. If provided, it should be a named list of
# booleans, e.g. list(Gamma=TRUE, Tau=TRUE, Rho=FALSE). The order of
# phenos_all determines the processing order.
no_corr_map <- snakemake@params[["no_correlation_filter"]]
phenos_to_process <- character(0)
if (!is.null(no_corr_map) && length(no_corr_map) > 0) {
  # If user provided a character vector, treat it as an ordered list of
  # phenotypes to filter (preserves user order).
  if (is.character(no_corr_map)) {
    phenos_to_process <- no_corr_map[no_corr_map %in% phenos_all]
  } 
} else {
  # If no mapping provided, default to processing all phenos in phenos_all
  phenos_to_process <- phenos_all
}

# Determine replicates (if provided) to identify expected OI columns (R1, R2, ...)
replicates <- snakemake@params[["replicates"]]
if (is.null(replicates)) {
  # default fallback if not provided
  replicates <- c("R1", "R2")
}

# Build the list of orientation-independent columns for ALL phenotypes (not just filtered ones)
# This ensures Rho OI columns are joined even if Rho isn't in NO_CORRELATION_FILTER
expected_oi_cols <- c("GuideCombinationID")
for (p in phenos_all) {
  for (r in replicates) expected_oi_cols <- c(expected_oi_cols, paste0(p, ".OI.", r))
  expected_oi_cols <- c(expected_oi_cols, paste0(p, ".OI.Avg"))
}

# Fail fast if expected columns are missing
missing_cols <- setdiff(expected_oi_cols, colnames(orind_phenotypes))
if (length(missing_cols) > 0) {
  stop(sprintf("apply_correlation_filter: missing expected orientation-independent columns: %s",
               paste(missing_cols, collapse = ", ")))
}

# Join only the required OI columns into the raw phenotype table
raw_phenotypes <- left_join(raw_phenotypes,
                            orind_phenotypes[, expected_oi_cols],
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
message(sprintf("[%s] Starting per-phenotype correlation filtering (threshold=%s)", 
                Sys.time(), snakemake@params[["no_correlation_threshold"]]))

# Get filter mode from params (default to avg_only)
filter_mode <- snakemake@params[["correlation_filter_mode"]]
if (is.null(filter_mode)) {
  filter_mode <- "avg_only"
}
message(sprintf("[%s] Correlation filter mode: %s", Sys.time(), filter_mode))

# We'll collect correlation tables and summaries for all phenotypes
all_corrs <- list()
all_summaries <- list()
pheno_filtered <- raw_phenotypes
single_pheno_filtered <- single_phenotypes

for (p in phenos_to_process) {
  message(sprintf("[%s] Processing phenotype: %s", Sys.time(), p))

  res <- filt_nocorrelation(combphenos = pheno_filtered,
                            singlephenos = single_pheno_filtered,
                            phenotype = p,
                            filter_mode = filter_mode,
                            filterthresh = snakemake@params[["no_correlation_threshold"]])

  cors_tbl <- res[[1]]
  failing_sgrnas <- res[[2]]
  summary_stats <- res[[3]]

  message(sprintf("[%s] %s: correlation table rows=%d; failing sgRNAs=%d", 
                  Sys.time(), p, nrow(cors_tbl), length(failing_sgrnas)))
  if (length(failing_sgrnas) > 0) {
    message(sprintf("[%s] Example %s failures: %s", 
                    Sys.time(), p, paste(head(failing_sgrnas, 5), collapse = ", ")))
  }
  
  # Log summary statistics
  message(sprintf("[%s] %s summary by column:", Sys.time(), p))
  for (i in seq_len(nrow(summary_stats))) {
    message(sprintf("[%s]   %s: %d/%d failed (%.1f percent)", 
                    Sys.time(), summary_stats$Column[i], 
                    summary_stats$N_Failed[i], summary_stats$N_Tested[i],
                    summary_stats$Pct_Failed[i]))
  }

  # Attach phenotype label and store
  # Rename correlation columns to generic names to allow rbind across phenotypes
  cor_cols <- grep("^.+\\.OI\\.", colnames(cors_tbl), value = TRUE)
  for (j in seq_along(cor_cols)) {
    old_name <- cor_cols[j]
    new_name <- paste0("Correlation_", j)
    colnames(cors_tbl)[colnames(cors_tbl) == old_name] <- new_name
  }
  cors_tbl$Phenotype <- p
  cors_tbl$TestedColumns <- paste(cor_cols, collapse = ", ")
  summary_stats$Phenotype <- p
  all_corrs[[p]] <- cors_tbl
  all_summaries[[p]] <- summary_stats

  # Pull out filtered phenotypes for interaction scores for this phenotype
  pheno_filtered <- pheno_filtered[!(pheno_filtered$FirstPosition %in% failing_sgrnas) &
                                     !(pheno_filtered$SecondPosition %in% failing_sgrnas), ]
  single_pheno_filtered <- single_pheno_filtered[!(single_pheno_filtered$sgRNA.ID %in% failing_sgrnas), ]

  # Update global filter flags (preserve existing non-NA flags)
  flag_msg <- sprintf("No correlation - %s phenotype - no %sGI scores", p, tolower(p))
  to_filt$Flag[(to_filt$FirstPosition %in% failing_sgrnas |
                to_filt$SecondPosition %in% failing_sgrnas) &
               is.na(to_filt$Flag)] <- flag_msg

  # Write per-phenotype outputs if the output keys exist (they will if rule was generated dynamically)
  out_constructs_key <- paste0("output_filtered_", tolower(p), "_phenotypes")
  out_single_key <- paste0("output_filtered_", tolower(p), "_single_sgRNA_phenotypes")

  if (!is.null(snakemake@output[[out_constructs_key]])) {
    write_tsv(pheno_filtered, snakemake@output[[out_constructs_key]])
    message(sprintf("[%s] Wrote %s filtered constructs to %s", Sys.time(), p, snakemake@output[[out_constructs_key]]))
  } else {
    message(sprintf("[%s] No snakemake output configured for %s constructs (key=%s); skipping write", 
                    Sys.time(), p, out_constructs_key))
  }

  if (!is.null(snakemake@output[[out_single_key]])) {
    write_tsv(single_pheno_filtered, snakemake@output[[out_single_key]])
    message(sprintf("[%s] Wrote %s filtered single-sgRNA phenotypes to %s", 
                    Sys.time(), p, snakemake@output[[out_single_key]]))
  } else {
    message(sprintf("[%s] No snakemake output configured for %s single sgRNA (key=%s); skipping write", 
                     Sys.time(), p, out_single_key))
  }
}

# Combine all per-phenotype correlation tables and summaries
combined_corrs <- do.call(rbind, all_corrs)
combined_summaries <- do.call(rbind, all_summaries)
message(sprintf("[%s] Combined correlation table rows=%d", Sys.time(), nrow(combined_corrs)))
message(sprintf("[%s] Combined summary table rows=%d", Sys.time(), nrow(combined_summaries)))

# Summarize flags after correlation filter
flag_summary <- to_filt %>% 
                  mutate(Flag = ifelse(is.na(Flag), "NOT_FLAGGED", Flag)) %>%
                  group_by(Flag) %>% 
                  summarise(n = n()) %>% 
                  arrange(desc(n))
message(sprintf("[%s] Flag summary after correlation filter:\n%s", 
                Sys.time(), paste(apply(flag_summary, 1, function(r) paste0(r[1], ": ", r[2])), collapse = "\n")))

# Save flags to file
write_tsv(to_filt, snakemake@output[["output_full_filter_flags"]])

# Save combined correlation results to file
write_tsv(combined_corrs, snakemake@output[["output_correlation_results"]])
message(sprintf("[%s] Wrote combined correlation results to %s", 
                Sys.time(), snakemake@output[["output_correlation_results"]]))

# Save combined summary statistics to file (if output exists)
if (!is.null(snakemake@output[["output_correlation_summary"]])) {
  write_tsv(combined_summaries, snakemake@output[["output_correlation_summary"]])
  message(sprintf("[%s] Wrote correlation summary statistics to %s", 
                  Sys.time(), snakemake@output[["output_correlation_summary"]]))
}