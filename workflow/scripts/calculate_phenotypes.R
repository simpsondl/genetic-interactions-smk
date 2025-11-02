library(readr)

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
reps <- snakemake@params[["replicates"]]

res <- resolve_count_prefixes(colnames(counts), snakemake@params[["counts_prefixes"]])
cols_for_cond <- res$indices
prefixes <- res$prefixes

message(sprintf("[%s] calculate_phenotypes.R starting; pseudocount=%s, normalize=%s", Sys.time(), pc, norm))
message(sprintf("[%s] doublings provided: %s", Sys.time(), paste(doubs, collapse = ", ")))
message(sprintf("[%s] columns used for conditions: %s", 
                Sys.time(), paste(colnames(counts)[cols_for_cond], collapse = ", ")))

# Define conditions in each arm
cond_df <- data.frame(Colname = colnames(counts)[cols_for_cond],
                      Samplename = gsub("\\.R.*", "", colnames(counts)[cols_for_cond]),
                      Replicate = gsub(".*\\.", "", colnames(counts)[cols_for_cond]))

# Calculate phenotypes for non-flagged constructs
keep_constructs <- to_filt$ConstructID[is.na(to_filt$Flag)]
message(sprintf("[%s] Starting calculate_phenotypes on %d constructs", Sys.time(), length(keep_constructs)))
raw_phenotypes <- calculate_phenotypes(counts = counts[counts$ConstructID %in% keep_constructs, ],
                                       conds = cond_df,
                                       pseudocount = pc,
                                       # do not normalize here; normalization moved to its own rule
                                       normalize = FALSE,
                                       doublings = doubs,
                                       phenotype_prefix_map = snakemake@params[["phenotype_prefix_map"]],
                                       replicates = reps,
                                       # do not compute replicate averages here; that is done downstream
                                       compute_avgs = FALSE)
message(sprintf("[%s] calculate_phenotypes returned %d rows", Sys.time(), nrow(raw_phenotypes)))
message(sprintf("[%s] Writing raw (per-replicate) phenotypes to %s", Sys.time(), snakemake@output[["output_phenotypes_raw"]]))
write_tsv(raw_phenotypes, snakemake@output[["output_phenotypes_raw"]])