library(readr)
library(dplyr)
library(broom)

source("scripts/helper_functions.R")

# Redirect R output/messages to Snakemake log if provided
if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
    log_file <- snakemake@log[[1]]
    log_con <- file(log_file, open = "wt")
    sink(log_con, type = "output")
    sink(log_con, type = "message")
    options(warn = 1)
    on.exit({
        sink(type = "message")
        sink(type = "output")
        close(log_con)
    }, add = TRUE)
}

orientation_indep_phenotypes <- read_tsv(snakemake@input[["input_orientation_indep_phenotypes"]])
single_phenotypes <- read_tsv(snakemake@input[["input_single_sgRNA_phenotypes"]])

# single score for this Snakemake job
score <- snakemake@params[["score"]]

# output directory for this score
out_dir <- snakemake@output[["output_dir"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message(sprintf("[%s] calculate_gi_scores.R starting; score=%s", Sys.time(), score))
message(sprintf("[%s] Inputs: orientation_indep=%s (%d rows), single_phenotypes=%s (%d rows)",
                Sys.time(), snakemake@input[["input_orientation_indep_phenotypes"]], nrow(orientation_indep_phenotypes),
                snakemake@input[["input_single_sgRNA_phenotypes"]], nrow(single_phenotypes)))
message(sprintf("[%s] Writing individual per-sgRNA files into: %s", Sys.time(), out_dir))

# Pre-allocate lists to collect results for this score
n <- length(single_phenotypes$sgRNA.ID)
all_gis_list <- vector("list", n)
ests_list <- vector("list", n)
stats_list <- vector("list", n)

for (idx in seq_along(single_phenotypes$sgRNA.ID)) {
    i <- single_phenotypes$sgRNA.ID[idx]

    # Log progress every 100 sgRNAs
    if (idx %% 100 == 0) {
        message(sprintf("[%s] processing sgRNA %d/%d (id=%s)", score, idx, n, i))
    }

    gi_scores <- compute_gis(i, single_phenotypes, orientation_indep_phenotypes, score)

    # Handle cases where compute_gis returns NA (no data)
    if (length(gi_scores) == 1 && is.na(gi_scores)) {
        all_gis_list[[idx]] <- tibble()
        ests_list[[idx]] <- tibble()
        stats_list[[idx]] <- tibble()
        next
    }

    # Safe extraction of model objects
    model_obj <- gi_scores[[3]]
    all_gis_list[[idx]] <- gi_scores[[1]]
    if (!is.null(model_obj)) {
        ests_list[[idx]] <- data.frame(sgRNA.ID = i, tidy(model_obj))
        stats_list[[idx]] <- data.frame(sgRNA.ID = i, glance(model_obj))
    } else {
        ests_list[[idx]] <- tibble()
        stats_list[[idx]] <- tibble()
    }

    # Write individual results tables for this score
    write_tsv(all_gis_list[[idx]], file.path(out_dir, paste0(i, ".txt")))
}

# Combine and write out aggregated tables for this score
combined_all <- bind_rows(all_gis_list)
combined_ests <- bind_rows(ests_list)
combined_stats <- bind_rows(stats_list)

message(sprintf("[%s] Combined results: constructs=%d, model_est_rows=%d, model_stat_rows=%d", Sys.time(),
                                nrow(combined_all), nrow(combined_ests), nrow(combined_stats)))

message(sprintf("[%s] Writing combined_all to %s", Sys.time(), snakemake@output[["output_all_scores"]]))
write_tsv(combined_all, snakemake@output[["output_all_scores"]])
message(sprintf("[%s] Writing combined model estimates to %s", 
                Sys.time(), snakemake@output[["output_model_estimates"]]))
write_tsv(combined_ests, snakemake@output[["output_model_estimates"]])
message(sprintf("[%s] Writing combined model stats to %s", 
                Sys.time(), snakemake@output[["output_model_stats"]]))
write_tsv(combined_stats, snakemake@output[["output_model_stats"]])