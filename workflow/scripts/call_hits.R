library(readr)
library(dplyr)

# Redirect output/messages to Snakemake log if provided
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

infile <- snakemake@input[["input_gi_scores"]]
outfile <- snakemake@output[["output_hits"]]

message(sprintf("[%s] call_hits.R starting; infile=%s, outfile=%s", 
                Sys.time(), infile, outfile))

scores <- read_tsv(infile, show_col_types = FALSE)

discr_col <- "Discriminant.score"
message(sprintf("[%s] Using discriminant column: %s", Sys.time(), discr_col))

# threshold parameter (quantile)
# Accept either 'hit_threshold' or 'threshold' param from the rule
hit_threshold <- as.numeric(snakemake@params[["threshold"]])
if(is.na(hit_threshold) || hit_threshold <= 0 || hit_threshold >= 1){
  stop("hit_threshold must be a numeric between 0 and 1 (exclusive), e.g. 0.995)")
}
message(sprintf("[%s] Using quantile threshold: %s", Sys.time(), hit_threshold))

# Baseline: discriminant scores for categories X+NT or NT+NT
baseline_vals <- scores %>% filter(Category %in% c("X+NT", "NT+NT")) %>% pull(!!sym(discr_col))

# compute quantile
qval <- as.numeric(quantile(baseline_vals, probs = hit_threshold, na.rm = TRUE))
message(sprintf("[%s] Computed discriminant quantile (prob=%s) = %s", Sys.time(), hit_threshold, qval))

# Add Hit logical column
scores$Hit <- scores[[discr_col]] > qval

# Report counts
n_hits <- sum(scores$Hit, na.rm = TRUE)
n_total <- nrow(scores)
message(sprintf("[%s] Hit summary: %d/%d (%.2f%%)", Sys.time(), n_hits, n_total, 100 * n_hits / n_total))

# Write out
message(sprintf("[%s] Writing output to %s", Sys.time(), outfile))
write_tsv(scores, outfile)
message(sprintf("[%s] call_hits.R completed", Sys.time()))
