library(readr)
library(dplyr)

source("scripts/r_precise_io.R")

# Dual logging to both console and log file when running under Snakemake
if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
  source("scripts/dual_logging.R")
  .dual_cleanup <- setup_dual_logging(snakemake@log[[1]])
  on.exit({ 
    .dual_cleanup() 
  }, add = TRUE)
}

infile <- snakemake@input[["input_scores"]]
outfile <- snakemake@output[["output_hits"]]

message(sprintf("[%s] call_hits.R starting; infile=%s, outfile=%s", 
                Sys.time(), infile, outfile))

scores <- NULL
if (!is.null(snakemake@input[["input_discriminant_workspace"]])) {
  message(sprintf("[%s] Loading discriminant workspace: %s", 
                  Sys.time(), snakemake@input[["input_discriminant_workspace"]]))
  disc_ws <- load_workspace(snakemake@input[["input_discriminant_workspace"]])
  validate_workspace(disc_ws, c("gene_gis_var"))
  scores <- disc_ws$gene_gis_var
} else {
  scores <- read_tsv(infile, show_col_types = FALSE)
}

discr_col <- "Discriminant"
message(sprintf("[%s] Using discriminant column: %s", Sys.time(), discr_col))

# threshold parameter (quantile)
hit_threshold <- as.numeric(snakemake@params[["threshold"]])
if (is.na(hit_threshold) || hit_threshold < 0 || hit_threshold > 1) {
  stop("hit_threshold must be a numeric between 0 and 1 (inclusive), e.g. 0.995)")
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
message(sprintf("[%s] Hit summary: %d/%d (%.2f percent)", Sys.time(), n_hits, n_total, 100 * n_hits / n_total))

# Write out
message(sprintf("[%s] Writing output to %s", Sys.time(), outfile))
write_tsv(scores, outfile)
message(sprintf("[%s] call_hits.R completed", Sys.time()))

# Save high-precision hits workspace
hits_ws <- list(
  meta = list(
    screen = snakemake@params[["screen"]],
    score = snakemake@params[["score"]],
    created = Sys.time(),
    rule = "call_hits",
    threshold = hit_threshold,
    quantile = qval
  ),
  hits_with_calls = scores
)
save_workspace(hits_ws, snakemake@output[["output_hits_workspace"]])
message(sprintf("[%s] Saved hits workspace to %s", Sys.time(), snakemake@output[["output_hits_workspace"]]))
