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

# Helper to get a named input or fallback to first input
get_input <- function(name=NULL){
  if(!is.null(name) && name %in% names(snakemake@input)) return(snakemake@input[[name]])
  return(snakemake@input[[1]])
}
# Helper to get a named output or fallback to first output
get_output <- function(name=NULL){
  if(!is.null(name) && name %in% names(snakemake@output)) return(snakemake@output[[name]])
  return(snakemake@output[[1]])
}

infile <- get_input("input_gi_scores")
outfile <- get_output("output_hits")

message(sprintf("[%s] call_hits.R starting; infile=%s, outfile=%s", Sys.time(), infile, outfile))

tab <- read_tsv(infile, show_col_types = FALSE)

# Determine discriminant column name (case-insensitive match for 'discriminant')
discr_cols <- grep("discriminant", colnames(tab), ignore.case = TRUE, value = TRUE)
if(length(discr_cols) == 0){
  stop("No discriminant column found in input (expected name containing 'discriminant')")
}
# prefer exact 'Discriminant.score' or 'Discriminant' if present
if("Discriminant.score" %in% discr_cols) discr_col <- "Discriminant.score"
else if("Discriminant" %in% discr_cols) discr_col <- "Discriminant"
else discr_col <- discr_cols[1]

message(sprintf("[%s] Using discriminant column: %s", Sys.time(), discr_col))

# threshold parameter (quantile)
# Accept either 'hit_threshold' or 'threshold' param from the rule
if("hit_threshold" %in% names(snakemake@params)){
  hit_threshold <- as.numeric(snakemake@params[["hit_threshold"]])
} else if("threshold" %in% names(snakemake@params)){
  hit_threshold <- as.numeric(snakemake@params[["threshold"]])
} else {
  stop("No hit threshold provided in snakemake params (expected 'hit_threshold' or 'threshold')")
}
if(is.na(hit_threshold) || hit_threshold <= 0 || hit_threshold >= 1){
  stop("hit_threshold must be a numeric between 0 and 1 (exclusive), e.g. 0.995)")
}
message(sprintf("[%s] Using quantile threshold: %s", Sys.time(), hit_threshold))

# Baseline: discriminant scores for categories X+NT or NT+NT
if(!"Category" %in% colnames(tab)){
  warning("No 'Category' column found; using all rows to compute quantile")
  baseline_vals <- tab[[discr_col]]
} else {
  baseline_vals <- tab %>% filter(Category %in% c("X+NT", "NT+NT")) %>% pull(!!sym(discr_col))
  if(length(baseline_vals) == 0){
    warning("No rows with Category in (X+NT, NT+NT) found; falling back to all rows for quantile calculation")
    baseline_vals <- tab[[discr_col]]
  }
}

# compute quantile
qval <- as.numeric(quantile(baseline_vals, probs = hit_threshold, na.rm = TRUE))
message(sprintf("[%s] Computed discriminant quantile (prob=%s) = %s", Sys.time(), hit_threshold, qval))

# Add Hit logical column
tab$Hit <- tab[[discr_col]] > qval

# Report counts
n_hits <- sum(tab$Hit, na.rm = TRUE)
n_total <- nrow(tab)
message(sprintf("[%s] Hit summary: %d/%d (%.2f%%)", Sys.time(), n_hits, n_total, 100 * n_hits / n_total))

# Write out
message(sprintf("[%s] Writing output to %s", Sys.time(), outfile))
write_tsv(tab, outfile)
message(sprintf("[%s] call_hits.R completed", Sys.time()))
