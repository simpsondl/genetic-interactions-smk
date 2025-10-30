library(readr)
library(dplyr)

source("scripts/helper_functions.R")
source("scripts/r_precise_io.R")

# Dual logging to both console and log file when running under Snakemake
if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
    source("scripts/dual_logging.R")
    .dual_cleanup <- setup_dual_logging(snakemake@log[[1]])
    on.exit({ .dual_cleanup() }, add = TRUE)
}

gi_tbl <- NULL
if (!is.null(snakemake@input[["input_gi_workspace"]])) {
    message(sprintf("[%s] Loading GI workspace: %s", Sys.time(), snakemake@input[["input_gi_workspace"]]))
    gi_ws <- load_workspace(snakemake@input[["input_gi_workspace"]])
    validate_workspace(gi_ws, c("combined_all"))
    gi_tbl <- gi_ws$combined_all
} else {
    message(sprintf("[%s] Loading GI scores TSV: %s", Sys.time(), snakemake@input[["input_gi_scores"]]))
    gi_tbl <- read_tsv(snakemake@input[["input_gi_scores"]])
}
idmap <- read_tsv(snakemake@input[["input_idmap"]])

message(sprintf("[%s] calculate_gene_level_scores.R starting", Sys.time()))
message(sprintf("[%s] Input tables: gi_tbl_rows=%d, idmap_rows=%d",
                Sys.time(), nrow(gi_tbl), nrow(idmap)))
message(sprintf("[%s] Computing gene-level GI scores using score column 'GI.z'", Sys.time()))

gene_gis <- compute_gene_interaction_scores(gi_tbl, "GI.z")

message(sprintf("[%s] compute_gene_interaction_scores returned %d rows and %d unique PseudogeneCombinationID(s)",
                Sys.time(), nrow(gene_gis), length(unique(gene_gis$PseudogeneCombinationID))))

gene_gis <- left_join(gene_gis, idmap) %>%
    relocate(Pseudogene1:PseudogeneCombinationName)

write_tsv(gene_gis, snakemake@output[["output_gene_level_scores"]])
message(sprintf("[%s] Wrote gene-level GI scores to %s", 
                Sys.time(), snakemake@output[["output_gene_level_scores"]]))

# Save high-precision workspace for downstream rules
gene_level_ws <- list(
    meta = list(
        screen = snakemake@params[["screen"]],
        score = snakemake@params[["score"]],
        created = Sys.time(),
        rule = "calculate_gene_level_scores"
    ),
    gene_level_scores = gene_gis,
    combined_all = gi_tbl
)
save_workspace(gene_level_ws, snakemake@output[["output_gene_level_workspace"]])
message(sprintf("[%s] Saved gene-level workspace to %s", Sys.time(), snakemake@output[["output_gene_level_workspace"]]))