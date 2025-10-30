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

gi_scores <- read_tsv(snakemake@input[["input_gi_scores"]])
gene_scores <- NULL
if (!is.null(snakemake@input[["input_gene_level_workspace"]])) {
    message(sprintf("[%s] Loading gene-level workspace: %s", 
                    Sys.time(), snakemake@input[["input_gene_level_workspace"]]))
    gl_ws <- load_workspace(snakemake@input[["input_gene_level_workspace"]])
    validate_workspace(gl_ws, c("gene_level_scores"))
    gene_scores <- gl_ws$gene_level_scores
} else {
    gene_scores <- read_tsv(snakemake@input[["input_gene_level_scores"]])
}

message(sprintf("[%s] calculate_discriminant_scores.R starting", Sys.time()))
message(sprintf("[%s] Inputs: gi_scores=%s (%d rows), gene_scores=%s (%d rows)",
                Sys.time(), snakemake@input[["input_gi_scores"]], nrow(gi_scores),
                snakemake@input[["input_gene_level_scores"]], nrow(gene_scores)))

message(sprintf("[%s] Assessing variance of SGC scores (assess_sgcscore_variance)", Sys.time()))
gene_gis_var <- assess_sgcscore_variance(gi_scores, gene_scores)

message(sprintf("[%s] assess_sgcscore_variance returned %d rows", Sys.time(), nrow(gene_gis_var)))

message(sprintf("[%s] Calculating Discriminant = -log10(Variance.p) * abs(InteractionScore)", Sys.time()))
gene_gis_var$Discriminant <- -log10(gene_gis_var$Variance.p) * abs(gene_gis_var$InteractionScore)

message(sprintf("[%s] Writing discriminant scores to %s", 
                Sys.time(), snakemake@output[["output_discriminant_scores"]]))
write_tsv(gene_gis_var, snakemake@output[["output_discriminant_scores"]])

# Save high-precision discriminant workspace for hit calling
disc_ws <- list(
    meta = list(
        screen = snakemake@params[["screen"]],
        score = snakemake@params[["score"]],
        created = Sys.time(),
        rule = "calculate_discriminant_scores"
    ),
    gene_gis_var = gene_gis_var
)
save_workspace(disc_ws, snakemake@output[["output_discriminant_workspace"]])
message(sprintf("[%s] Saved discriminant workspace to %s", 
                Sys.time(), snakemake@output[["output_discriminant_workspace"]]))