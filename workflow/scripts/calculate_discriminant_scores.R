library(readr)
library(dplyr)

source("scripts/helper_functions.R")

gi_scores <- read_tsv(snakemake@input[["input_gi_scores"]])
gene_scores <- read_tsv(snakemake@input[["input_gene_level_scores"]])

message(sprintf("[%s] calculate_discriminant_scores.R starting", Sys.time()))
message(sprintf("[%s] Inputs: gi_scores=%s (%d rows), gene_scores=%s (%d rows)",
                Sys.time(), snakemake@input[["input_gi_scores"]], nrow(gi_scores),
                snakemake@input[["input_gene_level_scores"]], nrow(gene_scores)))

message(sprintf("[%s] Assessing variance of SGC scores (assess_sgcscore_variance)", Sys.time()))
gene_gis_var <- assess_sgcscore_variance(gi_scores, gene_scores)

message(sprintf("[%s] assess_sgcscore_variance returned %d rows", Sys.time(), nrow(gene_gis_var)))

message(sprintf("[%s] Calculating Discriminant.score = -log10(Variance.p) * abs(InteractionScore)", Sys.time()))
gene_gis_var$Discriminant.score <- -log10(gene_gis_var$Variance.p) * abs(gene_gis_var$InteractionScore)

message(sprintf("[%s] Writing discriminant scores to %s", 
                Sys.time(), snakemake@output[["output_discriminant_scores"]]))
write_tsv(gene_gis_var, snakemake@output[["output_discriminant_scores"]])