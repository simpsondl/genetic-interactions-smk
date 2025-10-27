library(readr)
library(dplyr)

source("scripts/helper_functions.R")

gi_scores <- read_tsv(snakemake@input[["input_gi_scores"]])
gene_scores <- read_tsv(snakemake@input[["input_gene_level_scores"]])

gene_gis_var <- assess_sgcscore_variance(gi_scores, gene_scores)

gene_gis_var$Discriminant.score <- -log10(gene_gis_var$Variance.p) * abs(gene_gis_var$InteractionScore)

write_tsv(gene_gis_var, snakemake@output[["output_discriminant_scores"]])