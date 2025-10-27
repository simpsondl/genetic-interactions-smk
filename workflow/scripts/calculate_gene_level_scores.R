library(readr)
library(dplyr)

source("scripts/helper_functions.R")

gi_scores <- read_tsv(snakemake@input[["input_gi_scores"]])
idmap <- read_tsv(snakemake@input[["input_idmap"]])

gene_gis <- compute_gene_interaction_scores(gi_scores, "GI.z")

gene_gis <- left_join(gene_gis, idmap) %>%
    relocate(Pseudogene1:PseudogeneCombinationName)

write_tsv(gene_gis, snakemake@output[["output_gene_level_scores"]])