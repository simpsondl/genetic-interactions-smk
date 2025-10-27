library(readr)
library(dplyr)

source("scripts/helper_functions.R")

gi_scores <- read_tsv(snakemake@input[["input_gi_scores"]])
idmap <- read_tsv(snakemake@input[["input_idmap"]])

message(sprintf("[%s] calculate_gene_level_scores.R starting", Sys.time()))
message(sprintf("[%s] Input files: gi_scores=%s (%d rows), idmap=%s (%d rows)",
                Sys.time(), snakemake@input[["input_gi_scores"]], nrow(gi_scores),
                snakemake@input[["input_idmap"]], nrow(idmap)))
message(sprintf("[%s] Computing gene-level GI scores using score column 'GI.z'", Sys.time()))

gene_gis <- compute_gene_interaction_scores(gi_scores, "GI.z")

message(sprintf("[%s] compute_gene_interaction_scores returned %d rows and %d unique GeneCombinationID(s)",
                Sys.time(), nrow(gene_gis), length(unique(gene_gis$GeneCombinationID))))

gene_gis <- left_join(gene_gis, idmap) %>%
    relocate(Pseudogene1:PseudogeneCombinationName)

write_tsv(gene_gis, snakemake@output[["output_gene_level_scores"]])
message(sprintf("[%s] Wrote gene-level GI scores to %s", 
                Sys.time(), snakemake@output[["output_gene_level_scores"]]))