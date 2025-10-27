library(readr)
library(dplyr)

source("scripts/helper_functions.R")

gamma <- read_tsv(snakemake@input[["input_gamma_gi_scores"]])
tau <- read_tsv(snakemake@input[["input_tau_gi_scores"]])
idmap <- read_tsv(snakemake@input[["input_idmap"]])

nu <- compute_construct_differential_scores(gamma, tau)

nu_genes <- compute_gene_interaction_scores(nu, "GI.z")
nu_genes <- left_join(nu_genes, idmap) %>% relocate(Pseudogene1:PseudogeneCombinationName)

nu_var <- assess_sgcscore_variance(nu, nu_genes)
nu_var$Discriminant <- -log10(nu_var$Variance.p) * abs(nu_var$InteractionScore)

write_tsv(nu, snakemake@output[["output_differential_scores"]])
write_tsv(nu_genes, snakemake@output[["output_gene_differential_scores"]])
write_tsv(nu_var, snakemake@output[["output_discriminant_differential_scores"]])