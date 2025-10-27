library(readr)
library(dplyr)

source("scripts/helper_functions.R")

gamma <- read_tsv(snakemake@input[["input_gamma_gi_scores"]])
tau <- read_tsv(snakemake@input[["input_tau_gi_scores"]])
idmap <- read_tsv(snakemake@input[["input_idmap"]])
message(sprintf("[%s] calculate_differential_scores.R starting", Sys.time()))
message(sprintf("[%s] Inputs: gamma=%s (%d rows), tau=%s (%d rows), idmap=%s (%d rows)",
                Sys.time(), snakemake@input[["input_gamma_gi_scores"]], nrow(gamma),
                snakemake@input[["input_tau_gi_scores"]], nrow(tau),
                snakemake@input[["input_idmap"]], nrow(idmap)))

message(sprintf("[%s] Computing construct-level differential scores (compute_construct_differential_scores)", Sys.time()))
nu <- compute_construct_differential_scores(gamma, tau)
message(sprintf("[%s] compute_construct_differential_scores returned %d rows", Sys.time(), nrow(nu)))

message(sprintf("[%s] Aggregating to gene-level scores (compute_gene_interaction_scores)", Sys.time()))
nu_genes <- compute_gene_interaction_scores(nu, "GI.z")
message(sprintf("[%s] compute_gene_interaction_scores returned %d rows", Sys.time(), nrow(nu_genes)))

nu_genes <- left_join(nu_genes, idmap) %>% relocate(Pseudogene1:PseudogeneCombinationName)
message(sprintf("[%s] After joining idmap, gene-level table has %d rows and columns: %s", Sys.time(),
                nrow(nu_genes), paste(colnames(nu_genes), collapse=", ")))

message(sprintf("[%s] Assessing variance for gene-level differential scores (assess_sgcscore_variance)", Sys.time()))
nu_var <- assess_sgcscore_variance(nu, nu_genes)
message(sprintf("[%s] assess_sgcscore_variance returned %d rows", Sys.time(), nrow(nu_var)))

message(sprintf("[%s] Calculating Discriminant = -log10(Variance.p) * abs(InteractionScore)", Sys.time()))
nu_var$Discriminant <- -log10(nu_var$Variance.p) * abs(nu_var$InteractionScore)

message(sprintf("[%s] Writing outputs: differential=%s, gene_differential=%s, discriminant=%s", Sys.time(),
                snakemake@output[["output_differential_scores"]], snakemake@output[["output_gene_differential_scores"]],
                snakemake@output[["output_discriminant_differential_scores"]]))

write_tsv(nu, snakemake@output[["output_differential_scores"]])
write_tsv(nu_genes, snakemake@output[["output_gene_differential_scores"]])
write_tsv(nu_var, snakemake@output[["output_discriminant_differential_scores"]])