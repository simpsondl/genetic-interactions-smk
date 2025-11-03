library(readr)
library(dplyr)

source("scripts/helper_functions.R")
source("scripts/r_precise_io.R")

# Dual logging to both console and log file when running under Snakemake
if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
    source("scripts/dual_logging.R")
    .dual_cleanup <- setup_dual_logging(snakemake@log[[1]])
    on.exit({ 
        .dual_cleanup() 
    }, add = TRUE)
}

reference <- NULL
treated <- NULL
# pheno names passed by Snakemake params (reference_pheno, treated_pheno)
reference_name <- if (!is.null(snakemake@params[["reference_pheno"]])) {
                    snakemake@params[["reference_pheno"]]
                  } else {
                    "Gamma"
                  }
treated_name <- if (!is.null(snakemake@params[["treated_pheno"]])) {
                    snakemake@params[["treated_pheno"]]
                } else {
                    "Tau"
                }   

if (!is.null(snakemake@input[["input_reference_workspace"]])) {
    message(sprintf("[%s] Loading %s workspace: %s", 
                    Sys.time(), reference_name, snakemake@input[["input_reference_workspace"]]))
    reference_ws <- load_workspace(snakemake@input[["input_reference_workspace"]])
    validate_workspace(reference_ws, c("combined_all"))
    reference <- reference_ws$combined_all
} else {
    reference <- read_tsv(snakemake@input[["input_reference_gi_scores"]])
}
if (!is.null(snakemake@input[["input_treated_workspace"]])) {
    message(sprintf("[%s] Loading %s workspace: %s", 
                    Sys.time(), treated_name, snakemake@input[["input_treated_workspace"]]))
    treated_ws <- load_workspace(snakemake@input[["input_treated_workspace"]])
    validate_workspace(treated_ws, c("combined_all"))
    treated <- treated_ws$combined_all
} else {
    treated <- read_tsv(snakemake@input[["input_treated_gi_scores"]])
}


idmap <- read_tsv(snakemake@input[["input_idmap"]])
message(sprintf("[%s] calculate_differential_scores.R starting", Sys.time()))
message(sprintf("[%s] Inputs: %s=%s (%d rows), %s=%s (%d rows), idmap=%s (%d rows)",
                Sys.time(), reference_name, snakemake@input[["input_reference_gi_scores"]], nrow(reference),
                treated_name, snakemake@input[["input_treated_gi_scores"]], nrow(treated),
                snakemake@input[["input_idmap"]], nrow(idmap)))

message(sprintf("[%s] Computing construct-level differential scores (compute_construct_diff_scores)", Sys.time()))
nu <- compute_construct_diff_scores(reference, treated, reference_name, treated_name)
message(sprintf("[%s] compute_construct_diff_scores returned %d rows", Sys.time(), nrow(nu)))

message(sprintf("[%s] Aggregating to gene-level scores (compute_gene_interaction_scores)", Sys.time()))
nu_genes <- compute_gene_interaction_scores(nu, "GI.z")
message(sprintf("[%s] compute_gene_interaction_scores returned %d rows", Sys.time(), nrow(nu_genes)))

nu_genes <- left_join(nu_genes, idmap) %>% relocate(PseudogeneA:PseudogeneCombinationName)
message(sprintf("[%s] After joining idmap, gene-level table has %d rows and columns: %s", Sys.time(),
                nrow(nu_genes), paste(colnames(nu_genes), collapse = ", ")))

message(sprintf("[%s] Assessing variance for gene-level differential scores (assess_sgcscore_variance)", 
                Sys.time()))
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

# Save high-precision differential workspace for downstream hits calling
diff_ws <- list(
    meta = list(
        screen = snakemake@params[["screen"]],
        rep = snakemake@params[["rep"]],
        created = Sys.time(),
        rule = "calculate_differential_scores"
    ),
    differential_constructs = nu,
    differential_genes = nu_genes,
    differential_discriminant = nu_var
)

save_workspace(diff_ws, snakemake@output[["output_diff_workspace"]])
message(sprintf("[%s] Saved differential workspace to %s", 
                Sys.time(), snakemake@output[["output_diff_workspace"]]))