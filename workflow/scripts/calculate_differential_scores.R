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

gamma <- NULL
tau <- NULL
if (!is.null(snakemake@input[["input_gamma_workspace"]])) {
    message(sprintf("[%s] Loading Gamma workspace: %s", Sys.time(), snakemake@input[["input_gamma_workspace"]]))
    gamma_ws <- load_workspace(snakemake@input[["input_gamma_workspace"]])
    validate_workspace(gamma_ws, c("combined_all"))
    gamma <- gamma_ws$combined_all
} else {
    gamma <- read_tsv(snakemake@input[["input_gamma_gi_scores"]])
}
if (!is.null(snakemake@input[["input_tau_workspace"]])) {
    message(sprintf("[%s] Loading Tau workspace: %s", Sys.time(), snakemake@input[["input_tau_workspace"]]))
    tau_ws <- load_workspace(snakemake@input[["input_tau_workspace"]])
    validate_workspace(tau_ws, c("combined_all"))
    tau <- tau_ws$combined_all
} else {
    tau <- read_tsv(snakemake@input[["input_tau_gi_scores"]])
}
idmap <- read_tsv(snakemake@input[["input_idmap"]])
message(sprintf("[%s] calculate_differential_scores.R starting", Sys.time()))
message(sprintf("[%s] Inputs: gamma=%s (%d rows), tau=%s (%d rows), idmap=%s (%d rows)",
                Sys.time(), snakemake@input[["input_gamma_gi_scores"]], nrow(gamma),
                snakemake@input[["input_tau_gi_scores"]], nrow(tau),
                snakemake@input[["input_idmap"]], nrow(idmap)))

message(sprintf("[%s] Computing construct-level differential scores (compute_construct_diff_scores)", 
                Sys.time()))
nu <- compute_construct_diff_scores(gamma, tau)
message(sprintf("[%s] compute_construct_diff_scores returned %d rows", Sys.time(), nrow(nu)))

message(sprintf("[%s] Aggregating to gene-level scores (compute_gene_interaction_scores)", Sys.time()))
nu_genes <- compute_gene_interaction_scores(nu, "GI.z")
message(sprintf("[%s] compute_gene_interaction_scores returned %d rows", Sys.time(), nrow(nu_genes)))

nu_genes <- left_join(nu_genes, idmap) %>% relocate(Pseudogene1:PseudogeneCombinationName)
message(sprintf("[%s] After joining idmap, gene-level table has %d rows and columns: %s", Sys.time(),
                nrow(nu_genes), paste(colnames(nu_genes), collapse = ", ")))

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
message(sprintf("[%s] Saved differential workspace to %s", Sys.time(), snakemake@output[["output_diff_workspace"]]))