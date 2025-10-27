library(readr)
library(dplyr)

source("scripts/helper_functions.R")

orientation_indep_phenotypes <- read_tsv(snakemake@input[["input_orientation_indep_phenotypes"]])
single_sgRNA_phenotypes <- read_tsv(snakemake@input[["input_single_sgRNA_phenotypes"]])

scores <- snakemake@params[["scores"]]

# Initialize empty lists to store results
all_gis_list <- vector("list", length(single.pheno.filt$sgRNA.ID))
ests_list <- vector("list", length(single.pheno.filt$sgRNA.ID))
stats_list <- vector("list", length(single.pheno.filt$sgRNA.ID))

for(idx in seq_along(single.pheno.filt$sgRNA.ID)) {
    i <- single.pheno.filt$sgRNA.ID[idx]
    
    # Log progress every 100 sgRNAs
    if(idx %% 100 == 0){
        print(idx)
    }

    gi_scores <- compute_gis(i, single.pheno.filt, phenos.filt, "Gamma.R1")
    
    # Store results in lists
    all_gis_list[[idx]] <- gi_scores[[1]]
    ests_list[[idx]] <- data.frame(sgRNA.ID = i, tidy(gi_scores[[3]]))
    stats_list[[idx]] <- data.frame(sgRNA.ID = i, glance(gi_scores[[3]]))
    
    # Write individual results tables
    write_tsv(gi_scores[[1]], paste("Interaction_Scores/Gamma.R1/", i, ".txt", sep = ""))
}

# Combine all results at once
res <- list(
    bind_rows(all_gis_list),
    bind_rows(ests_list),
    bind_rows(stats_list))

