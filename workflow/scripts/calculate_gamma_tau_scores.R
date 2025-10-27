# Calculate gamma and tau interaction scores, per replicate and averaged
# Loop over all sgRNAs, modeling each individually to get scores
for(i in single.pheno.filt$sgRNA.ID){
  
  if(which(i == single.pheno.filt$sgRNA.ID) %% 100 == 0){
    print(which(i == single.pheno.filt$sgRNA.ID))
  }
  
  # Compute interactions for sgRNA i
  gamma.r1.gi <- compute_gis(i, single.pheno.filt, phenos.filt, "Gamma.R1")
  gamma.r2.gi <- compute_gis(i, single.pheno.filt, phenos.filt, "Gamma.R2")
  gamma.avg.gi <- compute_gis(i, single.pheno.filt, phenos.filt, "Gamma.Avg")
  gamma.oi.avg.gi <- compute_gis(i, single.pheno.filt, phenos.filt, "Gamma.OI.Avg")
  
  tau.r1.gi <- compute_gis(i, single.pheno.filt, phenos.filt, "Tau.R1")
  tau.r2.gi <- compute_gis(i, single.pheno.filt, phenos.filt, "Tau.R2")
  tau.avg.gi <- compute_gis(i, single.pheno.filt, phenos.filt, "Tau.Avg")
  tau.oi.avg.gi <- compute_gis(i, single.pheno.filt, phenos.filt, "Tau.OI.Avg")
  # Save all results together
  if(i == single.pheno.filt$sgRNA.ID[1]){
    gamma.r1.res <- create_interaction_result_df(gamma.r1.gi)
    gamma.r2.res <- create_interaction_result_df(gamma.r2.gi)
    gamma.avg.res <- create_interaction_result_df(gamma.avg.gi)
    gamma.oi.avg.res <- create_interaction_result_df(gamma.oi.avg.gi)
    tau.r1.res <- create_interaction_result_df(tau.r1.gi)
    tau.r2.res <- create_interaction_result_df(tau.r2.gi)
    tau.avg.res <- create_interaction_result_df(tau.avg.gi)
    tau.oi.avg.res <- create_interaction_result_df(tau.oi.avg.gi)
  } else {
    gamma.r1.res <- update_interaction_result_df(gamma.r1.gi, gamma.r1.res)
    gamma.r2.res <- update_interaction_result_df(gamma.r2.gi, gamma.r2.res)
    gamma.avg.res <- update_interaction_result_df(gamma.avg.gi, gamma.avg.res)
    gamma.oi.avg.res <- update_interaction_result_df(gamma.oi.avg.gi, gamma.oi.avg.res)
    tau.r1.res <- update_interaction_result_df(tau.r1.gi, tau.r1.res)
    tau.r2.res <- update_interaction_result_df(tau.r2.gi, tau.r2.res)
    tau.avg.res <- update_interaction_result_df(tau.avg.gi, tau.avg.res)
    tau.oi.avg.res <- update_interaction_result_df(tau.oi.avg.gi, tau.oi.avg.res)
  }
  # Write individual results tables
  write_tsv(gamma.r1.gi[[1]], paste("Interaction_Scores/Gamma.R1/", i, ".txt", sep = ""))
  write_tsv(gamma.r2.gi[[1]], paste("Interaction_Scores/Gamma.R2/", i, ".txt", sep = ""))
  write_tsv(gamma.avg.gi[[1]], paste("Interaction_Scores/Gamma.Avg/", i, ".txt", sep = ""))
  write_tsv(gamma.oi.avg.gi[[1]], paste("Interaction_Scores/Gamma.OI.Avg/", i, ".txt", sep = ""))
  write_tsv(tau.r1.gi[[1]], paste("Interaction_Scores/Tau.R1/", i, ".txt", sep = ""))
  write_tsv(tau.r2.gi[[1]], paste("Interaction_Scores/Tau.R2/", i, ".txt", sep = ""))
  write_tsv(tau.avg.gi[[1]], paste("Interaction_Scores/Tau.Avg/", i, ".txt", sep = ""))
  write_tsv(tau.oi.avg.gi[[1]], paste("Interaction_Scores/Tau.OI.Avg/", i, ".txt", sep = ""))
}

# Order table with all construct interaction scores by strength
gamma.r1.res[[1]] <- gamma.r1.res[[1]][order(gamma.r1.res[[1]]$GI.z),]
gamma.r2.res[[1]] <- gamma.r2.res[[1]][order(gamma.r2.res[[1]]$GI.z),]
gamma.avg.res[[1]] <- gamma.avg.res[[1]][order(gamma.avg.res[[1]]$GI.z),]
gamma.oi.avg.res[[1]] <- gamma.oi.avg.res[[1]][order(gamma.oi.avg.res[[1]]$GI.z),]
tau.r1.res[[1]] <- tau.r1.res[[1]][order(tau.r1.res[[1]]$GI.z),]
tau.r2.res[[1]] <- tau.r2.res[[1]][order(tau.r2.res[[1]]$GI.z),]
tau.avg.res[[1]] <- tau.avg.res[[1]][order(tau.avg.res[[1]]$GI.z),]
tau.oi.avg.res[[1]] <- tau.oi.avg.res[[1]][order(tau.oi.avg.res[[1]]$GI.z),]
# Save individual construct scores
con_path_pfx <- "Interaction_Scores/Compiled/Construct_Scores/"
write_tsv(gamma.r1.res[[1]], paste(con_path_pfx, "all_interaction_scores_gamma_r1.txt", sep = ""))
write_tsv(gamma.r2.res[[1]], paste(con_path_pfx, "all_interaction_scores_gamma_r2.txt", sep = ""))
write_tsv(gamma.avg.res[[1]], paste(con_path_pfx, "all_interaction_scores_gamma_avg.txt", sep = ""))
write_tsv(gamma.oi.avg.res[[1]], paste(con_path_pfx, "all_interaction_scores_gamma_oi_avg.txt", sep = ""))
write_tsv(tau.r1.res[[1]], paste(con_path_pfx, "all_interaction_scores_tau_r1.txt", sep = ""))
write_tsv(tau.r2.res[[1]], paste(con_path_pfx, "all_interaction_scores_tau_r2.txt", sep = ""))
write_tsv(tau.avg.res[[1]], paste(con_path_pfx, "all_interaction_scores_tau_avg.txt", sep = ""))
write_tsv(tau.oi.avg.res[[1]], paste(con_path_pfx, "all_interaction_scores_tau_oi_avg.txt", sep = ""))
# Save model estimate information
mod_path_pfx <- "Interaction_Scores/Compiled/Model_Results/"
write_tsv(gamma.r1.res[[2]], paste(mod_path_pfx, "model_estimates_gamma_r1.txt", sep = ""))
write_tsv(gamma.r2.res[[2]], paste(mod_path_pfx, "model_estimates_gamma_r2.txt", sep = ""))
write_tsv(gamma.avg.res[[2]], paste(mod_path_pfx, "model_estimates_gamma_avg.txt", sep = ""))
write_tsv(gamma.oi.avg.res[[2]], paste(mod_path_pfx, "model_estimates_gamma_oi_avg.txt", sep = ""))
write_tsv(tau.r1.res[[2]], paste(mod_path_pfx, "model_estimates_tau_r1.txt", sep = ""))
write_tsv(tau.r2.res[[2]], paste(mod_path_pfx, "model_estimates_tau_r2.txt", sep = ""))
write_tsv(tau.avg.res[[2]], paste(mod_path_pfx, "model_estimates_tau_avg.txt", sep = ""))
write_tsv(tau.oi.avg.res[[2]], paste(mod_path_pfx, "model_estimates_tau_oi_avg.txt", sep = ""))
# Save model statistic information
write_tsv(gamma.r1.res[[3]], paste(mod_path_pfx, "model_statistics_gamma_r1.txt", sep = ""))
write_tsv(gamma.r2.res[[3]], paste(mod_path_pfx, "model_statistics_gamma_r2.txt", sep = ""))
write_tsv(gamma.avg.res[[3]], paste(mod_path_pfx, "model_statistics_gamma_avg.txt", sep = ""))
write_tsv(gamma.oi.avg.res[[3]], paste(mod_path_pfx, "model_statistics_gamma_oi_avg.txt", sep = ""))
write_tsv(tau.r1.res[[3]], paste(mod_path_pfx, "model_statistics_tau_r1.txt", sep = ""))
write_tsv(tau.r2.res[[3]], paste(mod_path_pfx, "model_statistics_tau_r2.txt", sep = ""))
write_tsv(tau.avg.res[[3]], paste(mod_path_pfx, "model_statistics_tau_avg.txt", sep = ""))
write_tsv(tau.oi.avg.res[[3]], paste(mod_path_pfx, "model_statistics_tau_oi_avg.txt", sep = ""))

# Calculate guide combination (orientation-independent) interaction scores
gamma.r1.all.gis.orient <- compute_sgc_interaction_scores(gamma.r1.res[[1]])
gamma.r2.all.gis.orient <- compute_sgc_interaction_scores(gamma.r2.res[[1]])
gamma.avg.all.gis.orient <- compute_sgc_interaction_scores(gamma.avg.res[[1]])
gamma.oi.avg.all.gis.orient <- compute_sgc_interaction_scores(gamma.oi.avg.res[[1]])
tau.r1.all.gis.orient <- compute_sgc_interaction_scores(tau.r1.res[[1]])
tau.r2.all.gis.orient <- compute_sgc_interaction_scores(tau.r2.res[[1]])
tau.avg.all.gis.orient <- compute_sgc_interaction_scores(tau.avg.res[[1]])
tau.oi.avg.all.gis.orient <- compute_sgc_interaction_scores(tau.oi.avg.res[[1]])
# Calculate gene combination interaction scores
gamma.r1.gene.gis <- compute_gene_interaction_scores(gamma.r1.res[[1]], "GI.z")
gamma.r2.gene.gis <- compute_gene_interaction_scores(gamma.r2.res[[1]], "GI.z")
gamma.avg.gene.gis <- compute_gene_interaction_scores(gamma.avg.res[[1]], "GI.z")
gamma.oi.avg.gene.gis <- compute_gene_interaction_scores(gamma.oi.avg.res[[1]], "GI.z")
tau.r1.gene.gis <- compute_gene_interaction_scores(tau.r1.res[[1]], "GI.z")
tau.r2.gene.gis <- compute_gene_interaction_scores(tau.r2.res[[1]], "GI.z")
tau.avg.gene.gis <- compute_gene_interaction_scores(tau.avg.res[[1]], "GI.z")
tau.oi.avg.gene.gis <- compute_gene_interaction_scores(tau.oi.avg.res[[1]], "GI.z")

# Save guide combination (orientation-independent) interaction scores
sgc_path_pfx <- "Interaction_Scores/Compiled/GuideCombination_Scores/"
write_tsv(gamma.r1.all.gis.orient, paste(sgc_path_pfx, "guide_combination_interaction_scores_gamma_r1.txt", sep = ""))
write_tsv(gamma.r2.all.gis.orient, paste(sgc_path_pfx, "guide_combination_interaction_scores_gamma_r2.txt", sep = ""))
write_tsv(gamma.avg.all.gis.orient, paste(sgc_path_pfx, "guide_combination_interaction_scores_gamma_avg.txt", sep = ""))
write_tsv(gamma.oi.avg.all.gis.orient, paste(sgc_path_pfx, "guide_combination_interaction_scores_gamma_oi_avg.txt", sep = ""))
write_tsv(tau.r1.all.gis.orient, paste(sgc_path_pfx, "guide_combination_interaction_scores_tau_r1.txt", sep = ""))
write_tsv(tau.r2.all.gis.orient, paste(sgc_path_pfx, "guide_combination_interaction_scores_tau_r2.txt", sep = ""))
write_tsv(tau.avg.all.gis.orient, paste(sgc_path_pfx, "guide_combination_interaction_scores_tau_avg.txt", sep = ""))
write_tsv(tau.oi.avg.all.gis.orient, paste(sgc_path_pfx, "guide_combination_interaction_scores_tau_oi_avg.txt", sep = ""))

# Save gene combination interaction scores
gc_path_pfx <- "Interaction_Scores/Compiled/GeneCombination_Scores/"
write_tsv(gamma.r1.gene.gis, paste(gc_path_pfx, "gene_combination_interaction_scores_gamma_r1.txt", sep = ""))
write_tsv(gamma.r2.gene.gis, paste(gc_path_pfx, "gene_combination_interaction_scores_gamma_r2.txt", sep = ""))
write_tsv(gamma.avg.gene.gis, paste(gc_path_pfx, "gene_combination_interaction_scores_gamma_avg.txt", sep = ""))
write_tsv(gamma.oi.avg.gene.gis, paste(gc_path_pfx, "gene_combination_interaction_scores_gamma_oi_avg.txt", sep = ""))
write_tsv(tau.r1.gene.gis, paste(gc_path_pfx, "gene_combination_interaction_scores_tau_r1.txt", sep = ""))
write_tsv(tau.r2.gene.gis, paste(gc_path_pfx, "gene_combination_interaction_scores_tau_r2.txt", sep = ""))
write_tsv(tau.avg.gene.gis, paste(gc_path_pfx, "gene_combination_interaction_scores_tau_avg.txt", sep = ""))
write_tsv(tau.oi.avg.gene.gis, paste(gc_path_pfx, "gene_combination_interaction_scores_tau_oi_avg.txt", sep = ""))