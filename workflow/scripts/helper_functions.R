library(dplyr)
library(broom)
library(data.table)

### Function for identifying individual sgRNAs which are poorly represented
### in either position A or position B at T0
### Inputs:
### ### counts - raw count df, import from file
### ### conds - condition df, maps sample and replicate information to column names
### ### filtersamp - which arm of experiment to filter
### ### filterthresh - minimum median representation across all constructs/replicates
### Outputs:
### ### vector with sgRNA ids to with mean median representation below threshold
filt_low_representation <- function(counts, conds, filtersamp, filterthresh = 30) {
  # For each sgRNA in each position, get median count across all matching constructs
  # Expected output: 2 df with rows equal to number of sgRNAs and and 6 measurements across arms/replicates
  med_pos1 <- counts %>% group_by(FirstPosition) %>% summarise(across(conds$Colname, median))
  med_pos2 <- counts %>% group_by(SecondPosition) %>% summarise(across(conds$Colname, median))
  # For each sgRNA in each position, calculate mean median coverage per arm
  # Expected output: 2 df with rows equal to number of sgRNAs and 3 measurements across arms
  mean_pos1 <- as.data.frame(lapply(unique(conds$Samplename), 
                                    function(i) {
                                      apply(med_pos1[, colnames(med_pos1) %in% 
                                                        conds$Colname[conds$Samplename == i]], 
                                            1, mean)
                                    }))
  mean_pos2 <- as.data.frame(lapply(unique(conds$Samplename), 
                                    function(i) {
                                      apply(med_pos2[, colnames(med_pos2) %in% 
                                                        conds$Colname[conds$Samplename == i]], 
                                            1, mean)
                                    }))
  # Fix names from lapply
  colnames(mean_pos1) <- unique(conds$Samplename)
  colnames(mean_pos2) <- unique(conds$Samplename)
  # Set for return, only need pos1 since indexing is identical
  rownames(mean_pos1) <- med_pos1$FirstPosition
  # Return sgRNA ids which have mean median representation below threshold
  return(rownames(mean_pos1[mean_pos1[, filtersamp] <= filterthresh | 
                              mean_pos2[, filtersamp] <= filterthresh, ]))
}

### Function for identifying individual sgRNA combinations which are poorly represented
### Inputs:
### ### counts - raw count df, import from file
### ### conds - condition df, maps sample and replicate information to column names
### ### filtersamp - which arm of experiment to filter
### ### filterthresh - minimum median representation across all constructs/replicates
### Outputs:
### ### vector with sgRNA combination ids to with representation below threshold
filt_combinations <- function(counts, conds, filtersamp, filterthresh) {
  # Get column indices corresponding to arms to filter
  idx <- which(colnames(counts) %in% conds$Colname[conds$Samplename == filtersamp])
  # Create vector where nonzero entry corresponds to a count below threshold in a relevant arm
  tmp <- rowSums(counts[, idx] <= filterthresh)
  # Return guidecombinationids for sgRNA combinations to be filtered
  return(counts$ConstructID[which(tmp != 0)])
}

### Function for calculating individual sgRNA combination phenotypes
### Inputs:
### ### counts - raw count df, import from file
### ### conds - condition df, maps sample and replicate information to column names
### ### pseudocount - pseudocount to add to all raw count values, default 10
### ### normalize - normalize phenotypes to non-targeting controls and population doublings, default TRUE
### ### doublings - vector containing total population doublings in the following order:
### ### ### DMSO Replicate 1
### ### ### DMSO Replicate 2
### ### ### Niraparib Replicate 1
### ### ### Niraparib Replicate 2
### Outputs:
### ### df containing gamma, tau, and rho phenotypes for all combinations present
calculate_phenotypes <- function(counts, conds, pseudocount = 10, normalize = TRUE, doublings) {
  pseudo <- counts[, colnames(counts) %in% conds$Colname] + pseudocount
  fracs <- sweep(pseudo, 2, colSums(pseudo), FUN = "/")
  phenos <- data.frame(Gamma.R1 = log2(fracs$DMSO.R1 / fracs$T0.R1),
                       Gamma.R2 = log2(fracs$DMSO.R2 / fracs$T0.R2),
                       Tau.R1 = log2(fracs$NIRAP.R1 / fracs$T0.R1),
                       Tau.R2 = log2(fracs$NIRAP.R2 / fracs$T0.R2),
                       Rho.R1 = log2(fracs$NIRAP.R1 / fracs$DMSO.R1),
                       Rho.R2 = log2(fracs$NIRAP.R2 / fracs$DMSO.R2))
  phenos <- cbind(counts[, 1:13], phenos)
  
  if (normalize) {
    nt_gamma_r1 <- median(phenos$Gamma.R1[phenos$Category == "NT+NT"])
    nt_gamma_r2 <- median(phenos$Gamma.R2[phenos$Category == "NT+NT"])
    nt_tau_r1 <- median(phenos$Tau.R1[phenos$Category == "NT+NT"])
    nt_tau_r2 <- median(phenos$Tau.R2[phenos$Category == "NT+NT"])
    nt_rho_r1 <- median(phenos$Rho.R1[phenos$Category == "NT+NT"])
    nt_rho_r2 <- median(phenos$Rho.R2[phenos$Category == "NT+NT"])

    phenos$Gamma.R1 <- (phenos$Gamma.R1 - nt_gamma_r1) / doublings[1]
    phenos$Gamma.R2 <- (phenos$Gamma.R2 - nt_gamma_r2) / doublings[2]
    phenos$Tau.R1 <- (phenos$Tau.R1 - nt_tau_r1) / doublings[3]
    phenos$Tau.R2 <- (phenos$Tau.R2 - nt_tau_r2) / doublings[4]
    phenos$Rho.R1 <- (phenos$Rho.R1 - nt_rho_r1) / (doublings[1] - doublings[3])
    phenos$Rho.R2 <- (phenos$Rho.R2 - nt_rho_r2) / (doublings[2] - doublings[4])
  }

  phenos$Gamma.Avg <- rowMeans(phenos[, c("Gamma.R1", "Gamma.R2")])
  phenos$Tau.Avg <- rowMeans(phenos[, c("Tau.R1", "Tau.R2")])
  phenos$Rho.Avg <- rowMeans(phenos[, c("Rho.R1", "Rho.R2")])
  # rearranges dataframe to have averages next to individual replicates
  return(phenos[, c(1:15, 20, 16, 17, 21, 18, 19, 22)])
}

### Function for calculating orientation independent, averaged sgRNA combination phenotypes
### Inputs:
### ### phenos - output from calculate_phenotypes function, OR a df containing columns:
### ### ### Gamma.R1, Gamma.R2
### ### ### Tau.R1, Tau.R2
### ### ### GuideCombinationID with IDs beginning with sgc_
### Outputs:
### ### df containing orientation-independent phenotypes for all replicates and avg replicate
calculate_averaged_phenotypes <- function(phenos) {
  orind <- phenos %>% 
            group_by(GuideCombinationID) %>% 
              summarise(Gamma.OI.R1 = mean(Gamma.R1),
                        Gamma.OI.R2 = mean(Gamma.R2),
                        Gamma.OI.Avg = mean(Gamma.Avg),
                        Tau.OI.R1 = mean(Tau.R1),
                        Tau.OI.R2 = mean(Tau.R2),
                        Tau.OI.Avg = mean(Tau.Avg),
                        Rho.OI.R1 = mean(Rho.R1),
                        Rho.OI.R2 = mean(Rho.R2),
                        Rho.OI.Avg = mean(Rho.Avg),
                        N = n())
  orind <- orind[order(as.numeric(gsub("sgc_", "", orind$GuideCombinationID))), ]
  return(orind)
}

### Function for calculating phenotypes for individual sgRNAs using non-targeting guides
### Inputs:
### ### phenos - output from calculate_phenotypes function
### Outputs:
### ### df containing single sgRNA phenotypes by replicate
calculate_single_sgrna_phenotypes <- function(phenos) {
  # Create dataframe for saving results
  single_pheno <- data.frame("sgRNA.ID" = unique(c(phenos$FirstPosition, phenos$SecondPosition)), 
                             "Gamma.OI.R1" = 0, 
                             "Gamma.OI.R2" = 0,
                             "Gamma.OI.Avg" = 0,
                             "Tau.OI.R1" = 0,
                             "Tau.OI.R2" = 0,
                             "Tau.OI.Avg" = 0,
                             "Rho.OI.R1" = 0,
                             "Rho.OI.R2" = 0,
                             "Rho.OI.Avg" = 0,
                             "N" = 0)
  
  for (i in seq_len(nrow(single_pheno))){
    if (i %% 100 == 0) {
            progress_pct <- round((i / nrow(single_pheno)) * 100, 1)
            message(sprintf("[%s] Processing sgRNA %d/%d (%s percent) - ID: %s", 
                            Sys.time(), i, nrow(single_pheno), progress_pct, single_pheno$sgRNA.ID[i]))
        }

    # Handle non-targeting case
    if (grepl("non-targeting", single_pheno$sgRNA.ID[i])) {
      # Extract all non-targeting guide combinations with desired non-targeting guide in position A or B
      # If this is not handled explicitly, will aggregate all NT guides into a glob
      tmp <- phenos[(phenos$FirstPosition == single_pheno$sgRNA.ID[i] &
                                phenos$SecondPosition %in% unique(phenos$FirstPosition[phenos$Category == "NT+NT"])) | 
                               (phenos$SecondPosition == single_pheno$sgRNA.ID[i] &
                                  phenos$FirstPosition %in% unique(phenos$FirstPosition[phenos$Category == "NT+NT"])), ]
      single_pheno$Gamma.OI.R1[i] <- mean(tmp$Gamma.R1)
      single_pheno$Gamma.OI.R2[i] <- mean(tmp$Gamma.R2)
      single_pheno$Gamma.OI.Avg[i] <- mean(tmp$Gamma.Avg)
      single_pheno$Tau.OI.R1[i] <- mean(tmp$Tau.R1)
      single_pheno$Tau.OI.R2[i] <- mean(tmp$Tau.R2)
      single_pheno$Tau.OI.Avg[i] <- mean(tmp$Tau.Avg)
      single_pheno$Rho.OI.R1[i] <- mean(tmp$Rho.R1)
      single_pheno$Rho.OI.R2[i] <- mean(tmp$Rho.R2)
      single_pheno$Rho.OI.Avg[i] <- mean(tmp$Rho.Avg)
      single_pheno$N[i] <- nrow(tmp)
    }  else {     
      # Handle targeting case
      # Extract all combinations with desired targeting guide in position A or B and non-targeting guide in other
      # Easier to do this since we can just grab all non-targeting at once
      tmp <- phenos[(phenos$FirstPosition == single_pheno$sgRNA.ID[i] | 
                      phenos$SecondPosition == single_pheno$sgRNA.ID[i]) & 
                    (phenos$FirstPosition %in% unique(phenos$FirstPosition[phenos$Category == "NT+NT"]) | 
                      phenos$SecondPosition %in% unique(phenos$FirstPosition[phenos$Category == "NT+NT"])), ]
      single_pheno$Gamma.OI.R1[i] <- mean(tmp$Gamma.R1)
      single_pheno$Gamma.OI.R2[i] <- mean(tmp$Gamma.R2)
      single_pheno$Gamma.OI.Avg[i] <- mean(tmp$Gamma.Avg)
      single_pheno$Tau.OI.R1[i] <- mean(tmp$Tau.R1)
      single_pheno$Tau.OI.R2[i] <- mean(tmp$Tau.R2)
      single_pheno$Tau.OI.Avg[i] <- mean(tmp$Tau.Avg)
      single_pheno$Rho.OI.R1[i] <- mean(tmp$Rho.R1)
      single_pheno$Rho.OI.R2[i] <- mean(tmp$Rho.R2)
      single_pheno$Rho.OI.Avg[i] <- mean(tmp$Rho.Avg)
      single_pheno$N[i] <- nrow(tmp)
    }
  }
  return(single_pheno)
}

### Function for calculating correlations between single and combinatorial phenotypes for all sgRNAs
### Inputs:
### ### combphenos - output from calculate_phenotypes function
### ### singlephenos - output from calculate_single_sgRNA_phenotypes
### ### filterthresh - the minimum correlation to not flag an sgRNA for filtering, default 0.25
### Outputs:
### ### list containing four elements, in order:
### ### ### df of all calculated correlations
### ### ### sgRNA ids below threshold, gamma phenotypes
### ### ### sgRNA ids below threshold, tau phenotypes
### ### ### sgRNA ids below threshold, rho phenotypes
filt_nocorrelation <- function(combphenos, singlephenos, filterthresh = 0.25) {
  cors <- data.frame(sgRNA.ID = singlephenos$sgRNA.ID,
                     Gamma.OI.Correlation = NA,
                     Tau.OI.Correlation = NA,
                     Rho.OI.Correlation = NA)
  
  for (i in singlephenos$sgRNA.ID){
    tmp <- combphenos[combphenos$SecondPosition == i, c(colnames(combphenos)[1:9],
                                                      "Gamma.OI.Avg", "Tau.OI.Avg", "Rho.OI.Avg")]
    tmp <- merge(tmp, singlephenos[, c("sgRNA.ID", "Gamma.OI.Avg", "Tau.OI.Avg", "Rho.OI.Avg")],
                 by.x = "FirstPosition", by.y = "sgRNA.ID", suffix = c(".Comb", ".Single"))
    cors$Gamma.OI.Correlation[cors$sgRNA.ID == i] <- cor(tmp[, "Gamma.OI.Avg.Comb"], tmp[, "Gamma.OI.Avg.Single"])
    cors$Tau.OI.Correlation[cors$sgRNA.ID == i] <- cor(tmp[, "Tau.OI.Avg.Comb"], tmp[, "Tau.OI.Avg.Single"])
    cors$Rho.OI.Correlation[cors$sgRNA.ID == i] <- cor(tmp[, "Rho.OI.Avg.Comb"], tmp[, "Rho.OI.Avg.Single"])
  }
  
  return(list(cors,
              cors$sgRNA.ID[cors$Gamma.OI.Correlation < filterthresh],
              cors$sgRNA.ID[cors$Tau.OI.Correlation < filterthresh],
              cors$sgRNA.ID[cors$Rho.OI.Correlation < filterthresh]))
}

### Function for calculating interaction scores
### Inputs:
### ### query - sgRNA ID to calculate GI scores for
### ### singlepheno_df - filtered output from calculate_single_sgRNA_phenotypes and filt_nocorrelation
### ### pairpheno_df - filtered output from calculate_phenotypes and filt_nocorrelation
### ### phenocol - column name that contains the phenotypes to calculate GI scores from
### Outputs:
### ### list containing three elements, in order:
### ### ### df with all called interaction scores
### ### ### a vector containing expected GI scores determined from the model
### ### ### the model that was built
compute_gis <- function(query, singlepheno_df, pairpheno_df, phenocol) {
  ### Get the single sgRNA phenotypes
  tmp_data <- data.frame(query = query, 
                         sgRNA.id = singlepheno_df$sgRNA.ID, 
                         single = singlepheno_df[, phenocol])
  colnames(tmp_data)[3] <- "single"
  ### Merge single with interaction phenotypes
  tmp_data_merge <- merge(tmp_data, pairpheno_df[pairpheno_df$SecondPosition == query, 
                                                 c("ConstructID",  "GuideCombinationID", "Identical",
                                                   "FirstPseudogene", "SecondPseudogene", "PseudogeneCombinationID",
                                                   "Orientation", "GeneCombinationID", "Category",
                                                   "FirstPosition", phenocol)], 
                          by.x = "sgRNA.id", by.y = "FirstPosition", all.x = TRUE)
  
  ### Identify control combinations
  tmp_data_merge$Control <- grepl("non-targeting", tmp_data_merge$sgRNA.id)
  
  ### Remove any missing information
  tmp_data_merge <- tmp_data_merge[!is.na(tmp_data_merge[, phenocol]), ]
  
  if (nrow(tmp_data_merge) == 0) {
    return(NA)
  }

  try({
    ### Fit the model
    lm_fit <- lm(as.formula(paste0(phenocol, " ~ single + I(single^2)")), data = tmp_data_merge)
    
    ### Generate points for plotting the fit
    single_pred <- seq(-.35, .05, .01)
    pair_pred <- predict(lm_fit, list(single = single_pred))
    
    ### Define interaction score as observed - expected
    tmp_data_merge$Expected <- predict(lm_fit, list(single = tmp_data_merge$single))
    tmp_data_merge$GI <- lm_fit$residuals
    ### z-standardize GIs to negative controls
    tmp_data_merge$GI.z <- (tmp_data_merge$GI) / sd(tmp_data_merge$GI[tmp_data_merge$Control])
    
    ### Order dataframe
    tmp_data_merge <- tmp_data_merge[order(tmp_data_merge$GI), c("ConstructID", "GuideCombinationID", 
                                                                 "PseudogeneCombinationID", "GeneCombinationID", 
                                                                 "FirstPseudogene", "SecondPseudogene",
                                                                 "Orientation", "Identical", "Category", "Control",
                                                                 "sgRNA.id", "query", "single", phenocol, "Expected",
                                                                 "GI", "GI.z")]
    
    ### Want to label strongest 10 interactions [both buffering and synthetic lethal]
    n <- nrow(tmp_data_merge)
    tmp_data_merge$lb <- c(gsub("_.*", "", tmp_data_merge$sgRNA.id[1:(10)]),
                           rep("", n - 20),
                           gsub("_.*", "", tmp_data_merge$sgRNA.id[(n - 9):n]))
    
  }, silent = TRUE)
  
  return(list(tmp_data_merge, pair_pred, lm_fit))
}

### Function for cleanly saving results from compute_gis - used to establish save df's
### Inputs:
### ### gis - output from calculate_gis function
### Outputs:
### ### list containing three elements, in order:
### ### ### df with gi scores
### ### ### df with model coefficients
### ### ### df with model statistics
create_interaction_result_df <- function(gis) {
  all_gis <- gis[[1]]
  ests <- data.frame(sgRNA.ID = i, tidy(gis[[3]]))
  stats <- data.frame(sgRNA.ID = i, glance(gis[[3]]))

  return(list(all_gis, ests, stats))
}

### Function for cleanly saving results from compute_gis - used after create_interaction_result_df
### ### gis - output from calculate_gis function
### ### prev - output from either create_interaction_result_df or update_interaction_result_df
### Outputs:
### ### list containing three elements, in order:
### ### ### df with gi scores
### ### ### df with model coefficients
### ### ### df with model statistics
update_interaction_result_df <- function(gis, prev) {
  all_gis <- rbind(prev[[1]], gis[[1]])
  ests <- rbind(prev[[2]], data.frame(sgRNA.ID = i, tidy(gis[[3]])))
  stats <- rbind(prev[[3]], data.frame(sgRNA.ID = i, glance(gis[[3]])))
  
  return(list(all_gis, ests, stats))
}

### Function for aggregating the results from compute_gis across constructs
### Inputs:
### ### gis - final full output from update_interaction_result_df
### ### scorecol - column name that contains the interaction scores to be aggregated to gene level
### Outputs:
### ### df containing gene-level GI scores
compute_gene_interaction_scores <- function(gis, scorecol) {
  info_cols <- c("GuideCombinationID", "PseudogeneCombinationID", "Category", "Control", scorecol)
  gene_gis <- gis[, colnames(gis) %in% info_cols] %>% 
                group_by(PseudogeneCombinationID) %>% 
                mutate(InteractionScore = mean(!!sym(scorecol)), N = n()) %>%
                select(c(PseudogeneCombinationID, Category, InteractionScore, N)) %>% 
                unique()
  return(gene_gis)
}

### Helper: safely run a Wilcoxon test with group-size checks
### Inputs:
### ### df - data.frame or data.table containing the data to test
### ### score_col - column name (string) with numeric scores (default "GI.z")
### ### group_col - column name (string) with grouping indicator (logical/factor) (default "testdist")
### Output: list(pval = numeric or NA_real_, n_test = integer, n_control = integer)
safe_wilcox_test <- function(df, score_col = "GI.z", group_col = "testdist") {
  # Make sure group column is present
  if (!(group_col %in% colnames(df)) || !(score_col %in% colnames(df))) {
    return(list(pval = NA_real_, n_test = 0L, n_control = 0L))
  }

  if (length(unique(df[[group_col]])) < 2) {
    return(list(pval = NA_real_, n_test = 0L, n_control = 0L))
  }
  
  grp <- df[[group_col]]
  # Ensure logical vector for counting
  grp_logical <- as.logical(grp)
  n_test <- as.integer(sum(grp_logical, na.rm = TRUE))
  n_control <- as.integer(sum(!grp_logical, na.rm = TRUE))

  # If either group is empty, return NA
  if (n_test < 1L || n_control < 1L) {
    return(list(pval = NA_real_, n_test = n_test, n_control = n_control))
  }

  # Try the test and return NA on error
  wt <- try(stats::wilcox.test(as.formula(paste0(score_col, " ~ ", group_col)), data = df), silent = TRUE)
  pval <- if (inherits(wt, "try-error")) NA_real_ else wt$p.value
  return(list(pval = pval, n_test = n_test, n_control = n_control))
}

### Helper: build selection index for assess_sgcscore_variance
### Inputs:
### ### tmp - data.table subset keyed by SecondPseudogene
### ### i - PseudogeneCombinationID under test
### ### gene1, gene2 - pseudogene names
### Output: logical vector sel_idx of length nrow(tmp)
build_selection_index <- function(tmp, i, gene1, gene2) {
  # Determine whether gene names match the NTPG pattern
  is_ntpg1 <- grepl("NTPG", gene1)
  is_ntpg2 <- grepl("NTPG", gene2)

  if (is_ntpg1 && is_ntpg2) {
    sel_idx <- (tmp$PseudogeneCombinationID == i) |
                  (tmp$is_NTNT & (tmp$SecondPseudogene == gene1 | tmp$SecondPseudogene == gene2))
  } else if (is_ntpg1) {
    sel_idx <- (tmp$PseudogeneCombinationID == i) |
                  (tmp$is_NTNT & tmp$SecondPseudogene == gene1) | (tmp$is_XNT & tmp$SecondPseudogene == gene2)
  } else if (is_ntpg2) {
    sel_idx <- (tmp$PseudogeneCombinationID == i) |
                  (tmp$is_NTNT & tmp$SecondPseudogene == gene2) | (tmp$is_XNT & tmp$SecondPseudogene == gene1)
  } else {
    sel_idx <- (tmp$PseudogeneCombinationID == i) | tmp$is_XNT
  }

  # Ensure a logical vector is always returned
  if (length(sel_idx) == 0) {
    return(logical(0))
  }
  return(sel_idx)
}

### Function for assessing the variance of sgRNA-level GI scores that contribute to gene-level GI scores
### Inputs:
### ### congis - final full output from update_interaction_result_df
### ### genegis - output from compute_gene_interaction_scores, must have columns Gene1 and Gene2 included
### Outputs:
### ### df containing GI scores averaged across sgRNA IDs
assess_sgcscore_variance <- function(congis, genegis) {
  # Purpose: compute a Wilcoxon test comparing GI.z for constructs belonging to the
  # test pseudogene-combination vs relevant control constructs. This refactor
  # uses data.table for faster subsetting and clearer intent.
  
  # Convert to data.table (non-destructive to original inputs)
  dt_con <- data.table::as.data.table(congis)
  dt_gen <- data.table::as.data.table(genegis)

  # Precompute category flags to speed repeated checks
  dt_con[, `:=`(is_NTNT = Category == "NT+NT", is_XNT = Category == "X+NT")]

  # Key by SecondPseudogene so lookups by gene are fast
  data.table::setkey(dt_con, SecondPseudogene)

  # Prepare result columns on dt_gen (we'll fill vectors and assign once for speed)
  n_gen <- nrow(dt_gen)
  var_p <- rep(NA_real_, n_gen)
  var_n_nt <- integer(n_gen)
  var_n_test <- integer(n_gen)

  ids <- dt_gen$PseudogeneCombinationID
  total_ids <- length(ids)

  for (k in seq_along(ids)){
    i <- ids[k]
    if (k %% 5000 == 0) {
      message(sprintf("[%s] assess_sgcscore_variance: processed %d/%d combinations; current PseudogeneCombinationID=%s",
                      Sys.time(), k, total_ids, i))
    }

    row_gen <- dt_gen[PseudogeneCombinationID == i]
    if (nrow(row_gen) == 0) next
    gene1 <- row_gen$Pseudogene1[1]
    gene2 <- row_gen$Pseudogene2[1]

    # Fast subset: get all rows where SecondPseudogene is gene1 or gene2
    tmp <- dt_con[J(c(gene1, gene2)), nomatch = 0]
    if (nrow(tmp) == 0) next

    # Build selection index using helper to reduce cyclomatic complexity
    sel_idx <- build_selection_index(tmp, i, gene1, gene2)

    if (!any(sel_idx)) next
    tmp2 <- tmp[sel_idx, ]
    tmp2 <- tmp2[tmp2$Identical == FALSE, ]
    if (nrow(tmp2) == 0) next

    tmp2$testdist <- tmp2$PseudogeneCombinationID == i

    # Run Wilcoxon test safely via helper (returns p-value and group sizes)
    wres <- safe_wilcox_test(tmp2, score_col = "GI.z", group_col = "testdist")
    var_p[k] <- wres$pval
    var_n_nt[k] <- wres$n_control
    var_n_test[k] <- wres$n_test
  }

  # Attach computed results back onto dt_gen using data.table assignment for speed
  # Use clear column names so downstream code can read: Variance.p, Variance.N.NT, Variance.N.Test
  dt_gen[, `:=`(Variance.p = var_p,
                Variance.N.NT = var_n_nt,
                Variance.N.Test = var_n_test)]

  n_computed <- sum(!is.na(dt_gen$Variance.p))
  message(sprintf("[%s] assess_sgcscore_variance: completed; computed p-values for %d/%d combinations",
                  Sys.time(), n_computed, total_ids))

  # Return the same type as input (data.frame)
  return(as.data.frame(dt_gen))
}

compute_construct_diff_scores <- function(gamma, tau) {
  info_cols <- c(colnames(gamma)[1:12], "GI.z") #meta info columns + GI scores
  pm <- inner_join(gamma[, colnames(gamma) %in% info_cols], 
                   tau[, colnames(tau) %in% info_cols],
                   by = info_cols[1:12],
                   suffix = c(".Gamma", ".Tau"))
  pm$GI.z <- pm$GI.z.Tau - pm$GI.z.Gamma
  return(pm)
}

### Function to create initial filter dataframe from counts data
### Used in apply_filters.R logic
create_filter_df <- function(counts, sgrna_filt, combination_filt) {
  # Create initial filter dataframe
  to_filt <- counts[, 1:13]
  to_filt$Flag <- NA
  
  # Flag combinations involving identified sgRNAs that have low median representation
  to_filt$Flag[to_filt$FirstPosition %in% sgrna_filt | to_filt$SecondPosition %in% sgrna_filt] <- 
    "Low median representation at T0"
  
  # Flag individual combinations which have low representation at T0
  to_filt$Flag[to_filt$ConstructID %in% combination_filt & is.na(to_filt$Flag)] <- "Low representation at T0"
  
  return(to_filt)
}

### Function to apply correlation filter and return all results
### Consolidates logic from apply_correlation_filter.R
apply_correlation_filter <- function(raw_phenotypes, orind_phenotypes, single_phenotypes, filter_flags, threshold) {
  
  # Join orientation-independent phenotypes with raw phenotypes
  raw_phenotypes <- left_join(raw_phenotypes,
                             orind_phenotypes[, c("GuideCombinationID", 
                                                 "Gamma.OI.R1", "Gamma.OI.R2", "Gamma.OI.Avg", 
                                                 "Tau.OI.R1", "Tau.OI.R2", "Tau.OI.Avg", 
                                                 "Rho.OI.R1", "Rho.OI.R2", "Rho.OI.Avg")],
                             by = "GuideCombinationID")
  
  # Apply correlation filter
  nocorr_filt <- filt_nocorrelation(combphenos = raw_phenotypes,
                                   singlephenos = single_phenotypes,
                                   filterthresh = threshold)
  
  # Update filter flags with correlation results
  updated_filter_flags <- filter_flags
  updated_filter_flags$Flag[updated_filter_flags$FirstPosition %in% nocorr_filt[[2]] | 
                           updated_filter_flags$SecondPosition %in% nocorr_filt[[2]]] <- "Gamma no correlation"
  updated_filter_flags$Flag[updated_filter_flags$FirstPosition %in% nocorr_filt[[3]] | 
                           updated_filter_flags$SecondPosition %in% nocorr_filt[[3]]] <- "Tau no correlation"
  
  # Apply filters to phenotype data
  phenos_filt <- raw_phenotypes[is.na(updated_filter_flags$Flag), ]
  single_pheno_filt <- single_phenotypes[!(single_phenotypes$sgRNA.ID %in% c(nocorr_filt[[2]], nocorr_filt[[3]])), ]
  
  return(list(
    filtered_phenotypes = phenos_filt,
    filtered_single_phenotypes = single_pheno_filt,
    updated_filter_flags = updated_filter_flags,
    correlation_results = nocorr_filt[[1]]
  ))
}