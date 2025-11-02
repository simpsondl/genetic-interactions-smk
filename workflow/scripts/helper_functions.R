library(dplyr)

### Helper: resolve column-name prefixes to column indices
### Inputs:
### ### colnames_vec - character vector of column names (e.g. colnames(counts))
### ### prefixes_param - either a character vector, a single comma-separated string,
### ###                  or a list (from Snakemake params) of prefixes to match
### Output: list(indices = integer vector of 1-based column indices, prefixes = character vector)
resolve_count_prefixes <- function(colnames_vec, prefixes_param) {
  if (is.null(prefixes_param)) {
    stop("counts_prefixes not provided")
  }
  # normalize to character vector
  if (!is.character(prefixes_param)) {
    prefixes <- unlist(prefixes_param)
  } else if (length(prefixes_param) == 1 && grepl(",", prefixes_param)) {
    prefixes <- strsplit(prefixes_param, ",")[[1]]
    prefixes <- trimws(prefixes)
  } else {
    prefixes <- prefixes_param
  }

  if (length(prefixes) == 1) {
    matches <- startsWith(colnames_vec, prefixes)
  } else {
    matches_list <- lapply(prefixes, function(p) startsWith(colnames_vec, p))
    matches <- Reduce(`|`, matches_list)
  }

  indices <- which(matches)
  if (length(indices) == 0) {
    stop(sprintf("No columns matched prefixes: %s", paste(prefixes, collapse = ",")))
  }

  message(sprintf("[%s] resolve_count_prefixes: matched columns for prefixes %s: %s", 
                  Sys.time(), paste(prefixes, collapse = ","), 
                  paste(colnames_vec[matches], collapse = ", ")))

  return(list(indices = indices, prefixes = prefixes))
}

### Helper: get doubling value for a given prefix and replicate suffix
get_doubling_value <- function(doublings, prefix, suf) {
  
  val <- NULL
  
  if (is.null(doublings)) return(NA)
  # Expected doubles provided as nested list
  if (is.list(doublings) && !is.null(doublings[[prefix]])) {
    val <- doublings[[prefix]][[suf]]
  }
  
  if (is.null(val)) {
    return(NA)
  } else {
    return(as.numeric(val))
  }
}

### Helper: compute non-targeting medians per phenotype/replicate
compute_nt_medians <- function(out, phenotype_prefix_map) {
  nt_medians <- list()
  for (ph in names(phenotype_prefix_map)) {
    nt_medians[[ph]] <- list()
    rep_cols <- grep(paste0("^", ph, "\\."), colnames(out), value = TRUE)
    rep_cols <- rep_cols[!grepl("\\.Avg$", rep_cols)]
    for (cname in rep_cols) {
      suf <- sub(paste0("^", ph, "\\."), "", cname)
      nt_medians[[ph]][[suf]] <- median(out[[cname]][out$Category == "NT+NT"], na.rm = TRUE)
    }
  }
  return(nt_medians)
}

### Helper: apply normalization per-phenotype/replicate using doublings mapping
normalize_phenotypes <- function(out, phenotype_prefix_map, doublings) {
  nt_medians <- compute_nt_medians(out, phenotype_prefix_map)
  for (ph in names(phenotype_prefix_map)) {
    prefs <- phenotype_prefix_map[[ph]]
    num_pref <- prefs[1]
    den_pref <- prefs[2]
    rep_cols <- grep(paste0("^", ph, "\\."), colnames(out), value = TRUE)
    rep_cols <- rep_cols[!grepl("\\.Avg$", rep_cols)]
    for (cname in rep_cols) {
      suf <- sub(paste0("^", ph, "\\."), "", cname)
      nt_val <- nt_medians[[ph]][[suf]]
      if (ph != "Rho") {
        dval <- get_doubling_value(doublings, num_pref, suf)
        if (is.na(dval)) {
          stop(sprintf("[%s] Doubling value not found for %s.%s - update DOUBLINGS to a mapping in config", 
                       Sys.time(), num_pref, suf))
        } 
        
        out[[cname]] <- (out[[cname]] - nt_val) / dval
        
      } else {
        d_treated <- get_doubling_value(doublings, den_pref, suf)
        d_basal <- get_doubling_value(doublings, num_pref, suf)
        if (is.na(d_treated) || is.na(d_basal)) {
          stop(sprintf("[%s] Doubling values not found for Rho replicate %s (need %s.%s and %s.%s)", 
                       Sys.time(), suf, den_pref, suf, num_pref, suf))
        }
        out[[cname]] <- (out[[cname]] - nt_val) / (d_treated - d_basal)
      }
    }
  }
  return(out)
}

### Helper: compute per-phenotype replicate averages (adds/updates *.Avg columns)
compute_phenotype_averages <- function(df, phenotype_prefix_map) {
  for (ph in names(phenotype_prefix_map)) {
    rep_cols <- grep(paste0("^", ph, "\\."), colnames(df), value = TRUE)
    rep_cols <- rep_cols[!grepl("\\.Avg$", rep_cols)]
    if (length(rep_cols) > 0) {
      df[[paste0(ph, ".Avg")]] <- rowMeans(df[, rep_cols, drop = FALSE])
    }
  }
  return(df)
}

### Helper: identify replicate suffixes and matched column names for numerator/denominator
### Inputs:
### ### colnames_vec - character vector of available column names (e.g., conds$Colname)
### ### num_pref - numerator prefix (e.g., "DMSO")
### ### den_pref - denominator prefix (e.g., "T0")
### ### replicates - optional character vector of replicate suffixes to look for (e.g., c("R1","R2"))
### Output: list(common_suf=character(), matched_num=list or NULL, matched_den=list or NULL)
identify_replicate_matches <- function(colnames_vec, num_pref, den_pref, replicates = NULL) {
  num_cols <- colnames_vec[startsWith(colnames_vec, num_pref)]
  den_cols <- colnames_vec[startsWith(colnames_vec, den_pref)]
  if (length(num_cols) == 0 || length(den_cols) == 0) {
    return(list(common_suf = character(0), matched_num = NULL, matched_den = NULL))
  }

  if (!is.null(replicates)) {
    matched_num <- list()
    matched_den <- list()
    common_suf <- character(0)
    for (suf in replicates) {
      # accept optional separators (., _, -) between prefix and suffix
      num_pattern <- paste0("^", num_pref, "(?:\\.|_|-)?", suf, "$")
      den_pattern <- paste0("^", den_pref, "(?:\\.|_|-)?", suf, "$")
      nm <- grep(num_pattern, num_cols, value = TRUE, perl = TRUE)
      dm <- grep(den_pattern, den_cols, value = TRUE, perl = TRUE)
      if (length(nm) > 0 && length(dm) > 0) {
        matched_num[[suf]] <- nm[1]
        matched_den[[suf]] <- dm[1]
        common_suf <- c(common_suf, suf)
      }
    }
    return(list(common_suf = common_suf, matched_num = matched_num, matched_den = matched_den))
  } else {
    # derive replicate suffixes by removing optional separator after prefix
    num_suf <- sub(paste0("^", num_pref, "(?:\\.|_|-)?"), "", num_cols)
    den_suf <- sub(paste0("^", den_pref, "(?:\\.|_|-)?"), "", den_cols)
    common_suf <- intersect(num_suf, den_suf)
    return(list(common_suf = common_suf, matched_num = NULL, matched_den = NULL))
  }
}

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
calculate_phenotypes <- function(counts, conds, pseudocount = 10, normalize = TRUE,
                                 doublings = NULL, phenotype_prefix_map = NULL, replicates = NULL,
                                 compute_avgs = TRUE) {
  # phenotype_prefix_map should be a named list, e.g.
  # list(Gamma = c("DMSO","T0"), Tau = c("NIRAP","T0"), Rho = c("NIRAP","DMSO"))
  if (is.null(phenotype_prefix_map)) {
    phenotype_prefix_map <- list(Gamma = c("DMSO", "T0"),
                                 Tau = c("TREATED", "T0"),
                                 Rho = c("TREATED", "DMSO"))
  }

  message(sprintf("[%s] calculate_phenotypes: starting on %d rows; cond columns: %s", 
                  Sys.time(), nrow(counts), paste(conds$Colname, collapse = ", ")))

  pseudo <- counts[, colnames(counts) %in% conds$Colname] + pseudocount
  fracs <- sweep(pseudo, 2, colSums(pseudo), FUN = "/")

  # initialize phenos with same number of rows as counts
  phenos <- data.frame(matrix(nrow = nrow(counts), ncol = 0))
  # For each phenotype (Gamma, Tau, Rho), compute replicate columns based on common suffixes
  for (ph in names(phenotype_prefix_map)) {
    prefs <- phenotype_prefix_map[[ph]]
    if (length(prefs) < 2) stop(sprintf("phenotype_prefix_map entry for %s must contain two prefixes (num,den)", ph))
    num_pref <- prefs[1]
    den_pref <- prefs[2]

    # find matching column names within the selected conds columns
    num_cols <- conds$Colname[startsWith(conds$Colname, num_pref)]
    den_cols <- conds$Colname[startsWith(conds$Colname, den_pref)]
    if (length(num_cols) == 0 || length(den_cols) == 0) {
      stop(sprintf("No columns found for phenotype %s prefixes: %s / %s", ph, num_pref, den_pref))
    }

    # identify replicate suffixes and matched column names (delegated to helper)
    match_res <- identify_replicate_matches(conds$Colname, num_pref, den_pref, replicates)
    common_suf <- match_res$common_suf
    matched_num <- match_res$matched_num
    matched_den <- match_res$matched_den
    if (length(common_suf) == 0) {
      stop(sprintf("No common replicate suffixes for phenotype %s between %s and %s", ph, num_pref, den_pref))
    }

    # for each common suffix, compute phenotype replicate column
    for (suf in common_suf) {
      # use matched column names when available (handles different separators)
      if (!is.null(matched_num) && !is.null(matched_num[[suf]])) {
        num_colname <- matched_num[[suf]]
      } else {
        num_colname <- paste0(num_pref, ".", suf)
      }
      if (!is.null(matched_den) && !is.null(matched_den[[suf]])) {
        den_colname <- matched_den[[suf]]
      } else {
        den_colname <- paste0(den_pref, ".", suf)
      }
      ph_colname <- paste0(ph, ".", suf)
      phenos[[ph_colname]] <- log2(fracs[[num_colname]] / fracs[[den_colname]])
    }
  }

  # combine metadata and computed phenotype columns
  out <- cbind(counts[, 1:13], phenos)

  if (normalize) {
    message(sprintf("[%s] Starting normalization step", Sys.time()))
    out <- normalize_phenotypes(out, phenotype_prefix_map, doublings)
    message(sprintf("[%s] Normalization complete", Sys.time()))
  }

  # compute averages optionally (replicate averaging / adding *.Avg columns)
  if (isTRUE(compute_avgs)) {
    out <- compute_phenotype_averages(out, phenotype_prefix_map)
  }

  return(out)
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

### Function for aggregating the results from compute_gis across constructs
### Inputs:
### ### gis - final full output from update_interaction_result_df
### ### scorecol - column name that contains the interaction scores to be aggregated to gene level
### Outputs:
### ### df containing gene-level GI scores
compute_gene_interaction_scores <- function(gis, scorecol) {
  info_cols <- c("GuideCombinationID", "PseudogeneCombinationID", "Identical", "Category", "Control", scorecol)
  gene_gis <- gis[, colnames(gis) %in% info_cols] %>%
                filter(Identical == FALSE) %>% 
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
  grp <- df[[group_col]]
  # Ensure logical vector for counting
  grp_logical <- as.logical(grp)
  n_test <- as.integer(sum(grp_logical, na.rm = TRUE))
  n_control <- as.integer(sum(!grp_logical, na.rm = TRUE))

  # Try the test and return NA on error
  wt <- try(wilcox.test(as.formula(paste0(score_col, " ~ ", group_col)), data = df), silent = TRUE)
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
    gene1 <- row_gen$Pseudogene1[1]
    gene2 <- row_gen$Pseudogene2[1]

    # Fast subset: get all rows where SecondPseudogene is gene1 or gene2
    tmp <- dt_con[J(c(gene1, gene2)), nomatch = 0]
    
    # Build selection index
    sel_idx <- build_selection_index(tmp, i, gene1, gene2)

    tmp2 <- tmp[sel_idx, ]
    tmp2 <- tmp2[tmp2$Identical == FALSE, ]

    tmp2$testdist <- tmp2$PseudogeneCombinationID == i
    tmp2 <- tmp2 %>% distinct()

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
