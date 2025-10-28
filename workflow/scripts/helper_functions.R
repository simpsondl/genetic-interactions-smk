library(dplyr)
library(broom)

### Function for identifying individual sgRNAs which are poorly represented
### in either position A or position B at T0
### Inputs:
### ### counts - raw count df, import from file
### ### conds - condition df, maps sample and replicate information to column names
### ### filtersamp - which arm of experiment to filter
### ### filterthresh - minimum median representation across all constructs/replicates
### Outputs:
### ### vector with sgRNA ids to with mean median representation below threshold
filt_low_representation <- function(counts, conds, filtersamp, filterthresh = 30){
  # For each sgRNA in each position, get median count across all matching constructs
  # Expected output: 2 df with rows equal to number of sgRNAs and and 6 measurements across arms/replicates
  med.pos1 <- counts %>% group_by(FirstPosition) %>% summarise(across(conds$Colname, median))
  med.pos2 <- counts %>% group_by(SecondPosition) %>% summarise(across(conds$Colname, median))
  # For each sgRNA in each position, calculate mean median coverage per arm
  # Expected output: 2 df with rows equal to number of sgRNAs and 3 measurements across arms
  mean.pos1 <- as.data.frame(lapply(unique(conds$Samplename), 
                                    function(i) apply(med.pos1[,colnames(med.pos1) %in% 
                                                                 conds$Colname[conds$Samplename == i]], 
                                                      1, mean)))
  mean.pos2 <- as.data.frame(lapply(unique(conds$Samplename), 
                                    function(i) apply(med.pos2[,colnames(med.pos2) %in% 
                                                                 conds$Colname[conds$Samplename == i]], 
                                                      1, mean)))
  # Fix names from lapply
  colnames(mean.pos1) <- unique(conds$Samplename)
  colnames(mean.pos2) <- unique(conds$Samplename)
  # Set for return, only need pos1 since indexing is identical
  rownames(mean.pos1) <- med.pos1$FirstPosition
  # Return sgRNA ids which have mean median representation below threshold
  return(rownames(mean.pos1[mean.pos1[,filtersamp] <= filterthresh | mean.pos2[,filtersamp] <= filterthresh,]))
}

### Function for identifying individual sgRNA combinations which are poorly represented
### Inputs:
### ### counts - raw count df, import from file
### ### conds - condition df, maps sample and replicate information to column names
### ### filtersamp - which arm of experiment to filter
### ### filterthresh - minimum median representation across all constructs/replicates
### Outputs:
### ### vector with sgRNA combination ids to with representation below threshold
filt_combinations <- function(counts, conds, filtersamp, filterthresh){
  # Get column indices corresponding to arms to filter
  idx <- which(colnames(counts) %in% conds$Colname[conds$Samplename == filtersamp])
  # Create vector where nonzero entry corresponds to a count below threshold in a relevant arm
  tmp <- rowSums(counts[,idx] <= filterthresh)
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
calculate_phenotypes <- function(counts, conds, pseudocount = 10, normalize = TRUE, doublings){
  pseudo <- counts[,colnames(counts) %in% conds$Colname] + pseudocount
  fracs <- sweep(pseudo, 2, colSums(pseudo), FUN = "/")
  phenos <- data.frame(Gamma.R1 = log2(fracs$DMSO.R1/fracs$T0.R1),
                       Gamma.R2 = log2(fracs$DMSO.R2/fracs$T0.R2),
                       Tau.R1 = log2(fracs$NIRAP.R1/fracs$T0.R1),
                       Tau.R2 = log2(fracs$NIRAP.R2/fracs$T0.R2),
                       Rho.R1 = log2(fracs$NIRAP.R1/fracs$DMSO.R1),
                       Rho.R2 = log2(fracs$NIRAP.R2/fracs$DMSO.R2))
  phenos <- cbind(counts[,1:13], phenos)
  
  if(normalize){
    nt.gamma.r1 <- median(phenos$Gamma.R1[phenos$Category == "NT+NT"])
    nt.gamma.r2 <- median(phenos$Gamma.R2[phenos$Category == "NT+NT"])
    nt.tau.r1 <- median(phenos$Tau.R1[phenos$Category == "NT+NT"])
    nt.tau.r2 <- median(phenos$Tau.R2[phenos$Category == "NT+NT"])
    nt.rho.r1 <- median(phenos$Rho.R1[phenos$Category == "NT+NT"])
    nt.rho.r2 <- median(phenos$Rho.R2[phenos$Category == "NT+NT"])
    
    phenos$Gamma.R1 <- (phenos$Gamma.R1 - nt.gamma.r1)/doublings[1]
    phenos$Gamma.R2 <- (phenos$Gamma.R2 - nt.gamma.r2)/doublings[2]
    phenos$Tau.R1 <- (phenos$Tau.R1 - nt.tau.r1)/doublings[3]
    phenos$Tau.R2 <- (phenos$Tau.R2 - nt.tau.r2)/doublings[4]
    phenos$Rho.R1 <- (phenos$Rho.R1 - nt.rho.r1)/(doublings[1]-doublings[3])
    phenos$Rho.R2 <- (phenos$Rho.R2 - nt.rho.r2)/(doublings[2]-doublings[4])
  }
  
  phenos$Gamma.Avg <- rowMeans(phenos[,c("Gamma.R1", "Gamma.R2")])
  phenos$Tau.Avg <- rowMeans(phenos[,c("Tau.R1", "Tau.R2")])
  phenos$Rho.Avg <- rowMeans(phenos[,c("Rho.R1", "Rho.R2")])
  # rearranges dataframe to have averages next to individual replicates
  return(phenos[,c(1:15,20,16,17,21,18,19,22)])
}

### Function for calculating orientation independent, averaged sgRNA combination phenotypes
### Inputs:
### ### phenos - output from calculate_phenotypes function, OR a df containing columns:
### ### ### Gamma.R1, Gamma.R2
### ### ### Tau.R1, Tau.R2
### ### ### GuideCombinationID with IDs beginning with sgc_
### Outputs:
### ### df containing orientation-independent phenotypes for all replicates and avg replicate
calculate_averaged_phenotypes <- function(phenos){
  orind <- phenos %>% group_by(GuideCombinationID) %>% summarise(Gamma.OI.R1 = mean(Gamma.R1),
                                                                 Gamma.OI.R2 = mean(Gamma.R2),
                                                                 Gamma.OI.Avg = mean(Gamma.Avg),
                                                                 Tau.OI.R1 = mean(Tau.R1),
                                                                 Tau.OI.R2 = mean(Tau.R2),
                                                                 Tau.OI.Avg = mean(Tau.Avg),
                                                                 Rho.OI.R1 = mean(Rho.R1),
                                                                 Rho.OI.R2 = mean(Rho.R2),
                                                                 Rho.OI.Avg = mean(Rho.Avg),
                                                                 N = n())
  orind <- orind[order(as.numeric(gsub("sgc_","",orind$GuideCombinationID))),]
  return(orind)
}

### Function for calculating phenotypes for individual sgRNAs using non-targeting guides
### Inputs:
### ### phenos - output from calculate_phenotypes function
### Outputs:
### ### df containing single sgRNA phenotypes by replicate
calculate_single_sgRNA_phenotypes <- function(phenos){
  # Create dataframe for saving results
  single.pheno <- data.frame("sgRNA.ID" = unique(c(phenos$FirstPosition, phenos$SecondPosition)), 
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
  
  for(i in 1:nrow(single.pheno)){
    # Handle non-targeting case
    if(grepl("non-targeting",single.pheno$sgRNA.ID[i])){
      # Extract all non-targeting guide combinations with desired non-targeting guide in position A or B
      # If this is not handled explicitly, will aggregate all NT guides into a glob
      tmp <- phenos[(phenos$FirstPosition == single.pheno$sgRNA.ID[i] &
                                phenos$SecondPosition %in% unique(phenos$FirstPosition[phenos$Category == "NT+NT"])) | 
                               (phenos$SecondPosition == single.pheno$sgRNA.ID[i] &
                                  phenos$FirstPosition %in% unique(phenos$FirstPosition[phenos$Category == "NT+NT"])),]
      single.pheno$Gamma.OI.R1[i] <- mean(tmp$Gamma.R1)
      single.pheno$Gamma.OI.R2[i] <- mean(tmp$Gamma.R2)
      single.pheno$Gamma.OI.Avg[i] <- mean(tmp$Gamma.Avg)
      single.pheno$Tau.OI.R1[i] <- mean(tmp$Tau.R1)
      single.pheno$Tau.OI.R2[i] <- mean(tmp$Tau.R2)
      single.pheno$Tau.OI.Avg[i] <- mean(tmp$Tau.Avg)
      single.pheno$Rho.OI.R1[i] <- mean(tmp$Rho.R1)
      single.pheno$Rho.OI.R2[i] <- mean(tmp$Rho.R2)
      single.pheno$Rho.OI.Avg[i] <- mean(tmp$Rho.Avg)
      single.pheno$N[i] <- nrow(tmp)
    }
    # Handle targeting case
    else {
      # Extract all combinations with desired targeting guide in position A or B and non-targeting guide in other
      # Easier to do this since we can just grab all non-targeting at once
      tmp <- phenos[(phenos$FirstPosition == single.pheno$sgRNA.ID[i] | phenos$SecondPosition == single.pheno$sgRNA.ID[i]) & 
                (phenos$FirstPosition %in% unique(phenos$FirstPosition[phenos$Category == "NT+NT"]) | 
                   phenos$SecondPosition %in% unique(phenos$FirstPosition[phenos$Category == "NT+NT"])),]
      single.pheno$Gamma.OI.R1[i] <- mean(tmp$Gamma.R1)
      single.pheno$Gamma.OI.R2[i] <- mean(tmp$Gamma.R2)
      single.pheno$Gamma.OI.Avg[i] <- mean(tmp$Gamma.Avg)
      single.pheno$Tau.OI.R1[i] <- mean(tmp$Tau.R1)
      single.pheno$Tau.OI.R2[i] <- mean(tmp$Tau.R2)
      single.pheno$Tau.OI.Avg[i] <- mean(tmp$Tau.Avg)
      single.pheno$Rho.OI.R1[i] <- mean(tmp$Rho.R1)
      single.pheno$Rho.OI.R2[i] <- mean(tmp$Rho.R2)
      single.pheno$Rho.OI.Avg[i] <- mean(tmp$Rho.Avg)
      single.pheno$N[i] <- nrow(tmp)
    }
  }
  return(single.pheno)
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
filt_nocorrelation <- function(combphenos, singlephenos, filterthresh = 0.25){
  cors <- data.frame(sgRNA.ID = singlephenos$sgRNA.ID,
                     Gamma.OI.Correlation = NA,
                     Tau.OI.Correlation = NA,
                     Rho.OI.Correlation = NA)
  
  for(i in singlephenos$sgRNA.ID){
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
### ### singlepheno.df - filtered output from calculate_single_sgRNA_phenotypes and filt_nocorrelation
### ### pairpheno.df - filtered output from calculate_phenotypes and filt_nocorrelation
### ### phenocol - column name that contains the phenotypes to calculate GI scores from
### Outputs:
### ### list containing three elements, in order:
### ### ### df with all called interaction scores
### ### ### a vector containing expected GI scores determined from the model
### ### ### the model that was built
compute_gis <- function(query, singlepheno.df, pairpheno.df, phenocol){
  ### Get the single sgRNA phenotypes
  tmp.data <- data.frame(query = query, 
                         sgRNA.id = singlepheno.df$sgRNA.ID, 
                         single = singlepheno.df[,phenocol])
  colnames(tmp.data)[3] <- "single"
  ### Merge single with interaction phenotypes
  tmp.data.merge <- merge(tmp.data, pairpheno.df[pairpheno.df$SecondPosition == query, 
                                                 c("ConstructID",  "GuideCombinationID", "Identical",
                                                   "FirstPseudogene", "SecondPseudogene", "PseudogeneCombinationID",
                                                   "Orientation", "GeneCombinationID", "Category",
                                                   "FirstPosition", phenocol)], 
                          by.x = "sgRNA.id", by.y = "FirstPosition", all.x = TRUE)
  
  ### Identify control combinations
  tmp.data.merge$Control <- grepl("non-targeting", tmp.data.merge$sgRNA.id)
  
  ### Remove any missing information
  tmp.data.merge <- tmp.data.merge[!is.na(tmp.data.merge[,phenocol]),]
  
  if(nrow(tmp.data.merge) == 0){
    return(NA)
  }
  try({
    ### Fit the model
    lm.fit <- lm(as.formula(paste0(phenocol," ~ single + I(single^2)")), data = tmp.data.merge)
    
    ### Generate points for plotting the fit
    single.pred <- seq(-.35, .05, .01)
    pair.pred <- predict(lm.fit, list(single = single.pred))
    
    ### Define interaction score as observed - expected
    tmp.data.merge$Expected <- predict(lm.fit, list(single = tmp.data.merge$single))
    tmp.data.merge$GI <- lm.fit$residuals
    ### z-standardize GIs to negative controls
    tmp.data.merge$GI.z <- (tmp.data.merge$GI)/sd(tmp.data.merge$GI[tmp.data.merge$Control])
    
    ### Order dataframe
    tmp.data.merge <- tmp.data.merge[order(tmp.data.merge$GI),c("ConstructID", "GuideCombinationID", 
                                                                "PseudogeneCombinationID", "GeneCombinationID", 
                                                                "FirstPseudogene", "SecondPseudogene",
                                                                "Orientation", "Identical", "Category", "Control",
                                                                "sgRNA.id", "query", "single", phenocol, "Expected",
                                                                "GI", "GI.z")]
    
    ### Want to label strongest 10 interactions [both buffering and synthetic lethal]
    n <- nrow(tmp.data.merge)
    tmp.data.merge$lb <- c(gsub("_.*", "", tmp.data.merge$sgRNA.id[1:(10)]),
                           rep("", n-20),
                           gsub("_.*", "", tmp.data.merge$sgRNA.id[(n - 9):n]))
    
  }, silent = TRUE)
  
  return(list(tmp.data.merge, pair.pred, lm.fit))
}

### Function for cleanly saving results from compute_gis - used to establish save df's
### Inputs:
### ### gis - output from calculate_gis function
### Outputs:
### ### list containing three elements, in order:
### ### ### df with gi scores
### ### ### df with model coefficients
### ### ### df with model statistics
create_interaction_result_df <- function(gis){
  all_gis <- gis[[1]]
  ests <- data.frame(sgRNA.ID = i, tidy(gis[[3]]))
  stats <- data.frame(sgRNA.ID = i,glance(gis[[3]]))
  
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
update_interaction_result_df <- function(gis, prev){
  all_gis <- rbind(prev[[1]], gis[[1]])
  ests <- rbind(prev[[2]], data.frame(sgRNA.ID = i,tidy(gis[[3]])))
  stats <- rbind(prev[[3]], data.frame(sgRNA.ID = i,glance(gis[[3]])))
  
  return(list(all_gis, ests, stats))
}

### Function for aggregating the results from compute_gis across constructs
### Inputs:
### ### gis - final full output from update_interaction_result_df
### ### scorecol - column name that contains the interaction scores to be aggregated to gene level
### Outputs:
### ### df containing gene-level GI scores
compute_gene_interaction_scores <- function(gis, scorecol){
  info.cols <- c("GuideCombinationID", "PseudogeneCombinationID", "Category", "Control", scorecol)
  gene.gis <- gis[,colnames(gis) %in% info.cols] %>% 
                group_by(PseudogeneCombinationID) %>% 
                mutate(InteractionScore = mean(!!sym(scorecol)), N = n()) %>%
                select(c(PseudogeneCombinationID, Category, InteractionScore, N)) %>% unique()
  return(gene.gis)
}

### Function for assessing the variance of sgRNA-level GI scores that contribute to gene-level GI scores
### Inputs:
### ### congis - final full output from update_interaction_result_df
### ### genegis - output from compute_gene_interaction_scores, must have columns Gene1 and Gene2 included
### Outputs:
### ### df containing GI scores averaged across sgRNA IDs
assess_sgcscore_variance <- function(congis, genegis){
  
  genegis$Variance.p <- NA
  genegis$Variance.N.NT <- NA
  genegis$Variance.N.Test <- NA
  ids <- unique(genegis$PseudogeneCombinationID)
  total_ids <- length(ids)
  k <- 1
  for (i in ids){
    
    if(k %% 5000 == 0){
      message(sprintf("[%s] assess_sgcscore_variance: processed %d/%d combinations; current PseudogeneCombinationID=%s",
                      Sys.time(), k, total_ids, i))
    }
    
    gene1 <- genegis$Pseudogene1[genegis$PseudogeneCombinationID == i]
    gene2 <- genegis$Pseudogene2[genegis$PseudogeneCombinationID == i]
  
    tmp <- congis[congis$SecondPseudogene == gene1 | congis$SecondPseudogene == gene2,]
    
    if(grepl("NTPG", gene1) & grepl("NTPG", gene2)) {
      tmp <- tmp[tmp$PseudogeneCombinationID == i | 
                    (tmp$Category == "NT+NT" & tmp$SecondPseudogene == gene1) |
                    (tmp$Category == "NT+NT" & tmp$SecondPseudogene == gene2),]
    } else if(grepl("NTPG", gene1)){
      tmp <- tmp[tmp$PseudogeneCombinationID == i | 
                    (tmp$Category == "NT+NT" & tmp$SecondPseudogene == gene1) |
                    (tmp$Category == "X+NT" & tmp$SecondPseudogene == gene2),]
    } else if(grepl("NTPG", gene2)) {
      tmp <- tmp[tmp$PseudogeneCombinationID == i | 
                    (tmp$Category == "NT+NT" & tmp$SecondPseudogene == gene2) |
                    (tmp$Category == "X+NT" & tmp$SecondPseudogene == gene1),]
    } else {
      tmp <- tmp[tmp$PseudogeneCombinationID == i | tmp$Category == "X+NT",]
    }
    
    tmp$TestDist <- tmp$PseudogeneCombinationID == i
    tmp <- tmp[!tmp$Identical,]

    if(length(unique(tmp$TestDist)) == 1){
      next
    }
    
    genegis$Variance.p[genegis$PseudogeneCombinationID == i] <- wilcox.test(GI.z ~ TestDist, data = tmp)$p.value
    genegis$Variance.N.NT[genegis$PseudogeneCombinationID == i] <- sum(!tmp$TestDist)
    genegis$Variance.N.Test[genegis$PseudogeneCombinationID == i] <- sum(tmp$TestDist)
    
    k <- k+1
  }
  
  # Final log: how many p-values were computed
  n_computed <- sum(!is.na(genegis$Variance.p))
  message(sprintf("[%s] assess_sgcscore_variance: completed; computed p-values for %d/%d combinations",
                  Sys.time(), n_computed, total_ids))
  return(genegis)
}

compute_construct_diff_scores <- function(gamma, tau){
  info_cols <- c(colnames(gamma)[1:12], "GI.z") #meta info columns + GI scores
  pm <- inner_join(gamma[, colnames(gamma) %in% info_cols], 
                   tau[, colnames(tau) %in% info_cols],
                   by = info_cols[1:12],
                   suffix = c(".Gamma", ".Tau"))
  pm$GI.z <- pm$GI.z.Tau - pm$GI.z.Gamma
  return(pm)
}