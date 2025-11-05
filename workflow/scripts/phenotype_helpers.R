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
  for (ph in names(phenotype_prefix_map)) {
    prefs <- phenotype_prefix_map[[ph]]
    num_pref <- prefs[1]
    den_pref <- prefs[2]
    rep_cols <- grep(paste0("^", ph, "\\."), colnames(out), value = TRUE)
    rep_cols <- rep_cols[!grepl("\\.Avg$", rep_cols)]
    for (cname in rep_cols) {
      suf <- sub(paste0("^", ph, "\\."), "", cname)
      if (ph != "Rho") {
        dval <- get_doubling_value(doublings, num_pref, suf)
        if (is.na(dval)) {
          stop(sprintf("[%s] Doubling value not found for %s.%s - update DOUBLINGS to a mapping in config", 
                       Sys.time(), num_pref, suf))
        } 
        
        out[[cname]] <- out[[cname]] / dval
        
      } else {
        d_treated <- get_doubling_value(doublings, den_pref, suf)
        d_basal <- get_doubling_value(doublings, num_pref, suf)
        if (is.na(d_treated) || is.na(d_basal)) {
          stop(sprintf("[%s] Doubling values not found for Rho replicate %s (need %s.%s and %s.%s)", 
                       Sys.time(), suf, den_pref, suf, num_pref, suf))
        }
        out[[cname]] <- out[[cname]] / (d_treated - d_basal)
      }
    }
  }
  return(out)
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

### Helper: resolve column name for a given prefix and replicate suffix
### Tries matched_list first (if provided), then searches available colnames for
### common separators (., _, -) and a no-separator form. Falls back to prefix.suf
resolve_pref_suf_colname <- function(prefix, suf, matched_list = NULL, colnames_vec = NULL,
                                     seps = c(".", "_", "-"), strict = FALSE) {
  
  if (!is.null(matched_list) && !is.null(matched_list[[suf]])) {
    return(matched_list[[suf]])
  }

  # 2) if we have available column names, build candidate list in preferred order
  if (!is.null(colnames_vec)) {
    # candidates from separators, then no-sep, then dot as fallback
    sep_candidates <- paste0(prefix, seps, suf)
    candidates <- c(sep_candidates, paste0(prefix, suf), paste0(prefix, ".", suf))
    found <- candidates[candidates %in% colnames_vec]
    if (length(found) > 0) {
      return(found[1])
    }
  }

  # 3) nothing matched: either error (strict) or return legacy dot form
  if (isTRUE(strict)) {
    stop(sprintf("No column matched for prefix='%s' suf='%s'", prefix, suf))
  }
  return(paste0(prefix, ".", suf))
}

### Helper: build fractional abundances matrix from raw counts
### - counts: original counts data.frame
### - cond_colnames: character vector of column names to use (must exist in counts)
### - pseudocount: numeric to add to raw counts
build_fracs <- function(counts, cond_colnames, pseudocount = 10) {
  if (is.null(cond_colnames) || length(cond_colnames) == 0) stop("No condition column names provided to build_fracs")
  # Subset counts by the provided condition columns (preserve order of cond_colnames)
  missing_cols <- setdiff(cond_colnames, colnames(counts))
  if (length(missing_cols) > 0) {
    stop(sprintf("Requested condition columns not found in counts: %s", 
                 paste(missing_cols, collapse = ", ")))
  }

  mat <- as.data.frame(counts[, cond_colnames, drop = FALSE])
  # Ensure numeric
  mat[] <- lapply(mat, function(x) as.numeric(x))
  # Add pseudocount
  mat <- mat + pseudocount

  col_sums <- colSums(mat, na.rm = TRUE)
  if (any(col_sums == 0)) warning(sprintf("One or more condition columns have zero total after pseudocount: %s",
                                        paste(names(col_sums)[col_sums == 0], collapse = ", ")))

  fracs <- sweep(mat, 2, col_sums, FUN = "/")
  return(as.data.frame(fracs))
}

### Helper: precompute mapping from phenotype -> replicate suffixes -> actual column names
### Returns a named list where each element is a list with elements:
###   common_suf (character vector), num_map (named character vector suf->colname), den_map
build_prefix_column_map <- function(conds_colnames, phenotype_prefix_map, replicates = NULL) {
  if (is.null(phenotype_prefix_map) || length(phenotype_prefix_map) == 0) stop("phenotype_prefix_map is required")
  res <- list()
  for (ph in names(phenotype_prefix_map)) {
    prefs <- phenotype_prefix_map[[ph]]
    if (length(prefs) < 2) stop(sprintf("phenotype_prefix_map entry for %s must contain two prefixes (num,den)", ph))
    num_pref <- prefs[1]
    den_pref <- prefs[2]

    # verify there are candidate columns for each prefix
    num_cols <- conds_colnames[startsWith(conds_colnames, num_pref)]
    den_cols <- conds_colnames[startsWith(conds_colnames, den_pref)]
    if (length(num_cols) == 0 || length(den_cols) == 0) {
      stop(sprintf("No columns found for phenotype %s prefixes: %s / %s", ph, num_pref, den_pref))
    }

    match_res <- identify_replicate_matches(conds_colnames, num_pref, den_pref, replicates)
    common_suf <- match_res$common_suf
    matched_num <- match_res$matched_num
    matched_den <- match_res$matched_den
    if (length(common_suf) == 0) {
      stop(sprintf("No common replicate suffixes for phenotype %s between %s and %s", 
                   ph, num_pref, den_pref))
    }

    num_map <- character(0)
    den_map <- character(0)
    for (suf in common_suf) {
      ncol <- resolve_pref_suf_colname(num_pref, suf, matched_num, conds_colnames)
      dcol <- resolve_pref_suf_colname(den_pref, suf, matched_den, conds_colnames)
      # Ensure the resolved names actually exist in conds_colnames
      if (!(ncol %in% conds_colnames)) {
        stop(sprintf("Resolved numerator column '%s' not found for phenotype %s (prefix %s, suf %s)", 
                     ncol, ph, num_pref, suf))
      }
      if (!(dcol %in% conds_colnames)) {
        stop(sprintf("Resolved denominator column '%s' not found for phenotype %s (prefix %s, suf %s)", 
                     dcol, ph, den_pref, suf))
      }
      num_map[suf] <- ncol
      den_map[suf] <- dcol
    }

    res[[ph]] <- list(common_suf = common_suf, num_map = num_map, den_map = den_map,
                      num_pref = num_pref, den_pref = den_pref)
  }
  return(res)
}

### Helper: compute phenotype replicate columns (vectorized-ish)
### - fracs: data.frame/matrix of fractional abundances (columns named with sample colnames)
### - prefix_map: output of build_prefix_column_map
compute_phenotype_matrix <- function(fracs, prefix_map) {
  phenos <- data.frame(matrix(nrow = nrow(fracs), ncol = 0))
  for (ph in names(prefix_map)) {
    entry <- prefix_map[[ph]]
    for (suf in entry$common_suf) {
      ncol <- entry$num_map[[suf]]
      dcol <- entry$den_map[[suf]]
      nidx <- match(ncol, colnames(fracs))
      didx <- match(dcol, colnames(fracs))
      if (is.na(nidx) || is.na(didx)) {
        stop(sprintf("Column index resolution failed for %s.%s (%s vs %s)", ph, suf, ncol, dcol))
      }
      vec <- log2(fracs[[nidx]] / fracs[[didx]])
      phenos[[paste0(ph, ".", suf)]] <- vec
    }
  }
  return(phenos)
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
calculate_phenotypes <- function(counts, conds, pseudocount = 10, doublings = NULL, 
                                 phenotype_prefix_map = NULL, replicates = NULL) {
  # phenotype_prefix_map should be a named list, e.g.
  # list(Gamma = c("DMSO","T0"), Tau = c("NIRAP","T0"), Rho = c("NIRAP","DMSO"))
  if (is.null(phenotype_prefix_map)) {
    phenotype_prefix_map <- list(Gamma = c("DMSO", "T0"),
                                 Tau = c("TREATED", "T0"),
                                 Rho = c("TREATED", "DMSO"))
  }

  message(sprintf("[%s] calculate_phenotypes: starting on %d rows; cond columns: %s", 
          Sys.time(), nrow(counts), paste(conds$Colname, collapse = ", ")))

  # build fractional abundances (pseudocount and column normalization)
  fracs <- build_fracs(counts, conds$Colname, pseudocount = pseudocount)

  # build prefix -> suffix -> column maps once
  prefix_map <- build_prefix_column_map(conds$Colname, phenotype_prefix_map, replicates)

  # compute phenotype replicate columns using the precomputed maps
  phenos <- compute_phenotype_matrix(fracs, prefix_map)

  # combine metadata and computed phenotype columns (preserve original metadata cols)
  out <- cbind(counts[, 1:13], phenos)
  
  # Adjust by NT medians (subtract median of NT+NT controls from each phenotype column)
  message(sprintf("[%s] Adjusting phenotypes by NT medians", Sys.time()))
  nt_medians <- compute_nt_medians(out, phenotype_prefix_map)
  for (ph in names(phenotype_prefix_map)) {
    rep_cols <- grep(paste0("^", ph, "\\."), colnames(out), value = TRUE)
    rep_cols <- rep_cols[!grepl("\\.Avg$", rep_cols)]
    for (cname in rep_cols) {
      suf <- sub(paste0("^", ph, "\\."), "", cname)
      nt_val <- nt_medians[[ph]][[suf]]
      out[[cname]] <- out[[cname]] - nt_val
      message(sprintf("[%s]   %s: NT median = %.4f", Sys.time(), cname, nt_val))
    }
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
calculate_averaged_phenotypes <- function(phenos, phenotype_prefix_map = NULL, replicates = NULL) {
  # Allow caller to provide the phenotype names (Gamma/Tau/Rho etc.)
  if (is.null(phenotype_prefix_map)) {
    phenotype_prefix_map <- list(Gamma = c("DMSO", "T0"),
                                 Tau = c("TREATED", "T0"),
                                 Rho = c("TREATED", "DMSO"))
  }

  phs <- names(phenotype_prefix_map)
  if (length(phs) == 0) stop("phenotype_prefix_map must contain at least one phenotype name")

  # Default replicate suffixes if not provided
  if (is.null(replicates)) replicates <- c("R1", "R2", "Avg")

  conn_pat <- "[._-]?"
  
  # Build a regex to match per-phenotype replicate columns using the provided replicates
  ph_pattern <- paste0("^(", 
                       paste(phs, collapse = "|"), 
                       ")", conn_pat, "(", 
                       paste(replicates, collapse = "|"), 
                       ")$")

  # Compute per-GuideCombinationID means for all matched phenotype replicate columns
  orind <- phenos %>%
    group_by(GuideCombinationID) %>%
    summarise(across(.cols = matches(ph_pattern), .fns = ~ mean(.x, na.rm = TRUE)),
              N = n())

  # Rename phenotype columns to include orientation-independent marker ".OI"
  # e.g. Gamma.R1 -> Gamma.OI.R1
  orind <- orind %>%
    rename_with(.fn = function(x) {
      # Build a replacement pattern that mirrors the matching ph_pattern used
      # above so we insert ".OI." between the phenotype name and suffix
      pat <- paste0("^(", paste(phs, collapse = "|"), ")", conn_pat, "(.*)$")
      sub(pat, "\\1.OI.\\2", x, perl = TRUE)
    }, .cols = -GuideCombinationID)

  # Preserve the original ordering used previously (numeric part of sgc_...)
  orind <- orind[order(as.numeric(gsub("sgc_", "", orind$GuideCombinationID))), ]
  return(orind)
}

### Function for calculating phenotypes for individual sgRNAs using non-targeting guides
### Inputs:
### ### phenos - output from calculate_phenotypes function
### Outputs:
### ### df containing single sgRNA phenotypes by replicate
calculate_single_sgrna_phenotypes <- function(phenos, phenotype_prefix_map = NULL, replicates = NULL) {
  # Allow caller to provide phenotype names and replicates; fall back to defaults
  if (is.null(phenotype_prefix_map)) {
    phenotype_prefix_map <- list(Gamma = c("DMSO", "T0"),
                                 Tau = c("TREATED", "T0"),
                                 Rho = c("TREATED", "DMSO"))
  }
  phs <- names(phenotype_prefix_map)
  if (is.null(phs) || length(phs) == 0) stop("phenotype_prefix_map must contain at least one phenotype")

  if (is.null(replicates)) replicates <- c("R1", "R2", "Avg")

  conn_pat <- "[._-]?"
  
  # Prepare output columns dynamically: <Phenotype>.OI.<Replicate>
  oi_cols <- unlist(lapply(phs, function(ph) paste0(ph, ".OI.", replicates)))

  # Unique sgRNA IDs (from either position column)
  sg_ids <- unique(c(phenos$FirstPosition, phenos$SecondPosition))

  # Initialize output dataframe
  single_pheno <- data.frame("sgRNA.ID" = sg_ids, stringsAsFactors = FALSE)
  for (cname in oi_cols) single_pheno[[cname]] <- 0
  single_pheno$N <- 0

  for (i in seq_len(nrow(single_pheno))) {
    if (i %% 100 == 0) {
      progress_pct <- round((i / nrow(single_pheno)) * 100, 1)
      message(sprintf("[%s] Processing sgRNA %d/%d (%s percent) - ID: %s", 
                      Sys.time(), i, nrow(single_pheno), progress_pct, single_pheno$sgRNA.ID[i]))
    }

    sg <- single_pheno$sgRNA.ID[i]

    # construct tmp: combinations where sgRNA is in A or B and the other is non-targeting (or both NT for NT case)
    if (grepl("non-targeting", sg)) {
      tmp <- phenos[(phenos$FirstPosition == sg & 
                      phenos$SecondPosition %in% unique(phenos$FirstPosition[phenos$Category == "NT+NT"])) |
                    (phenos$SecondPosition == sg & 
                      phenos$FirstPosition %in% unique(phenos$FirstPosition[phenos$Category == "NT+NT"])), ]
    } else {
      tmp <- phenos[(phenos$FirstPosition == sg | phenos$SecondPosition == sg) & 
                    (phenos$FirstPosition %in% unique(phenos$FirstPosition[phenos$Category == "NT+NT"]) |
                     phenos$SecondPosition %in% unique(phenos$FirstPosition[phenos$Category == "NT+NT"])), ]
    }

    # For each phenotype and replicate, compute mean across matching combos (preserve NA if column missing)
    for (ph in phs) {
      for (suf in replicates) {
        src_col <- paste0(ph, ".", suf)
        out_col <- paste0(ph, ".OI.", suf)
        if (src_col %in% colnames(tmp)) {
          single_pheno[[out_col]][i] <- mean(tmp[[src_col]], na.rm = TRUE)
        } else {
          single_pheno[[out_col]][i] <- NA
        }
      }
    }
    single_pheno$N[i] <- nrow(tmp)
  }
  return(single_pheno)
}