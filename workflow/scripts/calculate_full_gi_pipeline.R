#!/usr/bin/env Rscript

# Full GI Pipeline Script - Maintains floating point precision throughout
# This script runs the entire pipeline from raw counts to differential hits
# while maintaining full precision by avoiding intermediate file I/O

library(readr)
library(dplyr)
library(broom)
library(data.table)

source("scripts/helper_functions.R")

# Set high precision for numerical output
options(digits = 17, scipen = 999)

# Checkpoint system for recovery from partial runs
checkpoint_dir <- function(outputs) {
    # Scope checkpoints per output gi_scores file (i.e., per score) to avoid collisions
    gi_dir <- dirname(outputs$gi_scores)
    gi_base <- basename(outputs$gi_scores)
    gi_stem <- sub("\\.[^.]*$", "", gi_base)
    checkpoint_path <- file.path(gi_dir, "checkpoints", gi_stem)
    dir.create(checkpoint_path, recursive = TRUE, showWarnings = FALSE)
    return(checkpoint_path)
}

save_checkpoint <- function(step_name, data, checkpoint_path) {
    checkpoint_file <- file.path(checkpoint_path, paste0(step_name, "_checkpoint.rds"))
    message(sprintf("[%s] Saving checkpoint: %s", Sys.time(), step_name))
    
    # Save with maximum precision and no compression to avoid artifacts
    # Use version 3 format for best precision preservation
    old_digits <- getOption("digits")
    old_scipen <- getOption("scipen")
    options(digits = 17, scipen = 999)
    
    saveRDS(data, checkpoint_file, compress = FALSE, version = 3)
    
    # Restore original options
    options(digits = old_digits, scipen = old_scipen)
    
    message(sprintf("[%s] Checkpoint saved: %s (%.2f MB)", 
                   Sys.time(), step_name, file.size(checkpoint_file) / 1024^2))
    return(checkpoint_file)
}

load_checkpoint <- function(step_name, checkpoint_path) {
    checkpoint_file <- file.path(checkpoint_path, paste0(step_name, "_checkpoint.rds"))
    if (file.exists(checkpoint_file)) {
        message(sprintf("[%s] Loading checkpoint: %s", Sys.time(), step_name))
        
        # Ensure high precision options are set before loading
        old_digits <- getOption("digits")
        old_scipen <- getOption("scipen")
        options(digits = 17, scipen = 999)
        
        data <- readRDS(checkpoint_file)
        
        # Keep high precision options (don't restore old ones)
        # since we want high precision throughout the pipeline
        
        message(sprintf("[%s] Checkpoint loaded: %s (%.2f MB)", 
                       Sys.time(), step_name, file.size(checkpoint_file) / 1024^2))
        return(data)
    }
    return(NULL)
}

# Test precision preservation through checkpoint save/load cycle
test_checkpoint_precision <- function(checkpoint_path) {
    # Create test data with high precision numbers
    test_data <- list(
        high_precision_number = 1.234567890123456789,
        scientific_notation = 1.23e-15,
        very_small = .Machine$double.eps,
        data_frame = data.frame(
            x = c(1.111111111111111, 2.222222222222222, 3.333333333333333),
            y = c(4.444444444444444e-10, 5.555555555555555e-10, 6.666666666666666e-10)
        )
    )
    
    # Save and load
    test_file <- save_checkpoint("precision_test", test_data, checkpoint_path)
    loaded_data <- load_checkpoint("precision_test", checkpoint_path)
    
    # Check precision preservation
    precision_preserved <- TRUE
    
    if (!identical(test_data$high_precision_number, loaded_data$high_precision_number)) {
        message("WARNING: High precision number not preserved!")
        precision_preserved <- FALSE
    }
    
    if (!identical(test_data$scientific_notation, loaded_data$scientific_notation)) {
        message("WARNING: Scientific notation not preserved!")
        precision_preserved <- FALSE
    }
    
    if (!identical(test_data$very_small, loaded_data$very_small)) {
        message("WARNING: Very small number not preserved!")
        precision_preserved <- FALSE
    }
    
    if (!identical(test_data$data_frame, loaded_data$data_frame)) {
        message("WARNING: Data frame precision not preserved!")
        precision_preserved <- FALSE
    }
    
    # Clean up test file
    if (file.exists(test_file)) {
        file.remove(test_file)
    }
    
    if (precision_preserved) {
        message(sprintf("[%s] Checkpoint precision verification PASSED", Sys.time()))
    } else {
        message(sprintf("[%s] Checkpoint precision verification FAILED", Sys.time()))
    }
    
    return(precision_preserved)
}

## (removed) check_outputs_exist: upstream resume logic simplification made this helper unnecessary

get_resume_point <- function(outputs, checkpoint_path) {
    # Determine the NEXT step to run based on which steps are complete.
    # Precision pipeline runs only GI -> gene_level -> discriminant -> hits.
    steps <- list(
        list(name = "gi_scores", 
             outputs = c("gi_scores", "model_estimates", "model_stats"), 
             checkpoint = "step4_gi_scores"),
        list(name = "gene_level", 
             outputs = c("gene_level_scores"), 
             checkpoint = "step5_gene_level"),
        list(name = "discriminant", 
             outputs = c("discriminant_scores"), 
             checkpoint = "step6_discriminant"),
        list(name = "hits", 
             outputs = c("hits"), 
             checkpoint = "step7_hits")
    )

    for (step in steps) {
        # Determine if this step is applicable to this rule by checking if any outputs are declared
        step_output_paths <- unlist(lapply(step$outputs, function(nm) outputs[[nm]]))
        step_output_paths <- step_output_paths[!vapply(step_output_paths, is.null, logical(1))]
        if (length(step_output_paths) == 0) {
            # This step is not part of this rule; skip it without requiring a checkpoint
            next
        }
        outputs_ok <- all(file.exists(step_output_paths))
        checkpoint_ok <- file.exists(file.path(checkpoint_path, paste0(step$checkpoint, "_checkpoint.rds")))
        if (outputs_ok && checkpoint_ok) {
            message(sprintf("Found complete step: %s (outputs and checkpoint exist)", step$name))
            next
        } else {
            if (!outputs_ok) message(sprintf("Step %s incomplete (missing outputs) - will resume here", step$name))
            if (outputs_ok && !checkpoint_ok){
                message(sprintf("Step %s has outputs but no checkpoint - will resume here", step$name))
            }
            return(list(resume_from = step$name, checkpoint_data = NULL))
        }
    }

    # If we got here, all steps are complete
    return(list(resume_from = "done", checkpoint_data = NULL))
}

# Function to set up logging for each step with dual output (console + log file)
setup_step_logging <- function(step_name, log_path) {
    if (!is.null(log_path) && length(log_path) > 0) {
        log_file <- log_path
        log_con <- file(log_file, open = "wt")
        
        # Create a custom message function that writes to both console and log
        assign("dual_message", function(...) {
            msg <- sprintf(...)
            # Write to console (stdout)
            cat(msg, "\n", file = stdout())
            # Write to log file
            cat(msg, "\n", file = log_con)
            flush(log_con)
        }, envir = .GlobalEnv)
        
        # Override the message function temporarily for this step
        assign("original_message", message, envir = .GlobalEnv)
        assign("message", dual_message, envir = .GlobalEnv)
        
        # Set warning level
        options(warn = 1)
        
        dual_message("=== %s - %s ===", step_name, Sys.time())
        
        # Return cleanup function
        return(function() {
            # Restore original message function
            if (exists("original_message", envir = .GlobalEnv)) {
                assign("message", get("original_message", envir = .GlobalEnv), envir = .GlobalEnv)
                rm("original_message", envir = .GlobalEnv)
            }
            # Clean up dual_message function
            if (exists("dual_message", envir = .GlobalEnv)) {
                rm("dual_message", envir = .GlobalEnv)
            }
            close(log_con)
        })
    } else {
        # If no log file, just use regular console output
        message(sprintf("=== %s - %s ===", step_name, Sys.time()))
        return(function() {}) # No-op cleanup function
    }
}

# High-precision write function
write_tsv_precision <- function(data, path, precision = 17) {
    # Save current options
    old_digits <- getOption("digits")
    old_scipen <- getOption("scipen")
    
    # Set high precision options
    options(digits = precision, scipen = 999)
    
    # Use data.table for output with high precision settings
    data.table::fwrite(data, path, sep = "\t", quote = FALSE)
    
    # Restore original options
    options(digits = old_digits, scipen = old_scipen)
}

# Main pipeline function
run_full_gi_pipeline <- function(snakemake_obj) {
    
    start_time <- Sys.time()

    # Extract inputs from snakemake object (precision pipeline consumes filtered inputs only)
    input_filtered_phenotypes_path <- snakemake_obj@input[["input_filtered_phenotypes"]]
    input_filtered_single_path <- snakemake_obj@input[["input_filtered_single_sgRNA_phenotypes"]]
    
    # Get all output paths (precision pipeline focuses on GI -> hits only)
    outputs <- list(
        # From GI score calculation
        gi_scores = snakemake_obj@output[["output_all_scores"]],
        model_estimates = snakemake_obj@output[["output_model_estimates"]],
        model_stats = snakemake_obj@output[["output_model_stats"]],
        individual_scores_dir = snakemake_obj@output[["output_dir"]],
        
        # Additional outputs for gene-level, discriminant scores, and hits
        gene_level_scores = snakemake_obj@output[["output_gene_level_scores"]],
        discriminant_scores = snakemake_obj@output[["output_discriminant_scores"]],
        hits = snakemake_obj@output[["output_hits"]]
    )
    
    # Get all log paths (remove upstream step logs)
    logs <- list(
        gi_scores = snakemake_obj@log[["gi_scores"]],
        gene_level = snakemake_obj@log[["gene_level"]],
        discriminant = snakemake_obj@log[["discriminant"]],
        hit_calling = snakemake_obj@log[["hit_calling"]]
    )
    
    # Get all parameters
    params <- list(
        score = snakemake_obj@params[["score"]],
        screen = snakemake_obj@params[["screen"]],
        hit_threshold = snakemake_obj@params[["hit_threshold"]],
        idmap_path = snakemake_obj@input[["input_idmap"]]
    )
    
    
    # ================================
    # CHECKPOINT RECOVERY SYSTEM
    # ================================
    checkpoint_path <- checkpoint_dir(outputs)
    
    # Test checkpoint precision preservation
    if (!test_checkpoint_precision(checkpoint_path)) {
        stop("Checkpoint precision test failed! Pipeline aborted to prevent precision loss.")
    }
    
    recovery_info <- get_resume_point(outputs, checkpoint_path)
    resume_from <- recovery_info$resume_from
    checkpoint_data <- recovery_info$checkpoint_data
    
    message(sprintf("========================================"))
    message(sprintf("HIGH-PRECISION GI PIPELINE"))
    message(sprintf("Screen: %s, Score: %s", params$screen, params$score))
    message(sprintf("Resume from: %s", resume_from))
    message(sprintf("Time: %s", Sys.time()))
    message(sprintf("========================================"))

    # If upstream steps are not part of this rule but filtered inputs were provided,
    # preload them so we can start directly at GI.
    preloaded_filtered_inputs <- FALSE
    if (!is.null(input_filtered_phenotypes_path)) {
        message(sprintf("[%s] Loading filtered inputs from standard pipeline outputs", Sys.time()))
        phenos_filt <- read_tsv(input_filtered_phenotypes_path, show_col_types = FALSE)
        if (!is.null(input_filtered_single_path)) {
            single_pheno_filt <- read_tsv(input_filtered_single_path, show_col_types = FALSE)
        }
        preloaded_filtered_inputs <- TRUE
    }

    # If everything is already complete, exit early
    if (identical(resume_from, "done")) {
        message("All steps complete. Nothing to do. Exiting.")
        return(invisible(TRUE))
    }
    
    # Initialize variables for potential resumption
    # upstream intermediates removed for precision path
    phenos_filt <- NULL
    single_pheno_filt <- NULL
    combined_all <- NULL
    combined_ests <- NULL
    combined_stats <- NULL
    gene_level_scores <- NULL
    
    
    # Steps 1–3 (upstream) removed: precision pipeline starts at GI with filtered inputs provided by standard pipeline.
    
    
    
    
    
    
    
    
    # ================================
    # Calculate GI Scores
    # ================================
    if (resume_from %in% c("gi_scores")) {
    cleanup_gi <- setup_step_logging("Calculate GI Scores", logs$gi_scores)
    
    message(sprintf("[%s] Starting GI score calculation for score: %s", Sys.time(), params$score))
    
    # Create output directory for individual scores
    if (!is.null(outputs$individual_scores_dir)) {
        dir.create(outputs$individual_scores_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Ensure filtered inputs are present; if missing, try to reload from provided paths
    if (is.null(phenos_filt) && !is.null(input_filtered_phenotypes_path)) {
        message(sprintf("[%s] Reloading filtered phenotypes from %s", Sys.time(), input_filtered_phenotypes_path))
        phenos_filt <- read_tsv(input_filtered_phenotypes_path, show_col_types = FALSE)
    }
    if (is.null(single_pheno_filt) && !is.null(input_filtered_single_path)) {
        message(sprintf("[%s] Reloading filtered single sgRNA phenotypes from %s", 
                        Sys.time(), input_filtered_single_path))
        single_pheno_filt <- read_tsv(input_filtered_single_path, show_col_types = FALSE)
    }
    # Validate preloaded inputs and expected columns
    if (is.null(single_pheno_filt) || is.null(phenos_filt)) {
        stop("Filtered phenotypes and filtered single sgRNA phenotypes are required but missing.")
    }
    if (!("sgRNA.ID" %in% colnames(single_pheno_filt))) {
        stop("Expected column 'sgRNA.ID' not found in filtered single sgRNA phenotypes.")
    }
    if (!(params$score %in% colnames(single_pheno_filt))) {
        stop(sprintf("Expected score column '%s' not found in filtered single sgRNA phenotypes.", params$score))
    }
    if (!(params$score %in% colnames(phenos_filt))) {
        stop(sprintf("Expected score column '%s' not found in filtered phenotypes.", params$score))
    }

    # Quick diagnostics
    message(sprintf("[%s] Filtered singles: %d rows, %d cols; Filtered combos: %d rows, %d cols", 
                    Sys.time(), nrow(single_pheno_filt), ncol(single_pheno_filt), nrow(phenos_filt), ncol(phenos_filt)))

    # Pre-allocate lists to collect results
    n <- length(single_pheno_filt$sgRNA.ID)
    all_gis_list <- vector("list", n)
    ests_list <- vector("list", n)
    stats_list <- vector("list", n)
    
    message(sprintf("[%s] Processing %d sgRNAs for score %s", Sys.time(), n, params$score))
    
    for (idx in seq_along(single_pheno_filt$sgRNA.ID)) {
        i <- single_pheno_filt$sgRNA.ID[idx]
        
        # Log progress every 100 sgRNAs and at key milestones
        if (idx %% 100 == 0) {
            progress_pct <- round((idx / n) * 100, 1)
            message(sprintf("[%s] Processing sgRNA %d/%d (%s percent) - ID: %s", 
                           Sys.time(), idx, n, progress_pct, i))
        } 
        
        # Calculate GI scores for this sgRNA
        gi_scores <- compute_gis(i, single_pheno_filt, phenos_filt, params$score)
        
        # Handle cases where compute_gis returns NA (no data)
        if (length(gi_scores) == 1 && is.na(gi_scores)) {
            all_gis_list[[idx]] <- tibble()
            ests_list[[idx]] <- tibble()
            stats_list[[idx]] <- tibble()
            next
        }
        
        # Extract results maintaining full precision
        model_obj <- gi_scores[[3]]
        all_gis_list[[idx]] <- gi_scores[[1]]
        if (!is.null(model_obj)) {
            ests_list[[idx]] <- data.frame(sgRNA.ID = i, tidy(model_obj))
            stats_list[[idx]] <- data.frame(sgRNA.ID = i, glance(model_obj))
        } else {
            ests_list[[idx]] <- tibble()
            stats_list[[idx]] <- tibble()
        }
        
        # Write individual results if directory is specified
        if (!is.null(outputs$individual_scores_dir)) {
            write_tsv_precision(all_gis_list[[idx]], 
                               file.path(outputs$individual_scores_dir, paste0(i, ".txt")))
        }
    }
    
    # Clear the progress line
    cat("\n")
    
    # Combine results maintaining precision
    combined_all <- bind_rows(all_gis_list)
    combined_ests <- bind_rows(ests_list)
    combined_stats <- bind_rows(stats_list)
    
    message(sprintf("[%s] Combined GI results: constructs=%d, model_est_rows=%d, model_stat_rows=%d", 
                   Sys.time(), nrow(combined_all), nrow(combined_ests), nrow(combined_stats)))
    
    # Write GI score outputs with high precision
    write_tsv_precision(combined_all, outputs$gi_scores)
    write_tsv_precision(combined_ests, outputs$model_estimates)
    write_tsv_precision(combined_stats, outputs$model_stats)
    
    message(sprintf("[%s] Wrote GI score outputs", Sys.time()))
    
    cleanup_gi()
    
    # Save checkpoint for GI scores
    step4_data <- list(
        phenos_filt = phenos_filt,
        single_pheno_filt = single_pheno_filt,
        combined_all = combined_all,
        combined_ests = combined_ests,
        combined_stats = combined_stats
    )
    save_checkpoint("step4_gi_scores", step4_data, checkpoint_path)

    message(sprintf(">>> GI SCORES COMPLETE: Calculated GI scores (%d constructs processed) <<<", 
                   nrow(combined_all)))
    
    } else if (resume_from %in% c("gene_level", "discriminant", "hits", "done")) {
    # Load checkpoint data for GI scores
    message(sprintf(">>> GI SCORES SKIPPED: Loading from checkpoint <<<"))
        step4_data <- load_checkpoint("step4_gi_scores", checkpoint_path)
        phenos_filt <- step4_data$phenos_filt
        single_pheno_filt <- step4_data$single_pheno_filt
        combined_all <- step4_data$combined_all
        combined_ests <- step4_data$combined_ests
        combined_stats <- step4_data$combined_stats
    }
    
    
    # ================================
    # Calculate Gene-Level Scores (if needed)
    # ================================
    if (!is.null(outputs$gene_level_scores)) {
        if (resume_from %in% c("gene_level")) {
        cleanup_gene <- setup_step_logging("Calculate Gene Level Scores", logs$gene_level)
        
        message(sprintf("[%s] Starting gene-level score calculation", Sys.time()))
        
        # Read ID mapping
        idmap <- read_tsv(params$idmap_path, show_col_types = FALSE)
        
        # Calculate gene-level scores maintaining precision
        gene_level_scores <- compute_gene_interaction_scores(combined_all, params$score)

        message(sprintf("[%s] Calculated %d gene-level scores", Sys.time(), nrow(gene_level_scores)))

        # Join with ID mapping to add Pseudogene1/Pseudogene2 columns needed by assess_sgcscore_variance
        gene_level_scores <- left_join(gene_level_scores, idmap) %>% 
                            relocate(Pseudogene1:PseudogeneCombinationName)

        # Write with high precision
        write_tsv_precision(gene_level_scores, outputs$gene_level_scores)
        
        cleanup_gene()
        
    # Save checkpoint for gene-level
        step5_data <- list(
            combined_all = combined_all,
            gene_level_scores = gene_level_scores
        )
        save_checkpoint("step5_gene_level", step5_data, checkpoint_path)
        
        } else if (resume_from %in% c("discriminant", "hits", "done")) {
            # Load checkpoint data for gene-level
            message(sprintf(">>> GENE-LEVEL SKIPPED: Loading from checkpoint <<<"))
            step5_data <- load_checkpoint("step5_gene_level", checkpoint_path)
            combined_all <- step5_data$combined_all
            gene_level_scores <- step5_data$gene_level_scores
        }
    }
    
    
    # ================================
    # Calculate Discriminant Scores (if needed)
    # ================================
    if (!is.null(outputs$discriminant_scores)) {
        if (resume_from %in% c("discriminant")) {
        cleanup_disc <- setup_step_logging("Calculate Discriminant Scores", logs$discriminant)
        
        message(sprintf("[%s] Starting discriminant score calculation", Sys.time()))
        
        message(sprintf("[%s] Assessing variance of SGC scores (assess_sgcscore_variance)", Sys.time()))
        # Calculate discriminant scores using gene-level scores
        gene_gis_var <- assess_sgcscore_variance(combined_all, gene_level_scores)
        
        message(sprintf("[%s] Calculated %d discriminant scores", Sys.time(), nrow(gene_gis_var)))
        
        message(sprintf("[%s] Calculating Discriminant = -log10(Variance.p) * abs(InteractionScore)", Sys.time()))
        gene_gis_var$Discriminant <- -log10(gene_gis_var$Variance.p) * abs(gene_gis_var$InteractionScore)


        # Write with high precision
        write_tsv_precision(gene_gis_var, outputs$discriminant_scores)
        
        cleanup_disc()
        
    # Save checkpoint for discriminant
        step6_data <- list(
            combined_all = combined_all,
            gene_level_scores = gene_level_scores,
            gene_gis_var = gene_gis_var
        )
        save_checkpoint("step6_discriminant", step6_data, checkpoint_path)
        
        }
    }
    
    
    # ================================
    # Call Hits (if needed)
    # ================================
    if (!is.null(outputs$hits)) {
        if (resume_from %in% c("hits")) {
        cleanup_hits <- setup_step_logging("Call Hits", logs$hit_calling)
        
        message(sprintf("[%s] Starting hit calling", Sys.time()))
        
        # Use discriminant scores for hit calling
        if (!exists("gene_gis_var")) {
            # Load from checkpoint if not already in memory
            if (file.exists(file.path(checkpoint_path, "step6_discriminant_checkpoint.rds"))) {
                step6_data <- load_checkpoint("step6_discriminant", checkpoint_path)
                gene_gis_var <- step6_data$gene_gis_var
            } else {
                stop("Discriminant scores not available for hit calling")
            }
        }
        
        # Apply hit threshold
        hit_threshold <- params$hit_threshold
        message(sprintf("[%s] Applying hit threshold: %s", Sys.time(), hit_threshold))
        
        # Identify baseline categories for quantile calculation
        baseline_vals <- gene_gis_var %>% 
            filter(Category %in% c("X+NT", "NT+NT")) %>% 
            pull(Discriminant)
        
        # Compute quantile threshold
        qval <- as.numeric(quantile(baseline_vals, probs = hit_threshold, na.rm = TRUE))
        message(sprintf("[%s] Computed discriminant quantile (prob=%s) = %s", 
                       Sys.time(), hit_threshold, qval))
        
        # Add Hit logical column
        hits_with_calls <- gene_gis_var %>%
            mutate(Hit = Discriminant > qval)
        
        # Report hit counts
        n_hits <- sum(hits_with_calls$Hit, na.rm = TRUE)
        n_total <- nrow(hits_with_calls)
        message(sprintf("[%s] Hit summary: %d/%d (%.2f percent)", 
                       Sys.time(), n_hits, n_total, 100 * n_hits / n_total))
        
        # Write hits with high precision
        write_tsv_precision(hits_with_calls, outputs$hits)
        message(sprintf("[%s] Wrote hits to %s", Sys.time(), outputs$hits))
        
        cleanup_hits()
        
    # Save checkpoint for hits
        step7_data <- list(
            gene_gis_var = gene_gis_var,
            hits_with_calls = hits_with_calls
        )
        save_checkpoint("step7_hits", step7_data, checkpoint_path)
        
    message(sprintf(">>> HITS COMPLETE: Called hits (%d hits identified) <<<", n_hits))
        
        } else {
            # Load checkpoint data for hits if available
            if (file.exists(file.path(checkpoint_path, "step7_hits_checkpoint.rds"))) {
                message(sprintf(">>> HITS SKIPPED: Loading from checkpoint <<<"))
                step7_data <- load_checkpoint("step7_hits", checkpoint_path)
                hits_with_calls <- step7_data$hits_with_calls
            }
        }
    }
    
    
    # ================================
    # STEP 8: Create Precision Handoff for Differential Analysis (if needed)
    # ================================
    # Create a minimal workspace containing only construct scores for differential analysis
    if (length(grep("(Gamma|Tau)\\.(OI|noOI)\\.(R1|R2|Avg)$", params$score)) > 0) {
        message(sprintf("[%s] Creating precision handoff for differential analysis", Sys.time()))
        
        # Create differential analysis handoff directory
        handoff_dir <- file.path(dirname(checkpoint_path), "differential_handoff")
        dir.create(handoff_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Extract just the essential data for differential analysis
        # Only keep: construct scores, basic identifiers, no intermediate processing data
        minimal_construct_data <- combined_all %>%
            select(any_of(c("ConstructID", "sgRNA.ID", "PseudogeneCombinationName", 
                           "Gene1", "Gene2", "Pseudogene1", "Pseudogene2",
                           "GI.z", "EffectSize", "Category")))
        
        # Save minimal handoff workspace
        handoff_file <- file.path(handoff_dir, sprintf("%s_%s_construct_scores.rds", 
                                                      params$screen, params$score))
        saveRDS(minimal_construct_data, handoff_file, compress = FALSE, version = 3)
        
        handoff_size_mb <- file.size(handoff_file) / 1024^2
        message(sprintf("[%s] Created precision handoff: %s (%.2f MB)", 
                       Sys.time(), handoff_file, handoff_size_mb))
        
        # Clean up memory by removing large intermediate objects not needed for differential analysis
        rm(list = setdiff(ls(), c("minimal_construct_data", "params", "outputs", 
                                 "checkpoint_path", "handoff_file", "handoff_size_mb")))
        gc()
        
        message(sprintf("[%s] Memory cleaned for differential analysis handoff", Sys.time()))
    }
    
    
    message(sprintf("========================================"))
    message(sprintf("HIGH-PRECISION GI PIPELINE COMPLETED"))
    message(sprintf("Screen: %s, Score: %s", params$screen, params$score))
    total_time <- Sys.time() - start_time
    message(sprintf("Total time: %s", total_time))
    message(sprintf("========================================"))
}

# Run the pipeline if called via Snakemake
if (exists("snakemake")) {
    run_full_gi_pipeline(snakemake)
}