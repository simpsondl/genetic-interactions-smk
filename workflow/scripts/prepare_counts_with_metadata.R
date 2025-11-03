library(readr)
library(dplyr)

source("scripts/helper_functions.R")

# Dual logging to both console and log file when running under Snakemake
if (exists("snakemake") && !is.null(snakemake@log) && length(snakemake@log) > 0) {
  source("scripts/dual_logging.R")
  .dual_cleanup <- setup_dual_logging(snakemake@log[[1]])
  on.exit({ 
    .dual_cleanup() 
  }, add = TRUE)
}

message(sprintf("[%s] prepare_counts_with_metadata.R starting", Sys.time()))

# Read input counts file (supports both TSV and ZIP containing TSV)
input_counts_path <- snakemake@input[["input_counts"]]
message(sprintf("[%s] Input counts file: %s", Sys.time(), input_counts_path))

# Support both a plain TSV or a ZIP containing the TSV (expected to contain *_raw_counts.tsv or *_counts.tsv)
if (grepl("\\.zip$", input_counts_path, ignore.case = TRUE)) {
  message(sprintf("[%s] Detected ZIP file, extracting counts TSV", Sys.time()))
  zlist <- utils::unzip(input_counts_path, list = TRUE)
  # prefer a file that ends with _raw_counts.tsv or _counts.tsv
  target_name <- zlist$Name[grepl("_counts(_only|_with_metadata)?\\.tsv$", zlist$Name, ignore.case = TRUE)]
  if (length(target_name) == 0) {
    stop(sprintf("No file matching '*_counts.tsv' or '*_raw_counts.tsv' found inside zip: %s", input_counts_path))
  }
  message(sprintf("[%s] Found %d matching TSV file(s) in ZIP, using: %s", 
                  Sys.time(), length(target_name), target_name[1]))
  # read the first matching entry
  counts <- read_tsv(unz(input_counts_path, target_name[1]), show_col_types = FALSE)
} else {
  counts <- read_tsv(input_counts_path, show_col_types = FALSE)
}

message(sprintf("[%s] Number of rows: %d", Sys.time(), nrow(counts)))
message(sprintf("[%s] Number of columns: %d", Sys.time(), ncol(counts)))
message(sprintf("[%s] Column names: %s", Sys.time(), paste(colnames(counts), collapse = ", ")))

# Get expected columns from config
metadata_cols <- snakemake@params[["required_meta_columns"]]
if (!is.null(metadata_cols) && !is.character(metadata_cols)) {
  metadata_cols <- unlist(metadata_cols)
}

cols_lower <- tolower(colnames(counts))
metadata_cols_lower <- tolower(metadata_cols)
present_metadata <- metadata_cols[metadata_cols_lower %in% cols_lower]
missing_metadata <- metadata_cols[!metadata_cols_lower %in% cols_lower]

message(sprintf("[%s] Metadata status: %d/%d columns present", 
                Sys.time(), length(present_metadata), length(metadata_cols)))
if (length(present_metadata) > 0) {
  message(sprintf("[%s]   Present: %s", Sys.time(), paste(present_metadata, collapse = ", ")))
}
if (length(missing_metadata) > 0) {
  message(sprintf("[%s]   Missing: %s", Sys.time(), paste(missing_metadata, collapse = ", ")))
}

# Determine processing mode
needs_generation <- length(missing_metadata) > 0

if (!needs_generation) {
  message(sprintf("[%s] All metadata columns present - validating and copying to output", Sys.time()))
} else if (length(present_metadata) == 0 || 
           (!("FirstPosition" %in% present_metadata) && !("SecondPosition" %in% present_metadata))) {
  message(sprintf("[%s] No metadata or missing sgRNA IDs - generating all metadata from scratch", Sys.time()))
} else {
  message(sprintf("[%s] Partial metadata present - will generate missing columns", Sys.time()))
}

if (needs_generation) {
  # Need to generate at least some metadata
  message(sprintf("[%s] Generating missing metadata columns", Sys.time()))
  
  # Ensure FirstPosition and SecondPosition exist
  if (!("FirstPosition" %in% colnames(counts)) && !("SecondPosition" %in% colnames(counts))) {
    # Assumes that the first two columns are the sgRNA identity columns
    message(sprintf("[%s] Renaming first two columns to FirstPosition and SecondPosition", Sys.time()))
    colnames(counts)[1:2] <- c("FirstPosition", "SecondPosition")
  }
  
  ######################################################################
  # Generate FirstPseudogene and SecondPseudogene if missing
  ######################################################################
  
  needs_first_pseudogene <- !("FirstPseudogene" %in% colnames(counts))
  needs_second_pseudogene <- !("SecondPseudogene" %in% colnames(counts))
  
  if (needs_first_pseudogene || needs_second_pseudogene) {
    message(sprintf("[%s] Generating Pseudogene mappings (First: %s, Second: %s)", 
                    Sys.time(),
                    ifelse(needs_first_pseudogene, "missing", "present"),
                    ifelse(needs_second_pseudogene, "missing", "present")))
    
    # If one pseudogene column exists, extract the mapping from it
    # Otherwise, generate the mapping from scratch
    if (!needs_first_pseudogene || !needs_second_pseudogene) {
      # At least one pseudogene column exists - extract the mapping
      existing_col <- if (!needs_first_pseudogene) "FirstPseudogene" else "SecondPseudogene"
      position_col <- if (!needs_first_pseudogene) "FirstPosition" else "SecondPosition"
      
      message(sprintf("[%s] Extracting existing mapping from %s", Sys.time(), existing_col))
      
      id_map <- counts %>%
        select(sgRNA_ID = !!sym(position_col), Pseudogene_ID = !!sym(existing_col)) %>%
        distinct()
      
      message(sprintf("[%s] Extracted %d unique sgRNA-to-Pseudogene mappings", 
                      Sys.time(), nrow(id_map)))
      
    } else {
      # Neither column exists - generate mapping from scratch
      message(sprintf("[%s] Generating new Pseudogene mappings from scratch", Sys.time()))
      
      # Generate sgRNA to Pseudogene ID mapping for non-targeting controls
      # Identify non-targeting control ids
      ntc_ids <- unique(c(
        counts$FirstPosition[grepl("^non-targeting_", counts$FirstPosition)],
        counts$SecondPosition[grepl("^non-targeting_", counts$SecondPosition)]
      ))
      
      message(sprintf("[%s] Found %d non-targeting control sgRNA ids", Sys.time(), length(ntc_ids)))
      
      # Randomly assign non-targeting controls to pseudogenes
      set.seed(31)
      # Identify how many pseudogenes are needed
      num_pseudogenes <- ceiling(length(ntc_ids) / 2)
      pseudogene_ids <- paste0("NTPG_", seq(1, num_pseudogenes))
      
      # Ensure that each pseudogene is used at most twice
      pseudogene_pool <- rep(pseudogene_ids, each = 2)
      pseudogene_assignment <- sample(pseudogene_pool, length(ntc_ids), replace = FALSE)
      pseudogene_id_map <- data.frame(
        sgRNA_ID = ntc_ids,
        Pseudogene_ID = pseudogene_assignment
      )
      
      message(sprintf("[%s] Assigned %d pseudogene ids (unique pseudogenes: %d)", 
              Sys.time(), nrow(pseudogene_id_map), length(unique(pseudogene_id_map$Pseudogene_ID))))
      
      # Identify non-control sgRNA ids
      targeting_ids <- unique(c(
        counts$FirstPosition[!counts$FirstPosition %in% ntc_ids],
        counts$SecondPosition[!counts$SecondPosition %in% ntc_ids]
      ))
      
      message(sprintf("[%s] Found %d targeting sgRNA ids", Sys.time(), length(targeting_ids)))
      
      # Extract gene names from targeting sgRNA ids
      # Assumes sgRNA ids are in the format GENE_STRAND_POSITION.VERSION-TRANSCRIPT
      gene_names <- gsub("_.*", "", targeting_ids)
      targeting_id_map <- data.frame(
        sgRNA_ID = targeting_ids,
        Pseudogene_ID = gene_names
      )
      
      message(sprintf("[%s] Built targeting id map with %d rows", Sys.time(), nrow(targeting_id_map)))
      
      # Combine both mappings
      id_map <- rbind(pseudogene_id_map, targeting_id_map)
      
      message(sprintf("[%s] Combined id_map has %d rows (pseudogene mappings: %d, targeting mappings: %d)",
                      Sys.time(), nrow(id_map), nrow(pseudogene_id_map), nrow(targeting_id_map)))
    }
    
    # Add gene names to counts table (only for missing columns)
    if (needs_first_pseudogene) {
      counts <- counts %>%
        left_join(id_map, by = c("FirstPosition" = "sgRNA_ID")) %>%
        rename(FirstPseudogene = Pseudogene_ID)
      message(sprintf("[%s] Added FirstPseudogene column", Sys.time()))
    }
    
    if (needs_second_pseudogene) {
      counts <- counts %>%
        left_join(id_map, by = c("SecondPosition" = "sgRNA_ID")) %>%
        rename(SecondPseudogene = Pseudogene_ID)
      message(sprintf("[%s] Added SecondPseudogene column", Sys.time()))
    }
    
    message(sprintf("[%s] Completed Pseudogene mapping", Sys.time()))
  }
  
  ######################################################################
  # Generate GuideCombinationID if missing
  ######################################################################
  
  if (!("GuideCombinationID" %in% colnames(counts))) {
    message(sprintf("[%s] Generating GuideCombinationID", Sys.time()))
    
    sgrna_combinations <- counts[, c("FirstPosition", "SecondPosition")]
    sgrna_combinations$FirstAlpha <- apply(sgrna_combinations, 1, min)
    sgrna_combinations$SecondAlpha <- apply(sgrna_combinations, 1, max)
    sgrna_combinations_map <- sgrna_combinations %>%
      select(FirstAlpha, SecondAlpha) %>%
      distinct() %>%
      mutate(GuideCombinationID = paste0("sgc_", row_number()))
    sgrna_combinations_map <- sgrna_combinations %>%
      left_join(sgrna_combinations_map, by = c("FirstAlpha", "SecondAlpha")) %>%
      select(FirstPosition, SecondPosition, GuideCombinationID)
    
    message(sprintf("[%s] Generated %d unique sgRNA combination mappings", Sys.time(), nrow(sgrna_combinations_map)))
    
    counts <- counts %>%
      left_join(sgrna_combinations_map, by = c("FirstPosition", "SecondPosition"))
    
    message(sprintf("[%s] Added GuideCombinationID column", Sys.time()))
  }
  
  ######################################################################
  # Generate ConstructID if missing
  ######################################################################
  
  if (!("ConstructID" %in% colnames(counts))) {
    message(sprintf("[%s] Generating ConstructID", Sys.time()))
    counts <- counts %>%
      mutate(ConstructID = paste0("con_", row_number()))
    message(sprintf("[%s] Added ConstructID column", Sys.time()))
  }
  
  ######################################################################
  # Generate PseudogeneCombinationID and PseudogeneCombinationName if missing
  ######################################################################
  
  # Check which pseudogene combination columns are missing
  needs_pgc_id <- !("PseudogeneCombinationID" %in% colnames(counts))
  needs_pgc_name <- !("PseudogeneCombinationName" %in% colnames(counts))
  
  if (needs_pgc_id || needs_pgc_name) {
    message(sprintf("[%s] Generating Pseudogene combination mappings (ID: %s, Name: %s)", 
                    Sys.time(), 
                    ifelse(needs_pgc_id, "missing", "present"),
                    ifelse(needs_pgc_name, "missing", "present")))
    
    gene_combinations_map <- counts %>%
      select(FirstPseudogene, SecondPseudogene) %>%
      mutate(PseudogeneA = apply(.[, c("FirstPseudogene", "SecondPseudogene")], 1, min),
             PseudogeneB = apply(.[, c("FirstPseudogene", "SecondPseudogene")], 1, max)) %>%
      select(PseudogeneA, PseudogeneB) %>%
      distinct() %>%
      mutate(PseudogeneCombinationID = paste0("pgc_", row_number()),
             PseudogeneCombinationName = paste0(PseudogeneA, ":", PseudogeneB)) %>%
      rename(FirstPseudogene = PseudogeneA,
             SecondPseudogene = PseudogeneB)
    
    message(sprintf("[%s] Generated %d unique pseudogene combination mappings", 
                    Sys.time(), nrow(gene_combinations_map)))
    
    # Merge in Pseudogene combination ids
    counts <- counts %>%
      mutate(PseudogeneA = apply(.[, c("FirstPseudogene", "SecondPseudogene")], 1, min),
             PseudogeneB = apply(.[, c("FirstPseudogene", "SecondPseudogene")], 1, max)) %>%
      left_join(gene_combinations_map, by = c("PseudogeneA" = "FirstPseudogene", 
                                               "PseudogeneB" = "SecondPseudogene"),
                suffix = c("", ".new")) %>%
      select(-PseudogeneA, -PseudogeneB)
    
    # Only keep the new columns if they were missing
    if (!needs_pgc_id && "PseudogeneCombinationID.new" %in% colnames(counts)) {
      counts <- counts %>% select(-PseudogeneCombinationID.new)
    } else if ("PseudogeneCombinationID.new" %in% colnames(counts)) {
      counts <- counts %>% 
        select(-PseudogeneCombinationID) %>%
        rename(PseudogeneCombinationID = PseudogeneCombinationID.new)
    }
    
    if (!needs_pgc_name && "PseudogeneCombinationName.new" %in% colnames(counts)) {
      counts <- counts %>% select(-PseudogeneCombinationName.new)
    } else if ("PseudogeneCombinationName.new" %in% colnames(counts)) {
      counts <- counts %>% 
        select(-PseudogeneCombinationName) %>%
        rename(PseudogeneCombinationName = PseudogeneCombinationName.new)
    }
    
    message(sprintf("[%s] Added missing Pseudogene combination columns", Sys.time()))
  }
  
  ######################################################################
  # Generate Category, Control, Identical, Orientation if missing
  ######################################################################
  
  missing_category_cols <- c("Category", "Control", "Identical", "Orientation")[
    !c("Category", "Control", "Identical", "Orientation") %in% colnames(counts)]
  
  if (length(missing_category_cols) > 0) {
    message(sprintf("[%s] Generating category columns: %s", Sys.time(), 
                    paste(missing_category_cols, collapse = ", ")))
    
    counts <- counts %>%
      mutate(
        Category = if ("Category" %in% missing_category_cols) {
          case_when(
            grepl("^non-targeting_", FirstPosition) & grepl("^non-targeting_", SecondPosition) ~ "NT+NT",
            grepl("^non-targeting_", FirstPosition) | grepl("^non-targeting_", SecondPosition) ~ "X+NT",
            FirstPseudogene == SecondPseudogene ~ "X+X",
            TRUE ~ "X+Y"
          )
        } else Category,
        Control = if ("Control" %in% missing_category_cols) {
          grepl("^non-targeting_", FirstPosition) | grepl("^non-targeting_", SecondPosition)
        } else Control,
        Identical = if ("Identical" %in% missing_category_cols) {
          FirstPosition == SecondPosition
        } else Identical,
        Orientation = if ("Orientation" %in% missing_category_cols) {
          ifelse(FirstPosition <= SecondPosition, "AB", "BA")
        } else Orientation
      )
    
    message(sprintf("[%s] Added category columns", Sys.time()))
  }
  
  ######################################################################
  # Arrange columns in standard order
  ######################################################################
  
  # Rearrange columns to put metadata first
  counts_mapped <- counts %>%
    select(FirstPosition, any_of("FirstPseudogene"), 
           SecondPosition, any_of("SecondPseudogene"),
           any_of("ConstructID"), any_of("GuideCombinationID"), 
           any_of("PseudogeneCombinationID"), any_of("PseudogeneCombinationName"),
           any_of("Category"), any_of("Control"), any_of("Identical"), any_of("Orientation"),
           everything())
  
  # Write out counts with metadata
  write_tsv(counts_mapped, snakemake@output[["output_counts_w_metadata"]])
  message(sprintf("[%s] Wrote counts with metadata to %s (%d rows, %d cols)", 
                  Sys.time(), snakemake@output[["output_counts_w_metadata"]], 
                  nrow(counts_mapped), ncol(counts_mapped)))
  
} else {
  # All metadata columns present - validate and copy
  message(sprintf("[%s] Validating metadata columns", Sys.time()))
  
  if (!is.null(metadata_cols) && length(metadata_cols) > 0) {
    # Check for expected columns (case-insensitive)
    cols_lower <- tolower(colnames(counts))
    expected_lower <- tolower(metadata_cols)
    
    missing_cols <- metadata_cols[!expected_lower %in% cols_lower]
    
    if (length(missing_cols) > 0) {
      stop(sprintf("[%s] Input file is missing expected metadata columns: %s", 
                      Sys.time(), paste(missing_cols, collapse = ", ")))
    }
  }
  
  # Write out (this essentially copies the file to the output location)
  write_tsv(counts, snakemake@output[["output_counts_w_metadata"]])
  message(sprintf("[%s] Wrote counts with metadata to %s (%d rows, %d cols)", 
                  Sys.time(), snakemake@output[["output_counts_w_metadata"]], 
                  nrow(counts), ncol(counts)))
}

message(sprintf("[%s] prepare_counts_with_metadata.R completed successfully", Sys.time()))
