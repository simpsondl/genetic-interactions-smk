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

message(sprintf("[%s] generate_id_maps_from_counts_metadata.R starting", Sys.time()))

# Read the counts with metadata file
counts_mapped <- read_tsv(snakemake@input[["input_counts_w_metadata"]], show_col_types = FALSE)

message(sprintf("[%s] Input counts with metadata file: %s", Sys.time(), snakemake@input[["input_counts_w_metadata"]]))
message(sprintf("[%s] Number of rows: %d", Sys.time(), nrow(counts_mapped)))
message(sprintf("[%s] Number of columns: %d", Sys.time(), ncol(counts_mapped)))

######################################################################
# Extract construct map
######################################################################

construct_map <- counts_mapped %>%
  select(FirstPosition, SecondPosition, ConstructID) %>%
  distinct()

message(sprintf("[%s] Extracted construct map with %d unique constructs", 
                Sys.time(), nrow(construct_map)))

######################################################################
# Extract guide combination map
######################################################################

guide_combination_map <- counts_mapped %>%
  select(FirstPosition, SecondPosition, GuideCombinationID) %>%
  distinct()

message(sprintf("[%s] Extracted guide combination map with %d unique combinations", 
                Sys.time(), nrow(guide_combination_map)))

######################################################################
# Extract gene/pseudogene combination map
######################################################################

gene_combination_map <- counts_mapped %>%
  select(FirstPseudogene, SecondPseudogene, PseudogeneCombinationID, PseudogeneCombinationName) %>%
  distinct() %>%
  # Ensure consistent ordering (alphabetical)
  mutate(PseudogeneA = apply(.[, c("FirstPseudogene", "SecondPseudogene")], 1, min),
         PseudogeneB = apply(.[, c("FirstPseudogene", "SecondPseudogene")], 1, max)) %>%
  select(PseudogeneA, PseudogeneB, PseudogeneCombinationID, PseudogeneCombinationName) %>%
  distinct()

message(sprintf("[%s] Extracted gene combination map with %d unique combinations", 
                Sys.time(), nrow(gene_combination_map)))

######################################################################
# Extract pseudogene map (sgRNA to Pseudogene mapping)
######################################################################

# Get mappings from FirstPosition
first_map <- counts_mapped %>%
  select(sgRNA_ID = FirstPosition, Pseudogene_ID = FirstPseudogene) %>%
  distinct()

# Get mappings from SecondPosition
second_map <- counts_mapped %>%
  select(sgRNA_ID = SecondPosition, Pseudogene_ID = SecondPseudogene) %>%
  distinct()

# Combine and deduplicate
pseudogene_map <- bind_rows(first_map, second_map) %>%
  distinct()

message(sprintf("[%s] Extracted pseudogene map with %d unique sgRNA to pseudogene mappings", 
                Sys.time(), nrow(pseudogene_map)))

######################################################################
# Write output files
######################################################################

write_tsv(construct_map, snakemake@output[["output_construct_map"]])
message(sprintf("[%s] Wrote construct map to %s (%d rows)", 
                Sys.time(), snakemake@output[["output_construct_map"]], nrow(construct_map)))

write_tsv(guide_combination_map, snakemake@output[["output_guidecombination_map"]])
message(sprintf("[%s] Wrote guide combination map to %s (%d rows)", 
                Sys.time(), snakemake@output[["output_guidecombination_map"]], nrow(guide_combination_map)))

write_tsv(gene_combination_map, snakemake@output[["output_genecombination_map"]])
message(sprintf("[%s] Wrote gene combination map to %s (%d rows)", 
                Sys.time(), snakemake@output[["output_genecombination_map"]], nrow(gene_combination_map)))

write_tsv(pseudogene_map, snakemake@output[["output_pseudogene_map"]])
message(sprintf("[%s] Wrote pseudogene map to %s (%d rows)", 
                Sys.time(), snakemake@output[["output_pseudogene_map"]], nrow(pseudogene_map)))

message(sprintf("[%s] generate_id_maps_from_counts_metadata.R completed successfully", Sys.time()))
