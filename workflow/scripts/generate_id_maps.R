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

message(sprintf("[%s] generate_id_maps.R starting", Sys.time()))

counts <- read_tsv(snakemake@input[["input_counts"]], show_col_types = FALSE)

message(sprintf("[%s] Input counts file: %s", Sys.time(), snakemake@input[["input_counts"]]))
message(sprintf("[%s] Number of rows: %d", Sys.time(), nrow(counts)))
message(sprintf("[%s] Number of columns: %d", Sys.time(), ncol(counts)))
message(sprintf("[%s] Column names: %s", Sys.time(), paste(colnames(counts), collapse = ", ")))

# Assumes that the first two columns are the sgRNA identity columns
colnames(counts)[1:2] <- c("FirstPosition", "SecondPosition")

######################################################################
# Generate sgRNA to Pseudogene ID mapping for non-targeting controls
######################################################################

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

######################################################################
# Generate sgRNA combination ids
######################################################################

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

######################################################################
# Merge in id mappings
######################################################################

# Add gene names and construct id to counts table
counts_mapped <- counts %>%
  left_join(id_map, by = c("FirstPosition" = "sgRNA_ID")) %>%
  rename(FirstPseudogene = Pseudogene_ID) %>%
  left_join(id_map, by = c("SecondPosition" = "sgRNA_ID")) %>%
  rename(SecondPseudogene = Pseudogene_ID) %>%
    left_join(sgrna_combinations_map, 
                by = c("FirstPosition", "SecondPosition")) %>%
  mutate(ConstructID = paste0("con_", row_number()))

message(sprintf("[%s] counts_mapped has %d rows and %d columns after merging id maps and guide combinations", 
        Sys.time(), nrow(counts_mapped), ncol(counts_mapped)))

######################################################################
# Generate Pseudogene combination ids
######################################################################

gene_combinations_map <- counts_mapped %>%
  select(FirstPseudogene, SecondPseudogene) %>%
  mutate(PseudogeneA = apply(.[, c("FirstPseudogene", "SecondPseudogene")], 1, min),
         PseudogeneB = apply(.[, c("FirstPseudogene", "SecondPseudogene")], 1, max)) %>%
  select(PseudogeneA, PseudogeneB) %>%
  distinct() %>%
  mutate(PseudogeneCombinationID = paste0("pgc_", row_number()),
         PseudogeneCombinationName = paste0(PseudogeneA, ":", PseudogeneB))

message(sprintf("[%s] Generated %d unique pseudogene combination mappings", 
                Sys.time(), nrow(gene_combinations_map)))

# Merge in Pseudogene combination ids
counts_mapped <- counts_mapped %>%
  mutate(PseudogeneA = apply(.[, c("FirstPseudogene", "SecondPseudogene")], 1, min),
         PseudogeneB = apply(.[, c("FirstPseudogene", "SecondPseudogene")], 1, max)) %>%
  left_join(gene_combinations_map, by = c("PseudogeneA", "PseudogeneB")) %>%
  select(-PseudogeneA, -PseudogeneB)

######################################################################
# Add in additional category columns
######################################################################

counts_mapped <- counts_mapped %>%
  mutate(
    Category = case_when(
      grepl("^non-targeting_", FirstPosition) & grepl("^non-targeting_", SecondPosition) ~ "NT+NT",
      grepl("^non-targeting_", FirstPosition) | grepl("^non-targeting_", SecondPosition) ~ "X+NT",
      FirstPseudogene == SecondPseudogene ~ "X+X",
      TRUE ~ "X+Y"
    ),
    Control = grepl("^non-targeting_", FirstPosition) | grepl("^non-targeting_", SecondPosition),
    Identical = FirstPosition == SecondPosition,
    Orientation = ifelse(FirstPosition <= SecondPosition, "AB", "BA")
  )

# Rearrange columns
counts_mapped <- counts_mapped %>%
  select(FirstPosition, FirstPseudogene, SecondPosition, SecondPseudogene,
         ConstructID, GuideCombinationID, PseudogeneCombinationID, PseudogeneCombinationName,
         Category, Control, Identical, Orientation,
         everything())

# Write out
write_tsv(counts_mapped[, c("FirstPosition", "SecondPosition", "ConstructID")], 
          snakemake@output[["output_construct_map"]])
message(sprintf("[%s] Wrote construct map to %s (%d rows)", 
                Sys.time(), snakemake@output[["output_construct_map"]], nrow(counts_mapped)))
                
write_tsv(sgrna_combinations_map, snakemake@output[["output_guidecombination_map"]])
message(sprintf("[%s] Wrote guide combination map to %s (%d rows)", 
        Sys.time(), snakemake@output[["output_guidecombination_map"]], nrow(sgrna_combinations_map)))

write_tsv(gene_combinations_map, snakemake@output[["output_genecombination_map"]])
message(sprintf("[%s] Wrote gene combination map to %s (%d rows)", 
                Sys.time(), snakemake@output[["output_genecombination_map"]], nrow(gene_combinations_map)))

write_tsv(id_map, snakemake@output[["output_pseudogene_map"]])
message(sprintf("[%s] Wrote pseudogene id map to %s (%d rows)", 
        Sys.time(), snakemake@output[["output_pseudogene_map"]], nrow(id_map)))

write_tsv(counts_mapped, snakemake@output[["output_counts_w_metadata"]])
message(sprintf("[%s] Wrote counts with metadata to %s (%d rows, %d cols)", 
                Sys.time(), snakemake@output[["output_counts_w_metadata"]], 
                nrow(counts_mapped), ncol(counts_mapped)))
