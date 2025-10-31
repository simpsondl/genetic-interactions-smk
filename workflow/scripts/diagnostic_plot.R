library(readr)
library(WGCNA)
library(dplyr)
library(impute)

gene_gis <- read_tsv(snakemake@input[["input_gene_level_scores"]])
idmap <- read_tsv(snakemake@input[["input_idmap"]])

# Rearrange data for heatmap creation
# Remove interactions involving non-targeting guides
noncontrol_gene_gis <- gene_gis[!grepl("non-targeting", gene_gis$GeneCombinationName),  ]

# Make a grid with gene names
gene_grid <- expand.grid(unique(noncontrol_gene_gis$Gene1), unique(noncontrol_gene_gis$Gene1))

# Merge in interaction scores
gene_grid <- left_join(gene_grid, 
                       noncontrol_gene_gis, 
                       by = c("Var1" = "Gene2", "Var2" = "Gene1"))
gene_grid <- left_join(gene_grid, 
                       noncontrol_gene_gis, 
                       by = c("Var1" = "Gene1", "Var2" = "Gene2"))
# Combine the columns
gene_grid$Gene.GI <- ifelse(is.na(gene_grid$InteractionScore.x), 
                            gene_grid$InteractionScore.y, gene_grid$InteractionScore.x)
gene_grid$N <- ifelse(is.na(gene_grid$N.x), gene_grid$N.y, gene_grid$N.x)

# Remove extraneous columns
gene_grid <- gene_grid[, c("Var1", "Var2", "Gene.GI", "N")]
gene_grid$Gene.GI <- as.numeric(gene_grid$Gene.GI)
gene_grid$N <- as.numeric(gene_grid$N)

# Go from long to wide
gene_mtx <- pivot_wider(gene_grid[, 1:3], names_from = Var2, values_from = Gene.GI)
gene_mtx <- as.data.frame(gene_mtx)
### Rename rows
rownames(gene_mtx) <- gene_mtx$Var1
### Remove redundant column
gene_mtx <- gene_mtx[, 2:ncol(gene_mtx)]

# Impute missing values
gene_mtx <- impute.knn(as.matrix(gene_mtx))$data

# identify soft thresholds
sft <- pickSoftThreshold(gene_mtx, powerVector = 1:20)

# save plot
svg(snakemake@output[["output_diagnostic_plot"]], width = 5, height = 5)
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2])
dev.off()