library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(impute)
library(colorBlindness)
library(ComplexHeatmap)
library(seriation)
library(circlize)

# Load data
nu.gene.gis <- read_tsv("gene_combination_interaction_scores_discriminant_significant_nu_oi_avg.txt")
id.map <- read_tsv("pseudogenecombination_id_map.txt")
clusters <- read_tsv("nu_clusters.txt")

# Rearrange data for heatmap creation
# Remove interactions involving non-targeting guides
noncontrol.nu.gene.gis <- nu.gene.gis[!grepl("NTPG_", nu.gene.gis$PseudogeneCombinationName),]


# Make a grid with gene names
gene.grid <- expand.grid(unique(noncontrol.nu.gene.gis$Pseudogene1), unique(noncontrol.nu.gene.gis$Pseudogene1))
# Merge in interaction scores
gene.grid <- left_join(gene.grid, 
                       noncontrol.nu.gene.gis, 
                       by = c("Var1" = "Pseudogene2", "Var2" = "Pseudogene1"))
gene.grid <- left_join(gene.grid, 
                       noncontrol.nu.gene.gis, 
                       by = c("Var1" = "Pseudogene1", "Var2" = "Pseudogene2"))
# Combine the columns
gene.grid$Gene.GI <- ifelse(is.na(gene.grid$InteractionScore.x), gene.grid$InteractionScore.y, gene.grid$InteractionScore.x)
gene.grid$N <- ifelse(is.na(gene.grid$N.x), gene.grid$N.y, gene.grid$N.x)


# Remove extraneous columns
gene.grid <- gene.grid[,c("Var1", "Var2", "Gene.GI", "N")]
gene.grid$Gene.GI <- as.numeric(gene.grid$Gene.GI)
gene.grid$N <- as.numeric(gene.grid$N)


# Go from long to wide
gene.mtx <- pivot_wider(gene.grid[,1:3], names_from = Var2, values_from = Gene.GI)
gene.mtx <- as.data.frame(gene.mtx)
### Rename rows
rownames(gene.mtx) <- gene.mtx$Var1
### Remove redundant column
gene.mtx <- gene.mtx[,2:ncol(gene.mtx)]


# Impute missing values
gene.mtx <- impute.knn(as.matrix(gene.mtx))$data


# Get genes in clusters
to.keep <- unique(c(clusters$gene[clusters$Cluster > 0]))


tmp.mtx <- gene.mtx[to.keep,to.keep]
d <- as.dist(1 - cor(as.matrix(tmp.mtx)))
dt <- as.dist(t(1 - cor(as.matrix(tmp.mtx))))
# Rearrange using OLO algorithm from seriation
o1 <- seriate(d, method = "OLO_average")
o2 <- seriate(dt, method = "GW_average")
maxmag <- 3.5
# Define color spectrum ranges
col_fun <- colorRamp2(c(-1 * maxmag, 0, maxmag), c("#33716B", "white", "#D81B60"))

clust <- data.frame(gene = rownames(tmp.mtx))
clust <- left_join(clust, clusters[c("gene", "Cluster")])

clust$Cluster[clust$Cluster == 0] <- NA
clust$gene <- factor(clust$gene, levels = to.keep)
clust <- clust[order(clust$gene),]
clust$Cluster <- factor(as.character(clust$Cluster), 
                        levels = as.character(1:12))

#row
c.both <- rowAnnotation(NuClust = clust$Cluster,
                        show_annotation_name = FALSE, 
                        col = list(NuClust = setNames(c("#E5E5E5", "#F06FAA", "#4D2D89",
                                                        "#9673B3", "#376DB5", "#70BF44", 
                                                        "#BA2C32", "#96D2B0", "#924C21",
                                                        "#DA6F27", "#009292", "#7DB2E0"), as.character(1:12))), 
                        na_col = "white")


c.both.col <- columnAnnotation(NuClust = clust$Cluster,
                               show_annotation_name = FALSE, 
                               show_legend = FALSE,
                               col = list(NuClust = setNames(c("#E5E5E5", "#F06FAA", "#4D2D89",
                                                               "#9673B3", "#376DB5", "#70BF44", 
                                                               "#BA2C32", "#96D2B0", "#924C21",
                                                               "#DA6F27", "#009292", "#7DB2E0"), as.character(1:12))), 
                               na_col = "white")

hm <- Heatmap(t(tmp.mtx), 
              rect_gp = gpar(type = "none"), 
              cluster_rows = as.dendrogram(o1[[1]]), 
              cluster_columns = as.dendrogram(o1[[1]]), 
              show_row_dend = TRUE,
              show_column_dend = FALSE,
              row_dend_side = "left",
              heatmap_legend_param = list(at = c(-3, -1.5, 0, 1.5, 3)), 
              row_names_side = "left",
              cell_fun = function(j, i, x, y, w, h, fill) {
                if(i == j){
                  grid.rect(x,y,w,h,gp = gpar(fill = "#bbbbbb", col = "#bbbbbb"))
                } else if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
                  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                }
              },
              col = col_fun,
              name = "GI",
              left_annotation = c.both, 
              bottom_annotation = c.both.col,
              show_column_names = FALSE,
              show_row_names = FALSE,
              column_names_side = "top",
              row_names_gp = gpar(fontsize = 4),
              column_names_gp = gpar(fontsize = 4),
              column_dend_height = unit(2, "cm"))

pdf("figure_2d_heatmap.pdf", width = 10, height = 10)
draw(hm)
dev.off()
svg("figure_2d_heatmap.svg", width = 10, height = 10)
draw(hm)
dev.off()
