library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggExtra)
library(reshape2)
library(impute)

# Load data
nu_r1 <- read_tsv(snakemake@input[["input_nu_r1"]])
nu_r2 <- read_tsv(snakemake@input[["input_nu_r2"]])
nu_gene <- read_tsv(snakemake@input[["input_nu_gene_level"]])
id_map <- read_tsv(snakemake@input[["input_id_map"]])
#clusts <- read_tsv(snakemake@input[["input_clusters"]])

# identify FA genes
fancs <- clusts$gene[clusts$Cluster == 2]
# Remove these genes for consistency with 2E
fancs <- fancs[!fancs %in% c("CEP89", "RIF1", "SLX4", "MUS81")]
nu_gene$fanc <- nu_gene$Pseudogene1 %in% fancs & nu_gene$Pseudogene2 %in% fancs

# Rearrange data for heatmap creation
# Remove interactions involving non-targeting guides
noncontrol_nu_r1 <- nu_r1[nu_r1$Category == "X+Y", ]
noncontrol_nu_r2 <- nu_r2[nu_r2$Category == "X+Y", ]

noncontrol_nu_r1 <- inner_join(noncontrol_nu_r1, id_map)
noncontrol_nu_r2 <- inner_join(noncontrol_nu_r2, id_map)

# Make a grid with gene names
gene_grid_r1 <- expand.grid(unique(noncontrol_nu_r1$Pseudogene1), unique(noncontrol_nu_r1$Pseudogene1))
gene_grid_r2 <- expand.grid(unique(noncontrol_nu_r2$Pseudogene1), unique(noncontrol_nu_r2$Pseudogene1))
# Merge in interaction scores
gene_grid_r1 <- left_join(gene_grid_r1, 
                          noncontrol_nu_r1, 
                          by = c("Var1"= "Pseudogene2", "Var2" = "Pseudogene1"))
gene_grid_r1 <- left_join(gene_grid_r1, 
                          noncontrol_nu_r1, 
                          by = c("Var1"_= "Pseudogene1", "Var2" = "Pseudogene2"))
gene.grid.r2 <- left_join(gene.grid.r2, 
                          noncontrol.nu.r2, 
                          by = c("Var1" = "Pseudogene2", "Var2" = "Pseudogene1"))
gene.grid.r2 <- left_join(gene.grid.r2, 
                          noncontrol.nu.r2, 
                          by = c("Var1" = "Pseudogene1", "Var2" = "Pseudogene2"))
# Combine the columns
gene_grid_r1$Gene.GI <- ifelse(is.na(gene_grid_r1$InteractionScore.x), 
                               gene_grid_r1$InteractionScore.y, gene_grid_r1$InteractionScore.x)
gene_grid_r1$N <- ifelse(is.na(gene_grid_r1$N.x), 
                         gene_grid_r1$N.y, gene_grid_r1$N.x)
gene_grid_r2$Gene.GI <- ifelse(is.na(gene.grid.r2$InteractionScore.x), 
                               gene.grid.r2$InteractionScore.y, gene.grid.r2$InteractionScore.x)
gene.grid.r2$N <- ifelse(is.na(gene.grid.r2$N.x), 
                         gene.grid.r2$N.y, gene.grid.r2$N.x)


# Remove extraneous columns
gene_grid_r1 <- gene_grid_r1[,c("Var1", "Var2", "Gene.GI", "N")]
gene_grid_r1$Gene.GI <- as.numeric(gene_grid_r1$Gene.GI)
gene_grid_r1$N <- as.numeric(gene_grid_r1$N)
gene_grid_r2 <- gene.grid.r2[,c("Var1", "Var2", "Gene.GI", "N")]
gene.grid.r2$Gene.GI <- as.numeric(gene.grid.r2$Gene.GI)
gene.grid.r2$N <- as.numeric(gene.grid.r2$N)

# Go from long to wide
gene.mtx.r1 <- pivot_wider(gene_grid_r1[,1:3], names_from = Var2, values_from = Gene.GI)
gene.mtx.r1 <- as.data.frame(ge_e.mt_.r1)
gene.mtx.r2 <- pivot_wider(gene.grid.r2[,1:3], names_from = Var2, values_from = Gene.GI)
gene.mtx.r2 <- as.data.frame(gene.mtx.r2)
# Rename rows
rownames(gene.mtx.r1) <- gene.mtx.r1$Var1
rownames(gene.mtx.r2) <- gene.mtx.r2$Var1
# Remove redundant column
gene.mtx.r1 <- gene.mtx.r1[,2:ncol(gene.mtx.r1)]
gene.mtx.r2 <- gene.mtx.r2[,2:ncol(gene.mtx.r2)]


# Impute missing values
gene.mtx.r1 <- impute.knn(as.matrix(gene.mtx.r1))$data
gene.mtx.r2 <- impute.knn(as.matrix(gene.mtx.r2))$data


# Calculate correlations
gene.cor.r1 <- cor(as.matrix(gene.mtx.r1))
gene.cor.r2 <- cor(as.matrix(gene.mtx.r2))
# Wide to long
cors.r1 <- melt(gene.cor.r1)
cors.r2 <- melt(gene.cor.r2)
cors <- inner_join(cors.r1, cors.r2,
                   by = c("Var1", "Var2"),
                   suffix = c(".R1",".R2"))
cors2 <- cors[cors$Var1 != cors$Var2,]
# Remove duplicates
cors2$First <- apply(cors2[,c(1,2)], 1, min)
cors2$Second <- apply(cors2[,c(1,2)], 1, max)
cors2$PseudogeneCombinationName <- paste(cors2$First, cors2$Second, sep = ":")
cors3 <- inner_join(cors2, id.map)
cors3$fanc <- cors3$Pseudogene1 %in% fancs & cors3$Pseudogene2 %in% fancs
cors4 <- unique(cors3[,c("PseudogeneCombinationID", "value.R1", "value.R2", "fanc")])

# Add interaction significance
cors4$Significance <- "NS"
cors4$Significance[cors4$PseudogeneCombinationID %in% 
                     nu.gene$PseudogeneCombinationID[nu.gene$Discriminant >
                                                       quantile(nu.gene$Discriminant[nu.gene$Category %in% 
                                                                                       c("X+NT", "NT+NT")], 1)[[1]]]] <- 
  "Significant"
cors4$Significance[cors4$fanc] <- "FA pathway"

# Rearrange for plotting
cors.tmp <- rbind(cors4[cors4$Significance == "NS",],
                  cors4[cors4$Significance == "Significant",],
                  cors4[cors4$Significance == "FA pathway",])


cors.tmp$Significance <- factor(cors.tmp$Significance, levels = c("NS", "Significant", "FA pathway"))

cors.samp <- cors4[sample(1:nrow(cors4), 10000),]
cors.samp$Significance <- factor(cors.samp$Significance, levels = c("NS", "Significant", "FA pathway"))


# Label-free theme
nolabel_theme <-   theme(axis.text.x = element_blank(),
                         axis.title.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.title.y = element_blank(),
                         panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         panel.grid = element_blank(),
                         panel.border = element_blank(),
                         axis.line = element_line(),
                         axis.ticks.length = unit(.2, "cm"),
                         legend.position = "none")
nolabel_axes <- list(xlab(""), ylab(""), ggtitle("")) 


# Plot
# No labels for flat PNG
plt <- ggplot(cors.tmp, aes(value.R1, value.R2, col = Significance)) +
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(alpha = .7, size = .25, pch = 16) +
  geom_point(data = cors.tmp[cors.tmp$Significance == "Significant",], alpha = .7, size = .3, pch = 16) +
  geom_point(data = cors.tmp[cors.tmp$Significance == "FA pathway",], alpha = 1, size = .5, pch = 16) +
  theme_bw() + 
  nolabel_axes +
  nolabel_theme +
  scale_color_manual(values = c("#999999", "#BB5566", "#33bbee")) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(-.6, .6, .3)) + 
  scale_y_continuous(breaks = seq(-.6, .6, .3))


plt.lab <- ggplot(cors.samp, aes(value.R1, value.R2, col = Significance)) +
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(alpha = .7, size = .25, pch = 16) +
  geom_point(data = cors.tmp[cors.tmp$Significance == "Significant",], alpha = .7, size = .3, pch = 16) +
  geom_point(data = cors.tmp[cors.tmp$Significance == "FA pathway",], alpha = 1, size = .5, pch = 16) +
  theme_bw() + 
  removeGrid() +
  xlab("Replicate 1") + ylab("Replicate 2") +
  scale_color_manual(values = c("#999999", "#BB5566", "#33bbee")) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(-.6, .6, .3)) + 
  scale_y_continuous(breaks = seq(-.6, .6, .3)) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  annotate(geom = "text", x = .3, y = -.25, label = paste0("r = ", 
                                                           round(cor(cors4$value.R1, cors4$value.R2), 3))) +
  annotate(geom = "text", x = .3, y = -.3, 
            label = paste0("r = ", 
                           round(cor(cors4$value.R1[cors4$Significance == "Significant"], 
                           cors4$value.R2[cors4$Significance == "Significant"]), 3)))


# Save plots
ggsave("nu_profile_replicate_correlation_scatterplot.png", plt, 
       device = "png", width = 3, height = 3, dpi = 300)
ggsave("nu_profile_replicate_correlation_scatterplot_labels.svg", plt.lab, 
       device = "svg", width = 8, height = 8, dpi = 300)
