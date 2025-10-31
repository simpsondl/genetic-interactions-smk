library(readr)
library(circlize)
library(dplyr)

# Import data
all.clusts <- read_tsv("nu_clusters.txt")
gamma.gene.gi <- read_tsv("gene_combination_interaction_scores_discriminant_significant_gamma_oi_avg.txt")
tau.gene.gi <- read_tsv("gene_combination_interaction_scores_discriminant_significant_tau_oi_avg.txt")
nu.gene.gi <- read_tsv("gene_combination_interaction_scores_discriminant_significant_nu_oi_avg.txt")
gene.id.map <- read_tsv("pseudogenecombination_id_map.txt")

# Filter to positive nu
gamma.gene.filt <- gamma.gene.gi[gamma.gene.gi$PseudogeneCombinationID %in%
                                   nu.gene.gi$PseudogeneCombinationID[nu.gene.gi$Sig & 
                                                                        nu.gene.gi$InteractionScore > 0],] %>%
  filter(Category == "X+Y")


tau.gene.filt <- tau.gene.gi[tau.gene.gi$PseudogeneCombinationID %in%
                               nu.gene.gi$PseudogeneCombinationID[nu.gene.gi$Sig & 
                                                                    nu.gene.gi$InteractionScore > 0],] %>%
  filter(Category == "X+Y")


nu.gene.filt <- nu.gene.gi[nu.gene.gi$Sig & nu.gene.gi$InteractionScore > 0,] %>%
  filter(Category == "X+Y")

clust.names <- c("PARP1 interactors", 
                 "FA pathway",
                 "CST pathway",
                 "911 complex",
                 "CIA complex",
                 "RAD54L",
                 "EGFR signaling",
                 "BCDX2",
                 "COP9 signalosome",
                 "BRCA1-A complex",
                 "Protein phosphatase 2A",
                 "MRN complex")


gis <- expand.grid(all.clusts$gene[all.clusts$Cluster > 0], all.clusts$gene[all.clusts$Cluster > 0])
gis$First <- apply(gis[,1:2], 1, min)
gis$Second <- apply(gis[,1:2], 1, max)
gis$name <- paste(gis$First, gis$Second, sep = ":")
gis$Same <- 0


for(i in 1:nrow(gis)){
  j <- all.clusts$Cluster[all.clusts$gene == gis$First[i]]
  k <- all.clusts$Cluster[all.clusts$gene == gis$Second[i]]
  gis$Same[i] <- j == k
}


to.keep <- nu.gene.filt$PseudogeneCombinationName[nu.gene.filt$Sig & nu.gene.filt$InteractionScore > 0 & nu.gene.filt$Category == "X+Y"]
gis <- gis[gis$name %in% to.keep & gis$Same == 0,]
gis <- gis[!duplicated(gis$name),]
gis <- gis[!(gis$First %in% c("PARP1", "PARP2")) & !(gis$Second %in% c("PARP1", "PARP2")),]
gis$gamma <- ifelse(gis$name %in% gamma.gene.gi$PseudogeneCombinationName[gamma.gene.gi$Sig], 1, 0)
gis$tau <- ifelse(gis$name %in% tau.gene.gi$PseudogeneCombinationName[tau.gene.gi$Sig], 1, 0)
gis$nu <- ifelse(gis$name %in% nu.gene.gi$PseudogeneCombinationName[nu.gene.gi$Sig], 1, 0)


gis2 <- left_join(gis[,3:9], gamma.gene.gi[,c(3,6)], by = c("name" = "PseudogeneCombinationName"))
gis2 <- left_join(gis2, tau.gene.gi[,c(3,6)], c("name" = "PseudogeneCombinationName"))
gis2 <- left_join(gis2, nu.gene.gi[,c(3,6)], c("name" = "PseudogeneCombinationName"))


colnames(gis2)[1:2] <- c("from", "to")
grouping <- structure(all.clusts$Cluster[all.clusts$gene %in% unique(c(gis2$from, gis2$to))], names = all.clusts$gene[all.clusts$gene %in% unique(c(gis2$from, gis2$to))])
gene.info <- all.clusts[all.clusts$gene %in% names(grouping),]
gene.info <- gene.info[order(gene.info$Cluster),]
# gene.info$color <- "gray20"
# gene.info$color[gene.info$Cluster %% 2 == 0] <- "gray80"
gene.info$color <- sapply(gene.info$Cluster, 
                          function(x) c("#E5E5E5", "#F06FAA", "#4D2D89",
                                        "#9673B3", "#376DB5", "#70BF44", 
                                        "#BA2C32", "#96D2B0", "#924C21",
                                        "#DA6F27", "#009292", "#7DB2E0")[as.numeric(x)])


gis2$gammacolor <- "#FFFFFF00"
gis2$gammacolor[gis2$gamma == 1 & gis2$InteractionScore.x > 0] <- "#ff7f0030"
gis2$gammacolor[gis2$gamma == 1 & gis2$InteractionScore.x < 0] <- "#1e90ff30"
gis2$tau.color <- "#FFFFFF00"
gis2$nu.color <- "#FFFFFF00"
gis2$tau.color[gis2$tau == 1 & gis2$InteractionScore.y > 0] <- "#ff7f0030"
gis2$tau.color[gis2$tau == 1 & gis2$InteractionScore.y < 0] <- "#1e90ff30"
gis2$nu.color[gis2$nu == 1 & gis2$InteractionScore > 0] <- "#D81B6030"


gis2$gammacolor[(gis2$from == "PARP1" | gis2$to == "PARP1") & gis2$gammacolor == "#1e90ff30"] <- "#1e90ff90"
gis2$tau.color[(gis2$from == "PARP1" | gis2$to == "PARP1") & gis2$tau.color == "#ff7f0030"] <- "#ff7f0090"
gis2$nu.color[(gis2$from == "PARP1" | gis2$to == "PARP1") & gis2$nu.color == "#D81B6030"] <- "#D81B6090"


coloring <- structure(gene.info$color, names = gene.info$gene)


pdf("figure_3b_positivenu_gamma_cluster.pdf", width = 8, height = 8)
chordDiagram(gis2[,c(1,2,6)], 
             group = grouping, 
             grid.col = coloring, 
             small.gap = 0, 
             col = gis2$gammacolor, 
             annotationTrack = "grid", 
             annotationTrackHeight = c(.05), 
             preAllocateTracks = list(track.height = .05))
dev.off()


pdf("figure_3b_positivenu_tau_cluster.pdf", width = 8, height = 8)
chordDiagram(gis2[,c(1,2,6)], 
             group = grouping, 
             grid.col = coloring, 
             small.gap = 0, 
             col = gis2$tau.color, 
             annotationTrack = "grid", 
             annotationTrackHeight = c(.05), 
             preAllocateTracks = list(track.height = .05))
dev.off()


pdf("figure_3b_positivenu_nu_cluster.pdf", width = 8, height = 8)
chordDiagram(gis2[,c(1,2,6)], 
             group = grouping, 
             grid.col = coloring, 
             small.gap = 0, 
             col = gis2$nu.color, 
             annotationTrack = "grid", 
             annotationTrackHeight = c(.05), 
             preAllocateTracks = list(track.height = .05))

dev.off()


svg("nu_cluster_grouplabels.svg", width = 8, height = 8)
chordDiagram(gis2[,c(1,2,6)], 
             group = grouping, 
             grid.col = coloring, 
             small.gap = 0, 
             col = gis2$nu.color, 
             annotationTrack = "grid", 
             annotationTrackHeight = c(.01), 
             preAllocateTracks = list(track.height = .05))
for(i in 1:12){highlight.sector(names(grouping)[grouping == i], 
                                track.index = 1, 
                                col = "white", 
                                cex = .5, 
                                text = clust.names[i], 
                                text.col = "black", 
                                niceFacing = TRUE)}
dev.off()
