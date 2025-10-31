library(readr)
library(dplyr)
library(ggplot2)
library(ggExtra)

nu_gene_gi_r1 <- read_tsv(snakemake@input[["input_nu_r1"]])
nu_gene_gi_r2 <- read_tsv(snakemake@input[["input_nu_r2"]])
nu_sig <- read_tsv(snakemake@input[["input_nu_gene_level"]])

# Merge replicate datasets
gc_nu <- full_join(nu_gene_gi_r1[, c("PseudogeneCombinationID", "Category", "InteractionScore")], 
                   nu_gene_gi_r2[, c("PseudogeneCombinationID", "Category", "InteractionScore")],
                   by = c("PseudogeneCombinationID", "Category"),
                   suffix = c(".R1", ".R2"))


gc_nu$Sig <- gc_nu$PseudogeneCombinationID %in% nu_sig$PseudogeneCombinationID[nu_sig$Hit]

gc_nu$CatSig <- gc_nu$Category
gc_nu$CatSig[gc_nu$Sig & gc_nu$Category != "X+X"] <- gc_nu$Sig[gc_nu$Sig & gc_nu$Category != "X+X"]

# Rearrange data for plotting
tmp_df <- rbind(gc_nu[gc_nu$Category == "X+Y", ],
                gc_nu[gc_nu$Category == "X+NT", ],
                gc_nu[gc_nu$Category == "NT+NT", ],
                gc_nu[gc_nu$CatSig == TRUE, ])

tmp_df$CatSig <- factor(tmp_df$CatSig, levels = c("X+Y", "X+NT", "NT+NT", "TRUE"))

# Sample subset for plotting labeled pngs
tmp_df2 <- tmp_df[sample(nrow(tmp_df), 10000), ] 
tmp_df2 <- rbind(tmp_df2,
                 tmp_df[tmp_df$InteractionScore.R1 == min(tmp_df$InteractionScore.R1), ])
tmp_df2 <- rbind(tmp_df2,
                 tmp_df[tmp_df$InteractionScore.R2 == min(tmp_df$InteractionScore.R2), ])
tmp_df2 <- rbind(tmp_df2,
                 tmp_df[tmp_df$InteractionScore.R1 == max(tmp_df$InteractionScore.R1), ])
tmp_df2 <- rbind(tmp_df2,
                 tmp_df[tmp_df$InteractionScore.R2 == max(tmp_df$InteractionScore.R2), ])
# Check that the sampling chose both categories
stopifnot(length(table(tmp_df2$CatSig)) == 4)

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


# No labels for flat PNG - FIGURE 2C
p2c <- ggplot(tmp_df, aes(InteractionScore.R1, InteractionScore.R2, color = CatSig)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size = .35, alpha = .6, pch = 16) +
  scale_color_manual(values = c("#999999", "#abddde", "#046c9a", "#BB5566")) +
  theme_bw() + 
  nolabel_theme +
  nolabel_axes +
  coord_fixed() +
  scale_x_continuous(breaks = seq(-8, 12, 4)) +
  scale_y_continuous(breaks = seq(-8, 12, 4))


# Add marginal histograms - FIGURE 2C
p2ch <- ggMarginal(p2c, groupColour = TRUE, groupFill = FALSE)


# labels for SVG - FIGURE 2C
p2cl <- ggplot(tmp_df, aes(InteractionScore.R1, InteractionScore.R2, color = CatSig)) + 
  geom_hline(yintercept = 0, alpha = .5) +
  geom_vline(xintercept = 0, alpha = .5) +
  geom_point(size = 1, alpha = .3, pch = 16) +
  scale_color_manual(values = c("#999999", "#abddde", "#046c9a", "#BB5566")) +
  annotate("text", hjust = 0, x = -6, y = 10, 
           label = paste("r[a] ==", round(cor(tmp_df$InteractionScore.R1, 
                                              tmp_df$InteractionScore.R2), 3)), parse = TRUE) +
  annotate("text", hjust = 0, x = -6, y = 9.2, 
           label = paste("r[s] ==", round(cor(tmp_df$InteractionScore.R1[tmp_df$Sig], 
                                              tmp_df$InteractionScore.R2[tmp_df$Sig]), 3)), parse = TRUE) +
  annotate("text", hjust = 0, x = -6, y = 8.4, 
           label = paste("r[s > 0] ==", round(cor(tmp_df$InteractionScore.R1[tmp_df$CatSig %in% c("TRUE") & 
                                                                               tmp_df$InteractionScore.R1 > 0], 
                                                  tmp_df$InteractionScore.R2[tmp_df$CatSig %in% c("TRUE") & 
                                                                               tmp_df$InteractionScore.R1 > 0]), 3)), 
                                                                               parse = TRUE) +
  annotate("text", hjust = 0, x = -6, y = 7.6, 
           label = paste("r[s < 0] ==", round(cor(tmp_df$InteractionScore.R1[tmp_df$CatSig %in% c("TRUE") & 
                                                                               tmp_df$InteractionScore.R1 < 0], 
                                                  tmp_df$InteractionScore.R2[tmp_df$CatSig %in% c("TRUE") & 
                                                                               tmp_df$InteractionScore.R1 < 0]), 3)), 
                                                                               parse = TRUE) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks.length = unit(.2, "cm")) +
  xlab("Nu Gene-Gene Interaction Score Replicate 1") +
  ylab("Nu Gene-Gene Interaction Score Replicate 2") +
  labs(color = "") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  coord_fixed() +
  scale_x_continuous(breaks = seq(-8, 12, 4)) +
  scale_y_continuous(breaks = seq(-8, 12, 4))

# Save plots
ggsave(snakemake@output[["output_figure_2a"]], 
       p2ch, device = "png", width = 5, height = 5, dpi = 300)
ggsave(snakemake@output[["output_figure_2a_labels"]], 
       p2cl, device = "png", width = 6, height = 6, dpi = 300)
