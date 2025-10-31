library(readr)
library(ggplot2)
library(dplyr)

# Load functions
source("helper_functions.R")

gamma.gene.gi <- read_tsv("gene_combination_interaction_scores_discriminant_significant_gamma_oi_avg.txt")
tau.gene.gi <- read_tsv("gene_combination_interaction_scores_discriminant_significant_tau_oi_avg.txt")
nu.gene.gi <- read_tsv("gene_combination_interaction_scores_discriminant_significant_nu_oi_avg.txt")
gene.id.map <- read_tsv("pseudogenecombination_id_map.txt")

gi.cmp <-inner_join(gamma.gene.gi[,c(4:6,12)], tau.gene.gi[,c(4:6,12)], by = c("PseudogeneCombinationID", "Category"), 
                    suffix = c(".Gamma", ".Tau"))
gi.cmp <-inner_join(gi.cmp, nu.gene.gi[,c(4:6, 12)], by = c("PseudogeneCombinationID", "Category"))
colnames(gi.cmp)[7:8] <- c("InteractionScore.Nu", "Sig.Nu")

gi.cmp2 <- inner_join(gi.cmp, gene.id.map) %>%
  filter(Category == "X+Y")

gi.cmp2$Highlight <- NA
gi.cmp2$Highlight[gi.cmp2$InteractionScore.Nu > 0 & gi.cmp2$Sig.Nu] <- "Positive"
gi.cmp2$Highlight[gi.cmp2$InteractionScore.Nu < 0 & gi.cmp2$Sig.Nu] <- "Negative"
gi.cmp2$Highlight[is.na(gi.cmp2$Highlight)] <- "None"
gi.cmp2$Highlight <- factor(gi.cmp2$Highlight, levels = c("None",
                                                          "Positive", "Negative"))


gi.cmp2 <- gi.cmp2[order(gi.cmp2$Highlight),]
gi.samp <- gi.cmp2[sample(1:nrow(gi.cmp2),10000),]


plt <- ggplot(gi.cmp2, aes(InteractionScore.Gamma, InteractionScore.Tau)) +
  geom_abline(alpha = .7) +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_point(alpha = .5, size = .8) +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight != "None",], 
             aes(color = Highlight), alpha = .5, size = 1.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        legend.position = "None",
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()) + coord_fixed() +
  scale_x_continuous(limits = c(-12,8), breaks = seq(-12, 8, 4)) +
  scale_y_continuous(limits = c(-11,16.5), breaks = seq(-8, 16, 4)) +
  xlab("") + ylab("") +
  scale_color_manual(values = c("#D81B60", "#33716B"))


plt_lab <- ggplot(gi.samp[gi.samp$Category == "X+Y",], aes(InteractionScore.Gamma, InteractionScore.Tau)) +
  geom_abline(alpha = .7) +
  geom_hline(yintercept = 0, alpha = 0.7) +
  geom_vline(xintercept = 0, alpha = 0.7) +
  geom_point(alpha = .5, size = .8) +
  geom_point(data = gi.cmp2[gi.cmp2$Highlight != "None",], aes(color = Highlight), alpha = .5, size = 1.2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()) + coord_fixed() +
  scale_x_continuous(limits = c(-12,8), breaks = seq(-12, 8, 4)) +
  scale_y_continuous(limits = c(-11,16.5), breaks = seq(-8, 16, 4)) +
  scale_color_manual(values = c("#D81B60", "#33716B")) +
  xlab("Gamma IS") + ylab("Tau IS") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))


ggsave("figure_3a_gamma_tau_interaction_category_scatterplot.png", plt,
       device = "png", height = 5, width = 3)
ggsave("figure_3a_gamma_tau_interaction_category_scatterplot_labels.svg", plt_lab,
       device = "svg", height = 10, width = 6)
