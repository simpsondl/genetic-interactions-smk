library(readr)
library(ggplot2)
library(ggrepel)

# read the consolidated construct scores file and subset to the PARP1 per-sgRNA table
parp1_gamma_all <- read_tsv(snakemake@input[["input_gamma"]], 
                            show_col_types = FALSE)
parp1_gamma_model <- parp1_gamma_all[parp1_gamma_all$query == "PARP1_+_226595744.23-P1P2", ]

cap <- 3.5

parp1_gamma_model$GI.capped <- parp1_gamma_model$GI.z
parp1_gamma_model$GI.capped[parp1_gamma_model$GI.capped >= cap] <- cap
parp1_gamma_model$GI.capped[parp1_gamma_model$GI.capped <= -1 * cap] <- -1 * cap

parp1_gamma_model$lb <- NA
parp1_gamma_model$lb[parp1_gamma_model$sgRNA.id == "PARP2_-_20811828.23-P1P2"] <- "PARP2_1"
parp1_gamma_model$lb[parp1_gamma_model$sgRNA.id == "PARP2_-_20811834.23-P1P2"] <- "PARP2_2"
parp1_gamma_model$lb[parp1_gamma_model$sgRNA.id == "BRCA2_-_32889750.23-P1P2"] <- "BRCA2_1"
parp1_gamma_model$lb[parp1_gamma_model$sgRNA.id == "BRCA2_+_32889673.23-P1P2"] <- "BRCA2_2"
parp1_gamma_model$lb[parp1_gamma_model$sgRNA.id == "POLB_+_42196051.23-P1P2"] <- "POLB_1"
parp1_gamma_model$lb[parp1_gamma_model$sgRNA.id == "POLB_+_42196009.23-P1P2"] <- "POLB_2"
parp1_gamma_model$lb[parp1_gamma_model$sgRNA.id == "FANCA_-_89883018.23-P1P2"] <- "FANCA_1"
parp1_gamma_model$lb[parp1_gamma_model$sgRNA.id == "FANCA_+_89883007.23-P1P2"] <- "FANCA_2"

gamma_fit <- lm(Gamma.OI.Avg ~ single + I(single^2), data = parp1_gamma_model)

gamma_model_line <- data.frame(single = seq(-.37, .07, .01), pred = NA)
gamma_model_line$pred <- gamma_fit$coefficients[1] + 
  gamma_fit$coefficients[2] * gamma_model_line$single + 
  gamma_fit$coefficients[3] * gamma_model_line$single^2
gamma_model_line$lb <- NA

gamma_plt <- ggplot(parp1_gamma_model, aes(single, Gamma.OI.Avg, label = lb)) + 
  geom_line(data = gamma_model_line, aes(single, pred), 
            col = "darkred", linewidth = 2, alpha = 0.5) +
  geom_point(aes(color = GI.capped), 
             shape = 16) + 
  geom_point(data = parp1_gamma_model[!is.na(parp1_gamma_model$lb), ], 
             shape = 1, size = 3, col = "black") +
  geom_text_repel(data = parp1_gamma_model[!is.na(parp1_gamma_model$lb), ], 
                  size = 3, max.overlaps = Inf, 
                  xlim = c(0.05, NA), col = "black") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.length = unit(.2, "cm"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line()) +
  xlab("Single sgRNA Gamma") +
  ylab("Combination Gamma") +
  scale_color_gradient2(name = "IS",
                       low = "dodgerblue", 
                       mid = "#bbbbbb", 
                       high = "darkorange1", 
                       midpoint = mean(parp1_gamma_model$Gamma.OI.Avg[parp1_gamma_model$Control]),
                       breaks = c(-3, -1.5, 0, 1.5, 3))


ggsave(snakemake@output[["output_figure_1d"]], gamma_plt, 
       device = "svg", height = 6, width = 7)