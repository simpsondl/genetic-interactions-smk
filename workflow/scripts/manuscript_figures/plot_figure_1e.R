library(readr)
library(ggplot2)
library(ggrepel)

# read the consolidated construct scores file and subset to the PARP1 per-sgRNA table
parp1_tau_all <- read_tsv(snakemake@input[["input_tau"]], 
                          show_col_types = FALSE)
parp1_tau_model <- parp1_tau_all[parp1_tau_all$query == "PARP1_+_226595744.23-P1P2", ]

cap <- 3.5

parp1_tau_model$GI.capped <- parp1_tau_model$GI.z
parp1_tau_model$GI.capped[parp1_tau_model$GI.capped >= cap] <- cap
parp1_tau_model$GI.capped[parp1_tau_model$GI.capped <= -1 * cap] <- -1 * cap


parp1_tau_model$lb <- NA
parp1_tau_model$lb[parp1_tau_model$sgRNA.id == "PARP2_-_20811828.23-P1P2"] <- "PARP2_1"
parp1_tau_model$lb[parp1_tau_model$sgRNA.id == "PARP2_-_20811834.23-P1P2"] <- "PARP2_2"
parp1_tau_model$lb[parp1_tau_model$sgRNA.id == "BRCA2_-_32889750.23-P1P2"] <- "BRCA2_1"
parp1_tau_model$lb[parp1_tau_model$sgRNA.id == "BRCA2_+_32889673.23-P1P2"] <- "BRCA2_2"
parp1_tau_model$lb[parp1_tau_model$sgRNA.id == "POLB_+_42196051.23-P1P2"] <- "POLB_1"
parp1_tau_model$lb[parp1_tau_model$sgRNA.id == "POLB_+_42196009.23-P1P2"] <- "POLB_2"
parp1_tau_model$lb[parp1_tau_model$sgRNA.id == "FANCA_-_89883018.23-P1P2"] <- "FANCA_1"
parp1_tau_model$lb[parp1_tau_model$sgRNA.id == "FANCA_+_89883007.23-P1P2"] <- "FANCA_2"

tau_fit <- lm(Tau.OI.Avg ~ single + I(single^2), data = parp1_tau_model)

tau_model_line <- data.frame(single = seq(-.5, .25, .01), pred = NA)
tau_model_line$pred <- tau_fit$coefficients[1] + 
  tau_fit$coefficients[2] * tau_model_line$single + 
  tau_fit$coefficients[3] * tau_model_line$single^2
tau_model_line$lb <- NA

tau_plt <- ggplot(parp1_tau_model, aes(single, Tau.OI.Avg, label = lb)) + 
  geom_line(data = tau_model_line, aes(single, pred), 
            col = "darkred", linewidth = 2, alpha = 0.5) +
  geom_point(aes(color = GI.capped), 
             shape = 16) + 
  geom_point(data = parp1_tau_model[!is.na(parp1_tau_model$lb), ], 
             shape = 1, size = 3, col = "black") +
  geom_text_repel(data = parp1_tau_model[!is.na(parp1_tau_model$lb), ], 
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
  xlab("Single sgRNA Tau") +
  ylab("Combination Tau") +
  scale_color_gradient2(name = "IS",
                       low = "dodgerblue", 
                       mid = "#bbbbbb", 
                       high = "darkorange1", 
                       midpoint = mean(parp1_tau_model$Tau.OI.Avg[parp1_tau_model$Control]),
                       breaks = c(-3, -1.5, 0, 1.5, 3))

ggsave(snakemake@output[["output_figure_1e"]], tau_plt, 
       device = "svg", height = 6, width = 7)