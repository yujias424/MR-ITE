#' This is the code to check the CI and tau summary of different traits, focusing on DRIV.
#' 
#' @author Yujia Shi
#' @date 2022.08.15

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(cowplot)
})

# LDL CHD
tau.summary.LDL.CHD <- read.csv("~/Project/2021-11-10-individual_MR/res/data/LDL/heart_disease/treatment_effect/driv_train_te_ul.csv")
tau.summary.LDL.CHD.test <- read.csv("~/Project/2021-11-10-individual_MR/res/data/LDL/heart_disease/treatment_effect/driv_test_te_ul.csv")
tau.summary.LDL.CHD <- rbind(tau.summary.LDL.CHD, tau.summary.LDL.CHD.test)
tau.summary.LDL.CHD$sigP <- ifelse(tau.summary.LDL.CHD$lower_bound > 0, "Yes", "No")
tau.summary.LDL.CHD <- tau.summary.LDL.CHD[order(tau.summary.LDL.CHD$sigP, tau.summary.LDL.CHD$lower_bound),]
tau.summary.LDL.CHD$sigP <- factor(tau.summary.LDL.CHD$sigP)
row.names(tau.summary.LDL.CHD) <- NULL

ldl_chd <- ggplot(tau.summary.LDL.CHD, aes(x = as.numeric(rownames(tau.summary.LDL.CHD)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.LDL.CHD)),
                 y = point,
                 shape = sigP))  +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) + 
  labs(title = "CAD") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = expression(paste(hat(tau)))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))

tau.summary.LDL.stroke <- read.csv("~/Project/2021-11-10-individual_MR/res/data/LDL/stroke/treatment_effect/driv_train_te_ul.csv")
tau.summary.LDL.stroke.test <- read.csv("~/Project/2021-11-10-individual_MR/res/data/LDL/stroke/treatment_effect/driv_test_te_ul.csv")
tau.summary.LDL.stroke <- rbind(tau.summary.LDL.stroke, tau.summary.LDL.stroke.test)
tau.summary.LDL.stroke$sigP <- ifelse(tau.summary.LDL.stroke$lower_bound > 0, "Yes", "No")
tau.summary.LDL.stroke <- tau.summary.LDL.stroke[order(tau.summary.LDL.stroke$sigP, tau.summary.LDL.stroke$lower_bound),]
tau.summary.LDL.stroke$sigP <- factor(tau.summary.LDL.stroke$sigP)
row.names(tau.summary.LDL.stroke) <- NULL

ldl_stroke <- ggplot(tau.summary.LDL.stroke, aes(x = as.numeric(rownames(tau.summary.LDL.stroke)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.LDL.stroke)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(title = "Stroke") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = expression(paste(hat(tau)))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))

# HDL CHD
tau.summary.HDL.CHD <- read.csv("~/Project/2021-11-10-individual_MR/res/data/HDL/heart_disease/treatment_effect/driv_train_te_ul.csv")
tau.summary.HDL.CHD.test <- read.csv("~/Project/2021-11-10-individual_MR/res/data/HDL/heart_disease/treatment_effect/driv_test_te_ul.csv")
tau.summary.HDL.CHD <- rbind(tau.summary.HDL.CHD, tau.summary.HDL.CHD.test)
tau.summary.HDL.CHD$sigP <- ifelse(tau.summary.HDL.CHD$upper_bound < 0, "Yes", "No")
tau.summary.HDL.CHD <- tau.summary.HDL.CHD[order(tau.summary.HDL.CHD$sigP, tau.summary.HDL.CHD$lower_bound),]
tau.summary.HDL.CHD$sigP <- factor(tau.summary.HDL.CHD$sigP)
row.names(tau.summary.HDL.CHD) <- NULL

hdl_chd <- ggplot(tau.summary.HDL.CHD, aes(x = as.numeric(rownames(tau.summary.HDL.CHD)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.HDL.CHD)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(title = "CAD") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = expression(paste(hat(tau)))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))

tau.summary.HDL.stroke <- read.csv("~/Project/2021-11-10-individual_MR/res/data/HDL/stroke/treatment_effect/driv_train_te_ul.csv")
tau.summary.HDL.stroke.test <- read.csv("~/Project/2021-11-10-individual_MR/res/data/HDL/stroke/treatment_effect/driv_test_te_ul.csv")
tau.summary.HDL.stroke <- rbind(tau.summary.HDL.stroke, tau.summary.HDL.stroke.test)
tau.summary.HDL.stroke$sigP <- ifelse(tau.summary.HDL.stroke$upper_bound < 0, "Yes", "No")
tau.summary.HDL.stroke <- tau.summary.HDL.stroke[order(tau.summary.HDL.stroke$sigP, tau.summary.HDL.stroke$lower_bound),]
tau.summary.HDL.stroke$sigP <- factor(tau.summary.HDL.stroke$sigP)
row.names(tau.summary.HDL.stroke) <- NULL

hdl_stroke <- ggplot(tau.summary.HDL.stroke, aes(x = as.numeric(rownames(tau.summary.HDL.stroke)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.HDL.stroke)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(title = "Stroke") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = expression(paste(hat(tau)))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))

# TG CHD
tau.summary.TG.CHD <- read.csv("~/Project/2021-11-10-individual_MR/res/data/TG/heart_disease/treatment_effect/driv_train_te_ul.csv")
tau.summary.TG.CHD.test <- read.csv("~/Project/2021-11-10-individual_MR/res/data/TG/heart_disease/treatment_effect/driv_test_te_ul.csv")
tau.summary.TG.CHD <- rbind(tau.summary.TG.CHD, tau.summary.TG.CHD.test)
tau.summary.TG.CHD$sigP <- ifelse(tau.summary.TG.CHD$lower_bound > 0, "Yes", "No")
tau.summary.TG.CHD <- tau.summary.TG.CHD[order(tau.summary.TG.CHD$sigP, tau.summary.TG.CHD$lower_bound),]
tau.summary.TG.CHD$sigP <- factor(tau.summary.TG.CHD$sigP)
row.names(tau.summary.TG.CHD) <- NULL

tg_chd <- ggplot(tau.summary.TG.CHD, aes(x = as.numeric(rownames(tau.summary.TG.CHD)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.TG.CHD)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(title = "CAD") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = expression(paste(hat(tau)))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))

tau.summary.TG.stroke <- read.csv("~/Project/2021-11-10-individual_MR/res/data/TG/stroke/treatment_effect/driv_train_te_ul.csv")
tau.summary.TG.stroke.test <- read.csv("~/Project/2021-11-10-individual_MR/res/data/TG/stroke/treatment_effect/driv_test_te_ul.csv")
tau.summary.TG.stroke <- rbind(tau.summary.TG.stroke, tau.summary.TG.stroke.test)
tau.summary.TG.stroke$sigP <- ifelse(tau.summary.TG.stroke$lower_bound > 0, "Yes", "No")
tau.summary.TG.stroke <- tau.summary.TG.stroke[order(tau.summary.TG.stroke$sigP, tau.summary.TG.stroke$lower_bound),]
tau.summary.TG.stroke$sigP <- factor(tau.summary.TG.stroke$sigP)
row.names(tau.summary.TG.stroke) <- NULL

tg_stroke <- ggplot(tau.summary.TG.stroke, aes(x = as.numeric(rownames(tau.summary.TG.stroke)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.TG.stroke)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(y = expression(paste(hat(tau)))) + 
  labs(title = "Stroke") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))

# TC CHD
tau.summary.TC.CHD <- read.csv("~/Project/2021-11-10-individual_MR/res/data/TC/heart_disease/treatment_effect/driv_train_te_ul.csv")
tau.summary.TC.CHD.test <- read.csv("~/Project/2021-11-10-individual_MR/res/data/TC/heart_disease/treatment_effect/driv_test_te_ul.csv")
tau.summary.TC.CHD <- rbind(tau.summary.TC.CHD, tau.summary.TC.CHD.test)
tau.summary.TC.CHD$sigP <- ifelse(tau.summary.TC.CHD$lower_bound > 0, "Yes", "No")
tau.summary.TC.CHD <- tau.summary.TC.CHD[order(tau.summary.TC.CHD$sigP, tau.summary.TC.CHD$lower_bound),]
tau.summary.TC.CHD$sigP <- factor(tau.summary.TC.CHD$sigP)
row.names(tau.summary.TC.CHD) <- NULL

tc_chd <- ggplot(tau.summary.TC.CHD, aes(x = as.numeric(rownames(tau.summary.TC.CHD)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.TC.CHD)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(title = "CAD") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = expression(paste(hat(tau)))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))

tau.summary.TC.stroke <- read.csv("~/Project/2021-11-10-individual_MR/res/data/TC/stroke/treatment_effect/driv_train_te_ul.csv")
tau.summary.TC.stroke.test <- read.csv("~/Project/2021-11-10-individual_MR/res/data/TC/stroke/treatment_effect/driv_test_te_ul.csv")
tau.summary.TC.stroke <- rbind(tau.summary.TC.stroke, tau.summary.TC.stroke.test)
tau.summary.TC.stroke$sigP <- ifelse(tau.summary.TC.stroke$lower_bound > 0, 1, 0)
tau.summary.TC.stroke <- tau.summary.TC.stroke[order(tau.summary.TC.stroke$sigP, tau.summary.TC.stroke$lower_bound),]
tau.summary.TC.stroke$sigP <- factor(tau.summary.TC.stroke$sigP)
row.names(tau.summary.TC.stroke) <- NULL

tc_stroke <- ggplot(tau.summary.TC.stroke, aes(x = as.numeric(rownames(tau.summary.TC.stroke)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.TC.stroke)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(title = "Stroke") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = expression(paste(hat(tau)))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))

# Merge figure
legend <- get_legend(
  # create some space to the left of the legend
  ldl_chd + theme(legend.box.margin = margin(0, 0, 0, 12))
)
prow <- plot_grid(ldl_chd + theme(legend.position="none"),
                  ldl_stroke + theme(legend.position="none"),
                  hdl_chd + theme(legend.position="none"),
                  hdl_stroke + theme(legend.position="none"),
                  tg_chd + theme(legend.position="none"),
                  tg_stroke + theme(legend.position="none"),
                  tc_chd + theme(legend.position="none"),
                  tc_stroke + theme(legend.position="none"),
                  labels = "AUTO", label_size = 12, ncol = 2)
res <- plot_grid(prow, legend, rel_widths = c(3, .4))
ggsave("~/Project/2021-11-10-individual_MR/res/plot/Manuscript/tau_CI.png", plot = res, width = 20,
       height = 20, bg = "white")

table(tau.summary.LDL.CHD$sigP)[2]/22158
table(tau.summary.LDL.stroke$sigP)[2]/10938
table(tau.summary.HDL.CHD$sigP)[2]/22158
table(tau.summary.HDL.stroke$sigP)[2]/10938
table(tau.summary.TG.CHD$sigP)[2]/22158
table(tau.summary.TG.stroke$sigP)[2]/10938
table(tau.summary.TC.CHD$sigP)[2]/22158
table(tau.summary.TC.stroke$sigP)[2]/10938

length(which(tau.summary.LDL.CHD$point > 0))/22158
length(which(tau.summary.LDL.stroke$point > 0))/10938
length(which(tau.summary.HDL.CHD$point < 0))/22158
length(which(tau.summary.HDL.stroke$point < 0))/10938
length(which(tau.summary.TG.CHD$point > 0))/22158
length(which(tau.summary.TG.stroke$point > 0))/10938
length(which(tau.summary.TC.CHD$point > 0))/22158
length(which(tau.summary.TC.stroke$point > 0))/10938