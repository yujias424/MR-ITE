#' This is the code to check the CI and tau summary of different traits, focusing on DRIV.
#' 
#' @author Yujia Shi
#' @date 2022.08.15

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(cowplot)
})

# ==================================================== LDL ==================================================== 
# =======================
# LDL CAD
# =======================

# DRIV
# # tau.summary.LDL.CAD.train <- read.csv("~/Project/2022-09-01-individual_MR/res/table/LDL/heart_disease/treatment_effect/driv_train_binaryWZ_te_ul.csv")
# # tau.summary.LDL.CAD.test <- read.csv("~/Project/2022-09-01-individual_MR/res/table/LDL/heart_disease/treatment_effect/driv_test_binaryWZ_te_ul.csv")
# # tau.summary.LDL.CAD <- rbind(tau.summary.LDL.CAD.train, tau.summary.LDL.CAD.test)
# 
# tau.summary.LDL.CAD <- read.csv("~/Project/2022-09-01-individual_MR/res/table/LDL/heart_disease/treatment_effect/driv_full_binaryWZ_te_ul.csv")

# IV-GRF
tau.summary.LDL.CAD <- read.csv("~/Project/2022-09-01-individual_MR/res/table/LDL/heart_disease/treatment_effect/ivgrf_full_binaryWZ_var.csv")
colnames(tau.summary.LDL.CAD)[2] <- "point"
tau.summary.LDL.CAD$lower_bound <- tau.summary.LDL.CAD$point - 1.96 * sqrt(tau.summary.LDL.CAD$variance.estimates)
tau.summary.LDL.CAD$upper_bound <- tau.summary.LDL.CAD$point + 1.96 * sqrt(tau.summary.LDL.CAD$variance.estimates)

# get significant p-value element.
tau.summary.LDL.CAD$sigP <- ifelse(tau.summary.LDL.CAD$upper_bound < 0, "Yes", "No")
tau.summary.LDL.CAD <- tau.summary.LDL.CAD[order(tau.summary.LDL.CAD$sigP, tau.summary.LDL.CAD$upper_bound),]
tau.summary.LDL.CAD$sigP <- factor(tau.summary.LDL.CAD$sigP)
row.names(tau.summary.LDL.CAD) <- NULL

ldl_CAD <- ggplot(tau.summary.LDL.CAD, aes(x = as.numeric(rownames(tau.summary.LDL.CAD)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.LDL.CAD)),
                 y = point,
                 shape = sigP))  +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) + 
  labs(title = "Coronary Artery Disease") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = "CLATE") + # labs(y = expression(paste(hat(tau)))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))
# ldl_CAD

# =======================
# LDL IS
# =======================

# DRIV
# # tau.summary.LDL.IS.train <- read.csv("~/Project/2022-09-01-individual_MR/res/table/LDL/ischemic_stroke/treatment_effect/driv_train_binaryWZ_te_ul.csv")
# # tau.summary.LDL.IS.test <- read.csv("~/Project/2022-09-01-individual_MR/res/table/LDL/ischemic_stroke/treatment_effect/driv_test_binaryWZ_te_ul.csv")
# # tau.summary.LDL.IS <- rbind(tau.summary.LDL.IS.train, tau.summary.LDL.IS.test)
# 
# tau.summary.LDL.IS <- read.csv("~/Project/2022-09-01-individual_MR/res/table/LDL/ischemic_stroke/treatment_effect/driv_full_binaryWZ_te_ul.csv")

# IV-GRF
tau.summary.LDL.IS <- read.csv("~/Project/2022-09-01-individual_MR/res/table/LDL/ischemic_stroke/treatment_effect/ivgrf_full_binaryWZ_var.csv")
colnames(tau.summary.LDL.IS)[2] <- "point"
tau.summary.LDL.IS$lower_bound <- tau.summary.LDL.IS$point - 1.96 * sqrt(tau.summary.LDL.IS$variance.estimates)
tau.summary.LDL.IS$upper_bound <- tau.summary.LDL.IS$point + 1.96 * sqrt(tau.summary.LDL.IS$variance.estimates)

tau.summary.LDL.IS$sigP <- ifelse(tau.summary.LDL.IS$upper_bound < 0, "Yes", "No")
tau.summary.LDL.IS <- tau.summary.LDL.IS[order(tau.summary.LDL.IS$sigP, tau.summary.LDL.IS$upper_bound),]
tau.summary.LDL.IS$sigP <- factor(tau.summary.LDL.IS$sigP)
row.names(tau.summary.LDL.IS) <- NULL

ldl_IS <- ggplot(tau.summary.LDL.IS, aes(x = as.numeric(rownames(tau.summary.LDL.IS)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.LDL.IS)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(title = "Ischemic Stroke") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = "CLATE") + # labs(y = expression(paste(hat(tau)))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))
# ldl_IS

# ==================================================== HDL ==================================================== 
# =======================
# HDL CAD
# =======================

# DRIV
# tau.summary.HDL.CAD.train <- read.csv("~/Project/2022-09-01-individual_MR/res/table/HDL/heart_disease/treatment_effect/driv_train_binaryWZ_te_ul.csv")
# tau.summary.HDL.CAD.test <- read.csv("~/Project/2022-09-01-individual_MR/res/table/HDL/heart_disease/treatment_effect/driv_test_binaryWZ_te_ul.csv")
# tau.summary.HDL.CAD <- rbind(tau.summary.HDL.CAD.train, tau.summary.HDL.CAD.test)
# 
# # tau.summary.HDL.CAD <- read.csv("~/Project/2022-09-01-individual_MR/res/table/HDL/heart_disease/treatment_effect/driv_full_binaryWZ_te_ul.csv")

# IV-GRF
tau.summary.HDL.CAD <- read.csv("~/Project/2022-09-01-individual_MR/res/table/HDL/heart_disease/treatment_effect/ivgrf_full_binaryWZ_var.csv")
colnames(tau.summary.HDL.CAD)[2] <- "point"
tau.summary.HDL.CAD$lower_bound <- tau.summary.HDL.CAD$point - 1.96 * sqrt(tau.summary.HDL.CAD$variance.estimates)
tau.summary.HDL.CAD$upper_bound <- tau.summary.HDL.CAD$point + 1.96 * sqrt(tau.summary.HDL.CAD$variance.estimates)

tau.summary.HDL.CAD$sigP <- ifelse(tau.summary.HDL.CAD$upper_bound < 0, "Yes", "No")
tau.summary.HDL.CAD <- tau.summary.HDL.CAD[order(tau.summary.HDL.CAD$sigP, tau.summary.HDL.CAD$upper_bound),]
tau.summary.HDL.CAD$sigP <- factor(tau.summary.HDL.CAD$sigP)
row.names(tau.summary.HDL.CAD) <- NULL

hdl_CAD <- ggplot(tau.summary.HDL.CAD, aes(x = as.numeric(rownames(tau.summary.HDL.CAD)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.HDL.CAD)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(title = "Coronary Artery Disease") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = "CLATE") + # labs(y = expression(paste(hat(tau)))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))
# hdl_CAD

# =======================
# HDL IS
# =======================

# DRIV
# # tau.summary.HDL.IS.train <- read.csv("~/Project/2022-09-01-individual_MR/res/table/HDL/ischemic_stroke/treatment_effect/driv_train_binaryWZ_te_ul.csv")
# # tau.summary.HDL.IS.test <- read.csv("~/Project/2022-09-01-individual_MR/res/table/HDL/ischemic_stroke/treatment_effect/driv_test_binaryWZ_te_ul.csv")
# # tau.summary.HDL.IS <- rbind(tau.summary.HDL.IS.train, tau.summary.HDL.IS.test)
# 
# tau.summary.HDL.IS <- read.csv("~/Project/2022-09-01-individual_MR/res/table/HDL/ischemic_stroke/treatment_effect/driv_full_binaryWZ_te_ul.csv")

# IV-GRF
tau.summary.HDL.IS <- read.csv("~/Project/2022-09-01-individual_MR/res/table/HDL/ischemic_stroke/treatment_effect/ivgrf_full_binaryWZ_var.csv")
colnames(tau.summary.HDL.IS)[2] <- "point"
tau.summary.HDL.IS$lower_bound <- tau.summary.HDL.IS$point - 1.96 * sqrt(tau.summary.HDL.IS$variance.estimates)
tau.summary.HDL.IS$upper_bound <- tau.summary.HDL.IS$point + 1.96 * sqrt(tau.summary.HDL.IS$variance.estimates)

tau.summary.HDL.IS$sigP <- ifelse(tau.summary.HDL.IS$upper_bound < 0, "Yes", "No")
tau.summary.HDL.IS <- tau.summary.HDL.IS[order(tau.summary.HDL.IS$sigP, tau.summary.HDL.IS$upper_bound),]
tau.summary.HDL.IS$sigP <- factor(tau.summary.HDL.IS$sigP)
row.names(tau.summary.HDL.IS) <- NULL

hdl_IS <- ggplot(tau.summary.HDL.IS, aes(x = as.numeric(rownames(tau.summary.HDL.IS)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.HDL.IS)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(title = "Ischemic Stroke") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = "CLATE") + # labs(y = expression(paste(hat(tau)))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))
# hdl_IS

# ==================================================== TG ==================================================== 
# =======================
# TG CAD
# =======================

# DRIV
# tau.summary.TG.CAD.train <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TG/heart_disease/treatment_effect/driv_train_binaryWZ_te_ul.csv")
# tau.summary.TG.CAD.test <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TG/heart_disease/treatment_effect/driv_test_binaryWZ_te_ul.csv")
# tau.summary.TG.CAD <- rbind(tau.summary.TG.CAD.train, tau.summary.TG.CAD.test)
# 
tau.summary.TG.CAD <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TG/heart_disease/treatment_effect/driv_full_binaryWZ_te_ul.csv")

# IV-GRF
# tau.summary.TG.CAD <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TG/heart_disease/treatment_effect/ivgrf_full_binaryWZ_var.csv")
colnames(tau.summary.TG.CAD)[2] <- "point"
tau.summary.TG.CAD$lower_bound <- tau.summary.TG.CAD$point - 1.96 * sqrt(tau.summary.TG.CAD$variance.estimates)
tau.summary.TG.CAD$upper_bound <- tau.summary.TG.CAD$point + 1.96 * sqrt(tau.summary.TG.CAD$variance.estimates)

tau.summary.TG.CAD$sigP <- ifelse(tau.summary.TG.CAD$upper_bound < 0, "Yes", "No")
tau.summary.TG.CAD <- tau.summary.TG.CAD[order(tau.summary.TG.CAD$sigP, tau.summary.TG.CAD$upper_bound),]
tau.summary.TG.CAD$sigP <- factor(tau.summary.TG.CAD$sigP)
row.names(tau.summary.TG.CAD) <- NULL

tg_CAD <- ggplot(tau.summary.TG.CAD, aes(x = as.numeric(rownames(tau.summary.TG.CAD)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.TG.CAD)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(title = "Coronary Artery Disease") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = "CLATE") + # labs(y = expression(paste(hat(tau)))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))
# tg_CAD

# =======================
# TG IS
# =======================

# DRIV
# tau.summary.TG.IS.train <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TG/ischemic_stroke/treatment_effect/driv_train_binaryWZ_te_ul.csv")
# tau.summary.TG.IS.test <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TG/ischemic_stroke/treatment_effect/driv_test_binaryWZ_te_ul.csv")
# tau.summary.TG.IS <- rbind(tau.summary.TG.IS.train, tau.summary.TG.IS.test)
# 
tau.summary.TG.IS <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TG/ischemic_stroke/treatment_effect/driv_full_binaryWZ_te_ul.csv")

# IV-GRF
# tau.summary.TG.IS <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TG/ischemic_stroke/treatment_effect/ivgrf_full_binaryWZ_var.csv")
colnames(tau.summary.TG.IS)[2] <- "point"
tau.summary.TG.IS$lower_bound <- tau.summary.TG.IS$point - 1.96 * sqrt(tau.summary.TG.IS$variance.estimates)
tau.summary.TG.IS$upper_bound <- tau.summary.TG.IS$point + 1.96 * sqrt(tau.summary.TG.IS$variance.estimates)

tau.summary.TG.IS$sigP <- ifelse(tau.summary.TG.IS$upper_bound < 0, "Yes", "No")
tau.summary.TG.IS <- tau.summary.TG.IS[order(tau.summary.TG.IS$sigP, tau.summary.TG.IS$upper_bound),]
tau.summary.TG.IS$sigP <- factor(tau.summary.TG.IS$sigP)
row.names(tau.summary.TG.IS) <- NULL

tg_IS <- ggplot(tau.summary.TG.IS, aes(x = as.numeric(rownames(tau.summary.TG.IS)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.TG.IS)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(y = "CLATE") + # labs(y = expression(paste(hat(tau)))) + 
  labs(title = "Ischemic Stroke") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))
# tg_IS

# ==================================================== TC ==================================================== 
# =======================
# TC CAD
# =======================

# DRIV
# tau.summary.TC.CAD.train <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TC/heart_disease/treatment_effect/driv_train_binaryWZ_te_ul.csv")
# tau.summary.TC.CAD.test <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TC/heart_disease/treatment_effect/driv_test_binaryWZ_te_ul.csv")
# tau.summary.TC.CAD <- rbind(tau.summary.TC.CAD.train, tau.summary.TC.CAD.test)
# 
# # tau.summary.TC.CAD <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TC/heart_disease/treatment_effect/driv_full_binaryWZ_te_ul.csv")

# IV-GRF
tau.summary.TC.CAD <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TC/heart_disease/treatment_effect/ivgrf_full_binaryWZ_var.csv")
colnames(tau.summary.TC.CAD)[2] <- "point"
tau.summary.TC.CAD$lower_bound <- tau.summary.TC.CAD$point - 1.96 * sqrt(tau.summary.TC.CAD$variance.estimates)
tau.summary.TC.CAD$upper_bound <- tau.summary.TC.CAD$point + 1.96 * sqrt(tau.summary.TC.CAD$variance.estimates)

tau.summary.TC.CAD$sigP <- ifelse(tau.summary.TC.CAD$upper_bound < 0, "Yes", "No")
tau.summary.TC.CAD <- tau.summary.TC.CAD[order(tau.summary.TC.CAD$sigP, tau.summary.TC.CAD$upper_bound),]
tau.summary.TC.CAD$sigP <- factor(tau.summary.TC.CAD$sigP)
row.names(tau.summary.TC.CAD) <- NULL

tc_CAD <- ggplot(tau.summary.TC.CAD, aes(x = as.numeric(rownames(tau.summary.TC.CAD)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.TC.CAD)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(title = "Coronary Artery Disease") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = "CLATE") + # labs(y = expression(paste(hat(tau)))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))
# tc_CAD

# =======================
# TC IS
# =======================

# DRIV
# # tau.summary.TC.IS.train <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TC/ischemic_stroke/treatment_effect/driv_train_binaryWZ_te_ul.csv")
# # tau.summary.TC.IS.test <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TC/ischemic_stroke/treatment_effect/driv_test_binaryWZ_te_ul.csv")
# # tau.summary.TC.IS <- rbind(tau.summary.TC.IS.train, tau.summary.TC.IS.test)
# 
# tau.summary.TC.IS <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TC/ischemic_stroke/treatment_effect/driv_full_binaryWZ_te_ul.csv")

# IV-GRF
tau.summary.TC.IS <- read.csv("~/Project/2022-09-01-individual_MR/res/table/TC/ischemic_stroke/treatment_effect/ivgrf_full_binaryWZ_var.csv")
colnames(tau.summary.TC.IS)[2] <- "point"
tau.summary.TC.IS$lower_bound <- tau.summary.TC.IS$point - 1.96 * sqrt(tau.summary.TC.IS$variance.estimates)
tau.summary.TC.IS$upper_bound <- tau.summary.TC.IS$point + 1.96 * sqrt(tau.summary.TC.IS$variance.estimates)

tau.summary.TC.IS$sigP <- ifelse(tau.summary.TC.IS$upper_bound < 0, 1, 0)
tau.summary.TC.IS <- tau.summary.TC.IS[order(tau.summary.TC.IS$sigP, tau.summary.TC.IS$upper_bound),]
tau.summary.TC.IS$sigP <- factor(tau.summary.TC.IS$sigP)
row.names(tau.summary.TC.IS) <- NULL

tc_IS <- ggplot(tau.summary.TC.IS, aes(x = as.numeric(rownames(tau.summary.TC.IS)), y = point)) +
  geom_pointrange(aes(ymin = lower_bound,
                      ymax = upper_bound,
                      colour = sigP), fatten = .5) +
  geom_point(aes(x = as.numeric(rownames(tau.summary.TC.IS)),
                 y = point,
                 shape = sigP)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x=element_blank()) +
  labs(title = "Ischemic Stroke") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(y = "CLATE") + # labs(y = expression(paste(hat(tau)))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + guides(col = guide_legend(title= "P-value < 0.05"), shape = guide_legend(title= "P-value < 0.05"))
# tc_IS

# ==========================
# Make plot
# ==========================
# Merge figure

# WithTG
legend <- get_legend(
  # create some space to the left of the legend
  ldl_CAD + theme(legend.box.margin = margin(0, 0, 0, 12))
)
prow <- plot_grid(ldl_CAD + theme(legend.position="none"),
                  ldl_IS + theme(legend.position="none"),
                  hdl_CAD + theme(legend.position="none"),
                  hdl_IS + theme(legend.position="none"),
                  tc_CAD + theme(legend.position="none"),
                  tc_IS + theme(legend.position="none"), 
                  tg_CAD + theme(legend.position="none"),
                  tg_IS + theme(legend.position="none"),
                  labels = "AUTO", label_size = 12, ncol = 2)
res <- plot_grid(prow, legend, rel_widths = c(3, .4))
res
ggsave("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/withTG/Figure_sup1_binaryCI_ivgrf.png", plot = res, width = 20,
       height = 25, bg = "white")

# WithoutTG
legend <- get_legend(
  # create some space to the left of the legend
  ldl_CAD + theme(legend.box.margin = margin(0, 0, 0, 12))
)
prow <- plot_grid(ldl_CAD + theme(legend.position="none"),
                  ldl_IS + theme(legend.position="none"),
                  hdl_CAD + theme(legend.position="none"),
                  hdl_IS + theme(legend.position="none"),
                  tc_CAD + theme(legend.position="none"),
                  tc_IS + theme(legend.position="none"), # tg_CAD + theme(legend.position="none"), tg_IS + theme(legend.position="none"),
                  labels = "AUTO", label_size = 12, ncol = 2)
res <- plot_grid(prow, legend, rel_widths = c(3, .4))
res
ggsave("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/withoutTG/Figure_sup1_binaryCI_ivgrf.png", plot = res, width = 20,
       height = 20, bg = "white")

table(tau.summary.LDL.CAD$sigP)[2]/33237
table(tau.summary.LDL.IS$sigP)[2]/13308
table(tau.summary.HDL.CAD$sigP)[2]/33237
table(tau.summary.HDL.IS$sigP)[2]/13308
table(tau.summary.TC.CAD$sigP)[2]/33237
table(tau.summary.TC.IS$sigP)[2]/13308
table(tau.summary.TG.CAD$sigP)[2]/44316
table(tau.summary.TG.IS$sigP)[2]/16635

length(which(tau.summary.LDL.CAD$point < 0))/33237
length(which(tau.summary.LDL.IS$point < 0))/13308
length(which(tau.summary.HDL.CAD$point < 0))/33237
length(which(tau.summary.HDL.IS$point < 0))/13308
length(which(tau.summary.TC.CAD$point < 0))/33237
length(which(tau.summary.TC.IS$point < 0))/13308
length(which(tau.summary.TG.CAD$point < 0))/44316
length(which(tau.summary.TG.IS$point < 0))/16635

summary(tau.summary.LDL.CAD$point)
summary(tau.summary.LDL.IS$point)
summary(tau.summary.HDL.CAD$point)
summary(tau.summary.HDL.IS$point)
summary(tau.summary.TC.CAD$point)
summary(tau.summary.TC.IS$point)
summary(tau.summary.TG.CAD$point)
summary(tau.summary.TG.IS$point)