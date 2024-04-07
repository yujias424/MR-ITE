#' This code is to run a pilot study in using instrumental random forest and causal forest for analyzing the 
#' causal relationship between LDL and heart disease.
#' 
#' @author: Yujia Shi
#' @date: 2022.09.01

suppressPackageStartupMessages({
  library(grf)
  library(tidyverse)
  library(data.table)
  library(caret)
  library(ggplot2)
  library(xgboost)
  library(foreach)
  library(doParallel)
  library(DescTools)
  library(patchwork)
  library(latex2exp)
  library(AER)
  library(foreign)
  library(xlsx)
  library(ivreg)
  library(ggsignif) 
  library(dHSIC)
  library(GeneralisedCovarianceMeasure)
  library(chngpt)
  
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/psm.R")
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/test_heterogeneity.R")
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/find_best_tree.R")
})

# ==================
# Data Preprocessing
# ==================
# X, Y, T, Z from their original source file
X.set3 <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/LDL/CAD/ukbb.covariate.LDL.set3.gz", sep = "\t")
Z.set3 <- fread("~/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/LDL/CAD/set3/1e-08/LDL_prs.best", sep = " ") # score_constd
W <- fread("~/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/LDL/CAD/ukbb.phenotype.LDL.mgdL", sep = "\t")
Y.date3 <- fread("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date3_james2022.gz", sep = "\t")

# homonize the ID included in the analysis
selected_id_set3 <- Reduce(intersect, list(X.set3$IID, Z.set3$IID, W$IID, Y.date3$IID))

# model 3
X.model3 <- X.set3[X.set3$IID %in% selected_id_set3, ]
Y.model3 <- Y.date3[Y.date3$IID %in% selected_id_set3, ]
W.model3 <- W[W$IID %in% selected_id_set3, ]
Z.model3 <- Z.set3[Z.set3$IID %in% selected_id_set3, ]
X.model3 <- cbind(X.model3, Y.model3[, c("htn", "t2dm", "heart_failure", "hemorrhage_stroke", "ischemic_stroke")])
Y.model3 <- Y.model3[, c("IID", "CAD")]
message("Patient Info for model 3")
dim(X.model3); dim(Y.model3); dim(W.model3); dim(Z.model3);
print(table(Y.model3$CAD))

# generate mat file for model 3
W.model3.mat <- W.model3[, c("30780-0.0")]
X.model3.mat <- X.model3[, 3:dim(X.model3)[2]]
Y.model3.mat <- Y.model3[, c("CAD")]
Z.model3.mat <- Z.model3[, 4]
dim(W.model3.mat); dim(X.model3.mat); dim(Y.model3.mat); dim(Z.model3.mat)

# rename the mat colnames
X.model3.mat  <- X.model3.mat  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0",   
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0",
                         "BMI"="21001-0.0", "Body Fat Percentage"="23099-0.0", "Weight"="21002-0.0", # "Waist-hip-ratio"="whr", 
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "Apolipoprotein A"="30630-0.0", "HDL-C"="30760-0.0", "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", "Calcium"="30680-0.0",
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0", "Total protein"="30860-0.0", 
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")

# model 3
W.model3.vector <- as.vector(W.model3.mat$`30780-0.0`)
Y.model3.vector <- as.integer(as.vector(Y.model3.mat$CAD))
Z.model3.vector <- as.vector(Z.model3.mat$PRS)
Z.model3.vector.binary <- as.integer(ifelse(Z.model3.mat$PRS <= 0, 1, 0))
W.model3.vector.binary <- as.integer(ifelse(W.model3.vector <= 130, 1, 0)) # mimic usage of statin

# ===============
# Formal Analysis
# ===============

# Continuous W
taus <- read.csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/01_individual_treatment_effect/driv_full_continuousW_te_ul_model3.csv")

threshold_effects_list <- list("Gamma glutamyltransferase" = "segmented",
                               "Body Fat Percentage" = "segmented",
                               "Testosterone" = "segmented",
                               "C-reactive protein" = "segmented",
                               "Calcium" = "segmented",
                               "systolic blood pressure" = "segmented",
                               "Total protein" = "segmented",
                               "Age" = "segmented",
                               "Cystatin C" = "segmented",
                               "BMI" = "segmented")
                               
threshold_regress <- matrix(nrow = 10, ncol = 9)
rindex <- 1

for (i in names(threshold_effects_list)){
    
    message(paste0("Currently running covariate ", i))
    sample_data <- data.frame("taus" = taus$point, "covariate" = X.model3.mat[, ..i])
    colnames(sample_data)[2] <- "covariate"

    first_stage_fit <- chngptm(formula.1=taus~1, formula.2=~covariate, sample_data, type=threshold_effects_list[[i]], family="gaussian")
    second_stage_fit <- lincomb(first_stage_fit, comb=c(0,1,1), alpha=0.05)

    threshold_regress[rindex, 1:3] <- summary(first_stage_fit)$coefficients[2, c("est", "(lower", "upper)")]
    threshold_regress[rindex, 4:6] <- c(second_stage_fit[1],
                                        second_stage_fit[2],
                                        second_stage_fit[3])
    threshold_regress[rindex, 7:9] <- c(summary(first_stage_fit)$chngpt[1],
                                        summary(first_stage_fit)$chngpt[3],
                                        summary(first_stage_fit)$chngpt[4])
    rindex <- rindex + 1
    
}

threshold_regress <- as.data.frame(threshold_regress)
row.names(threshold_regress) <- names(threshold_effects_list)
colnames(threshold_regress) <- c("stage 1 est", "stage 1 lower", "stage 1 upper",
                                 "stage 2 est", "stage 2 lower", "stage 2 upper",
                                 "Threshold", "Threshold lower", "Threshold upper")

write.csv(threshold_regress, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/04_effct_modifiers/segmented_LS/LDL_CAD_continuousW.csv")

# run extra regression forest analysis and ordinary least-squares analysis
extra.forest <- grf::regression_forest(X.model3.mat, taus$point, num.threads = 10, num.trees = 5000, 
                                        sample.fraction = 0.025, min.node.size = 100, ci.group.size = 10)
message("/nRunning the regression forest analysis.")
print(test_calibration(extra.forest))

message("/nRunning the ordinary least-squares analysis.")
tau_vector <- taus$point
X_tmp <- as.matrix(X.model3.mat)
ols.analysis <- coeftest(lm(tau_vector ~ X_tmp), vcov = vcovHC)

ols.analysis.df <- as.data.frame(ols.analysis[,])
rownames(ols.analysis.df)[2:length(rownames(ols.analysis.df))] <- stringr::str_extract_all(rownames(ols.analysis.df)[2:length(rownames(ols.analysis.df))], "(?<=X_tmp).*")
sig.codes <- c()
for (i in 1:dim(ols.analysis.df)[1]){
    if (ols.analysis.df[i, 4]<0.001){
        sig.codes <- c(sig.codes, "***")
    } else if (ols.analysis.df[i, 4]<0.01 & ols.analysis.df[i, 4]>0.001) {
        sig.codes <- c(sig.codes, "**")
    } else if (ols.analysis.df[i, 4]<0.05 & ols.analysis.df[i, 4]>0.01) {
        sig.codes <- c(sig.codes, "*")
    } else if (ols.analysis.df[i, 4]<0.1 & ols.analysis.df[i, 4]>0.05) {
        sig.codes <- c(sig.codes, ".")
    } else {
        sig.codes <- c(sig.codes, " ")
    }
}
ols.analysis.df$sig_code <- sig.codes
ols.analysis.df <- ols.analysis.df[order(abs(ols.analysis.df$Estimate), decreasing = T), ]

write.csv(ols.analysis.df, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/04_effct_modifiers/OLS/LDL_CAD_continuousW.csv")

# Binary W
taus <- read.csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_model3.csv")

threshold_effects_list <- list("Body Fat Percentage" = "segmented",
                               "systolic blood pressure" = "segmented",
                               "Gamma glutamyltransferase" = "segmented",
                               "Calcium" = "segmented",
                               "Total protein" = "segmented",
                               "Vitamin D" = "segmented",
                               "Triglycerides" = "segmented",
                               "Cystatin C" = "segmented",
                               "Testosterone" = "segmented",
                               "HbA1c" = "segmented")
                               
threshold_regress <- matrix(nrow = 10, ncol = 9)
rindex <- 1

for (i in names(threshold_effects_list)){
    
    message(paste0("Currently running covariate ", i))
    sample_data <- data.frame("taus" = taus$point, "covariate" = X.model3.mat[, ..i])
    colnames(sample_data)[2] <- "covariate"

    first_stage_fit <- chngptm(formula.1=taus~1, formula.2=~covariate, sample_data, type=threshold_effects_list[[i]], family="gaussian")
    second_stage_fit <- lincomb(first_stage_fit, comb=c(0,1,1), alpha=0.05)

    threshold_regress[rindex, 1:3] <- summary(first_stage_fit)$coefficients[2, c("est", "(lower", "upper)")]
    threshold_regress[rindex, 4:6] <- c(second_stage_fit[1],
                                        second_stage_fit[2],
                                        second_stage_fit[3])
    threshold_regress[rindex, 7:9] <- c(summary(first_stage_fit)$chngpt[1],
                                        summary(first_stage_fit)$chngpt[3],
                                        summary(first_stage_fit)$chngpt[4])
    rindex <- rindex + 1
    
}

threshold_regress <- as.data.frame(threshold_regress)
row.names(threshold_regress) <- names(threshold_effects_list)
colnames(threshold_regress) <- c("stage 1 est", "stage 1 lower", "stage 1 upper",
                                 "stage 2 est", "stage 2 lower", "stage 2 upper",
                                 "Threshold", "Threshold lower", "Threshold upper")

write.csv(threshold_regress, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/04_effct_modifiers/segmented_LS/LDL_CAD_binaryW.csv")

# run extra regression forest analysis and ordinary least-squares analysis
extra.forest <- grf::regression_forest(X.model3.mat, taus$point, num.threads = 10, num.trees = 5000, 
                                        sample.fraction = 0.025, min.node.size = 100, ci.group.size = 10)
message("/nRunning the regression forest analysis.")
print(test_calibration(extra.forest))

message("/nRunning the ordinary least-squares analysis.")
tau_vector <- taus$point
X_tmp <- as.matrix(X.model3.mat)
ols.analysis <- coeftest(lm(tau_vector ~ X_tmp), vcov = vcovHC)

ols.analysis.df <- as.data.frame(ols.analysis[,])
rownames(ols.analysis.df)[2:length(rownames(ols.analysis.df))] <- stringr::str_extract_all(rownames(ols.analysis.df)[2:length(rownames(ols.analysis.df))], "(?<=X_tmp).*")
sig.codes <- c()
for (i in 1:dim(ols.analysis.df)[1]){
    if (ols.analysis.df[i, 4]<0.001){
        sig.codes <- c(sig.codes, "***")
    } else if (ols.analysis.df[i, 4]<0.01 & ols.analysis.df[i, 4]>0.001) {
        sig.codes <- c(sig.codes, "**")
    } else if (ols.analysis.df[i, 4]<0.05 & ols.analysis.df[i, 4]>0.01) {
        sig.codes <- c(sig.codes, "*")
    } else if (ols.analysis.df[i, 4]<0.1 & ols.analysis.df[i, 4]>0.05) {
        sig.codes <- c(sig.codes, ".")
    } else {
        sig.codes <- c(sig.codes, " ")
    }
}
ols.analysis.df$sig_code <- sig.codes
ols.analysis.df <- ols.analysis.df[order(abs(ols.analysis.df$Estimate), decreasing = T), ]

write.csv(ols.analysis.df, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/04_effct_modifiers/OLS/LDL_CAD_binaryW.csv")
