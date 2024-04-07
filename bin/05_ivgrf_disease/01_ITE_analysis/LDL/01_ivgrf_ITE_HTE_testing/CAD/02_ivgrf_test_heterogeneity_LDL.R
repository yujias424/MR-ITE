#' This code is to run a pilot study in using instrumental random forest and causal forest for analyzing the 
#' causal relationship between LDL and heart disease.
#' Testing the presence of heterogeneity.
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
  
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/psm.R")
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/test_heterogeneity.R")
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/find_best_tree.R")
})

#' =======================
#' Coronary Artery Disease
#' =======================
 
message("Starting analyze Coronary Artery Disease.\n")
set.seed(319)

# ==================
# Data Preprocessing
# ==================
# X, Y, T, Z from their original source file
X <- fread("~/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/LDL/CAD/ukbb.covariate.LDL.complete", sep = "\t")
Z <- fread("~/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/LDL/CAD/1e-07/LDL_prs.best", sep = " ") # score_constd
W <- fread("~/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/LDL/CAD/ukbb.phenotype.LDL.complete.mgdL", sep = "\t")
Y <- fread("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid", sep = "\t")

# homonize the ID included in the analysis
Z_id <- unique(Z$IID)
X_id <- unique(X$IID)
W_id <- unique(W$IID)
Y_id <- unique(Y$IID) 

selected_id <- Reduce(intersect, list(X_id, Y_id, W_id, Z_id))

X <- X[X$IID %in% selected_id, ]
Y <- Y[Y$IID %in% selected_id, ]
W <- W[W$IID %in% selected_id, ]
Z <- Z[Z$IID %in% selected_id, ]

# select covariates
selected.covariates <- c("FID", "IID", 
                         "22001-0.0", "4079-0.0", "4080-0.0", "189-0.0", "21022-0.0",                 
                         "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5",                 
                         "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10",
                         "whr", # "23099-0.0", "21001-0.0", 
                         "Blood_pressure_medication", "No_medication", "Insulin", "Cholesterol_lowering_medication", 
                         "Non_alcohol_drinker", "Previous_alcohol_drinker", "Current_alcohol_drinker",
                         "Non_smoker", "Previous_smoker", "Current_smoker",
                         "30630-0.0", "30760-0.0", "30870-0.0", # lipid-related covariates
                         "30700-0.0", "30710-0.0", "30720-0.0", "30730-0.0", "30740-0.0", "30750-0.0", 
                         "30650-0.0", "30660-0.0", "30670-0.0", "30770-0.0",
                         "30810-0.0", "30830-0.0", "30850-0.0", "30860-0.0", "30880-0.0", "30890-0.0") 
X <- dplyr::select(X, all_of(selected.covariates))
X$`22001-0.0` <- as.integer(X$`22001-0.0`)
X$`21022-0.0` <- as.numeric(X$`21022-0.0`)
X$`4079-0.0` <- as.numeric(X$`4079-0.0`)
X$`4080-0.0` <- as.numeric(X$`4080-0.0`)

X <- cbind(X, Y[, c("copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "stroke", "ischemic_stroke")])
Y <- Y[, c("IID", "CAD")]

samples <- X$IID

W.mat <- W[W$IID %in% samples, c("30780-0.0")]
X.mat <- X[X$IID %in% samples, 3:dim(X)[2]]
Y.mat <- Y[Y$IID %in% samples, 2]
Z.mat <- Z[Z$IID %in% samples, 4]

X.mat <- X.mat %>% 
           rename("Gender"="22001-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0", 
                  "Age"="21022-0.0",             
                  "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                  "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10",
                  "Waist-hip-ratio"="whr", # "Body Fat Percentage"="23099-0.0", "BMI"="21001-0.0", 
                  "Blood pressure medication"="Blood_pressure_medication", "No medication"="No_medication", "Insulin"="Insulin", "Cholesterol lowering medication"="Cholesterol_lowering_medication", 
                  "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                  "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                  "Apolipoprotein A"="30630-0.0", "HDL-C"="30760-0.0", "Triglycerides"="30870-0.0", # lipid-related covariates
                  "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", 
                  "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                  "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                  "Total protein"="30860-0.0", "Urate"="30880-0.0", "Vitamin D"="30890-0.0",
                  "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", 
                  "Venous thromboembolism history"="vte", "Chronic obstructive pulmonary disease history"="copd", "Atrial fibrillation history"="af", 
                  "Hemorrhage Stroke history"="hemorrhage_stroke", "Other stroke history"="stroke", "Ischemic Stroke history"="ischemic_stroke")

W.vector <- as.vector(W.mat$`30780-0.0`)
Y.vector <- as.vector(Y.mat$CAD)
Y.vector <- as.integer(Y.vector)
Z.vector <- as.vector(Z.mat$PRS)

W.vector.binary.cutoff <- as.integer(ifelse(W.vector <= 130, 1, 0)) # mimic usage of statin
W.vector.binary.cutoff <- as.integer(W.vector.binary.cutoff)

# ===============
# Formal Analysis
# ===============

# ==================================
# test the presence of heterogeneity
# ==================================
message("\n================\nStart testing the heterogeneity.\n")
# continuous W and Z
perm.pvalue.continuousW <- test_heterogeneity(X.mat, Y.vector, W.vector, Z.vector, 
                                              num.trees = 5000, seed = 309,
                                              sample.fraction = 0.05, min.node.size = 100, ci.group.size = 10, num.strap = 100, min.perm = 50)
print(perm.pvalue.continuousW) 

# binary W and Z
perm.pvalue.binaryW <- test_heterogeneity(X.mat, Y.vector, W.vector.binary.cutoff, Z.vector, 
                                          num.trees = 5000, seed = 309,
                                          sample.fraction = 0.05, min.node.size = 100, ci.group.size = 10, num.strap = 100, min.perm = 50)
print(perm.pvalue.binaryW)
    
permute.test.dat <- data.frame("Perm_Var" = c(perm.pvalue.continuousW[1], perm.pvalue.binaryW[1]),
                               "Perm_Risk" = c(perm.pvalue.continuousW[2], perm.pvalue.binaryW[2]))
row.names(permute.test.dat) <- c("continuous_W", "binary_W")
write.csv(permute.test.dat, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/02_hte_testing/permutationtest_LDL_CAD.csv")
