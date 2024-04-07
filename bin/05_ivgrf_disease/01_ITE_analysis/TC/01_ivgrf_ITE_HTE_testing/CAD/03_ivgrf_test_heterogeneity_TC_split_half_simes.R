#' This code is to run a pilot study in using instrumental random forest and causal forest for analyzing the 
#' causal relationship between TC and heart disease.
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
X <- fread("~/Project/2022-09-01-individual_MR/dat/04_pheno_covar_data/TC/ukbb.covariate.TC.complete", sep = "\t")
Z <- fread("~/Project/2022-09-01-individual_MR/dat/06_PRS_calculation/TC/CAD/1e-08/TC_prs.best", sep = " ") # score_constd
W <- fread("~/Project/2022-09-01-individual_MR/dat/04_pheno_covar_data/TC/ukbb.phenotype.TC.complete.mgdL", sep = "\t")
Y <- fread("~/Project/2022-09-01-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid", sep = "\t")

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
                         "30870-0.0", # lipid-related covariates
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

W.mat <- W[W$IID %in% samples, c("30690-0.0")]
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
                  "Triglycerides"="30870-0.0", # lipid-related covariates
                  "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", 
                  "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                  "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                  "Total protein"="30860-0.0", "Urate"="30880-0.0", "Vitamin D"="30890-0.0",
                  "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", 
                  "Venous thromboembolism history"="vte", "Chronic obstructive pulmonary disease history"="copd", "Atrial fibrillation history"="af", 
                  "Hemorrhage Stroke history"="hemorrhage_stroke", "Other stroke history"="stroke", "Ischemic Stroke history"="ischemic_stroke")

W.vector <- as.vector(W.mat$`30690-0.0`)
Y.vector <- as.vector(Y.mat$CAD)
Y.vector <- as.integer(Y.vector)
Z.vector <- as.vector(Z.mat$PRS)

W.vector.binary.cutoff <- as.integer(ifelse(W.vector <= 220, 1, 0)) # mimic usage of statin
W.vector.binary.cutoff <- as.integer(W.vector.binary.cutoff)

check.fstat <- function(X, W, Y, Z, name = NULL){
  dat <- cbind(Y, W, Z, X)
  dat$`Current alcohol drinker` <- NULL
  dat$`Current smoker` <- NULL

  colnames(dat)[1] <- "Y"
  colnames(dat)[2] <- "W"
  colnames(dat)[3] <- "Z"

  fn <- lm(W ~ . - Y - Z, data=dat)
  fs <- lm(W ~ . - Y, data = dat)

  message(paste0("F-test results of ", name, ": "))
  summary(fs)
  message(waldtest(fs, fn)$F[2]) 
  message(waldtest(fs, fn, vcov = vcovHC(fs, type="HC0"))$F[2]) 
}

# use all data to build the model
X.tmp <- X.mat
check.fstat(X.tmp, W.vector, Y.vector, Z.vector, "all")

# continuous W and continuous Z
# perm.pvalue.continuousW.all <- test_heterogeneity(X.tmp, Y.vector, W.vector, Z.vector,
#                                                   num.trees = 5000, seed = 309,
#                                                   sample.fraction = 0.15, min.node.size = 300, ci.group.size = 10, num.strap = 100, min.perm = 50)
# message("Permutation test on continuous treatment.")
# print(perm.pvalue.continuousW.all)

# # binary W and continuous Z
# perm.pvalue.binaryW.all <- test_heterogeneity(X.tmp, Y.vector, W.vector.binary.cutoff, Z.vector,
#                                               num.trees = 5000, seed = 309,
#                                               sample.fraction = 0.15, min.node.size = 300, ci.group.size = 10, num.strap = 100, min.perm = 50)
# message("Permutation test on binary treatment.")
# print(perm.pvalue.binaryW.all)
# message("\n\n=================================")

# initialize the matrix to store the permutation results
iteration_num <- 10
perm_var_mat <- matrix(nrow = iteration_num, ncol = 4)
perm_risk_mat <- matrix(nrow = iteration_num, ncol = 4)

for (inum in 1:iteration_num){

  # split training and testing dataset.
  trainIndex.continuous <- createDataPartition(W.vector, p = .5, 
                                               list = FALSE, 
                                               times = 1)
  
  check.fstat(X.tmp[trainIndex.continuous, ], W.vector[trainIndex.continuous], Y.vector[trainIndex.continuous], Z.vector[trainIndex.continuous], "training")
  check.fstat(X.tmp[-trainIndex.continuous, ], W.vector[-trainIndex.continuous], Y.vector[-trainIndex.continuous], Z.vector[-trainIndex.continuous], "testing")
  
  # continuous W and continuous Z
  perm.pvalue.continuousW.train <- test_heterogeneity(X.tmp[trainIndex.continuous, ], Y.vector[trainIndex.continuous], W.vector[trainIndex.continuous], Z.vector[trainIndex.continuous],
                                                      num.trees = 5000, seed = 309,
                                                      sample.fraction = 0.15, min.node.size = 150, ci.group.size = 10, num.strap = 100, min.perm = 50)
  perm.pvalue.continuousW.test <- test_heterogeneity(X.tmp[-trainIndex.continuous, ], Y.vector[-trainIndex.continuous], W.vector[-trainIndex.continuous], Z.vector[-trainIndex.continuous],
                                                      num.trees = 5000, seed = 309,
                                                      sample.fraction = 0.15, min.node.size = 150, ci.group.size = 10, num.strap = 100, min.perm = 50)
  
  perm_var_mat[inum, 1] <- perm.pvalue.continuousW.train[1]
  perm_var_mat[inum, 2] <- perm.pvalue.continuousW.test[1]
  perm_risk_mat[inum, 1] <- perm.pvalue.continuousW.train[2]
  perm_risk_mat[inum, 2] <- perm.pvalue.continuousW.test[2]
  
  message("Permutation test on continuous treatment.")
  print(perm.pvalue.continuousW.train)
  print(perm.pvalue.continuousW.test)

  # split training and testing dataset.
  trainIndex.binary <- createDataPartition(W.vector.binary.cutoff, p = .5, 
                                           list = FALSE, 
                                           times = 1)
  
  check.fstat(X.tmp[trainIndex.binary, ], W.vector[trainIndex.binary], Y.vector[trainIndex.binary], Z.vector[trainIndex.binary], "training")
  check.fstat(X.tmp[-trainIndex.binary, ], W.vector[-trainIndex.binary], Y.vector[-trainIndex.binary], Z.vector[-trainIndex.binary], "testing")
  
  # binary W and continuous Z
  perm.pvalue.binaryW.train <- test_heterogeneity(X.tmp[trainIndex.binary, ], Y.vector[trainIndex.binary], W.vector.binary.cutoff[trainIndex.binary], Z.vector[trainIndex.binary],
                                                  num.trees = 5000, seed = 309,
                                                  sample.fraction = 0.15, min.node.size = 150, ci.group.size = 10, num.strap = 100, min.perm = 50)
  perm.pvalue.binaryW.test <- test_heterogeneity(X.tmp[-trainIndex.binary, ], Y.vector[-trainIndex.binary], W.vector.binary.cutoff[-trainIndex.binary], Z.vector[-trainIndex.binary],
                                                  num.trees = 5000, seed = 309,
                                                  sample.fraction = 0.15, min.node.size = 150, ci.group.size = 10, num.strap = 100, min.perm = 50)
  
  perm_var_mat[inum, 3] <- perm.pvalue.binaryW.train[1]
  perm_var_mat[inum, 4] <- perm.pvalue.binaryW.test[1]
  perm_risk_mat[inum, 3] <- perm.pvalue.binaryW.train[2]
  perm_risk_mat[inum, 4] <- perm.pvalue.binaryW.test[2]
  
  message("Permutation test on binary treatment.")
  print(perm.pvalue.binaryW.train)
  print(perm.pvalue.binaryW.test)
}

perm_var_mat <- as.data.frame(perm_var_mat)
perm_risk_mat <- as.data.frame(perm_risk_mat)

colnames(perm_var_mat) <- c("continuous_train", "continuous_test", "binary_train", "binary_test")
colnames(perm_risk_mat) <- c("continuous_train", "continuous_test", "binary_train", "binary_test")

write.csv(perm_var_mat, "/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/06_validation_test_res/TC_CAD_var.csv", row.names = T)
write.csv(perm_risk_mat, "/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/06_validation_test_res/TC_CAD_risk.csv", row.names = T)

message("           ")
message("\n")
message("\n")

message("\n===================================================\n")
#' =============================================================================================================================================================