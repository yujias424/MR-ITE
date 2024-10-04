#' This code is to run a pilot study in using instrumental random forest and causal forest for analyzing the 
#' causal relationship between LDL and CAD.
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
#   source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/test_heterogeneity.R")
  source("/home/yujia/Project/2023-07-20-individual_MR/bin/05_ivgrf_disease/support_func/test_heterogeneity.R")
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/find_best_tree.R")
})

#' =======================
#' Coronary Artery Disease
#' =======================
 
message("Starting analyze Coronary Artery Disease.\n")
set.seed(1)

# ==================
# Data Preprocessing
# ==================
# X, Y, T, Z from their original source file
X.set1 <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.traits.set1.gz", sep = "\t")
X.set2 <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.traits.set2.gz", sep = "\t")
X.set3 <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/LDL/CAD/ukbb.covariate.LDL.set3.gz", sep = "\t")

Z.set1 <- fread("~/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/LDL/CAD/set1/1e-08/LDL_prs.best", sep = " ") # score_constd
Z.set2 <- fread("~/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/LDL/CAD/set2/1e-08/LDL_prs.best", sep = " ") # score_constd
Z.set3 <- fread("~/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/LDL/CAD/set3/1e-08/LDL_prs.best", sep = " ") # score_constd
W <- fread("~/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/LDL/CAD/ukbb.phenotype.LDL.mgdL", sep = "\t")

Y.date1 <- fread("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date1_james2022.gz", sep = "\t")
Y.date2 <- fread("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date2_james2022.gz", sep = "\t")
Y.date3 <- fread("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date3_james2022.gz", sep = "\t")

dim(X.set1); dim(X.set1); dim(X.set1); dim(X.set1); dim(Z.set1); dim(Z.set2); dim(Z.set3); dim(W); dim(Y.date1); dim(Y.date2); dim(Y.date3)

# homonize the ID included in the analysis
selected_id_set1 <- Reduce(intersect, list(X.set1$IID, Z.set1$IID, W$IID, Y.date1$IID))
selected_id_set2 <- Reduce(intersect, list(X.set2$IID, Z.set2$IID, W$IID, Y.date2$IID))
selected_id_set3 <- Reduce(intersect, list(X.set3$IID, Z.set3$IID, W$IID, Y.date3$IID))

selected_id_set1b <- Reduce(intersect, list(X.set1$IID, Z.set1$IID, W$IID, Y.date3$IID))
selected_id_set2b <- Reduce(intersect, list(X.set2$IID, Z.set2$IID, W$IID, Y.date3$IID))

# model 1
X.model1 <- X.set1[X.set1$IID %in% selected_id_set1, ]
Y.model1 <- Y.date1[Y.date1$IID %in% selected_id_set1, ]
W.model1 <- W[W$IID %in% selected_id_set1, ]
Z.model1 <- Z.set1[Z.set1$IID %in% selected_id_set1, ]
message("Patient Info for model 1")
dim(X.model1); dim(Y.model1); dim(W.model1); dim(Z.model1);
print(table(Y.model1$CAD))

# model 1b
X.model1b <- X.set1[X.set1$IID %in% selected_id_set1b, ]
Y.model1b <- Y.date3[Y.date3$IID %in% selected_id_set1b, ]
W.model1b <- W[W$IID %in% selected_id_set1b, ]
Z.model1b <- Z.set1[Z.set1$IID %in% selected_id_set1b, ]
message("Patient Info for model 1b")
dim(X.model1b); dim(X.model1b); dim(X.model1b); dim(X.model1b);
print(table(Y.model1b$CAD))

# model 2
X.model2 <- X.set2[X.set2$IID %in% selected_id_set2, ]
Y.model2 <- Y.date2[Y.date2$IID %in% selected_id_set2, ]
W.model2 <- W[W$IID %in% selected_id_set2, ]
Z.model2 <- Z.set2[Z.set2$IID %in% selected_id_set2, ]
X.model2 <- cbind(X.model2, Y.model2[, c("htn", "t2dm", "heart_failure", "hemorrhage_stroke", "ischemic_stroke")]) # "copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "stroke", "ischemic_stroke"
Y.model2 <- Y.model2[, c("IID", "CAD")]
message("Patient Info for model 2")
dim(X.model2); dim(Y.model2); dim(W.model2); dim(Z.model2);
print(table(Y.model2$CAD))

# model 2b
X.model2b <- X.set2[X.set2$IID %in% selected_id_set2b, ]
Y.model2b <- Y.date3[Y.date3$IID %in% selected_id_set2b, ]
W.model2b <- W[W$IID %in% selected_id_set2b, ]
Z.model2b <- Z.set2[Z.set2$IID %in% selected_id_set2b, ]
X.model2b <- cbind(X.model2b, Y.model2b[, c("htn", "t2dm", "heart_failure", "hemorrhage_stroke", "ischemic_stroke")]) # "copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "stroke", "ischemic_stroke"
Y.model2b <- Y.model2b[, c("IID", "CAD")]
message("Patient Info for model 2b")
dim(Y.model2b); dim(Y.model2b); dim(Y.model2b); dim(Y.model2b);
print(table(Y.model2b$CAD))

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

# generate mat file for three models
# model 1
W.model1.mat <- W.model1[, c("30780-0.0")]
X.model1.mat <- X.model1[, 3:dim(X.model1)[2]]
Y.model1.mat <- Y.model1[, c("CAD")]
Z.model1.mat <- Z.model1[, 4]
dim(W.model1.mat); dim(X.model1.mat); dim(Y.model1.mat); dim(Z.model1.mat)

# model 1b
W.model1b.mat <- W.model1b[, c("30780-0.0")]
X.model1b.mat <- X.model1b[, 3:dim(X.model1b)[2]]
Y.model1b.mat <- Y.model1b[, c("CAD")]
Z.model1b.mat <- Z.model1b[, 4]
dim(W.model1b.mat); dim(X.model1b.mat); dim(Y.model1b.mat); dim(Z.model1b.mat)

# model 2
W.model2.mat <- W.model2[, c("30780-0.0")]
X.model2.mat <- X.model2[, 3:dim(X.model2)[2]]
Y.model2.mat <- Y.model2[, c("CAD")]
Z.model2.mat <- Z.model2[, 4]
dim(W.model2.mat); dim(X.model2.mat); dim(Y.model2.mat); dim(Z.model2.mat)

# model 2b
W.model2b.mat <- W.model2b[, c("30780-0.0")]
X.model2b.mat <- X.model2b[, 3:dim(X.model2b)[2]]
Y.model2b.mat <- Y.model2b[, c("CAD")]
Z.model2b.mat <- Z.model2b[, 4]
dim(W.model2b.mat); dim(X.model2b.mat); dim(Y.model2b.mat); dim(Z.model2b.mat)

# model 3
W.model3.mat <- W.model3[, c("30780-0.0")]
X.model3.mat <- X.model3[, 3:dim(X.model3)[2]]
Y.model3.mat <- Y.model3[, c("CAD")]
Z.model3.mat <- Z.model3[, 4]
dim(W.model3.mat); dim(X.model3.mat); dim(Y.model3.mat); dim(Z.model3.mat)

# rename the mat colnames
X.model1.mat  <- X.model1.mat  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "Genotype Batch"="22000-0.0",
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10")

X.model1b.mat  <- X.model1b.mat  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "Genotype Batch"="22000-0.0",
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10")

X.model2.mat  <- X.model2.mat  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0", 
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0",
                         "BMI"="21001-0.0", "Body Fat Percentage"="23099-0.0",  "Weight"="21002-0.0", # "Waist-hip-ratio"="whr", 
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker")

X.model2b.mat  <- X.model2b.mat  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0", 
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0",
                         "BMI"="21001-0.0", "Body Fat Percentage"="23099-0.0", "Weight"="21002-0.0", # "Waist-hip-ratio"="whr",
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker")

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

# preparing vector for different model
# model 1
W.model1.vector <- as.vector(W.model1.mat$`30780-0.0`)
Y.model1.vector <- as.integer(as.vector(Y.model1.mat$CAD))
Z.model1.vector <- as.vector(Z.model1.mat$PRS)
W.model1.vector.binary <- as.integer(ifelse(W.model1.vector <= 130, 1, 0)) # mimic usage of statin

# model 1b
W.model1b.vector <- as.vector(W.model1b.mat$`30780-0.0`)
Y.model1b.vector <- as.integer(as.vector(Y.model1b.mat$CAD))
Z.model1b.vector <- as.vector(Z.model1b.mat$PRS)
W.model1b.vector.binary <- as.integer(ifelse(W.model1b.vector <= 130, 1, 0)) # mimic usage of statin

# model 2
W.model2.vector <- as.vector(W.model2.mat$`30780-0.0`)
Y.model2.vector <- as.integer(as.vector(Y.model2.mat$CAD))
Z.model2.vector <- as.vector(Z.model2.mat$PRS)
W.model2.vector.binary <- as.integer(ifelse(W.model2.vector <= 130, 1, 0)) # mimic usage of statin

# model 2b
W.model2b.vector <- as.vector(W.model2b.mat$`30780-0.0`)
Y.model2b.vector <- as.integer(as.vector(Y.model2b.mat$CAD))
Z.model2b.vector <- as.vector(Z.model2b.mat$PRS)
W.model2b.vector.binary <- as.integer(ifelse(W.model2b.vector <= 130, 1, 0)) # mimic usage of statin

# model 3
W.model3.vector <- as.vector(W.model3.mat$`30780-0.0`)
Y.model3.vector <- as.integer(as.vector(Y.model3.mat$CAD))
Z.model3.vector <- as.vector(Z.model3.mat$PRS)
W.model3.vector.binary <- as.integer(ifelse(W.model3.vector <= 130, 1, 0)) # mimic usage of statin

# Applying Wu-Hausman Test to exploring the validity of applying instrument in the analysis.
# # model 1
# WHtest.mat.model1 <- cbind(X.model1.mat, Y.model1.mat, W.model1.mat, Z.model1.mat)
# dim(WHtest.mat.model1)
# colnames(WHtest.mat.model1)[dim(WHtest.mat.model1)[2]-1] <- "W"
# res.model1 <- ivreg(CAD ~ . - W - PRS | W | PRS, data = WHtest.mat.model1)
# summary(res.model1, diagnostics = TRUE)

# model 1b
WHtest.mat.model1b <- cbind(X.model1b.mat, Y.model1b.mat, W.model1b.mat, Z.model1b.mat)
dim(WHtest.mat.model1b)
colnames(WHtest.mat.model1b)[dim(WHtest.mat.model1b)[2]-1] <- "W"
res.model1b <- ivreg(CAD ~ . - W - PRS | W | PRS, data = WHtest.mat.model1b)
summary(res.model1b, diagnostics = TRUE)

# # model 2
# WHtest.mat.model2 <- cbind(X.model2.mat, Y.model2.mat, W.model2.mat, Z.model2.mat)
# dim(WHtest.mat.model2)
# colnames(WHtest.mat.model2)[dim(WHtest.mat.model2)[2]-1] <- "W"
# res.model2 <- ivreg(CAD ~ . - W - PRS | W | PRS, data = WHtest.mat.model2)
# summary(res.model2, diagnostics = TRUE)

# model 2b
WHtest.mat.model2b <- cbind(X.model2b.mat, Y.model2b.mat, W.model2b.mat, Z.model2b.mat)
dim(WHtest.mat.model2b)
colnames(WHtest.mat.model2b)[dim(WHtest.mat.model2b)[2]-1] <- "W"
res.model2b <- ivreg(CAD ~ . - W - PRS | W | PRS, data = WHtest.mat.model2b)
summary(res.model2b, diagnostics = TRUE)

# model 3
WHtest.mat.model3 <- cbind(X.model3.mat, Y.model3.mat, W.model3.mat, Z.model3.mat)
dim(WHtest.mat.model3)
colnames(WHtest.mat.model3)[dim(WHtest.mat.model3)[2]-1] <- "W"
res.model3 <- ivreg(CAD ~ . - W - PRS | W | PRS, data = WHtest.mat.model3)
summary(res.model3, diagnostics = TRUE)

# # Check the F-statistics
# dat <- cbind(Y.vector, W.vector, Z.vector, X.mat)
# dat$`Current alcohol drinker` <- NULL
# dat$`Current smoker` <- NULL
# colnames(dat)[1] <- "Y"
# colnames(dat)[2] <- "W"
# colnames(dat)[3] <- "Z"

# fn <- lm(W ~ . - Y - Z, data=dat)
# fs <- lm(W ~ . - Y, data = dat)

# summary(fs)
# waldtest(fs, fn)$F[2]
# waldtest(fs, fn, vcov = vcovHC(fs, type="HC0"))$F[2]

# wald.test.res <- data.frame("Simple F-test" = c(waldtest(fs, fn)$F[2]), "Heteroscedasticity-consistent F-test" = waldtest(fs, fn, vcov = vcovHC(fs, type="HC0"))$F[2])
# write.csv(wald.test.res, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/01_individual_treatment_effect/wald_test.csv")

# ==================================
# test the presence of heterogeneity
# ==================================
message("\n================\nStart testing the heterogeneity.\n")
# continuous W and Z
# perm.pvalue.continuousW.model1 <- test_heterogeneity(X.model1.mat, Y.model1.vector, W.model1.vector, Z.model1.vector,
#                                                       num.trees = 1000, seed = 309,
#                                                       sample.fraction = 0.025, min.node.size = 100, ci.group.size = 2, num.strap = 100, min.perm = 50)
# message("Model 1")
# print(perm.pvalue.continuousW.model1) 

perm.pvalue.continuousW.model1b <- test_heterogeneity(X.model1b.mat, Y.model1b.vector, W.model1b.vector, Z.model1b.vector,
                                                      num.trees = 1000, seed = 309,
                                                      sample.fraction = 0.2, min.node.size = 1000, ci.group.size = 2, num.strap = 100, min.perm = 50)
message("Model 1b")
print(perm.pvalue.continuousW.model1b$p.value) 

# perm.pvalue.continuousW.model2 <- test_heterogeneity(X.model2.mat, Y.model2.vector, W.model2.vector, Z.model2.vector,
#                                                       num.trees = 1000, seed = 309,
#                                                       sample.fraction = 0.025, min.node.size = 100, ci.group.size = 2, num.strap = 100, min.perm = 50)
# message("Model 2")
# print(perm.pvalue.continuousW.model2) 

perm.pvalue.continuousW.model2b <- test_heterogeneity(X.model2b.mat, Y.model2b.vector, W.model2b.vector, Z.model2b.vector,
                                                      num.trees = 1000, seed = 309,
                                                      sample.fraction = 0.025, min.node.size = 100, ci.group.size = 2, num.strap = 100, min.perm = 50)
message("Model 2b")
print(perm.pvalue.continuousW.model2b) 

perm.pvalue.continuousW.model3 <- test_heterogeneity(X.model3.mat, Y.model3.vector, W.model3.vector, Z.model3.vector,
                                                      num.trees = 1000, seed = 309,
                                                      sample.fraction = 0.025, min.node.size = 100, ci.group.size = 2, num.strap = 100, min.perm = 50)
message("Model 3")
print(perm.pvalue.continuousW.model3) 

# ===========
iv.forest.continuousW.model3 <- instrumental_forest(X.model3.mat, Y.model3.vector, W.model3.vector, Z.model3.vector, 
                                             num.threads = 10, num.trees = 5000, 
                                             sample.fraction = 0.025, min.node.size = 100, ci.group.size = 100)
varimp <- variable_importance(iv.forest.continuousW.model3)
selected.idx <- which((varimp > mean(varimp)))
nonpc.idx <- colnames(X.model3.mat)[selected.idx][!startsWith(colnames(X.model3.mat)[selected.idx], "PC")]
selected.idx <- nonpc.idx

perm.pvalue.continuousW.model3 <- test_heterogeneity(X.model3.mat[, ..selected.idx], Y.model3.vector, W.model3.vector, Z.model3.vector,
                                                     Y.hat = iv.forest.continuousW.model3$Y.hat, 
                                                     W.hat = iv.forest.continuousW.model3$W.hat, 
                                                     Z.hat = iv.forest.continuousW.model3$Z.hat, 
                                                     num.trees = 1000, seed = 309,
                                                     sample.fraction = 0.025, min.node.size = 100, ci.group.size = 2, num.strap = 100, min.perm = 50)
message("Model 3")
print(perm.pvalue.continuousW.model3) 
# ===========

    
# permute.test.dat.continuousW <- data.frame("Perm_Var" = c(perm.pvalue.continuousW.model1[1], perm.pvalue.continuousW.model1b[1], perm.pvalue.continuousW.model2[1], perm.pvalue.continuousW.model2b[1], perm.pvalue.continuousW.model3[1]),
#                                            "Perm_Risk" = c(perm.pvalue.continuousW.model1[2], perm.pvalue.continuousW.model1b[2], perm.pvalue.continuousW.model2[2], perm.pvalue.continuousW.model2b[2], perm.pvalue.continuousW.model3[2]))
# row.names(permute.test.dat.continuousW) <- c("Model 1", "Model 1b", "Model 2", "Model 2b", "Model 3")
# write.csv(permute.test.dat.continuousW, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/02_hte_testing/permutationtest_LDL_CAD_continuousW.csv")

permute.test.dat.continuousW <- data.frame("Perm_Var" = c(perm.pvalue.continuousW.model1b[1], perm.pvalue.continuousW.model2b[1], perm.pvalue.continuousW.model3[1]),
                                           "Perm_Risk" = c(perm.pvalue.continuousW.model1b[2], perm.pvalue.continuousW.model2b[2], perm.pvalue.continuousW.model3[2]))
row.names(permute.test.dat.continuousW) <- c("Model 1b", "Model 2b", "Model 3")
write.csv(permute.test.dat.continuousW, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/02_hte_testing/permutationtest_LDL_CAD_continuousW.csv")

# binary W and Z
# perm.pvalue.binaryW.model1 <- test_heterogeneity(X.model1.mat, Y.model1.vector, W.model1.vector.binary, Z.model1.vector,
#                                           num.trees = 1000, seed = 309,
#                                           sample.fraction = 0.025, min.node.size = 100, ci.group.size = 2, num.strap = 100, min.perm = 50)
# message("Model 1")
# print(perm.pvalue.binaryW.model1)

perm.pvalue.binaryW.model1b <- test_heterogeneity(X.model1b.mat, Y.model1b.vector, W.model1b.vector.binary, Z.model1b.vector,
                                          num.trees = 1000, seed = 309,
                                          sample.fraction = 0.025, min.node.size = 100, ci.group.size = 2, num.strap = 100, min.perm = 50)
message("Model 1b")
print(perm.pvalue.binaryW.model1b)

# perm.pvalue.binaryW.model2 <- test_heterogeneity(X.model2.mat, Y.model2.vector, W.model2.vector.binary, Z.model2.vector,
#                                                       num.trees = 1000, seed = 309,
#                                                       sample.fraction = 0.025, min.node.size = 100, ci.group.size = 2, num.strap = 100, min.perm = 50)
# message("Model 2")
# print(perm.pvalue.binaryW.model2) 

perm.pvalue.binaryW.model2b <- test_heterogeneity(X.model2b.mat, Y.model2b.vector, W.model2b.vector.binary, Z.model2b.vector,
                                                      num.trees = 1000, seed = 309,
                                                      sample.fraction = 0.025, min.node.size = 100, ci.group.size = 2, num.strap = 100, min.perm = 50)
message("Model 2b")
print(perm.pvalue.binaryW.model2b) 

perm.pvalue.binaryW.model3 <- test_heterogeneity(X.model3.mat, Y.model3.vector, W.model3.vector.binary, Z.model3.vector,
                                                      num.trees = 1000, seed = 309,
                                                      sample.fraction = 0.025, min.node.size = 100, ci.group.size = 2, num.strap = 100, min.perm = 50)
message("Model 3")
print(perm.pvalue.binaryW.model3) 

# permute.test.dat.binaryW <- data.frame("Perm_Var" = c(perm.pvalue.binaryW.model1[1], perm.pvalue.binaryW.model1b[1], perm.pvalue.binaryW.model2[1], perm.pvalue.binaryW.model2b[1], perm.pvalue.binaryW.model3[1]),
#                                        "Perm_Risk" = c(perm.pvalue.binaryW.model1[2], perm.pvalue.binaryW.model1b[2], perm.pvalue.binaryW.model2[2], perm.pvalue.binaryW.model2b[2], perm.pvalue.binaryW.model3[2]))
# row.names(permute.test.dat.binaryW) <- c("Model 1", "Model 1b", "Model 2", "Model 2b", "Model 3")
# write.csv(permute.test.dat.binaryW, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/02_hte_testing/permutationtest_LDL_CAD_binaryW.csv")

permute.test.dat.binaryW <- data.frame("Perm_Var" = c(perm.pvalue.binaryW.model1b[1], perm.pvalue.binaryW.model2b[1], perm.pvalue.binaryW.model3[1]),
                                       "Perm_Risk" = c(perm.pvalue.binaryW.model1b[2], perm.pvalue.binaryW.model2b[2], perm.pvalue.binaryW.model3[2]))
row.names(permute.test.dat.binaryW) <- c("Model 1b", "Model 2b", "Model 3")
write.csv(permute.test.dat.binaryW, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/02_hte_testing/permutationtest_LDL_CAD_binaryW.csv")
