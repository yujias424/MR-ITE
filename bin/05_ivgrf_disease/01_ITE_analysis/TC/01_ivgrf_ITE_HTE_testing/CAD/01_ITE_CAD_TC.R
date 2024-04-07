#' This code is to run a pilot study in using instrumental random forest and causal forest for analyzing the 
#' causal relationship between TC and CAD.
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
set.seed(1)

# ==================
# Data Preprocessing
# ==================
# X, Y, T, Z from their original source file
X.set1 <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.traits.set1.gz", sep = "\t")
X.set2 <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.traits.set2.gz", sep = "\t")
X.set3 <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TC/CAD/ukbb.covariate.TC.set3.gz", sep = "\t")

Z.set1 <- fread("~/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/TC/CAD/set1/1e-08/TC_prs.best", sep = " ") # score_constd
Z.set2 <- fread("~/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/TC/CAD/set2/1e-08/TC_prs.best", sep = " ") # score_constd
Z.set3 <- fread("~/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/TC/CAD/set3/1e-08/TC_prs.best", sep = " ") # score_constd
W <- fread("~/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TC/CAD/ukbb.phenotype.TC.mgdL", sep = "\t")

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
W.model1.mat <- W.model1[, c("30690-0.0")]
X.model1.mat <- X.model1[, 3:dim(X.model1)[2]]
Y.model1.mat <- Y.model1[, c("CAD")]
Z.model1.mat <- Z.model1[, 4]
dim(W.model1.mat); dim(X.model1.mat); dim(Y.model1.mat); dim(Z.model1.mat)

# model 1b
W.model1b.mat <- W.model1b[, c("30690-0.0")]
X.model1b.mat <- X.model1b[, 3:dim(X.model1b)[2]]
Y.model1b.mat <- Y.model1b[, c("CAD")]
Z.model1b.mat <- Z.model1b[, 4]
dim(W.model1b.mat); dim(X.model1b.mat); dim(Y.model1b.mat); dim(Z.model1b.mat)

# model 2
W.model2.mat <- W.model2[, c("30690-0.0")]
X.model2.mat <- X.model2[, 3:dim(X.model2)[2]]
Y.model2.mat <- Y.model2[, c("CAD")]
Z.model2.mat <- Z.model2[, 4]
dim(W.model2.mat); dim(X.model2.mat); dim(Y.model2.mat); dim(Z.model2.mat)

# model 2b
W.model2b.mat <- W.model2b[, c("30690-0.0")]
X.model2b.mat <- X.model2b[, 3:dim(X.model2b)[2]]
Y.model2b.mat <- Y.model2b[, c("CAD")]
Z.model2b.mat <- Z.model2b[, 4]
dim(W.model2b.mat); dim(X.model2b.mat); dim(Y.model2b.mat); dim(Z.model2b.mat)

# model 3
W.model3.mat <- W.model3[, c("30690-0.0")]
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
                         "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", "Calcium"="30680-0.0",
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0", "Total protein"="30860-0.0", 
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")

# preparing vector for different model
# model 1
W.model1.vector <- as.vector(W.model1.mat$`30690-0.0`)
Y.model1.vector <- as.integer(as.vector(Y.model1.mat$CAD))
Z.model1.vector <- as.vector(Z.model1.mat$PRS)
Z.model1.vector.binary <- as.integer(ifelse(Z.model1.mat$PRS <= 0, 1, 0))
W.model1.vector.binary <- as.integer(ifelse(W.model1.vector <= 220, 1, 0)) # mimic usage of statin

# model 1b
W.model1b.vector <- as.vector(W.model1b.mat$`30690-0.0`)
Y.model1b.vector <- as.integer(as.vector(Y.model1b.mat$CAD))
Z.model1b.vector <- as.vector(Z.model1b.mat$PRS)
Z.model1b.vector.binary <- as.integer(ifelse(Z.model1b.mat$PRS <= 0, 1, 0))
W.model1b.vector.binary <- as.integer(ifelse(W.model1b.vector <= 220, 1, 0)) # mimic usage of statin

# model 2
W.model2.vector <- as.vector(W.model2.mat$`30690-0.0`)
Y.model2.vector <- as.integer(as.vector(Y.model2.mat$CAD))
Z.model2.vector <- as.vector(Z.model2.mat$PRS)
Z.model2.vector.binary <- as.integer(ifelse(Z.model2.mat$PRS <= 0, 1, 0))
W.model2.vector.binary <- as.integer(ifelse(W.model2.vector <= 220, 1, 0)) # mimic usage of statin

# model 2b
W.model2b.vector <- as.vector(W.model2b.mat$`30690-0.0`)
Y.model2b.vector <- as.integer(as.vector(Y.model2b.mat$CAD))
Z.model2b.vector <- as.vector(Z.model2b.mat$PRS)
Z.model2b.vector.binary <- as.integer(ifelse(Z.model2b.mat$PRS <= 0, 1, 0))
W.model2b.vector.binary <- as.integer(ifelse(W.model2b.vector <= 220, 1, 0)) # mimic usage of statin

# model 3
W.model3.vector <- as.vector(W.model3.mat$`30690-0.0`)
Y.model3.vector <- as.integer(as.vector(Y.model3.mat$CAD))
Z.model3.vector <- as.vector(Z.model3.mat$PRS)
Z.model3.vector.binary <- as.integer(ifelse(Z.model3.mat$PRS <= 0, 1, 0))
W.model3.vector.binary <- as.integer(ifelse(W.model3.vector <= 220, 1, 0)) # mimic usage of statin

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
# write.csv(wald.test.res, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/wald_test.csv")

# ===============
# Formal Analysis
# ===============

# ========================
####### continuous ####### 
# ======================== 
# # model 1
# # instrumental forest (continuous Z)
# message("\nContinuous W IV-grf model 1")
# iv.forest.continuousW.model1 <- instrumental_forest(X.model1.mat, Y.model1.vector, W.model1.vector, Z.model1.vector, 
#                                              num.threads = 10, num.trees = 5000, 
#                                              sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# iv.pred.continuousW.model1 <- predict(iv.forest.continuousW.model1, estimate.variance = TRUE)
# iv.pred.stats.continuousW.model1 <- compute_stats(iv.pred.continuousW.model1)
# summary(iv.pred.continuousW.model1)
# summary(iv.pred.stats.continuousW.model1)
# write.csv(iv.pred.continuousW.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model1.csv")
# write.csv(iv.pred.stats.continuousW.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model1.csv")

# # save variable importance
# iv.variable_importance.continuousW.model1 <- as.data.frame(variable_importance(iv.forest.continuousW.model1))
# iv.variable_importance.continuousW.model1$Variable <- colnames(X.model1.mat)
# colnames(iv.variable_importance.continuousW.model1) <- c("Importance", "Variable")
# iv.variable_importance.continuousW.model1 <- iv.variable_importance.continuousW.model1[, c("Variable", "Importance")]
# write.csv(iv.variable_importance.continuousW.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_continuousW_variable_importance_model1.csv")

# # causal forest
# message("\nContinuous W CF-grf model 1")
# causal.forest.model1 <- causal_forest(X.model1.mat, Y.model1.vector, W.model1.vector, 
#                                       num.threads = 10, num.trees = 5000, 
#                                       sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# cf.pred.model1 <- predict(causal.forest.model1, estimate.variance = TRUE)
# cf.pred.stats.model1 <- compute_stats(cf.pred.model1)
# summary(cf.pred.model1)
# summary(cf.pred.stats.model1)
# write.csv(cf.pred.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model1.csv")
# write.csv(cf.pred.stats.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model1.csv")

# model 1b
# instrumental forest (continuous Z)
message("\nContinuous W IV-grf model 1b")
iv.forest.continuousW.model1b <- instrumental_forest(X.model1b.mat, Y.model1b.vector, W.model1b.vector, Z.model1b.vector, 
                                             num.threads = 10, num.trees = 5000, 
                                             sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
iv.pred.continuousW.model1b <- predict(iv.forest.continuousW.model1b, estimate.variance = TRUE)
iv.pred.stats.continuousW.model1b <- compute_stats(iv.pred.continuousW.model1b)
summary(iv.pred.continuousW.model1b)
summary(iv.pred.stats.continuousW.model1b)
write.csv(iv.pred.continuousW.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model1b.csv")
write.csv(iv.pred.stats.continuousW.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model1b.csv")

# save variable importance
iv.variable_importance.continuousW.model1b <- as.data.frame(variable_importance(iv.forest.continuousW.model1b))
iv.variable_importance.continuousW.model1b$Variable <- colnames(X.model1b.mat)
colnames(iv.variable_importance.continuousW.model1b) <- c("Importance", "Variable")
iv.variable_importance.continuousW.model1b <- iv.variable_importance.continuousW.model1b[, c("Variable", "Importance")]
write.csv(iv.variable_importance.continuousW.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_continuousW_variable_importance_model1b.csv")

# causal forest
message("\nContinuous W CF-grf model 1b")
causal.forest.model1b <- causal_forest(X.model1b.mat, Y.model1b.vector, W.model1b.vector, 
                                      num.threads = 10, num.trees = 5000, 
                                      sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
cf.pred.model1b <- predict(causal.forest.model1b, estimate.variance = TRUE)
cf.pred.stats.model1b <- compute_stats(cf.pred.model1b)
summary(cf.pred.model1b)
summary(cf.pred.stats.model1b)
write.csv(cf.pred.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model1b.csv")
write.csv(cf.pred.stats.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model1b.csv")

# # model 2
# # instrumental forest (continuous Z)
# message("\nContinuous W IV-grf model 2")
# iv.forest.continuousW.model2 <- instrumental_forest(X.model2.mat, Y.model2.vector, W.model2.vector, Z.model2.vector, 
#                                              num.threads = 10, num.trees = 5000, 
#                                              sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# iv.pred.continuousW.model2 <- predict(iv.forest.continuousW.model2, estimate.variance = TRUE)
# iv.pred.stats.continuousW.model2 <- compute_stats(iv.pred.continuousW.model2)
# summary(iv.pred.continuousW.model2)
# summary(iv.pred.stats.continuousW.model2)
# write.csv(iv.pred.continuousW.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model2.csv")
# write.csv(iv.pred.stats.continuousW.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model2.csv")

# # save variable importance
# iv.variable_importance.continuousW.model2 <- as.data.frame(variable_importance(iv.forest.continuousW.model2))
# iv.variable_importance.continuousW.model2$Variable <- colnames(X.model2.mat)
# colnames(iv.variable_importance.continuousW.model2) <- c("Importance", "Variable")
# iv.variable_importance.continuousW.model2 <- iv.variable_importance.continuousW.model2[, c("Variable", "Importance")]
# write.csv(iv.variable_importance.continuousW.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_continuousW_variable_importance_model2.csv")

# # causal forest
# message("\nContinuous W CF-grf model 2")
# causal.forest.model2 <- causal_forest(X.model2.mat, Y.model2.vector, W.model2.vector, 
#                                       num.threads = 10, num.trees = 5000, 
#                                       sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# cf.pred.model2 <- predict(causal.forest.model2, estimate.variance = TRUE)
# cf.pred.stats.model2 <- compute_stats(cf.pred.model2)
# summary(cf.pred.model2)
# summary(cf.pred.stats.model2)
# write.csv(cf.pred.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model2.csv")
# write.csv(cf.pred.stats.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model2.csv")

# model 2b
# instrumental forest (continuous Z)
message("\nContinuous W IV-grf model 2b")
iv.forest.continuousW.model2b <- instrumental_forest(X.model2b.mat, Y.model2b.vector, W.model2b.vector, Z.model2b.vector, 
                                             num.threads = 10, num.trees = 5000, 
                                             sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
iv.pred.continuousW.model2b <- predict(iv.forest.continuousW.model2b, estimate.variance = TRUE)
iv.pred.stats.continuousW.model2b <- compute_stats(iv.pred.continuousW.model2b)
summary(iv.pred.continuousW.model2b)
summary(iv.pred.stats.continuousW.model2b)
write.csv(iv.pred.continuousW.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model2b.csv")
write.csv(iv.pred.stats.continuousW.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model2b.csv")

# save variable importance
iv.variable_importance.continuousW.model2b <- as.data.frame(variable_importance(iv.forest.continuousW.model2b))
iv.variable_importance.continuousW.model2b$Variable <- colnames(X.model2b.mat)
colnames(iv.variable_importance.continuousW.model2b) <- c("Importance", "Variable")
iv.variable_importance.continuousW.model2b <- iv.variable_importance.continuousW.model2b[, c("Variable", "Importance")]
write.csv(iv.variable_importance.continuousW.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_continuousW_variable_importance_model2b.csv")

# causal forest
message("\nContinuous W CF-grf model 2b")
causal.forest.model2b <- causal_forest(X.model2b.mat, Y.model2b.vector, W.model2b.vector, 
                                      num.threads = 10, num.trees = 5000, 
                                      sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
cf.pred.model2b <- predict(causal.forest.model2b, estimate.variance = TRUE)
cf.pred.stats.model2b <- compute_stats(cf.pred.model2b)
summary(cf.pred.model2b)
summary(cf.pred.stats.model2b)
write.csv(cf.pred.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model2b.csv")
write.csv(cf.pred.stats.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model2b.csv")

# model 3
# instrumental forest (continuous Z)
message("\nContinuous W IV-grf model 3")
iv.forest.continuousW.model3 <- instrumental_forest(X.model3.mat, Y.model3.vector, W.model3.vector, Z.model3.vector, 
                                             num.threads = 10, num.trees = 5000, 
                                             sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
iv.pred.continuousW.model3 <- predict(iv.forest.continuousW.model3, estimate.variance = TRUE)
iv.pred.stats.continuousW.model3<- compute_stats(iv.pred.continuousW.model3)
summary(iv.pred.continuousW.model3)
summary(iv.pred.stats.continuousW.model3)
write.csv(iv.pred.continuousW.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model3.csv")
write.csv(iv.pred.stats.continuousW.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model3.csv")

# save variable importance
iv.variable_importance.continuousW.model3 <- as.data.frame(variable_importance(iv.forest.continuousW.model3))
iv.variable_importance.continuousW.model3$Variable <- colnames(X.model3.mat)
colnames(iv.variable_importance.continuousW.model3) <- c("Importance", "Variable")
iv.variable_importance.continuousW.model3 <- iv.variable_importance.continuousW.model3[, c("Variable", "Importance")]
write.csv(iv.variable_importance.continuousW.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_continuousW_variable_importance_model3.csv")

# causal forest
message("\nContinuous W CF-grf model 3")
causal.forest.model3 <- causal_forest(X.model3.mat, Y.model3.vector, W.model3.vector, 
                                      num.threads = 10, num.trees = 5000, 
                                      sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
cf.pred.model3 <- predict(causal.forest.model3, estimate.variance = TRUE)
cf.pred.stats.model3 <- compute_stats(cf.pred.model3)
summary(cf.pred.model3)
summary(cf.pred.stats.model3)
write.csv(cf.pred.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model3.csv")
write.csv(cf.pred.stats.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model3.csv")

# ===================================
####### binary W continuous Z #######
# =================================== 
# # model 1
# # instrumental forest
# message("\nBinary W continuous Z IV-grf model 1")
# iv.forest.binaryW.model1 <- instrumental_forest(X.model1.mat, Y.model1.vector, W.model1.vector.binary, Z.model1.vector, 
#                                          num.threads = 10, num.trees = 5000,
#                                          sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# iv.pred.binaryW.model1 <- predict(iv.forest.binaryW.model1, estimate.variance = TRUE)
# iv.pred.binaryW.stats.model1 <- compute_stats(iv.pred.binaryW.model1)
# summary(iv.pred.binaryW.model1)
# summary(iv.pred.binaryW.stats.model1)
# write.csv(iv.pred.binaryW.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model1.csv")
# write.csv(iv.pred.binaryW.stats.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model1.csv")

# # save variable importance
# iv.variable_importance.binaryW.model1 <- as.data.frame(variable_importance(iv.forest.binaryW.model1))
# iv.variable_importance.binaryW.model1$Variable <- colnames(X.model1.mat)
# colnames(iv.variable_importance.binaryW.model1) <- c("Importance", "Variable")
# iv.variable_importance.binaryW.model1 <- iv.variable_importance.binaryW.model1[, c("Variable", "Importance")]
# write.csv(iv.variable_importance.binaryW.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_binaryW_continuousZ_variable_importance_model1.csv")

# # causal forest
# message("\nBinary W continuous Z CF-grf model 1")
# causal.forest.binaryW.model1 <- causal_forest(X.model1.mat, Y.model1.vector, W.model1.vector.binary,
#                                        num.threads = 10, num.trees = 5000,
# sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# cf.pred.binaryW.model1 <- predict(causal.forest.binaryW.model1, estimate.variance = TRUE)
# cf.pred.binaryW.stats.model1 <- compute_stats(cf.pred.binaryW.model1)
# summary(cf.pred.binaryW.model1)
# summary(cf.pred.binaryW.stats.model1)
# write.csv(cf.pred.binaryW.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model1.csv")
# write.csv(cf.pred.binaryW.stats.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model1.csv")

# model 1b
# instrumental forest (continuous Z)
message("\nBinary W continuous Z IV-grf model 1b")
iv.forest.binaryW.model1b <- instrumental_forest(X.model1b.mat, Y.model1b.vector, W.model1b.vector.binary, Z.model1b.vector, 
                                         num.threads = 10, num.trees = 5000,
                                         sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
iv.pred.binaryW.model1b <- predict(iv.forest.binaryW.model1b, estimate.variance = TRUE)
iv.pred.binaryW.stats.model1b <- compute_stats(iv.pred.binaryW.model1b)
summary(iv.pred.binaryW.model1b)
summary(iv.pred.binaryW.stats.model1b)
write.csv(iv.pred.binaryW.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model1b.csv")
write.csv(iv.pred.binaryW.stats.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model1b.csv")

# save variable importance
iv.variable_importance.binaryW.model1b <- as.data.frame(variable_importance(iv.forest.binaryW.model1b))
iv.variable_importance.binaryW.model1b$Variable <- colnames(X.model1b.mat)
colnames(iv.variable_importance.binaryW.model1b) <- c("Importance", "Variable")
iv.variable_importance.binaryW.model1b <- iv.variable_importance.binaryW.model1b[, c("Variable", "Importance")]
write.csv(iv.variable_importance.binaryW.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_binaryW_continuousZ_variable_importance_model1b.csv")

# causal forest
message("\nBinary W continuous Z CF-grf model 1b")
causal.forest.binaryW.model1b <- causal_forest(X.model1b.mat, Y.model1b.vector, W.model1b.vector.binary,
                                       num.threads = 10, num.trees = 5000,
                                       sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
cf.pred.binaryW.model1b <- predict(causal.forest.binaryW.model1b, estimate.variance = TRUE)
cf.pred.binaryW.stats.model1b <- compute_stats(cf.pred.binaryW.model1b)
summary(cf.pred.binaryW.model1b)
summary(cf.pred.binaryW.stats.model1b)
write.csv(cf.pred.binaryW.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model1b.csv")
write.csv(cf.pred.binaryW.stats.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model1b.csv")

# # model 2
# # instrumental forest (continuous Z)
# message("\nBinary W continuous Z IV-grf model 2")
# iv.forest.binaryW.model2 <- instrumental_forest(X.model2.mat, Y.model2.vector, W.model2.vector.binary, Z.model2.vector, 
#                                          num.threads = 10, num.trees = 5000,
#                                          sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# iv.pred.binaryW.model2 <- predict(iv.forest.binaryW.model2, estimate.variance = TRUE)
# iv.pred.binaryW.stats.model2 <- compute_stats(iv.pred.binaryW.model2)
# summary(iv.pred.binaryW.model2)
# summary(iv.pred.binaryW.stats.model2)
# write.csv(iv.pred.binaryW.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model2.csv")
# write.csv(iv.pred.binaryW.stats.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model2.csv")

# # save variable importance
# iv.variable_importance.binaryW.model2 <- as.data.frame(variable_importance(iv.forest.binaryW.model2))
# iv.variable_importance.binaryW.model2$Variable <- colnames(X.model2.mat)
# colnames(iv.variable_importance.binaryW.model2) <- c("Importance", "Variable")
# iv.variable_importance.binaryW.model2 <- iv.variable_importance.binaryW.model2[, c("Variable", "Importance")]
# write.csv(iv.variable_importance.binaryW.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_binaryW_continuousZ_variable_importance_model2.csv")

# # causal forest
# message("\nBinary W continuous Z CF-grf model 2")
# causal.forest.binaryW.model2 <- causal_forest(X.model2.mat, Y.model2.vector, W.model2.vector.binary,
#                                        num.threads = 10, num.trees = 5000,
#                                        sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# cf.pred.binaryW.model2 <- predict(causal.forest.binaryW.model2, estimate.variance = TRUE)
# cf.pred.binaryW.stats.model2 <- compute_stats(cf.pred.binaryW.model2)
# summary(cf.pred.binaryW.model2)
# summary(cf.pred.binaryW.stats.model2)
# write.csv(cf.pred.binaryW.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model2.csv")
# write.csv(cf.pred.binaryW.stats.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model2.csv")

# model 2b
# instrumental forest (continuous Z)
message("\nBinary W continuous Z IV-grf model 2b")
iv.forest.binaryW.model2b <- instrumental_forest(X.model2b.mat, Y.model2b.vector, W.model2b.vector.binary, Z.model2b.vector, 
                                         num.threads = 10, num.trees = 5000,
                                         sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
iv.pred.binaryW.model2b <- predict(iv.forest.binaryW.model2b, estimate.variance = TRUE)
iv.pred.binaryW.stats.model2b <- compute_stats(iv.pred.binaryW.model2b)
summary(iv.pred.binaryW.model2b)
summary(iv.pred.binaryW.stats.model2b)
write.csv(iv.pred.binaryW.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model2b.csv")
write.csv(iv.pred.binaryW.stats.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model2b.csv")

# save variable importance
iv.variable_importance.binaryW.model2b <- as.data.frame(variable_importance(iv.forest.binaryW.model2b))
iv.variable_importance.binaryW.model2b$Variable <- colnames(X.model2b.mat)
colnames(iv.variable_importance.binaryW.model2b) <- c("Importance", "Variable")
iv.variable_importance.binaryW.model2b <- iv.variable_importance.binaryW.model2b[, c("Variable", "Importance")]
write.csv(iv.variable_importance.binaryW.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_binaryW_continuousZ_variable_importance_model2b.csv")

# causal forest
message("\nBinary W continuous Z CF-grf model 2b")
causal.forest.binaryW.model2b <- causal_forest(X.model2b.mat, Y.model2b.vector, W.model2b.vector.binary,
                                       num.threads = 10, num.trees = 5000,
                                       sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
cf.pred.binaryW.model2b <- predict(causal.forest.binaryW.model2b, estimate.variance = TRUE)
cf.pred.binaryW.stats.model2b <- compute_stats(cf.pred.binaryW.model2b)
summary(cf.pred.binaryW.model2b)
summary(cf.pred.binaryW.stats.model2b)
write.csv(cf.pred.binaryW.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model2b.csv")
write.csv(cf.pred.binaryW.stats.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model2b.csv")

# model 3
# instrumental forest (continuous Z)
message("\nBinary W continuous Z IV-grf model 3")
iv.forest.binaryW.model3 <- instrumental_forest(X.model3.mat, Y.model3.vector, W.model3.vector.binary, Z.model3.vector, 
                                                num.threads = 10, num.trees = 5000,
                                                sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
iv.pred.binaryW.model3 <- predict(iv.forest.binaryW.model3, estimate.variance = TRUE)
iv.pred.binaryW.stats.model3 <- compute_stats(iv.pred.binaryW.model3)
summary(iv.pred.binaryW.model3)
summary(iv.pred.binaryW.stats.model3)
write.csv(iv.pred.binaryW.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model3.csv")
write.csv(iv.pred.binaryW.stats.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model3.csv")

# save variable importance
iv.variable_importance.binaryW.model3 <- as.data.frame(variable_importance(iv.forest.binaryW.model3))
iv.variable_importance.binaryW.model3$Variable <- colnames(X.model3.mat)
colnames(iv.variable_importance.binaryW.model3) <- c("Importance", "Variable")
iv.variable_importance.binaryW.model3 <- iv.variable_importance.binaryW.model3[, c("Variable", "Importance")]
write.csv(iv.variable_importance.binaryW.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_binaryW_continuousZ_variable_importance_model3.csv")

# causal forest
message("\nBinary W continuous Z CF-grf model 3")
causal.forest.binaryW.model3 <- causal_forest(X.model3.mat, Y.model3.vector, W.model3.vector.binary,
                                              num.threads = 10, num.trees = 5000,
                                              sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
cf.pred.binaryW.model3 <- predict(causal.forest.binaryW.model3, estimate.variance = TRUE)
cf.pred.binaryW.stats.model3 <- compute_stats(cf.pred.binaryW.model3)
summary(cf.pred.binaryW.model3)
summary(cf.pred.binaryW.stats.model3)
write.csv(cf.pred.binaryW.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model3.csv")
write.csv(cf.pred.binaryW.stats.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model3.csv")

# # ===================================
# ####### binary W binary Z #######
# # =================================== 
# # model 1
# # instrumental forest
# message("\nBinary W binary Z IV-grf model 1")
# iv.forest.binaryW.model1 <- instrumental_forest(X.model1.mat, Y.model1.vector, W.model1.vector.binary, Z.model1.vector.binary, 
#                                          num.threads = 10, num.trees = 5000,
#                                          sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# iv.pred.binaryW.model1 <- predict(iv.forest.binaryW.model1, estimate.variance = TRUE)
# iv.pred.binaryW.stats.model1 <- compute_stats(iv.pred.binaryW.model1)
# summary(iv.pred.binaryW.model1)
# summary(iv.pred.binaryW.stats.model1)
# write.csv(iv.pred.binaryW.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_binaryZ_var_model1.csv")
# write.csv(iv.pred.binaryW.stats.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_binaryZ_pvalue_model1.csv")

# # save variable importance
# iv.variable_importance.binaryW.model1 <- as.data.frame(variable_importance(iv.forest.binaryW.model1))
# iv.variable_importance.binaryW.model1$Variable <- colnames(X.model1.mat)
# colnames(iv.variable_importance.binaryW.model1) <- c("Importance", "Variable")
# iv.variable_importance.binaryW.model1 <- iv.variable_importance.binaryW.model1[, c("Variable", "Importance")]
# write.csv(iv.variable_importance.binaryW.model1, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_binaryW_binaryZ_variable_importance_model1.csv")

# # model 1b
# # instrumental forest (binary Z)
# message("\nBinary W binary Z IV-grf model 1b")
# iv.forest.binaryW.model1b <- instrumental_forest(X.model1b.mat, Y.model1b.vector, W.model1b.vector.binary, Z.model1b.vector.binary, 
#                                          num.threads = 10, num.trees = 5000,
#                                          sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# iv.pred.binaryW.model1b <- predict(iv.forest.binaryW.model1b, estimate.variance = TRUE)
# iv.pred.binaryW.stats.model1b <- compute_stats(iv.pred.binaryW.model1b)
# summary(iv.pred.binaryW.model1b)
# summary(iv.pred.binaryW.stats.model1b)
# write.csv(iv.pred.binaryW.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_binaryZ_var_model1b.csv")
# write.csv(iv.pred.binaryW.stats.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_binaryZ_pvalue_model1b.csv")

# # save variable importance
# iv.variable_importance.binaryW.model1b <- as.data.frame(variable_importance(iv.forest.binaryW.model1b))
# iv.variable_importance.binaryW.model1b$Variable <- colnames(X.model1b.mat)
# colnames(iv.variable_importance.binaryW.model1b) <- c("Importance", "Variable")
# iv.variable_importance.binaryW.model1b <- iv.variable_importance.binaryW.model1b[, c("Variable", "Importance")]
# write.csv(iv.variable_importance.binaryW.model1b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_binaryW_binaryZ_variable_importance_model1b.csv")

# # model 2
# # instrumental forest (binary Z)
# message("\nBinary W binary Z IV-grf model 2")
# iv.forest.binaryW.model2 <- instrumental_forest(X.model2.mat, Y.model2.vector, W.model2.vector.binary, Z.model2.vector.binary, 
#                                          num.threads = 10, num.trees = 5000,
#                                          sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# iv.pred.binaryW.model2 <- predict(iv.forest.binaryW.model2, estimate.variance = TRUE)
# iv.pred.binaryW.stats.model2 <- compute_stats(iv.pred.binaryW.model2)
# summary(iv.pred.binaryW.model2)
# summary(iv.pred.binaryW.stats.model2)
# write.csv(iv.pred.binaryW.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_binaryZ_var_model2.csv")
# write.csv(iv.pred.binaryW.stats.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_binaryZ_pvalue_model2.csv")

# # save variable importance
# iv.variable_importance.binaryW.model2 <- as.data.frame(variable_importance(iv.forest.binaryW.model2))
# iv.variable_importance.binaryW.model2$Variable <- colnames(X.model2.mat)
# colnames(iv.variable_importance.binaryW.model2) <- c("Importance", "Variable")
# iv.variable_importance.binaryW.model2 <- iv.variable_importance.binaryW.model2[, c("Variable", "Importance")]
# write.csv(iv.variable_importance.binaryW.model2, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_binaryW_binaryZ_variable_importance_model2.csv")

# # model 2b
# # instrumental forest (binary Z)
# message("\nBinary W binary Z IV-grf model 2b")
# iv.forest.binaryW.model2b <- instrumental_forest(X.model2b.mat, Y.model2b.vector, W.model2b.vector.binary, Z.model2b.vector.binary, 
#                                          num.threads = 10, num.trees = 5000,
#                                          sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# iv.pred.binaryW.model2b <- predict(iv.forest.binaryW.model2b, estimate.variance = TRUE)
# iv.pred.binaryW.stats.model2b <- compute_stats(iv.pred.binaryW.model2b)
# summary(iv.pred.binaryW.model2b)
# summary(iv.pred.binaryW.stats.model2b)
# write.csv(iv.pred.binaryW.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_binaryZ_var_model2b.csv")
# write.csv(iv.pred.binaryW.stats.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_binaryZ_pvalue_model2b.csv")

# # save variable importance
# iv.variable_importance.binaryW.model2b <- as.data.frame(variable_importance(iv.forest.binaryW.model2b))
# iv.variable_importance.binaryW.model2b$Variable <- colnames(X.model2b.mat)
# colnames(iv.variable_importance.binaryW.model2b) <- c("Importance", "Variable")
# iv.variable_importance.binaryW.model2b <- iv.variable_importance.binaryW.model2b[, c("Variable", "Importance")]
# write.csv(iv.variable_importance.binaryW.model2b, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_binaryW_binaryZ_variable_importance_model2b.csv")

# # model 3
# # instrumental forest (binary Z)
# message("\nBinary W binary Z IV-grf model 3")
# iv.forest.binaryW.model3 <- instrumental_forest(X.model3.mat, Y.model3.vector, W.model3.vector.binary, Z.model3.vector.binary, 
#                                          num.threads = 10, num.trees = 5000,
#                                          sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
# iv.pred.binaryW.model3 <- predict(iv.forest.binaryW.model3, estimate.variance = TRUE)
# iv.pred.binaryW.stats.model3 <- compute_stats(iv.pred.binaryW.model3)
# summary(iv.pred.binaryW.model3)
# summary(iv.pred.binaryW.stats.model3)
# write.csv(iv.pred.binaryW.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_binaryZ_var_model3.csv")
# write.csv(iv.pred.binaryW.stats.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_binaryZ_pvalue_model3.csv")

# # save variable importance
# iv.variable_importance.binaryW.model3 <- as.data.frame(variable_importance(iv.forest.binaryW.model3))
# iv.variable_importance.binaryW.model3$Variable <- colnames(X.model3.mat)
# colnames(iv.variable_importance.binaryW.model3) <- c("Importance", "Variable")
# iv.variable_importance.binaryW.model3 <- iv.variable_importance.binaryW.model3[, c("Variable", "Importance")]
# write.csv(iv.variable_importance.binaryW.model3, file = "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/ivgrf_binaryW_binaryZ_variable_importance_model3.csv")