#' This code is to run a pilot study in using instrumental random forest and causal forest for analyzing the 
#' causal relationship between HDL and CAD.
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
  #   library(xlsx)
  library(ivreg)
  library(ggsignif) 
  
  source("/mnt/md0/yujia/project/2023-07-20-individual_MR/bin/05_ivgrf_disease/support_func/psm.R")
  source("/mnt/md0/yujia/project/2023-07-20-individual_MR/bin/05_ivgrf_disease/support_func/test_heterogeneity.R")
  source("/mnt/md0/yujia/project/2023-07-20-individual_MR/bin/05_ivgrf_disease/support_func/find_best_tree.R")
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
X.set1 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/external_analysis_for_paper_revision/traits/ukbb.covariate.traits.set1.gz", sep = "\t")
X.set2 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/external_analysis_for_paper_revision/traits/ukbb.covariate.traits.set2.gz", sep = "\t")
X.set3 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/external_analysis_for_paper_revision/HDL/CAD/ukbb.covariate.HDL.set3.gz", sep = "\t")

Z.set1 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/external_analysis_for_paper_revision/PRSice2/HDL/CAD/set1/1e-08/HDL_prs.best", sep = " ") # score_constd
Z.set2 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/external_analysis_for_paper_revision/PRSice2/HDL/CAD/set2/1e-08/HDL_prs.best", sep = " ") # score_constd
Z.set3 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/external_analysis_for_paper_revision/PRSice2/HDL/CAD/set3/1e-08/HDL_prs.best", sep = " ") # score_constd
W <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/external_analysis_for_paper_revision/HDL/CAD/ukbb.phenotype.HDL.mgdL", sep = "\t")

Y.date1 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date1_james2022.gz", sep = "\t")
Y.date2 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date2_james2022.gz", sep = "\t")
Y.date3 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date3_james2022.gz", sep = "\t")

dim(X.set1); dim(X.set1); dim(X.set1); dim(X.set1); dim(Z.set1); dim(Z.set2); dim(Z.set3); dim(W); dim(Y.date1); dim(Y.date2); dim(Y.date3)

# homonize the ID included in the analysis
selected_id_set1b <- Reduce(intersect, list(X.set1$IID, Z.set1$IID, W$IID, Y.date3$IID))
selected_id_set3 <- Reduce(intersect, list(X.set3$IID, Z.set3$IID, W$IID, Y.date3$IID))

# model 1b
X.model1b <- X.set1[X.set1$IID %in% selected_id_set1b, ]
Y.model1b <- Y.date3[Y.date3$IID %in% selected_id_set1b, ]
W.model1b <- W[W$IID %in% selected_id_set1b, ]
Z.model1b <- Z.set1[Z.set1$IID %in% selected_id_set1b, ]
message("Patient Info for model 1b")
dim(X.model1b); dim(X.model1b); dim(X.model1b); dim(X.model1b);
print(table(Y.model1b$CAD))

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

# Get the ID of different ethnic background
white_ids <- X.model3[X.model3$ethnic_group == 1, ]$IID
asian_ids <- X.model3[X.model3$ethnic_group == 2, ]$IID
black_ids <- X.model3[X.model3$ethnic_group == 3, ]$IID
chinese_ids <- X.model3[X.model3$ethnic_group == 4, ]$IID

white_idx <- which(X.model3$IID %in% white_ids)
asian_idx <- which(X.model3$IID %in% asian_ids)
black_idx <- which(X.model3$IID %in% black_ids)
chinese_idx <- which(X.model3$IID %in% chinese_ids)

# generate mat file for three models
# model 1b
W.model1b.mat <- W.model1b[, c("30760-0.0")]
X.model1b.mat <- X.model1b[, 3:dim(X.model1b)[2]]
Y.model1b.mat <- Y.model1b[, c("CAD")]
Z.model1b.mat <- Z.model1b[, 4]
dim(W.model1b.mat); dim(X.model1b.mat); dim(Y.model1b.mat); dim(Z.model1b.mat)

# model 3
W.model3.mat <- W.model3[, c("30760-0.0")]
X.model3.mat <- X.model3[, 3:dim(X.model3)[2]]
Y.model3.mat <- Y.model3[, c("CAD")]
Z.model3.mat <- Z.model3[, 4]
dim(W.model3.mat); dim(X.model3.mat); dim(Y.model3.mat); dim(Z.model3.mat)

# rename the mat colnames
X.model1b.mat  <- X.model1b.mat  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "Genotype Batch"="22000-0.0",
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10")

X.model3.mat  <- X.model3.mat  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0",   
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0",
                         "BMI"="21001-0.0", "Body Fat Percentage"="23099-0.0", "Weight"="21002-0.0", # "Waist-hip-ratio"="whr", 
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "LDL-C"="30780-0.0", "Triglycerides"="30870-0.0", "Apolipoprotein B"="30640-0.0", "Lipoprotein A"="30790-0.0", # lipid-related covariates
                         "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", "Calcium"="30680-0.0",
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0", "Total protein"="30860-0.0", 
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")

# preparing vector for different model
# model 1b
W.model1b.vector <- as.vector(W.model1b.mat$`30760-0.0`)
Y.model1b.vector <- as.integer(as.vector(Y.model1b.mat$CAD))
Z.model1b.vector <- as.vector(Z.model1b.mat$PRS)
W.model1b.vector.binary <- as.integer(ifelse(W.model1b.vector <= 46, 1, 0)) # mimic usage of statin

# model 3
W.model3.vector <- as.vector(W.model3.mat$`30760-0.0`)
Y.model3.vector <- as.integer(as.vector(Y.model3.mat$CAD))
Z.model3.vector <- as.vector(Z.model3.mat$PRS)
W.model3.vector.binary <- as.integer(ifelse(W.model3.vector <= 46, 1, 0)) # mimic usage of statin

# Applying Wu-Hausman Test to exploring the validity of applying instrument in the analysis.
# model 1b
WHtest.mat.model1b <- cbind(X.model1b.mat, Y.model1b.mat, W.model1b.mat, Z.model1b.mat)
dim(WHtest.mat.model1b)
colnames(WHtest.mat.model1b)[dim(WHtest.mat.model1b)[2]-1] <- "W"
res.model1b <- ivreg(CAD ~ . - W - PRS | W | PRS, data = WHtest.mat.model1b)
summary(res.model1b, diagnostics = TRUE)

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
# write.csv(wald.test.res, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/wald_test.csv")

# ===============
# Formal Analysis
# ===============

dir.create("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect", showWarnings = F, recursive = T)
dir.create("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/03_variable_importance", showWarnings = F, recursive = T)

# ========================
####### continuous ####### 
# ======================== 
################################################ model 1b ################################################

##########################################
### instrumental forest (continuous Z) ###
##########################################
message("\nContinuous W IV-grf model 1b")
iv.forest.continuousW.model1b <- instrumental_forest(X.model1b.mat[white_idx, ], Y.model1b.vector[white_idx], W.model1b.vector[white_idx], Z.model1b.vector[white_idx], 
                                                     num.threads = 40, num.trees = 5000, 
                                                     sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100) # Train on white ancestry

# save variable importance
iv.variable_importance.continuousW.model1b <- as.data.frame(variable_importance(iv.forest.continuousW.model1b))
iv.variable_importance.continuousW.model1b$Variable <- colnames(X.model1b.mat)
colnames(iv.variable_importance.continuousW.model1b) <- c("Importance", "Variable")
iv.variable_importance.continuousW.model1b <- iv.variable_importance.continuousW.model1b[, c("Variable", "Importance")]
write.csv(iv.variable_importance.continuousW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/03_variable_importance/ivgrf_continuousW_variable_importance_model1b.csv")

# White Ancestry Results
message("\nWhite Ancestry")
iv.pred.continuousW.model1b <- predict(iv.forest.continuousW.model1b, estimate.variance = TRUE, num.threads = 40)
iv.pred.stats.continuousW.model1b <- compute_stats(iv.pred.continuousW.model1b)
iv.pred.continuousW.model1b$IID <- white_ids
iv.pred.stats.continuousW.model1b$IID <- white_ids
summary(iv.pred.continuousW.model1b)
summary(iv.pred.stats.continuousW.model1b)
write.csv(iv.pred.continuousW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model1b_white.csv")
write.csv(iv.pred.stats.continuousW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model1b_white.csv")

# Black Ancestry Results
message("\nBlack Ancestry")
iv.pred.continuousW.model1b <- predict(iv.forest.continuousW.model1b, newdata = X.model1b.mat[black_idx, ], estimate.variance = TRUE, num.threads = 40)
iv.pred.stats.continuousW.model1b <- compute_stats(iv.pred.continuousW.model1b)
iv.pred.continuousW.model1b$IID <- black_ids
iv.pred.stats.continuousW.model1b$IID <- black_ids
summary(iv.pred.continuousW.model1b)
summary(iv.pred.stats.continuousW.model1b)
write.csv(iv.pred.continuousW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model1b_black.csv")
write.csv(iv.pred.stats.continuousW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model1b_black.csv")

# Asian Ancestry Results
message("\nAsian Ancestry")
iv.pred.continuousW.model1b <- predict(iv.forest.continuousW.model1b, newdata = X.model1b.mat[asian_idx, ], estimate.variance = TRUE, num.threads = 40)
iv.pred.stats.continuousW.model1b <- compute_stats(iv.pred.continuousW.model1b)
iv.pred.continuousW.model1b$IID <- asian_ids
iv.pred.stats.continuousW.model1b$IID <- asian_ids
summary(iv.pred.continuousW.model1b)
summary(iv.pred.stats.continuousW.model1b)
write.csv(iv.pred.continuousW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model1b_asian.csv")
write.csv(iv.pred.stats.continuousW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model1b_asian.csv")

# Chinese Ancestry Results
message("\nChinese Ancestry")
iv.pred.continuousW.model1b <- predict(iv.forest.continuousW.model1b, newdata = X.model1b.mat[chinese_idx, ], estimate.variance = TRUE, num.threads = 40)
iv.pred.stats.continuousW.model1b <- compute_stats(iv.pred.continuousW.model1b)
iv.pred.continuousW.model1b$IID <- chinese_ids
iv.pred.stats.continuousW.model1b$IID <- chinese_ids
summary(iv.pred.continuousW.model1b)
summary(iv.pred.stats.continuousW.model1b)
write.csv(iv.pred.continuousW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model1b_chinese.csv")
write.csv(iv.pred.stats.continuousW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model1b_chinese.csv")

##########################################
############# causal forest ##############
##########################################
message("\nContinuous W CF-grf model 1b")
causal.forest.model1b <- causal_forest(X.model1b.mat[white_idx, ], Y.model1b.vector[white_idx], W.model1b.vector[white_idx], 
                                       num.threads = 40, num.trees = 5000, 
                                       sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)

# White Ancestry Results
message("\nWhite Ancestry")
cf.pred.model1b <- predict(causal.forest.model1b, estimate.variance = TRUE, num.threads = 40)
cf.pred.stats.model1b <- compute_stats(cf.pred.model1b)
cf.pred.model1b$IID <- white_ids
cf.pred.stats.model1b$IID <- white_ids
summary(cf.pred.model1b)
summary(cf.pred.stats.model1b)
write.csv(cf.pred.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model1b_white.csv")
write.csv(cf.pred.stats.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model1b_white.csv")

# Black Ancestry Results
message("\nBlack Ancestry")
cf.pred.model1b <- predict(causal.forest.model1b, newdata = X.model1b.mat[black_idx, ], estimate.variance = TRUE, num.threads = 40)
cf.pred.stats.model1b <- compute_stats(cf.pred.model1b)
cf.pred.model1b$IID <- black_ids
cf.pred.stats.model1b$IID <- black_ids
summary(cf.pred.model1b)
summary(cf.pred.stats.model1b)
write.csv(cf.pred.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model1b_black.csv")
write.csv(cf.pred.stats.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model1b_black.csv")

# Asian Ancestry Results
message("\nAsian Ancestry")
cf.pred.model1b <- predict(causal.forest.model1b, newdata = X.model1b.mat[asian_idx, ], estimate.variance = TRUE, num.threads = 40)
cf.pred.stats.model1b <- compute_stats(cf.pred.model1b)
cf.pred.model1b$IID <- asian_ids
cf.pred.stats.model1b$IID <- asian_ids
summary(cf.pred.model1b)
summary(cf.pred.stats.model1b)
write.csv(cf.pred.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model1b_asian.csv")
write.csv(cf.pred.stats.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model1b_asian.csv")

# Chinese Ancestry Results
message("\nChinese Ancestry")
cf.pred.model1b <- predict(causal.forest.model1b, newdata = X.model1b.mat[chinese_idx, ], estimate.variance = TRUE, num.threads = 40)
cf.pred.stats.model1b <- compute_stats(cf.pred.model1b)
cf.pred.model1b$IID <- chinese_ids
cf.pred.stats.model1b$IID <- chinese_ids
summary(cf.pred.model1b)
summary(cf.pred.stats.model1b)
write.csv(cf.pred.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model1b_chinese.csv")
write.csv(cf.pred.stats.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model1b_chinese.csv")

################################################ model 3 ################################################
##########################################
### instrumental forest (continuous Z) ###
##########################################
message("\nContinuous W IV-grf model 3")
iv.forest.continuousW.model3 <- instrumental_forest(X.model3.mat[white_idx, ], Y.model3.vector[white_idx], W.model3.vector[white_idx], Z.model3.vector[white_idx], 
                                                    num.threads = 40, num.trees = 5000, 
                                                    sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)

# save variable importance
iv.variable_importance.continuousW.model3 <- as.data.frame(variable_importance(iv.forest.continuousW.model3))
iv.variable_importance.continuousW.model3$Variable <- colnames(X.model3.mat)
colnames(iv.variable_importance.continuousW.model3) <- c("Importance", "Variable")
iv.variable_importance.continuousW.model3 <- iv.variable_importance.continuousW.model3[, c("Variable", "Importance")]
write.csv(iv.variable_importance.continuousW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/03_variable_importance/ivgrf_continuousW_variable_importance_model3.csv")

# White Ancestry Results
message("\nWhite Ancestry")
iv.pred.continuousW.model3 <- predict(iv.forest.continuousW.model3, estimate.variance = TRUE, num.threads = 40)
iv.pred.stats.continuousW.model3<- compute_stats(iv.pred.continuousW.model3)
iv.pred.continuousW.model3$IID <- white_ids
iv.pred.stats.continuousW.model3$IID <- white_ids
summary(iv.pred.continuousW.model3)
summary(iv.pred.stats.continuousW.model3)
write.csv(iv.pred.continuousW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model3_white.csv")
write.csv(iv.pred.stats.continuousW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model3_white.csv")

# Black Ancestry Results
message("\nBlack Ancestry")
iv.pred.continuousW.model3 <- predict(iv.forest.continuousW.model3, newdata = X.model3.mat[black_idx, ], estimate.variance = TRUE, num.threads = 40)
iv.pred.stats.continuousW.model3<- compute_stats(iv.pred.continuousW.model3)
iv.pred.continuousW.model3$IID <- black_ids
iv.pred.stats.continuousW.model3$IID <- black_ids
summary(iv.pred.continuousW.model3)
summary(iv.pred.stats.continuousW.model3)
write.csv(iv.pred.continuousW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model3_black.csv")
write.csv(iv.pred.stats.continuousW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model3_black.csv")

# Asian Ancestry Results
message("\nAsian Ancestry")
iv.pred.continuousW.model3 <- predict(iv.forest.continuousW.model3, newdata = X.model3.mat[asian_idx, ], estimate.variance = TRUE, num.threads = 40)
iv.pred.stats.continuousW.model3<- compute_stats(iv.pred.continuousW.model3)
iv.pred.continuousW.model3$IID <- asian_ids
iv.pred.stats.continuousW.model3$IID <- asian_ids
summary(iv.pred.continuousW.model3)
summary(iv.pred.stats.continuousW.model3)
write.csv(iv.pred.continuousW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model3_asian.csv")
write.csv(iv.pred.stats.continuousW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model3_asian.csv")

# Chinese Ancestry Results
message("\nChinese Ancestry")
iv.pred.continuousW.model3 <- predict(iv.forest.continuousW.model3, newdata = X.model3.mat[chinese_idx, ], estimate.variance = TRUE, num.threads = 40)
iv.pred.stats.continuousW.model3<- compute_stats(iv.pred.continuousW.model3)
iv.pred.continuousW.model3$IID <- chinese_ids
iv.pred.stats.continuousW.model3$IID <- chinese_ids
summary(iv.pred.continuousW.model3)
summary(iv.pred.stats.continuousW.model3)
write.csv(iv.pred.continuousW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_var_model3_chinese.csv")
write.csv(iv.pred.stats.continuousW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model3_chinese.csv")

##########################################
############# causal forest ##############
##########################################
message("\nContinuous W CF-grf model 3")
causal.forest.model3 <- causal_forest(X.model3.mat[white_idx, ], Y.model3.vector[white_idx], W.model3.vector[white_idx], 
                                      num.threads = 40, num.trees = 5000, 
                                      sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)

# White Ancestry Results
message("\nWhite Ancestry")
cf.pred.model3 <- predict(causal.forest.model3, estimate.variance = TRUE, num.threads = 40)
cf.pred.stats.model3 <- compute_stats(cf.pred.model3)
cf.pred.model3$IID <- white_ids
cf.pred.stats.model3$IID <- white_ids
summary(cf.pred.model3)
summary(cf.pred.stats.model3)
write.csv(cf.pred.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model3_white.csv")
write.csv(cf.pred.stats.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model3_white.csv")

# Black Ancestry Results
message("\nBlack Ancestry")
cf.pred.model3 <- predict(causal.forest.model3, newdata = X.model3.mat[black_idx, ], estimate.variance = TRUE, num.threads = 40)
cf.pred.stats.model3 <- compute_stats(cf.pred.model3)
cf.pred.model3$IID <- black_ids
cf.pred.stats.model3$IID <- black_ids
summary(cf.pred.model3)
summary(cf.pred.stats.model3)
write.csv(cf.pred.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model3_black.csv")
write.csv(cf.pred.stats.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model3_black.csv")

# Asian Ancestry Results
message("\nAsian Ancestry")
cf.pred.model3 <- predict(causal.forest.model3, newdata = X.model3.mat[asian_idx, ], estimate.variance = TRUE, num.threads = 40)
cf.pred.stats.model3 <- compute_stats(cf.pred.model3)
cf.pred.model3$IID <- asian_ids
cf.pred.stats.model3$IID <- asian_ids
summary(cf.pred.model3)
summary(cf.pred.stats.model3)
write.csv(cf.pred.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model3_asian.csv")
write.csv(cf.pred.stats.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model3_asian.csv")

# Chinese Ancestry Results
message("\nChinese Ancestry")
cf.pred.model3 <- predict(causal.forest.model3, newdata = X.model3.mat[chinese_idx, ], estimate.variance = TRUE, num.threads = 40)
cf.pred.stats.model3 <- compute_stats(cf.pred.model3)
cf.pred.model3$IID <- chinese_ids
cf.pred.stats.model3$IID <- chinese_ids
summary(cf.pred.model3)
summary(cf.pred.stats.model3)
write.csv(cf.pred.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_var_model3_chinese.csv")
write.csv(cf.pred.stats.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_continuousW_pvalue_model3_chinese.csv")

# ===================================
####### binary W continuous Z #######
# =================================== 
################################################ model 1b ################################################

##########################################
### instrumental forest (continuous Z) ###
##########################################
message("\nBinary W continuous Z IV-grf model 1b")
iv.forest.binaryW.model1b <- instrumental_forest(X.model1b.mat[white_idx, ], Y.model1b.vector[white_idx], W.model1b.vector.binary[white_idx], Z.model1b.vector[white_idx], 
                                                 num.threads = 40, num.trees = 5000,
                                                 sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100) # Train on white ancestry

# save variable importance
iv.variable_importance.binaryW.model1b <- as.data.frame(variable_importance(iv.forest.binaryW.model1b))
iv.variable_importance.binaryW.model1b$Variable <- colnames(X.model1b.mat)
colnames(iv.variable_importance.binaryW.model1b) <- c("Importance", "Variable")
iv.variable_importance.binaryW.model1b <- iv.variable_importance.binaryW.model1b[, c("Variable", "Importance")]
write.csv(iv.variable_importance.binaryW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/03_variable_importance/ivgrf_binaryW_continuousZ_variable_importance_model1b.csv")

# White Ancestry Results
message("\nWhite Ancestry")
iv.pred.binaryW.model1b <- predict(iv.forest.binaryW.model1b, estimate.variance = TRUE)
iv.pred.binaryW.stats.model1b <- compute_stats(iv.pred.binaryW.model1b)
iv.pred.binaryW.model1b$IID <- white_ids
iv.pred.binaryW.stats.model1b$IID <- white_ids
summary(iv.pred.binaryW.model1b)
summary(iv.pred.binaryW.stats.model1b)
write.csv(iv.pred.binaryW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model1b_white.csv")
write.csv(iv.pred.binaryW.stats.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model1b_white.csv")

# Black Ancestry Results
message("\nBlack Ancestry")
iv.pred.binaryW.model1b <- predict(iv.forest.binaryW.model1b, newdata = X.model1b.mat[black_idx, ],  estimate.variance = TRUE)
iv.pred.binaryW.stats.model1b <- compute_stats(iv.pred.binaryW.model1b)
iv.pred.binaryW.model1b$IID <- black_ids
iv.pred.binaryW.stats.model1b$IID <- black_ids
summary(iv.pred.binaryW.model1b)
summary(iv.pred.binaryW.stats.model1b)
write.csv(iv.pred.binaryW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model1b_black.csv")
write.csv(iv.pred.binaryW.stats.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model1b_black.csv")

# Asian Ancestry Results
message("\nAsian Ancestry")
iv.pred.binaryW.model1b <- predict(iv.forest.binaryW.model1b, newdata = X.model1b.mat[asian_idx, ],  estimate.variance = TRUE)
iv.pred.binaryW.stats.model1b <- compute_stats(iv.pred.binaryW.model1b)
iv.pred.binaryW.model1b$IID <- asian_ids
iv.pred.binaryW.stats.model1b$IID <- asian_ids
summary(iv.pred.binaryW.model1b)
summary(iv.pred.binaryW.stats.model1b)
write.csv(iv.pred.binaryW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model1b_asian.csv")
write.csv(iv.pred.binaryW.stats.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model1b_asian.csv")

# Chinese Ancestry Results
message("\nChinese Ancestry")
iv.pred.binaryW.model1b <- predict(iv.forest.binaryW.model1b, newdata = X.model1b.mat[chinese_idx, ],  estimate.variance = TRUE)
iv.pred.binaryW.stats.model1b <- compute_stats(iv.pred.binaryW.model1b)
iv.pred.binaryW.model1b$IID <- chinese_ids
iv.pred.binaryW.stats.model1b$IID <- chinese_ids
summary(iv.pred.binaryW.model1b)
summary(iv.pred.binaryW.stats.model1b)
write.csv(iv.pred.binaryW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model1b_chinese.csv")
write.csv(iv.pred.binaryW.stats.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model1b_chinese.csv")

##########################################
############# causal forest ##############
##########################################
message("\nBinary W continuous Z CF-grf model 1b")
causal.forest.binaryW.model1b <- causal_forest(X.model1b.mat[white_idx, ], Y.model1b.vector[white_idx], W.model1b.vector.binary[white_idx],
                                               num.threads = 10, num.trees = 5000,
                                               sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)

# White Ancestry Results                                        
message("\nWhite Ancestry")                                
cf.pred.binaryW.model1b <- predict(causal.forest.binaryW.model1b, estimate.variance = TRUE)
cf.pred.binaryW.stats.model1b <- compute_stats(cf.pred.binaryW.model1b)
cf.pred.binaryW.model1b$IID <- white_ids
cf.pred.binaryW.stats.model1b$IID <- white_ids
summary(cf.pred.binaryW.model1b)
summary(cf.pred.binaryW.stats.model1b)
write.csv(cf.pred.binaryW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model1b_white.csv")
write.csv(cf.pred.binaryW.stats.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model1b_white.csv")

# Black Ancestry Results                                        
message("\nBlack Ancestry")                                        
cf.pred.binaryW.model1b <- predict(causal.forest.binaryW.model1b, newdata = X.model1b.mat[black_idx, ], estimate.variance = TRUE)
cf.pred.binaryW.stats.model1b <- compute_stats(cf.pred.binaryW.model1b)
cf.pred.binaryW.model1b$IID <- black_ids
cf.pred.binaryW.stats.model1b$IID <- black_ids
summary(cf.pred.binaryW.model1b)
summary(cf.pred.binaryW.stats.model1b)
write.csv(cf.pred.binaryW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model1b_black.csv")
write.csv(cf.pred.binaryW.stats.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model1b_black.csv")

# Asian Ancestry Results                                        
message("\nAsian Ancestry")                                    
cf.pred.binaryW.model1b <- predict(causal.forest.binaryW.model1b, newdata = X.model1b.mat[asian_idx, ], estimate.variance = TRUE)
cf.pred.binaryW.stats.model1b <- compute_stats(cf.pred.binaryW.model1b)
cf.pred.binaryW.model1b$IID <- asian_ids
cf.pred.binaryW.stats.model1b$IID <- asian_ids
summary(cf.pred.binaryW.model1b)
summary(cf.pred.binaryW.stats.model1b)
write.csv(cf.pred.binaryW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model1b_asian.csv")
write.csv(cf.pred.binaryW.stats.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model1b_asian.csv")

# Chinese Ancestry Results                                        
message("\nChinese Ancestry")                                     
cf.pred.binaryW.model1b <- predict(causal.forest.binaryW.model1b, newdata = X.model1b.mat[chinese_idx, ], estimate.variance = TRUE)
cf.pred.binaryW.stats.model1b <- compute_stats(cf.pred.binaryW.model1b)
cf.pred.binaryW.model1b$IID <- chinese_ids
cf.pred.binaryW.stats.model1b$IID <- chinese_ids
summary(cf.pred.binaryW.model1b)
summary(cf.pred.binaryW.stats.model1b)
write.csv(cf.pred.binaryW.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model1b_chinese.csv")
write.csv(cf.pred.binaryW.stats.model1b, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model1b_chinese.csv")

################################################ model 3 ################################################
##########################################
### instrumental forest (continuous Z) ###
##########################################
message("\nBinary W continuous Z IV-grf model 3")
iv.forest.binaryW.model3 <- instrumental_forest(X.model3.mat[white_idx, ], Y.model3.vector[white_idx], W.model3.vector.binary[white_idx], Z.model3.vector[white_idx], 
                                                num.threads = 40, num.trees = 5000,
                                                sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100) # Train on white ancestry

# save variable importance
iv.variable_importance.binaryW.model3 <- as.data.frame(variable_importance(iv.forest.binaryW.model3))
iv.variable_importance.binaryW.model3$Variable <- colnames(X.model3.mat)
colnames(iv.variable_importance.binaryW.model3) <- c("Importance", "Variable")
iv.variable_importance.binaryW.model3 <- iv.variable_importance.binaryW.model3[, c("Variable", "Importance")]
write.csv(iv.variable_importance.binaryW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/03_variable_importance/ivgrf_binaryW_continuousZ_variable_importance_model3.csv")

# White Ancestry Results
message("\nWhite Ancestry")
iv.pred.binaryW.model3 <- predict(iv.forest.binaryW.model3, estimate.variance = TRUE)
iv.pred.binaryW.stats.model3 <- compute_stats(iv.pred.binaryW.model3)
iv.pred.binaryW.model3$IID <- white_ids
iv.pred.binaryW.stats.model3$IID <- white_ids
summary(iv.pred.binaryW.model3)
summary(iv.pred.binaryW.stats.model3)
write.csv(iv.pred.binaryW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model3_white.csv")
write.csv(iv.pred.binaryW.stats.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model3_white.csv")

# Black Ancestry Results
message("\nBlack Ancestry")
iv.pred.binaryW.model3 <- predict(iv.forest.binaryW.model3, newdata = X.model3.mat[black_idx, ],  estimate.variance = TRUE)
iv.pred.binaryW.stats.model3 <- compute_stats(iv.pred.binaryW.model3)
iv.pred.binaryW.model3$IID <- black_ids
iv.pred.binaryW.stats.model3$IID <- black_ids
summary(iv.pred.binaryW.model3)
summary(iv.pred.binaryW.stats.model3)
write.csv(iv.pred.binaryW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model3_black.csv")
write.csv(iv.pred.binaryW.stats.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model3_black.csv")

# Asian Ancestry Results
message("\nAsian Ancestry")
iv.pred.binaryW.model3 <- predict(iv.forest.binaryW.model3, newdata = X.model3.mat[asian_idx, ],  estimate.variance = TRUE)
iv.pred.binaryW.stats.model3 <- compute_stats(iv.pred.binaryW.model3)
iv.pred.binaryW.model3$IID <- asian_ids
iv.pred.binaryW.stats.model3$IID <- asian_ids
summary(iv.pred.binaryW.model3)
summary(iv.pred.binaryW.stats.model3)
write.csv(iv.pred.binaryW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model3_asian.csv")
write.csv(iv.pred.binaryW.stats.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model3_asian.csv")

# Chinese Ancestry Results
message("\nChinese Ancestry")
iv.pred.binaryW.model3 <- predict(iv.forest.binaryW.model3, newdata = X.model3.mat[chinese_idx, ],  estimate.variance = TRUE)
iv.pred.binaryW.stats.model3 <- compute_stats(iv.pred.binaryW.model3)
iv.pred.binaryW.model3$IID <- chinese_ids
iv.pred.binaryW.stats.model3$IID <- chinese_ids
summary(iv.pred.binaryW.model3)
summary(iv.pred.binaryW.stats.model3)
write.csv(iv.pred.binaryW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_var_model3_chinese.csv")
write.csv(iv.pred.binaryW.stats.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_model3_chinese.csv")

##########################################
############# causal forest ##############
##########################################
message("\nBinary W continuous Z CF-grf model 3")
causal.forest.binaryW.model3 <- causal_forest(X.model3.mat[white_idx, ], Y.model3.vector[white_idx], W.model3.vector.binary[white_idx],
                                              num.threads = 10, num.trees = 5000,
                                              sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)

# White Ancestry Results                                        
message("\nWhite Ancestry")                                      
cf.pred.binaryW.model3 <- predict(causal.forest.binaryW.model3, estimate.variance = TRUE)
cf.pred.binaryW.stats.model3 <- compute_stats(cf.pred.binaryW.model3)
cf.pred.binaryW.model3$IID <- white_ids
cf.pred.binaryW.stats.model3$IID <- white_ids
summary(cf.pred.binaryW.model3)
summary(cf.pred.binaryW.stats.model3)
write.csv(cf.pred.binaryW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model3_white.csv")
write.csv(cf.pred.binaryW.stats.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model3_white.csv")

# Black Ancestry Results                                        
message("\nBlack Ancestry")                                     
cf.pred.binaryW.model3 <- predict(causal.forest.binaryW.model3, newdata = X.model3.mat[black_idx, ], estimate.variance = TRUE)
cf.pred.binaryW.stats.model3 <- compute_stats(cf.pred.binaryW.model3)
cf.pred.binaryW.model3$IID <- black_ids
cf.pred.binaryW.stats.model3$IID <- black_ids
summary(cf.pred.binaryW.model3)
summary(cf.pred.binaryW.stats.model3)
write.csv(cf.pred.binaryW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model3_black.csv")
write.csv(cf.pred.binaryW.stats.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model3_black.csv")

# Asian Ancestry Results                                        
message("\nAsian Ancestry")                                      
cf.pred.binaryW.model3 <- predict(causal.forest.binaryW.model3, newdata = X.model3.mat[asian_idx, ], estimate.variance = TRUE)
cf.pred.binaryW.stats.model3 <- compute_stats(cf.pred.binaryW.model3)
cf.pred.binaryW.model3$IID <- asian_ids
cf.pred.binaryW.stats.model3$IID <- asian_ids
summary(cf.pred.binaryW.model3)
summary(cf.pred.binaryW.stats.model3)
write.csv(cf.pred.binaryW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model3_asian.csv")
write.csv(cf.pred.binaryW.stats.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model3_asian.csv")

# Chinese Ancestry Results                                        
message("\nChinese Ancestry")                                       
cf.pred.binaryW.model3 <- predict(causal.forest.binaryW.model3, newdata = X.model3.mat[chinese_idx, ], estimate.variance = TRUE)
cf.pred.binaryW.stats.model3 <- compute_stats(cf.pred.binaryW.model3)
cf.pred.binaryW.model3$IID <- chinese_ids
cf.pred.binaryW.stats.model3$IID <- chinese_ids
summary(cf.pred.binaryW.model3)
summary(cf.pred.binaryW.stats.model3)
write.csv(cf.pred.binaryW.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_var_model3_chinese.csv")
write.csv(cf.pred.binaryW.stats.model3, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/HDL/CAD/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_model3_chinese.csv")
