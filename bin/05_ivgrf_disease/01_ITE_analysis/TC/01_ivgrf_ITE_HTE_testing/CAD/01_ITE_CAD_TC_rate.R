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
                         "Waist-hip-ratio"="whr","Body Fat Percentage"="23099-0.0", "BMI"="21001-0.0", "Weight"="21002-0.0",
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "No medication"="No_medication", "Insulin"="Insulin", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker")

X.model2b.mat  <- X.model2b.mat  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0", 
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0",
                         "Waist-hip-ratio"="whr","Body Fat Percentage"="23099-0.0", "BMI"="21001-0.0", "Weight"="21002-0.0",
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "No medication"="No_medication", "Insulin"="Insulin", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker")

X.model3.mat  <- X.model3.mat  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0",   
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10",
                         "Waist-hip-ratio"="whr", "Body Fat Percentage"="23099-0.0", "BMI"="21001-0.0", "Weight"="21002-0.0", "Genotype Batch"="22000-0.0",
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "No medication"="No_medication", "Insulin"="Insulin", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Calcium"="30680-0.0","Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", 
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Total protein"="30860-0.0", "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0",
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")

# preparing vector for different model
# model 1
W.model1.vector <- as.vector(W.model1.mat$`30690-0.0`)
Y.model1.vector <- as.integer(as.vector(Y.model1.mat$CAD))
Z.model1.vector <- as.vector(Z.model1.mat$PRS)
W.model1.vector.binary <- as.integer(ifelse(W.model1.vector <= 220, 1, 0)) # mimic usage of statin

# model 1b
W.model1b.vector <- as.vector(W.model1b.mat$`30690-0.0`)
Y.model1b.vector <- as.integer(as.vector(Y.model1b.mat$CAD))
Z.model1b.vector <- as.vector(Z.model1b.mat$PRS)
Z.model1b.vector.binary <- as.integer(ifelse(Z.model1b.vector <= 0, 1, 0))
W.model1b.vector.binary <- as.integer(ifelse(W.model1b.vector <= 220, 1, 0)) # mimic usage of statin

# model 2
W.model2.vector <- as.vector(W.model2.mat$`30690-0.0`)
Y.model2.vector <- as.integer(as.vector(Y.model2.mat$CAD))
Z.model2.vector <- as.vector(Z.model2.mat$PRS)
W.model2.vector.binary <- as.integer(ifelse(W.model2.vector <= 220, 1, 0)) # mimic usage of statin

# model 2b
W.model2b.vector <- as.vector(W.model2b.mat$`30690-0.0`)
Y.model2b.vector <- as.integer(as.vector(Y.model2b.mat$CAD))
Z.model2b.vector <- as.vector(Z.model2b.mat$PRS)
W.model2b.vector.binary <- as.integer(ifelse(W.model2b.vector <= 220, 1, 0)) # mimic usage of statin

# model 3
W.model3.vector <- as.vector(W.model3.mat$`30690-0.0`)
Y.model3.vector <- as.integer(as.vector(Y.model3.mat$CAD))
Z.model3.vector <- as.vector(Z.model3.mat$PRS)
W.model3.vector.binary <- as.integer(ifelse(W.model3.vector <= 220, 1, 0)) # mimic usage of statin

# Applying Wu-Hausman Test to exploring the validity of applying instrument in the analysis.
# model 1
WHtest.mat.model1 <- cbind(X.model1.mat, Y.model1.mat, W.model1.mat, Z.model1.mat)
dim(WHtest.mat.model1)
colnames(WHtest.mat.model1)[dim(WHtest.mat.model1)[2]-1] <- "W"
res.model1 <- ivreg(CAD ~ . - W - PRS | W | PRS, data = WHtest.mat.model1)
summary(res.model1, diagnostics = TRUE)

# model 1b
WHtest.mat.model1b <- cbind(X.model1b.mat, Y.model1b.mat, W.model1b.mat, Z.model1b.mat)
dim(WHtest.mat.model1b)
colnames(WHtest.mat.model1b)[dim(WHtest.mat.model1b)[2]-1] <- "W"
res.model1b <- ivreg(CAD ~ . - W - PRS | W | PRS, data = WHtest.mat.model1b)
summary(res.model1b, diagnostics = TRUE)

# model 2
WHtest.mat.model2 <- cbind(X.model2.mat, Y.model2.mat, W.model2.mat, Z.model2.mat)
dim(WHtest.mat.model2)
colnames(WHtest.mat.model2)[dim(WHtest.mat.model2)[2]-1] <- "W"
res.model2 <- ivreg(CAD ~ . - W - PRS | W | PRS, data = WHtest.mat.model2)
summary(res.model2, diagnostics = TRUE)

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
# model 1b
# instrumental forest (continuous Z)
split1 <- createDataPartition(as.factor(Y.model1b.vector), p = 0.5)[[1]]
table(Y.model1b.vector[split1])
message("\nContinuous W IV-grf model 1b")
iv.forest.continuousW.model1b <- instrumental_forest(X.model1b.mat[split1, ], Y.model1b.vector[split1], W.model1b.vector.binary[split1], Z.model1b.vector.binary[split1], 
                                             num.threads = 10, num.trees = 5000, 
                                             sample.fraction = 0.02, min.node.size = 50, ci.group.size = 100)
iv.pred.continuousW.model1b <- predict(iv.forest.continuousW.model1b, estimate.variance = TRUE)
iv.pred.stats.continuousW.model1b <- compute_stats(iv.pred.continuousW.model1b)
summary(iv.pred.continuousW.model1b)
summary(iv.pred.stats.continuousW.model1b)

iv.pred.continuousW.model1b.split1 <- predict(iv.forest.continuousW.model1b, X = X.model1b.mat[-split1, ], estimate.variance = TRUE)

iv.forest.continuousW.model1b.full <- instrumental_forest(X.model1b.mat, Y.model1b.vector, W.model1b.vector.binary, Z.model1b.vector.binary, 
                                             num.threads = 10, num.trees = 5000, 
                                             sample.fraction = 0.02, min.node.size = 100, ci.group.size = 100)
average_treatment_effect(iv.forest.continuousW.model1b.full)

testrate <- rank_average_treatment_effect(iv.forest.continuousW.model1b, iv.pred.continuousW.model1b.split1$predictions)
plot(testrate)
-0.06946793 + 0.03005579 * 1.96
-0.04883862 + 0.01102224 * 1.96

mean(Z.model1b.vector[which(W.model1.vector.binary==1)], na.rm=T)
mean(Z.model1b.vector[which(W.model1.vector.binary==0)], na.rm=T)

# driv_model1b <- read.csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_continuousW_te_ul_model1b.csv")
# ivgrf_model1b <- read.csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model1b.csv")

# head(driv_model1b)
# driv_model1b <- driv_model1b[order(driv_model1b$IID), ]

# rate.res <- rank_average_treatment_effect.fit(
#   driv_model1b$point,
#   ivgrf_model1b$tau,
#   target = c("AUTOC", "QINI"),
#   q = seq(0.1, 1, by = 0.1),
#   R = 200,
#   sample.weights = NULL,
#   clusters = NULL
# )
# plot(rate.res)

# cor.test(ivgrf_model1b$tau, driv_model1b$point)
# plot(ivgrf_model1b$tau, driv_model1b$point)

# sort the ivreg.dat by tau
ivreg.dat <- cbind(X.model1b.mat[-split1, ],  Y.model1b.vector[-split1], W.model1b.vector[-split1], Z.model1b.vector[-split1], iv.pred.continuousW.model1b.split1$predictions)
ivreg.dat$ID <- selected_id_set1b[-split1]
colnames(ivreg.dat)[(length(colnames(ivreg.dat))-4): length(colnames(ivreg.dat))] <- c("Y", "W", "Z", "Tau", "ID")
ivreg.dat <- ivreg.dat[order(ivreg.dat$Tau, decreasing = T), ]

mean(ivgrf_model1b$tau)
test_iv <- ivreg(Y ~ . - W - Z - Tau - ID | W | Z, dat = ivreg.dat)
summary(test_iv)
ATE <- test_iv$coefficients['W']
print(ATE)
TOC_vector <- c()
for (i in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)){
       test_iv <- ivreg(Y ~ . - W - Z - Tau - ID | W | Z, dat = ivreg.dat[1:round(dim(ivreg.dat)[1] * i), ])
       assign(paste0("TOC_", i), test_iv$coefficients['W'] - ATE)
       TOC_vector <- c(TOC_vector, get(paste0("TOC_", i)))
}
mean(TOC_vector)
TOC_vector


priorities <- as.data.frame(ivgrf_model1b$tau, fix.empty.names = FALSE)

# dim(ress[["t"]]) <- c(5, dim(ress[["t0"]]))
# point.estimate <- ress[["t0"]]
# std.errors <- apply(ress[["t"]], c(2, 3), sd)
# point.estimate[abs(point.estimate) < 1e-15] <- 0
# std.errors[abs(std.errors) < 1e-15] <- 0
# output <- list()
# class(output) <- "rank_average_treatment_effect"
# output[["estimate"]] <- point.estimate[1, ]
# output[["std.err"]] <- std.errors[1, ]
# output

ate_ivreg <- ivreg(Y ~ . - W - Z - Tau - ID | W | Z, dat = ivreg.dat)
ATE <- ate_ivreg$coefficients['W']

rank_average_treatment_effect <- function(ivreg_dat,
                                          target = c("AUTOC", "QINI"),
                                          q = seq(0.1, 1, by = 0.1),
                                          R = 10,
                                          subset = NULL,
                                          debiasing.weights = NULL,
                                          compliance.score = NULL,
                                          num.trees.for.weights = 500) {
       
       ivreg_dat <- ivreg_dat[order(ivreg_dat$Tau), ]
       clusters <- 1:dim(ivreg_dat)[1]

       ate_ivreg <- ivreg(Y ~ . - W - Z - Tau - ID | W | Z, dat = ivreg_dat)
       ATE <- ate_ivreg$coefficients['W']
       
       boot.output <- boot_grf(
              data = ivreg_dat,
              # In case of two priorities do a paired bootstrap estimating both prios on same sample.
              statistic <- function(data, indices, q, ATE){
                     estimate_rate(data, indices, q, ATE)
              },
              R = R,
              clusters = clusters,
              half.sample = TRUE,
              q = q,
              ATE = ATE
       )

       boot.output
       # dim(boot.output[["t"]]) <- c(R, dim(boot.output[["t0"]]))
       # point.estimate <- boot.output[["t0"]]
       # std.errors <- apply(boot.output[["t"]], c(2, 3), sd)
       # point.estimate[abs(point.estimate) < 1e-15] <- 0
       # std.errors[abs(std.errors) < 1e-15] <- 0

       # return(c(point.estimate[1, ], std.errors[1, ]))
}

estimate_rate <- function(data, indices, q, ATE) {

       data <- data[order(data$Tau, decreasing = T), ]
       data <- data[indices, ]
       print(head(indices))
       print(head(data))
       TOC <- c()

       for (i in q){
              tmp_sample <- round(dim(data)[1]*i)
              print(tmp_sample)
              ate_q_ivreg <- ivreg(Y ~ . - W - Z - Tau - ID | W | Z, dat = data[1:tmp_sample, ])
              TOC <- c(TOC, (ate_q_ivreg$coefficients['W'] - ATE))
       }
       print(TOC)
       message("\n\n")
       RATE <- mean(TOC)
       RATE
}

boot_grf <- function(data, statistic, R, clusters, half.sample = TRUE, ...){

       samples.by.cluster <- split(seq_along(clusters), clusters)
       n <- length(samples.by.cluster)

       if (half.sample) {
              n.bs <- floor(n / 2)
              index.list <- replicate(R, unlist(samples.by.cluster[sample.int(n, n.bs, replace = FALSE)], use.names = FALSE), simplify = FALSE)
       } else {
              index.list <- replicate(R, unlist(samples.by.cluster[sample.int(n, replace = TRUE)], use.names = FALSE), simplify = FALSE)
       }

       t0 <- statistic(data, seq_len(NROW(data)), ...)
       t0 <- matrix(unlist(t0), ncol = length(t0))

       res <- lapply(seq_len(R), function(i) statistic(data, index.list[[i]][order(index.list[[i]])], ...))
       t <- matrix(, R, length(t0))
       for (r in seq_len(R)) {
              t[r, ] <- unlist(res[[r]])
       }

       list(t0 = t0, t = t)

}


ress <- rank_average_treatment_effect(ivreg.dat, q = c(0.1), R = 50, target = "AUTOC")
ress[["t0"]] - 1.96*sd(ress[["t"]])
sd(ress[["t"]])
ress[["t0"]]

dim(ress[["t"]]) <- c(10, dim(ress[["t0"]]))
apply(ress[["t"]], c(2, 3), sd)
sd(ress[["t"]])

priorities <- as.data.frame(ivgrf_model1b$tau, fix.empty.names = FALSE)
priorities <- priorities[order(priorities[,1], decreasing = T), ]
clusters <- 1:length(priorities)
samples.by.cluster <- split(seq_along(clusters), clusters)
n <- length(samples.by.cluster)
index.list <- replicate(5, unlist(samples.by.cluster[sample.int(n, floor(n/2), replace = FALSE)], use.names = FALSE), simplify = FALSE)
index.list[[1]]
head(index.list[[1]][order(index.list[[1]])])

aaa <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/03_outcome_data/First-diagnoses-before-20221219-included-GP-till20210930-and-hos-and-death-till20221219-tidied-at-20230915.csv")
head(aaa)


aaaa <- data.frame("W" =  W.model1b.vector.binary, "Z" = Z.model1b.vector)
split1 <- createDataPartition(as.factor(W.model1b.vector.binary), p = 0.1)[[1]]
ggplot(aaaa[split1, ], aes(Z, W)) + geom_point() + geom_smooth(method = "loess")

aaaa <- data.frame("Y" =  Y.model1b.vector, "Z" = Z.model1b.vector)
split1 <- createDataPartition(as.factor(Y.model1b.vector), p = 0.1)[[1]]
ggplot(aaaa[split1, ], aes(Z, Y)) + geom_point() + geom_smooth(method = "loess")
