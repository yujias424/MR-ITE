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
  library(xlsx)
  library(ivreg)
  library(ggsignif) 
  
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/psm.R")
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/test_heterogeneity.R")
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/find_best_tree.R")
})

# ==================
# Data Preprocessing
# ==================
# X, Y, T, Z from their original source file
X.set1 <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.traits.set1.gz", sep = "\t")
X.set2 <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.traits.set2.gz", sep = "\t")
X.set3 <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/HDL/CAD/ukbb.covariate.HDL.set3.gz", sep = "\t")

W <- fread("~/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/HDL/CAD/ukbb.phenotype.HDL.mgdL", sep = "\t")

Y.date1 <- fread("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date1_james2022.gz", sep = "\t")
Y.date2 <- fread("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date2_james2022.gz", sep = "\t")
Y.date3 <- fread("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date3_james2022.gz", sep = "\t")

dim(X.set1); dim(X.set1); dim(X.set1); dim(X.set1); dim(W); dim(Y.date1); dim(Y.date2); dim(Y.date3)

# homonize the ID included in the analysis
selected_id_set1 <- Reduce(intersect, list(X.set1$IID, W$IID, Y.date1$IID))
selected_id_set2 <- Reduce(intersect, list(X.set2$IID, W$IID, Y.date2$IID))
selected_id_set3 <- Reduce(intersect, list(X.set3$IID, W$IID, Y.date3$IID))

# model 1
X.model1 <- X.set1[X.set1$IID %in% selected_id_set1, ]
Y.model1 <- Y.date1[Y.date1$IID %in% selected_id_set1, ]
W.model1 <- W[W$IID %in% selected_id_set1, ]
X.model1 <- cbind(X.model1, Y.model1[, c("age_CAD")])
Y.model1 <- Y.model1[, c("IID", "CAD")]
X.model1$`21022-0.0` <- NULL
X.model1$`21022-0.0` <- X.model1$age_CAD
X.model1$age_CAD <- NULL
dim(X.model1); dim(Y.model1); dim(W.model1); 

# model 2
X.model2 <- X.set2[X.set2$IID %in% selected_id_set2, ]
Y.model2 <- Y.date2[Y.date2$IID %in% selected_id_set2, ]
W.model2 <- W[W$IID %in% selected_id_set2, ]
X.model2 <- cbind(X.model2, Y.model2[, c("htn", "t2dm", "heart_failure", "hemorrhage_stroke", "ischemic_stroke", "age_CAD")]) # "copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "stroke", "ischemic_stroke"
Y.model2 <- Y.model2[, c("IID", "CAD")]
X.model2$`21022-0.0` <- NULL
X.model2$`21022-0.0` <- X.model2$age_CAD
X.model2$age_CAD <- NULL
dim(X.model2); dim(Y.model2); dim(W.model2); 

# model 3
X.model3 <- X.set3[X.set3$IID %in% selected_id_set3, ]
Y.model3 <- Y.date3[Y.date3$IID %in% selected_id_set3, ]
W.model3 <- W[W$IID %in% selected_id_set3, ]
X.model3 <- cbind(X.model3, Y.model3[, c("htn", "t2dm", "heart_failure", "hemorrhage_stroke", "ischemic_stroke", "age_CAD")])
Y.model3 <- Y.model3[, c("IID", "CAD")]
X.model3$`21022-0.0` <- NULL
X.model3$`21022-0.0` <- X.model3$age_CAD
X.model3$age_CAD <- NULL
dim(X.model3); dim(Y.model3); dim(W.model3); 

# model 1
setnames(X.model1 , c("22001-0.0", "21022-0.0", "22000-0.0", 
                      "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
                      "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10"),
                    c("Gender", "Age", "Genotype_Batch", 
                      "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
fwrite(X.model1 , "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/HDL/CAD/ukbb.covariate.prs.HDL.set1", sep = "\t")

# model 2
setnames(X.model2 , c("22001-0.0", "21022-0.0", "22000-0.0", 
                      "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
                      "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10"),
                    c("Gender", "Age", "Genotype_Batch", 
                      "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
fwrite(X.model2 , "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/HDL/CAD/ukbb.covariate.prs.HDL.set2", sep = "\t")

# model 3
setnames(X.model3 , c("22001-0.0", "21022-0.0", "22000-0.0", 
                      "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
                      "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10"),
                    c("Gender", "Age", "Genotype_Batch", 
                      "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
fwrite(X.model3 , "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/HDL/CAD/ukbb.covariate.prs.HDL.set3", sep = "\t")