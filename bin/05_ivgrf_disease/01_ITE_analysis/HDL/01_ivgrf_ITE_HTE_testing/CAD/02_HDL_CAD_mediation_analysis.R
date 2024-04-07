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
  library(HIMA)
  
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/psm.R")
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/test_heterogeneity.R")
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/find_best_tree.R")
})

# ==================
# Data Preprocessing
# ==================
# X, Y, T, Z from their original source file
X <- fread("~/Project/2022-09-01-individual_MR/dat/04_pheno_covar_data/HDL/ukbb.covariate.HDL.complete", sep = "\t")
Z <- fread("~/Project/2022-09-01-individual_MR/dat/06_PRS_calculation/HDL/CAD/1e-05/HDL_prs.best", sep = " ") # score_constd
W <- fread("~/Project/2022-09-01-individual_MR/dat/04_pheno_covar_data/HDL/ukbb.phenotype.HDL.complete.mgdL", sep = "\t")
Y <- fread("~/Project/2022-09-01-individual_MR/dat/03_outcome_data/CAD_outcome", sep = "\t")

W <- W[W$`30760-0.0` <= 80, ]

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
                         "22001-0.0", "4079-0.0", "4080-0.0", "189-0.0", "21022-0.0", "21001-0.0", "23099-0.0",               
                         "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5",                 
                         "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10",
                         "whr", # "British", "Indian", "Caribbean", "African",
                         "Cholesterol_lowering_medication", "Blood_pressure_medication", "No_medication", "Insulin", 
                         "Non_alcohol_drinker", "Previous_alcohol_drinker", "Current_alcohol_drinker", 
                         "Non_smoker", "Previous_smoker", "Current_smoker",
                         "30780-0.0", "30870-0.0", "30640-0.0", "30790-0.0", # lipid-related covariates
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
W.mat <- W[W$IID %in% samples, c("30760-0.0")]
X.mat <- X[X$IID %in% samples, 3:dim(X)[2]]
Y.mat <- Y[Y$IID %in% samples, 2]
Z.mat <- Z[Z$IID %in% samples, 4]

W.vector <- as.vector(W.mat$`30760-0.0`)
Y.vector <- as.vector(Y.mat$CAD)
Y.vector <- as.integer(Y.vector)
Z.vector <- as.vector(Z.mat$PRS)

hima.res <- hima(W.vector, Y.vector, X.mat)
hima.res.selected.covar <- hima.res[hima.res$Bonferroni.p < 0.05, ]
tmp <- row.names(hima.res.selected.covar)
tmp <- ifelse(!startsWith(tmp, "`"), tmp, stringr::str_extract(tmp, "(?<=`).*(?=`)"))
tmp <- c(tmp)
X.mat <- as.data.frame(X.mat)
X.mat <- X.mat[, !(colnames(X.mat) %in% c(tmp))]