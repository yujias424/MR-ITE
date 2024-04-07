suppressPackageStartupMessages({
  library(grf)
  library(tidyverse)
  library(data.table)
  library(caret)
  library(ggplot2)
  library(xgboost)
  
  source("~/Project/2021-11-10-individual_MR/bin/propensity_score_matching/psm.R")
})

# X, Y, T, Z from their original source file
X <- fread("~/Project/2021-11-10-individual_MR/dat/analysis/pheno_data/ukbb.covariate.ite", sep = ",")
Z <- fread("~/Project/2021-11-10-individual_MR/dat/analysis/prs/ldl_prs.best", sep = " ")
W <- fread("~/Project/2021-11-10-individual_MR/dat/analysis/pheno_data/ukbb.LDL.complete", sep = "\t")
Y <- fread("~/Project/2021-11-10-individual_MR/dat/analysis/pheno_data/Heart_Disease_outcome_updated", sep = "\t")

# homonize the ID included in the analysis
Z_id = unique(Z$IID)
X_id = unique(X$IID)
W_id = unique(W$IID)
Y_id = unique(Y$IID) 

selected_id = Reduce(intersect, list(X_id, Y_id, W_id, Z_id))

X <- X[X$IID %in% selected_id, ]
Y <- Y[Y$IID %in% selected_id, ]
W <- W[W$IID %in% selected_id, ]
Z <- Z[Z$IID %in% selected_id, ]

samples <- read.csv("~/Project/2021-11-10-individual_MR/dat/analysis/psm_res/heart_disease_psm.csv")

Z <- Z[Z$IID %in% samples$sample]
W <- W[W$IID %in% samples$sample]

ctrl_samples <- Z$PRS[]

which(W$`30780-0.0`)

cor.test(W$`30780-0.0`, Z$PRS)

W.vector <- as.vector(W$`30780-0.0`)
W.vector.binary.cutoff1 <- as.integer(ifelse(W.vector <= 2.6, 0, 1))
W.vector.binary.cutoff1 <- as.integer(W.vector.binary.cutoff1)

W.vector.binary.cutoff2 <- as.integer(ifelse(W.vector <= 3.3, 0, 1))
W.vector.binary.cutoff2 <- as.integer(W.vector.binary.cutoff2)

ct1.ctrl <- Z$PRS[which(W.vector.binary.cutoff1==0)]
ct1.tret <- Z$PRS[which(W.vector.binary.cutoff1==1)]

t.test(ct1.ctrl, ct1.tret, alternative = c("two.sided", "less", "greater"))

ct2.ctrl <- Z$PRS[which(W.vector.binary.cutoff2==0)]
ct2.tret <- Z$PRS[which(W.vector.binary.cutoff2==1)]

t.test(ct2.ctrl, ct2.tret, alternative = c("two.sided", "less", "greater"))
