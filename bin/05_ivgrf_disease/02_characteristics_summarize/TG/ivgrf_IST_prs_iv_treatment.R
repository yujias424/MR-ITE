#' This code is to run a pilot study in using instrumental random forest and causal forest for analyzing the 
#' causal relationship between Triglycerides and ischemic stroke.
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
  
  source("~/Project/2021-11-10-individual_MR/bin/propensity_score_matching/psm.R")
  source("~/Project/2021-11-10-individual_MR/bin/simulation/support_functions/test_heterogeneity.R")
  source("~/Project/2021-11-10-individual_MR/bin/pilot_study/validation/find_best_tree.R")
})

# ==================
# Data Preprocessing
# ==================
# X, Y, T, Z from their original source file
X <- fread("~/Project/2022-09-01-individual_MR/dat/analysis/pheno_data/TG/ukbb.covariate.TG.complete", sep = "\t")
Z <- fread("~/Project/2022-09-01-individual_MR/dat/analysis/prs/mgdL/score_std/TG/ischemic_stroke/TG_prs.best", sep = " ")
W <- fread("~/Project/2022-09-01-individual_MR/dat/analysis/pheno_data/TG/ukbb.phenotype.TG.complete.mgdL", sep = "\t")
Y <- fread("~/Project/2022-09-01-individual_MR/dat/analysis/pheno_data/ischemic_stroke_outcome_updated", sep = "\t")

# remove patients got diabetes before.
X <- X[X$`2443-0.0` == 0 | X$`2443-0.0` == 1, ]

# convert drug to dummy variable 
medi_type_length <- length(table(X$medi_dia_bp_chole))
medi_mat <- matrix(0, nrow = dim(X)[1], ncol = medi_type_length)
for (i in 1:medi_type_length){
  medi_mat[which(X$medi_dia_bp_chole == i-1), i] <- 1
}
medi_mat <- as.data.frame(medi_mat)
medi_mat <- as.data.frame(lapply(medi_mat, as.integer))
colnames(medi_mat) <- c("No_medication", "Cholesterol_lowering_medication",
                        "Blood_pressure_medication", "Insulin", "Hormone_replacement_therapy",
                        "Oral_contraceptive_pill_or_minipill")
X$medi_dia_bp_chole <- NULL
X <- cbind(X, medi_mat)

# remove patients refused to report whether they drink alcohol before and set dummy variables.
X <- X[X$`20117-0.0` != -3, ]
alcohol_type_length <- length(table(X$`20117-0.0`))
alcohol_mat <- matrix(0, nrow = dim(X)[1], ncol = alcohol_type_length)
for (i in 1:alcohol_type_length){
  alcohol_mat[which(X$`20117-0.0` == names(table(X$`20117-0.0`))[i]), i] <- 1
}
alcohol_mat <- as.data.frame(alcohol_mat)
alcohol_mat <- as.data.frame(lapply(alcohol_mat, as.integer))
colnames(alcohol_mat) <- c("Non-alcohol drinker", "Previous alcohol drinker", "Current alcohol drinker")
X$`20117-0.0` <- NULL
X <- cbind(X, alcohol_mat)

# remove patients refused to report whether they smoke before and set dummy variables.
X <- X[X$`20116-0.0` != -3, ]
smoker_type_length <- length(table(X$`20116-0.0`))
smoker_mat <- matrix(0, nrow = dim(X)[1], ncol = smoker_type_length)
for (i in 1:smoker_type_length){
  smoker_mat[which(X$`20116-0.0` == names(table(X$`20116-0.0`))[i]), i] <- 1
}
smoker_mat <- as.data.frame(smoker_mat)
smoker_mat <- as.data.frame(lapply(smoker_mat, as.integer))
colnames(smoker_mat) <- c("Non-smoker", "Previous smoker", "Current smoker")
X$`20116-0.0` <- NULL
X <- cbind(X, smoker_mat)

# gender as integer
X$`22001-0.0` <- as.integer(X$`22001-0.0`)

# convert ethnic into dummy variable
X <- X[X$`21000-0.0` > 0]
table(X$`21000-0.0`)
selected_ethnic <- c(1001, 3001, 4001, 4002)
X <- X[X$`21000-0.0` %in% selected_ethnic, ]
ethnic_type_length <- length(table(X$`21000-0.0`))
ethnic_mat <- matrix(0, nrow = dim(X)[1], ncol = ethnic_type_length)
for (i in 1:ethnic_type_length){
  ethnic_mat[which(X$`21000-0.0` == names(table(X$`21000-0.0`))[i]), i] <- 1
}
ethnic_mat <- as.data.frame(ethnic_mat)
ethnic_mat <- as.data.frame(lapply(ethnic_mat, as.integer))
colnames(ethnic_mat) <- c("British", "Indian", "Caribbean", "African")
X$`21000-0.0` <- NULL
X <- cbind(X, ethnic_mat)

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
selected.covariates <- c("FID", "IID", "22001-0.0", "4079-0.0", "4080-0.0", "189-0.0", # "2443-0.0",              
                         "21001-0.0", "21022-0.0", "22009-0.1", "22009-0.2",                          
                         "22009-0.3", "22009-0.4", "22009-0.5", "22009-0.6", "22009-0.7", "22009-0.8",                          
                         "22009-0.9", "22009-0.10", 
                         "whr", "British", "Indian", "Caribbean", "African", 
                         "No_medication", "Cholesterol_lowering_medication", "Blood_pressure_medication", "Insulin", 
                         "Hormone_replacement_therapy", "Oral_contraceptive_pill_or_minipill",   
                         "Non-alcohol drinker", "Previous alcohol drinker", "Current alcohol drinker",
                         "Non-smoker", "Previous smoker", "Current smoker",
                         "30630-0.0", "30640-0.0", "30690-0.0", "30760-0.0", "30780-0.0", "30790-0.0", # lipid-related covariates
                         "30700-0.0", "30710-0.0", "30720-0.0", "30730-0.0", "30740-0.0", "30750-0.0", 
                         "30650-0.0", "30660-0.0", "30670-0.0", "30770-0.0",
                         "30810-0.0", "30830-0.0", "30850-0.0", "30860-0.0", "30880-0.0", "30890-0.0") 
X <- dplyr::select(X, all_of(selected.covariates))

X$`21022-0.0` <- as.numeric(X$`21022-0.0`)
X$`4079-0.0` <- as.numeric(X$`4079-0.0`)
X$`4080-0.0` <- as.numeric(X$`4080-0.0`)

X <- as.data.table(X)
setnames(X, c("FID", "IID", "22001-0.0", "4079-0.0", "4080-0.0", "189-0.0",                   
              "21001-0.0", "21022-0.0", "22009-0.1", "22009-0.2",                          
              "22009-0.3", "22009-0.4", "22009-0.5", "22009-0.6", "22009-0.7", "22009-0.8",                        
              "22009-0.9", "22009-0.10",
              "whr", "British", "Indian", "Caribbean", "African",
              "No_medication", "Cholesterol_lowering_medication", "Blood_pressure_medication", "Insulin",                            
              "Hormone_replacement_therapy", "Oral_contraceptive_pill_or_minipill",            
              "Non-alcohol drinker", "Previous alcohol drinker", "Current alcohol drinker",
              "Non-smoker", "Previous smoker", "Current smoker",
              "30630-0.0", "30640-0.0", "30690-0.0", "30760-0.0", "30780-0.0", "30790-0.0", # lipid-related covariates
              "30700-0.0", "30710-0.0", "30720-0.0", "30730-0.0", "30740-0.0", "30750-0.0", 
              "30650-0.0", "30660-0.0", "30670-0.0", "30770-0.0",
              "30810-0.0", "30830-0.0", "30850-0.0", "30860-0.0", "30880-0.0", "30890-0.0"),
            c("FID", "IID", "Gender", "diastolic blood pressure", "systolic blood pressure", "Townsend deprivation index",
              "BMI", "Age",  "PC1",  "PC2", 
              "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", 
              "PC9", "PC10",
              "whr", "British", "Indian", "Caribbean", "African",
              "No_medication", "Cholesterol_lowering_medication", "Blood_pressure_medication", "Insulin",                            
              "Hormone_replacement_therapy", "Oral_contraceptive_pill_or_minipill",            
              "Non-alcohol drinker", "Previous alcohol drinker", "Current alcohol drinker",
              "Non-smoker", "Previous smoker", "Current smoker",
              "Apolipoprotein A", "Apolipoprotein B", "Total Cholesterol", "HDL-C", "LDL-C", "Lipoprotein A", 
              "Creatinine", "C-reactive protein", "Cystatin C", "Gamma glutamyltransferase",  "Glucose", "HbA1c", 
              "Aspartate aminotransferase", "Direct bilirubin", "Urea", "IGF-1",
              "Phosphate", "SHBG", "Testosterone", "Total protein", "Urate", "Vitamin D"))
X <- as.data.frame(X)

X$Insulin <- NULL
X$Oral_contraceptive_pill_or_minipill <- NULL
X$No_medication <- NULL
X$Hormone_replacement_therapy <- NULL

X <- cbind(X, Y[, c("copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "stroke", "cad")])
Y <- Y[, c("IID", "ischemic_stroke")]

samples <- fread("~/Project/2022-09-01-individual_MR/dat/analysis/psm_res/TG/ischemic_stroke_psm_4.csv")
samples <- samples$sample

W <- W[W$IID %in% samples, c("30870-0.0")]
X <- X[X$IID %in% samples, 3:dim(X)[2]]
Y <- Y[Y$IID %in% samples, 2]
Z <- Z[Z$IID %in% samples, 4]

W.vector <- as.vector(W$`30870-0.0`)
Y.vector <- as.vector(Y$ischemic_stroke)
Y.vector <- as.integer(Y.vector)
Z.vector <- as.vector(Z$PRS)
Z.vector.binary <- ifelse(Z.vector < quantile(Z.vector, probs = 0.75), 0, 1)
Z.vector.binary <- as.integer(Z.vector.binary)

W.vector.binary.cutoff <- as.integer(ifelse(W.vector <= 150, 1, 0))
W.vector.binary.cutoff <- as.integer(W.vector.binary.cutoff)

# ===============
# Formal Analysis
# ===============
binary_var <- c("Gender", "British", "Indian", "Caribbean", "African", 
                "Cholesterol_lowering_medication", "Blood_pressure_medication", 
                "Non-alcohol drinker", "Previous alcohol drinker", "Current alcohol drinker", 
                "Non-smoker", "Previous smoker", "Current smoker",
                "copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "stroke", "cad")
continuous_summary <- as.data.frame(do.call(cbind, lapply(X[, !(colnames(X) %in% c(binary_var))], summary)))
binary_summary <- as.data.frame(do.call(cbind, lapply(X[, c(binary_var)], table)))

write.csv(continuous_summary, file = "~/Project/2022-09-01-individual_MR/res/table/summary_statistics/TG_IS_continuous.csv")
write.csv(binary_summary, file = "~/Project/2022-09-01-individual_MR/res/table/summary_statistics/TG_IS_binary.csv")
