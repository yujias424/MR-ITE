#' This is the code to format the covariates file and phenotype file to match PRS requirements
#' Mainly FID, IID, PC
#' HDL
#' 
#' @author: Yujia Shi
#' @date: 2022.09.01

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(dplyr)
  library(plyr)
})

covar <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.imputed.fullset.gz") # exclude non-genetic white and include all PC
covar_ethnic <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.gz")
selected_cov <- c("eid", "21000-0.0")
covar_ethnic <- covar_ethnic[, ..selected_cov]
colnames(covar_ethnic)[1] <- "IID"
covar <- covar[covar_ethnic, on=c("IID"), nomatch=FALSE]

table(covar$`21000-0.0`)

# We will focus on the highest level of the ethnic, defined by White = 1, Asian or Asian British = 3, Black or Black Britsh = 4, Chinese = 5)
covar <- covar[covar$`21000-0.0` %in% c("1", "1001", "1002", "1003", "3", "3001", "3002", "3003", "3004", "4", "4001", "4002", "4003", "5")]
table(covar$`21000-0.0`)
covar$ethnic_group <- mapvalues(covar$`21000-0.0`, 
                                from = c("1", "1001", "1002", "1003", "3", "3001", "3002", "3003", "3004", "4", "4001", "4002", "4003", "5"),
                                to = c("1", "1", "1", "1", "2", "2", "2", "2", "2", "3", "3", "3", "3", "4")) # White: 1; Asian: 2; Black: 3; Chinese: 4
covar$`21000-0.0` <- NULL

# # check whether NA is 0 across all covariates as we have done the imputation.
# covar %>% summarise(across(everything(), ~ sum(is.na(.))))

# remove patients refused to report whether they smoke before and set dummy variables.
covar <- covar[covar$`20116-0.0` != -3, ]
smoker_type_length <- length(table(covar$`20116-0.0`))
smoker_mat <- matrix(0, nrow = dim(covar)[1], ncol = smoker_type_length)
for (i in 1:smoker_type_length){
  smoker_mat[which(covar$`20116-0.0` == names(table(covar$`20116-0.0`))[i]), i] <- 1
}
smoker_mat <- as.data.frame(smoker_mat)
smoker_mat <- as.data.frame(lapply(smoker_mat, as.integer))
colnames(smoker_mat) <- c("Non_smoker", "Previous_smoker", "Current_smoker")
covar$`20116-0.0` <- NULL
covar <- cbind(covar, smoker_mat)

# remove patients refused to report whether they drink alcohol before and set dummy variables.
covar <- covar[covar$`20117-0.0` != -3, ]
alcohol_type_length <- length(table(covar$`20117-0.0`))
alcohol_mat <- matrix(0, nrow = dim(covar)[1], ncol = alcohol_type_length)
for (i in 1:alcohol_type_length){
  alcohol_mat[which(covar$`20117-0.0` == names(table(covar$`20117-0.0`))[i]), i] <- 1
}
alcohol_mat <- as.data.frame(alcohol_mat)
alcohol_mat <- as.data.frame(lapply(alcohol_mat, as.integer))
colnames(alcohol_mat) <- c("Non_alcohol_drinker", "Previous_alcohol_drinker", "Current_alcohol_drinker")
covar$`20117-0.0` <- NULL
covar <- cbind(covar, alcohol_mat)

# convert drug to dummy variable
covar <- covar[!covar$medi_dia_bp_chole %in% c(-3, -1), ]
covar$medi_dia_bp_chole[which(covar$medi_dia_bp_chole == -7)] <- 0
medi_type_length <- length(table(covar$medi_dia_bp_chole))
medi_mat <- matrix(0, nrow = dim(covar)[1], ncol = medi_type_length)
for (i in 1:medi_type_length){
  medi_mat[which(covar$medi_dia_bp_chole == i-1), i] <- 1
}
medi_mat <- as.data.frame(medi_mat)
medi_mat <- as.data.frame(lapply(medi_mat, as.integer))
colnames(medi_mat) <- c("No_medication", "Cholesterol_lowering_medication",
                        "Blood_pressure_medication", "Insulin", "Hormone_replacement_therapy",
                        "Oral_contraceptive_pill_or_minipill")
covar$medi_dia_bp_chole <- NULL
covar <- cbind(covar, medi_mat)

# reset id to char
covar$FID <- as.character(covar$FID)
covar$IID <- as.character(covar$IID)

# exclude several covariates that we are not interested at all.
covar <- covar[, .SD, .SDcols = !c('whr', 'No_medication')]

# Phenotype (CRP)
pheno.CRP.mgL <- covar[, c("FID", "IID", "30710-0.0")]

# convert mg/L to mg/dL
pheno.CRP.mgdL <- pheno.CRP.mgL
pheno.CRP.mgdL$`30710-0.0` <- pheno.CRP.mgdL$`30710-0.0` * 0.1

# Covariates set (three types of covariates sets)
covar.set.1 <- c("FID", "IID", "22001-0.0", "21022-0.0", 
                 "22000-0.0",
                 "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
                 "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10", "ethnic_group") # only gender and age + Genotype measurement batch 10 PCs
covar.set.2 <- c("FID", "IID", "22001-0.0", "21022-0.0", 
                 "22000-0.0",
                 "4079-0.0", "4080-0.0", "189-0.0", 
                 "Non_smoker", "Previous_smoker", "Current_smoker",
                 "Non_alcohol_drinker", "Previous_alcohol_drinker", "Current_alcohol_drinker",
                 "Cholesterol_lowering_medication", "Blood_pressure_medication", "Insulin", # "No_medication", 
                 "21001-0.0", "23099-0.0", "21002-0.0", # "whr", 
                 "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
                 "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10", "ethnic_group") # set 1 + covariates assessed at attending date.
covar.set.3 <- c("FID", "IID", "22001-0.0", "21022-0.0", 
                 "22000-0.0",
                 "4079-0.0", "4080-0.0", "189-0.0", 
                 "Non_smoker", "Previous_smoker", "Current_smoker",
                 "Non_alcohol_drinker", "Previous_alcohol_drinker", "Current_alcohol_drinker",
                 "Cholesterol_lowering_medication", "Blood_pressure_medication", "Insulin", # "No_medication", 
                 "21001-0.0", "23099-0.0", "21002-0.0", # "whr", 
                 "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
                 "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10",
                 "30630-0.0", "30640-0.0", "30650-0.0", "30660-0.0", "30670-0.0",
                 "30690-0.0", "30700-0.0", "30710-0.0", "30720-0.0", "30680-0.0", 
                 "30730-0.0", "30740-0.0", "30750-0.0", "30760-0.0", "30770-0.0", 
                 "30780-0.0", "30790-0.0", "30810-0.0", "30830-0.0", "30840-0.0", "30860-0.0",
                 "30850-0.0",  "30870-0.0", "30880-0.0", "30890-0.0", "ethnic_group") # set 2 + biological indexes 
covar.set.prs <- c("FID", "IID", "22001-0.0", "21022-0.0", "22000-0.0", "54-0.0",
                   "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
                   "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10")

# covars
covar.set.1.df <- covar[, ..covar.set.1]
covar.set.2.df <- covar[, ..covar.set.2]
covar.set.3.df <- covar[, ..covar.set.3]
covar.set.3.CRP.df <- covar.set.3.df[, .SD, .SDcols = !c("30710-0.0")] # CRP

# write covar data
dir.create("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/external_analysis_for_paper_revision/CRP/CAD/", showWarnings = FALSE, recursive = TRUE)

fwrite(covar.set.1.df, "/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/external_analysis_for_paper_revision/CRP/CAD/ukbb.covariate.traits.set1.gz", sep = "\t")
fwrite(covar.set.2.df, "/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/external_analysis_for_paper_revision/CRP/CAD/ukbb.covariate.traits.set2.gz", sep = "\t")
fwrite(covar.set.3.CRP.df, "/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/external_analysis_for_paper_revision/CRP/CAD/ukbb.covariate.CRP.set3.gz", sep = "\t")
fwrite(pheno.CRP.mgdL, "/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/external_analysis_for_paper_revision/CRP/CAD/ukbb.phenotype.CRP.mgdL", sep = "\t")


message("Finished the pipeline.")
