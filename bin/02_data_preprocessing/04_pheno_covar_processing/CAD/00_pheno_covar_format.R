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
})

covar <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.imputed.gz") # exclude non-genetic white and include all PC

# check whether NA is 0 across all covariates as we have done the imputation.
covar %>% summarise(across(everything(), ~ sum(is.na(.))))

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
# covar <- covar[, .SD, .SDcols = !c('30680-0.0', '30840-0.0', '30860-0.0')]
covar <- covar[, .SD, .SDcols = !c('whr', 'No_medication')]

# Phenotype (HDL, LDL, TC, TG)
pheno.HDL.mmolL <- covar[, c("FID", "IID", "30760-0.0")]
pheno.LDL.mmolL <- covar[, c("FID", "IID", "30780-0.0")]
pheno.TC.mmolL <- covar[, c("FID", "IID", "30690-0.0")]
pheno.TG.mmolL <- covar[, c("FID", "IID", "30870-0.0")]

# convert mmol/L to mg/dL
pheno.HDL.mgdL <- pheno.HDL.mmolL
pheno.HDL.mgdL$`30760-0.0` <- pheno.HDL.mgdL$`30760-0.0` * 38.67
pheno.LDL.mgdL <- pheno.LDL.mmolL
pheno.LDL.mgdL$`30780-0.0` <- pheno.LDL.mgdL$`30780-0.0` * 38.67
pheno.TC.mgdL <- pheno.TC.mmolL
pheno.TC.mgdL$`30690-0.0` <- pheno.TC.mgdL$`30690-0.0` * 38.67
pheno.TG.mgdL <- pheno.TG.mmolL
pheno.TG.mgdL$`30870-0.0` <- pheno.TG.mgdL$`30870-0.0` * 88.57

# Covariates set (three types of covariates sets)
covar.set.1 <- c("FID", "IID", "22001-0.0", "21022-0.0", 
                 "22000-0.0",
                 "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
                 "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10") # only gender and age + Genotype measurement batch 10 PCs
covar.set.2 <- c("FID", "IID", "22001-0.0", "21022-0.0", 
                 "22000-0.0",
                 "4079-0.0", "4080-0.0", "189-0.0", 
                 "Non_smoker", "Previous_smoker", "Current_smoker",
                 "Non_alcohol_drinker", "Previous_alcohol_drinker", "Current_alcohol_drinker",
                 "Cholesterol_lowering_medication", "Blood_pressure_medication", "Insulin", # "No_medication", 
                 "21001-0.0", "23099-0.0", "21002-0.0", # "whr", 
                 "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
                 "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10") # set 1 + covariates assessed at attending date.
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
                 "30850-0.0",  "30870-0.0", "30880-0.0", "30890-0.0") # set 2 + biological indexes 
covar.set.prs <- c("FID", "IID", "22001-0.0", "21022-0.0", "22000-0.0", "54-0.0",
                   "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
                   "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10")

# covars
covar.set.1.df <- covar[, ..covar.set.1]
covar.set.2.df <- covar[, ..covar.set.2]
covar.set.3.df <- covar[, ..covar.set.3]
covar.set.3.HDL.df <- covar.set.3.df[, .SD, .SDcols = !c("30760-0.0", "30690-0.0", "30630-0.0")] # HDL, TC, Apolipoprotein A
covar.set.3.LDL.df <- covar.set.3.df[, .SD, .SDcols = !c("30780-0.0", "30690-0.0", "30640-0.0", "30790-0.0")] # LDL, TC, Apolipoprotein B, Lipoprotein A
covar.set.3.TC.df <- covar.set.3.df[, .SD, .SDcols = !c("30780-0.0", "30760-0.0", "30690-0.0", "30630-0.0", "30640-0.0", "30790-0.0")] # LDL, HDL, TC, Apolipoprotein A, Apolipoprotein B, Lipoprotein A
covar.set.3.TG.df <- covar.set.3.df[, .SD, .SDcols = !c("30870-0.0")]
# covar.set.prs.df <- covar[, ..covar.set.prs]

# write covar data
fwrite(covar.set.1.df, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.traits.set1.gz", sep = "\t")
fwrite(covar.set.2.df, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.traits.set2.gz", sep = "\t")

fwrite(covar.set.3.HDL.df, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/HDL/CAD/ukbb.covariate.HDL.set3.gz", sep = "\t")
fwrite(covar.set.3.LDL.df, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/LDL/CAD/ukbb.covariate.LDL.set3.gz", sep = "\t")
fwrite(covar.set.3.TC.df, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TC/CAD/ukbb.covariate.TC.set3.gz", sep = "\t")
fwrite(covar.set.3.TG.df, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TG/CAD/ukbb.covariate.TG.set3.gz", sep = "\t")

fwrite(pheno.HDL.mmolL, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/HDL/CAD/ukbb.phenotype.HDL.mmolL", sep = "\t")
fwrite(pheno.HDL.mgdL, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/HDL/CAD/ukbb.phenotype.HDL.mgdL", sep = "\t")
fwrite(pheno.LDL.mmolL, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/LDL/CAD/ukbb.phenotype.LDL.mmolL", sep = "\t")
fwrite(pheno.LDL.mgdL, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/LDL/CAD/ukbb.phenotype.LDL.mgdL", sep = "\t")
fwrite(pheno.TC.mmolL, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TC/CAD/ukbb.phenotype.TC.mmolL", sep = "\t")
fwrite(pheno.TC.mgdL, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TC/CAD/ukbb.phenotype.TC.mgdL", sep = "\t")
fwrite(pheno.TG.mmolL, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TG/CAD/ukbb.phenotype.TG.mmolL", sep = "\t")
fwrite(pheno.TG.mgdL, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TG/CAD/ukbb.phenotype.TG.mgdL", sep = "\t")

# # rename PC covariate for PRS analysis
# setnames(covar.set.prs.df , c("22001-0.0", "21022-0.0", "22000-0.0", # "54-0.0", 
#                               "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
#                               "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10"),
#                             c("Gender", "Age", "Genotype_Batch", # "Assessment_Centre", 
#                               "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
# fwrite(covar.set.prs.df , "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.prs", sep = "\t")

# # rename PC covariate
# # model 1
# setnames(covar.set.1.df , c("22001-0.0", "21022-0.0", "22000-0.0", 
#                               "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
#                               "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10"),
#                             c("Gender", "Age", "Genotype_Batch", 
#                               "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
# fwrite(covar.set.1.df , "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.prs.set1", sep = "\t")

# # model 2
# setnames(covar.set.2.df , c("22001-0.0", "21022-0.0", "22000-0.0", 
#                               "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
#                               "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10"),
#                             c("Gender", "Age", "Genotype_Batch", 
#                               "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
# fwrite(covar.set.2.df , "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.prs.set2", sep = "\t")

# # model 3
# setnames(covar.set.3.HDL.df , c("22001-0.0", "21022-0.0", "22000-0.0", 
#                               "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
#                               "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10"),
#                             c("Gender", "Age", "Genotype_Batch", 
#                               "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
# fwrite(covar.set.3.HDL.df , "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.prs.HDL.set3", sep = "\t")

# setnames(covar.set.3.LDL.df , c("22001-0.0", "21022-0.0", "22000-0.0", 
#                               "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
#                               "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10"),
#                             c("Gender", "Age", "Genotype_Batch", 
#                               "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
# fwrite(covar.set.3.LDL.df , "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.prs.LDL.set3", sep = "\t")

# setnames(covar.set.3.TC.df , c("22001-0.0", "21022-0.0", "22000-0.0", 
#                               "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
#                               "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10"),
#                             c("Gender", "Age", "Genotype_Batch", 
#                               "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
# fwrite(covar.set.3.TC.df , "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.prs.TC.set3", sep = "\t")

# setnames(covar.set.3.TG.df , c("22001-0.0", "21022-0.0", "22000-0.0", 
#                               "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", 
#                               "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10"),
#                             c("Gender", "Age", "Genotype_Batch", 
#                               "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))
# fwrite(covar.set.3.TG.df , "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.prs.TG.set3", sep = "\t")

message("Finished the pipeline.")