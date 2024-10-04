#' This is the code to format the covariates file and phenotype file to match PRS requirements
#' Mainly FID, IID, PC
#' We will also perform imputation on the covariates
#' 
#' @author: Yujia Shi
#' @date: 2022.09.01

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(dplyr)
  library(doParallel)
  library(missRanger)
})

covar <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.gz") # exclude non-genetic white and include all PC

# only select patients with Genetic ethnic grouping = 1 and included in PCA calculation
covar <- covar[is.na(covar$`22006-0.0`)] # we only include Caucasian here, where 22006-0.0 = 1
# covar <- covar[covar$`22020-0.0` == 1] # we only include samples that are included in the PC computation.

# format the id and sex variable
# pheno.new <- data.table("FID" = covar$eid, "IID" = covar$eid, "30760-0.0" = covar$`30760-0.0`)
covar.new <- data.table("FID" = covar$eid, covar)
setnames(covar.new, "eid", "IID")

# convert the medicaiton to dummy variable
medication <- coalesce(covar.new$`6177-0.0`, covar.new$`6153-0.0`) # since 6153 is for female and 6177 for male
covar.new$medi_dia_bp_chole <- medication
covar.new$`6177-0.0` <- NULL
covar.new$`6153-0.0` <- NULL

covar.new.select <- covar.new[, c("FID", "IID", 
                                    "54-0.0", # UK Biobank assessment centre
                                    "189-0.0", # Townsend deprivation index at recruitment
                                    "4079-0.0", # Diastolic blood pressure, automated reading
                                    "4080-0.0", # Systolic blood pressure, automated reading
                                    "23099-0.0", # Body fat percentage
                                    "20116-0.0", # Smoking status
                                    "20117-0.0", # Alcohol drinker status
                                    "21001-0.0", # BMI
                                    "22001-0.0", # genetic sex
                                    "21002-0.0", # weight
                                    "21022-0.0", # age at recruitment
                                    "22000-0.0", # Genotype measurement batch
                                    "whr", # waist-hip-ratio
                                    "medi_dia_bp_chole", # medication subscription
                                    "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5",
                                    "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10",
                                    "30630-0.0", # Apolipoprotein A
                                    "30640-0.0", # Apolipoprotein B
                                    "30650-0.0", # Aspartate aminotransferase
                                    "30660-0.0", # Direct bilirubin
                                    "30670-0.0", # Uera
                                    "30680-0.0", # Calcium
                                    "30690-0.0", # Cholesterol
                                    "30700-0.0", # Creatinine
                                    "30710-0.0", # C-reactive protein
                                    "30720-0.0", # Cystatin C
                                    "30730-0.0", # Gamma glutamyltransferase
                                    "30740-0.0", # Glucose
                                    "30750-0.0", # Glycated haemoglobin (HbA1c)
                                    "30760-0.0", # HDL cholesterol
                                    "30770-0.0", # IGF-1
                                    "30780-0.0", # LDL direct
                                    "30790-0.0", # Lipoprotein A
                                    "30810-0.0", # Phosphate
                                    "30830-0.0", # SHBG
                                    "30840-0.0", # Total bilirubin
                                    "30850-0.0", # Testosterone
                                    "30860-0.0", # Total protein
                                    "30870-0.0", # Triglycerides
                                    "30880-0.0", # Urate
                                    "30890-0.0" # Vitamin D
                                  )]

print(covar.new.select %>% summarise(across(everything(), ~ sum(is.na(.)))))

covar.new.select <- as.data.frame(covar.new.select)

selected_cat_traits <- c("54-0.0", 
                         "22001-0.0", 
                         "20116-0.0", "22000-0.0",
                         "20117-0.0", "medi_dia_bp_chole")

# Convert the categorical features to factor
covar.new.select[c(selected_cat_traits)] <- lapply(covar.new.select[c(selected_cat_traits)], factor)

colnames(covar.new.select)[3:dim(covar.new.select)[2]] <- paste0("X", colnames(covar.new.select)[3:dim(covar.new.select)[2]])
colnames(covar.new.select)[3:dim(covar.new.select)[2]] <- stringr::str_replace(colnames(covar.new.select)[3:dim(covar.new.select)[2]], "-", "")

message("Start running missRanger")
covar.select.imp <- missRanger(data = covar.new.select, formula = . ~ . - FID - IID, returnOOB=T, pmm.k = 5)
message("Finished running missRanger")

save.image("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb_covariate_imputed_nonwhite.RData")

#' ===================
#' Generate the X data
#' ===================

load("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb_covariate_imputed_nonwhite.RData")

colnames(covar.select.imp) <- c("FID", "IID", 
                                "54-0.0", # UK Biobank assessment centre
                                "189-0.0", # Townsend deprivation index at recruitment
                                "4079-0.0", # Diastolic blood pressure, automated reading
                                "4080-0.0", # Systolic blood pressure, automated reading
                                "23099-0.0", # Body fat percentage
                                "20116-0.0", # Smoking status
                                "20117-0.0", # Alcohol drinker status
                                "21001-0.0", # BMI
                                "22001-0.0", # genetic sex
                                "21002-0.0", # weight
                                "21022-0.0", # age at recruitment
                                "22000-0.0", # Genotype measurement batch
                                "whr", # waist-hip-ratio
                                "medi_dia_bp_chole", # medication subscription
                                "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5",
                                "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10",
                                "30630-0.0", # Apolipoprotein A
                                "30640-0.0", # Apolipoprotein B
                                "30650-0.0", # Aspartate aminotransferase
                                "30660-0.0", # Direct bilirubin
                                "30670-0.0", # Uera
                                "30680-0.0", # Calcium
                                "30690-0.0", # Cholesterol
                                "30700-0.0", # Creatinine
                                "30710-0.0", # C-reactive protein
                                "30720-0.0", # Cystatin C
                                "30730-0.0", # Gamma glutamyltransferase
                                "30740-0.0", # Glucose
                                "30750-0.0", # Glycated haemoglobin (HbA1c)
                                "30760-0.0", # HDL cholesterol
                                "30770-0.0", # IGF-1
                                "30780-0.0", # LDL direct
                                "30790-0.0", # Lipoprotein A
                                "30810-0.0", # Phosphate
                                "30830-0.0", # SHBG
                                "30840-0.0", # Total bilirubin
                                "30850-0.0", # Testosterone
                                "30860-0.0", # Total protein
                                "30870-0.0", # Triglycerides
                                "30880-0.0", # Urate
                                "30890-0.0") # Vitamin D
fwrite(covar.select.imp, file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.imputed.nonwhite.gz")

