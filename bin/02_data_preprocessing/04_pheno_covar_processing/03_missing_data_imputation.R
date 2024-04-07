#' This is the code to perform imputation on dataset.
#' Mainly FID, IID, PC
#' 
#' @author: Yujia Shi
#' @date: 2022.08.12

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(dplyr)
  library(doParallel)
  library(missRanger)
})

#' ==================
#' Perform imputation
#' ==================
covar <- fread("/home/yujia/Project/2022-08-17-individual_MR_with_imputation/dat/analysis/pheno_data/ukbb.covariate.gz")

medication <- coalesce(covar$`6177-0.0`, covar$`6153-0.0`) # since 6153 is for female and 6177 for male
covar$medi_dia_bp_chole <- medication
covar$`6177-0.0` <- NULL
covar$`6153-0.0` <- NULL

covar.select <- covar[, c("eid", "22001-0.0", '2443-0.0', "4079-0.0", "4080-0.0", "20116-0.0", "189-0.0",
                            "20117-0.0", "21000-0.0", '21001-0.0',
                            "21002-0.0", "21022-0.0", "22009-0.1",
                            "22009-0.2", "22009-0.3", "22009-0.4",
                            "22009-0.5", "22009-0.6", "22009-0.7",
                            "22009-0.8", "22009-0.9", "22009-0.10", "22000-0.0",
                            "30630-0.0", "30640-0.0", "30650-0.0", "30660-0.0", "30670-0.0",
                            "30680-0.0", "30690-0.0", "30700-0.0", "30710-0.0", "30720-0.0",       
                            "30730-0.0", "30740-0.0", "30750-0.0", "30770-0.0", "30780-0.0",
                            "30760-0.0", "30790-0.0", "30810-0.0", "30830-0.0", "30850-0.0",
                            "30860-0.0", "30870-0.0", "30880-0.0", "30890-0.0", "whr", "medi_dia_bp_chole")]
covar.select <- as.data.frame(covar.select)

selected_cat_traits <- c("22001-0.0", '2443-0.0', "20116-0.0", "20117-0.0", "21000-0.0", "medi_dia_bp_chole")

# Convert the categorical features to factor
covar.select[c(selected_cat_traits)] <- lapply(covar.select[c(selected_cat_traits)], factor)

colnames(covar.select) <- c("eid", "X220010.0", 'X24430.0', "X40790.0", "X40800.0", "X201160.0", "X1890.0",
                            "X201170.0", "X210000.0", "X210010.0",
                            "X210020.0", "X210220.0", "X220090.1",
                            "X220090.2", "X220090.3", "X220090.4",
                            "X220090.5", "X220090.6", "X220090.7",
                            "X220090.8", "X220090.9", "X220090.10", "X220000.0",
                            "X306300.0", "X306400.0", "X306500.0", "X306600.0", "X306700.0",
                            "X306800.0", "X306900.0", "X307000.0", "X307100.0", "X307200.0",       
                            "X307300.0", "X307400.0", "X307500.0", "X307700.0", "X307800.0",
                            "X307600.0", "X307900.0", "X308100.0", "X308300.0", "X308500.0",
                            "X308600.0", "X308700.0", "X308800.0", "X308900.0", "whr", "medi_dia_bp_chole")

message("Start running missRanger")
covar.select.imp <- missRanger(data = covar.select, formula = . ~ . - eid, returnOOB=T, pmm.k = 5)
message("Finished running missRanger")

save.image("~/Project/2022-08-17-individual_MR_with_imputation/dat/analysis/pheno_data/ukbb_covariate_imputed.RData")

#' ===================
#' Generate the X data
#' ===================

load("~/Project/2022-08-17-individual_MR_with_imputation/dat/analysis/pheno_data/ukbb_covariate_imputed.RData")

colnames(covar.select.imp) <- c("eid", "22001-0.0", '2443-0.0', "4079-0.0", "4080-0.0", "20116-0.0", "189-0.0",
                                "20117-0.0", "21000-0.0", '21001-0.0',
                                "21002-0.0", "21022-0.0", "22009-0.1",
                                "22009-0.2", "22009-0.3", "22009-0.4",
                                "22009-0.5", "22009-0.6", "22009-0.7",
                                "22009-0.8", "22009-0.9", "22009-0.10", "22000-0.0",
                                "30630-0.0", "30640-0.0", "30650-0.0", "30660-0.0", "30670-0.0",
                                "30680-0.0", "30690-0.0", "30700-0.0", "30710-0.0", "30720-0.0",
                                "30730-0.0", "30740-0.0", "30750-0.0", "30770-0.0", "30780-0.0",
                                "30760-0.0", "30790-0.0", "30810-0.0", "30830-0.0", "30850-0.0",
                                "30860-0.0", "30870-0.0", "30880-0.0", "30890-0.0", "whr", "medi_dia_bp_chole")
fwrite(covar.select.imp, file = "~/Project/2022-08-17-individual_MR_with_imputation/dat/analysis/pheno_data/ukbb.covariate.imputed.gz")

