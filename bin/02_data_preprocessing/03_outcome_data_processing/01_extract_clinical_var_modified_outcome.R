#' This is the code to extract clinical variable from UKBB phenotype dataset
#' the extreacted clinical variable will further be used as covariates in the further analysis
#' 
#' The difference between this file and the extract_clinical_var.R is that we include the "53-0.0", # Date of attending assessment centre in the selected covariates, as we need this data
#' to ensure that the disease is happened after the initial LDL assessment. We will subsequently use the outcome_corrected_assess_date_earlier_than_instance_date.R file to obtain the correct
#' outcome estimation.
#' 
#' @author: Yujia Shi
#' @date: 2022.09.01

suppressPackageStartupMessages({
    library(data.table)
    library(tidyverse)
})

message("Started.")

clinic.var <- fread("/mnt/data/share/UKBB/phenotype/ukb37268.csv.gz", nThread = 20)
james.var <- fread("/home/yujia/Project/2022-09-01-individual_MR/dat/03_outcome_data/COVID-19-vaccine-effects-history-new-hospitalization-covariates-whole-UKBB-eid-Aug-30.csv.gz") # history data

#' ================================================================================
#' obtain the dataframe containing both selected covariates and the outcome (CAD_hx) 
#' ================================================================================

# covariates
eid.selected <- c("eid", 
                    "31-0.0", # sex
                    "34-0.0", # year of birth
                    "48-0.0", # waist circumference (Initial)
                    "49-0.0", # hip circumference
                    "50-0.0", # standing height
                    "53-0.0", # Date of attending assessment centre #############################
                    "102-0.0", # pulse rate, automated reading
                    "102-0.1", # pulse rate, automated reading
                    "1160-0.0", # sleep duration
                    "1239-0.0", # current tobacco smoking
                    "1249-0.0", # past tobacco smoking
                    "129-0.0", # smoking/smokers in household
                    "189-0.0", # Townsend deprivation index at recruitment
                    "1558-0.0", # alcohol intake frequency
                    "2443-0.0", # diabetes diagnosed by doctor
                    "3456-0.0", # number of cigarettes currently smoked daily (current cigarette smokers)
                    "4079-0.0", # diastolic blood pressure, automated reading
                    "4079-0.1", # diastolic blood pressure, automated reading
                    "4080-0.0", # systolic blood pressure, automated reading
                    "4080-0.1", # systolic blood pressure, automated reading
                    "6153-0.0", # Medication for cholesterol, blood pressure, diabetes, or take exogenous hormones
                    "6177-0.0", # medication for cholesterol, blood pressure or diabetes
                    "12143-2.0", # weight (pre-imaging)
                    "12144-2.0", # height
                    "20116-0.0", # smoking status
                    "20117-0.0", # alcohol drinker status
                    "21000-0.0", # ethnic background
                    "22006-0.0", # Genetic ethnic grouping ###
                    "22020-0.0", # used in genetic principal component ###
                    "21001-0.0", # body mass index (BMI)
                    "21002-0.0", # weight
                    "21022-0.0", # age at recruitment
                    "22001-0.0", # Genetic sex
                    "54-0.0",    # UK Biobank assessment centre (For PRS calculation.)
                    "22000-0.0", # Genotype measurement batch (For PRS calculation.)
                    "22009-0.1", # pc1
                    "22009-0.2", # pc2
                    "22009-0.3", # pc3
                    "22009-0.4", # pc4
                    "22009-0.5", # pc5
                    "22009-0.6", # pc6
                    "22009-0.7", # pc7
                    "22009-0.8", # pc8
                    "22009-0.9", # pc9
                    "22009-0.10", # pc10
                    "22126-0.0", # doctor diagnosed hayfever or allergic rhinitis
                    "22127-0.0", # doctor diagnosed asthma
                    "22128-0.0", # doctor diagnosed emphysema
                    "22129-0.0", # doctor diagnosed chronic bronchitis
                    "22130-0.0", # doctor diagnosed COPD
                    "22131-0.0", # Doctor diagnosed cystic fibrosis
                    "22132-0.0", # Doctor diagnosed cystic fibrosis
                    "22133-0.0", # Doctor diagnosed sarcoidosis
                    "22134-0.0", # Doctor diagnosed bronchiectasis
                    "22135-0.0", # Doctor diagnosed idiopathic pulmonary fibrosis
                    "22136-0.0", # Doctor diagnosed fibrosing alveolitis/unspecified alveolitis
                    "22137-0.0", # Doctor diagnosed tuberculosis
                    "22138-0.0", # Doctor diagnosed silicosis
                    "22139-0.0", # Doctor diagnosed asbestosis
                    "23099-0.0", # body fat percentage
                    "30000-0.0", # white blood cell (leukocyte) count
                    "30010-0.0", # red blood cell (erythrocyte) count
                    "30020-0.0", # haemoglobin concentration
                    "30080-0.0", # platelet count
                    "30120-0.0", # lymphocyte count
                    "30130-0.0", # monocyte count
                    "30140-0.0", # neutrophill count
                    "30150-0.0", # eosinophill count
                    "30160-0.0", # basophill count
                    "30170-0.0", # nucleated red blood cell count
                    "30250-0.0", # reticulocyte count
                    "30630-0.0", # Apolipoprotein A
                    "30631-0.0", 
                    "30640-0.0", # Apolipoprotein B
                    "30641-0.0",
                    "30650-0.0", # Aspartate aminotransferase
                    "30651-0.0",
                    "30660-0.0", # Direct bilirubin
                    "30661-0.0",
                    "30670-0.0", # Uera
                    "30671-0.0",
                    "30680-0.0", # Calcium
                    "30681-0.0",
                    "30690-0.0", # Cholesterol
                    "30691-0.0",
                    "30700-0.0", # Creatinine
                    "30701-0.0",
                    "30710-0.0", # C-reactive protein
                    "30711-0.0",
                    "30720-0.0", # Cystatin C
                    "30721-0.0",
                    "30730-0.0", # Gamma glutamyltransferase
                    "30731-0.0",
                    "30740-0.0", # Glucose
                    "30741-0.0",
                    "30750-0.0", # Glycated haemoglobin (HbA1c)
                    "30751-0.0",
                    "30760-0.0", # HDL cholesterol
                    "30761-0.0",
                    "30770-0.0", # IGF-1
                    "30771-0.0",
                    "30780-0.0", # LDL direct
                    "30781-0.0", 
                    "30790-0.0", # Lipoprotein A
                    "30791-0.0",
                    "30810-0.0", # Phosphate
                    "30811-0.0",
                    "30830-0.0", # SHBG
                    "30831-0.0",
                    "30840-0.0", # Total bilirubin
                    "30841-0.0",
                    "30850-0.0", # Testosterone
                    "30851-0.0",
                    "30860-0.0", # Total protein
                    "30861-0.0",
                    "30870-0.0", # Triglycerides
                    "30871-0.0",
                    "30880-0.0", # Urate
                    "30881-0.0",
                    "30890-0.0", # Vitamin D
                    "30891-0.0") 

# select covariates                    
clinic.var.selected <- clinic.var[, ..eid.selected]
eid.set <- clinic.var.selected[, c('eid')]$`eid`
james.eid.set <- james.var[, c('eid')]$`eid`
n_occur <- data.frame(table(james.eid.set))
duplicated.id <- as.array(n_occur[n_occur$Freq > 1, ]$`james.eid.set`)

# remove duplicated id from both clinic.var.selected and the james.var
clinic.var.selected.subset <- clinic.var.selected[!(clinic.var.selected$eid %in% duplicated.id), ]
james.var.selected.subset <- james.var[!(james.var$eid %in% duplicated.id), ]

# select disease outcome based on james's dataset
selected.var <- c("eid", "CAD_hx", "Stroke_hx", "COPD_hx", "HTN_hx", "T2DM_hx", "VTE_hx", "AF_hx", "Heart_Failure_hx") # include more outcome   
james.var.selected.subset <- james.var.selected.subset[, ..selected.var]

# join two table to get final results
join.table <- clinic.var.selected.subset[james.var.selected.subset, on = .(eid = eid)]
dim(clinic.var.selected.subset)
dim(join.table)

# calculate the waist-to-hip ratio (WHR)
join.table[, whr := (join.table$`48-0.0`/join.table$`49-0.0`)]
join.table[, c("48-0.0", "49-0.0"):=NULL]

# # only select patients with Genetic ethnic grouping = 1 and included in PCA calculation
# join.table <- join.table[`22006-0.0` == 1]
# join.table <- join.table[`22020-0.0` == 1]

fwrite(join.table, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.gz")

message("Finished.")

#' =======================================================================
#' Extract cov and pheno to calculate the PRS and for further hte analysis
#' =======================================================================

clinic.var.jt <- join.table

# LDL
pheno.table <- clinic.var.jt[, c('eid', '30780-0.0', '30781-0.0')]
fwrite(pheno.table, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/LDL/ukbb.LDL")

# HDL
pheno.table <- clinic.var.jt[, c('eid', '30760-0.0', "30761-0.0")]
fwrite(pheno.table, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/HDL/ukbb.HDL")

# TC
pheno.table <- clinic.var.jt[, c('eid', '30690-0.0', "30691-0.0")]
fwrite(pheno.table, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TC/ukbb.TC")

# TG
pheno.table <- clinic.var.jt[, c('eid', '30870-0.0', '30871-0.0')]
fwrite(pheno.table, "/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TG/ukbb.TG")

####################################

# DBP (used as outcome)
outcome.table <- clinic.var.jt[, c('eid', '4079-0.0', "53-0.0")]
fwrite(outcome.table, "/home/yujia/Project/2023-07-20-individual_MR/dat/03_outcome_data/ukbb.DBP")

# SBP (used as outcome)
outcome.table <- clinic.var.jt[, c('eid', '4080-0.0', "53-0.0")]
fwrite(outcome.table, "/home/yujia/Project/2023-07-20-individual_MR/dat/03_outcome_data/ukbb.SBP")

#' ==========================================
#' Extreact Outcome
#' ==========================================

outcome_write <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/ukbb.covariate.gz")
stroke.all <- read.csv("~/Project/2022-09-01-individual_MR/dat/03_outcome_data/COVID-19-vaccine-effects-history-new-hospitalization-covariates-hemorrhage-cerebral-infraction-whole-UKBB-eid-June-27.csv.gz")
stroke.all.histroke <- stroke.all[, c("eid", "hemorrhage_stroke.x", "cerebral_infarction.x")]
outcome.res <- merge(outcome_write, stroke.all.histroke, by = "eid")
outcome.res.write <- outcome.res[, c("eid", "34-0.0", "53-0.0", "CAD_hx", "Stroke_hx", "COPD_hx", "HTN_hx", "T2DM_hx", "VTE_hx", "AF_hx", "Heart_Failure_hx", "hemorrhage_stroke.x", "cerebral_infarction.x",
                                    "30631-0.0", 
                                    "30641-0.0",
                                    "30651-0.0",
                                    "30661-0.0",
                                    "30671-0.0",
                                    "30681-0.0",
                                    "30691-0.0",
                                    "30701-0.0",
                                    "30711-0.0",
                                    "30721-0.0",
                                    "30731-0.0",
                                    "30741-0.0",
                                    "30751-0.0",
                                    "30761-0.0",
                                    "30771-0.0",
                                    "30781-0.0", 
                                    "30791-0.0",
                                    "30811-0.0",
                                    "30831-0.0",
                                    "30841-0.0",
                                    "30851-0.0",
                                    "30861-0.0",
                                    "30871-0.0",
                                    "30881-0.0",
                                    "30891-0.0")]
outcome.res.write[outcome.res.write == ""] <- NA
setnames(outcome.res.write, c("eid", "34-0.0", "53-0.0", "CAD_hx", "Stroke_hx", "COPD_hx", "HTN_hx", "T2DM_hx", "VTE_hx", "AF_hx", "Heart_Failure_hx", "hemorrhage_stroke.x", "cerebral_infarction.x"), 
                            c("IID", "birth_date", "included_date", "cad", "stroke", "copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "ischemic_stroke"))
fwrite(outcome.res.write, "/home/yujia/Project/2023-07-20-individual_MR/dat/03_outcome_data/outcome_with_date", sep = "\t")

message("Finished the analysis.")