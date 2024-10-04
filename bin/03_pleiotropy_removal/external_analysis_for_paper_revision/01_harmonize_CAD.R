#' This code is to using twosampleMR to harmonize the data (including outcome and exposure) for further pleiotropy removal
#' 
#' @author Yujia Shi
#' @date 2022.09.01

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

# read the outcome summary data
outcome <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/QC_SAS/CAD.QC.gz"), sep = "\t")
outcome <- as.data.frame(outcome)
outcome <- outcome[, c("SNP", "A1", "A2", "N", "SE", "P", "BETA", "EAF")]
colnames(outcome) <- c("SNP", "effect_allele", "other_allele", "samplesize", 
                          "se", "pval", "beta", "eaf")
outcome$Phenotype <- rep("CAD", dim(outcome)[1])
outcome <- outcome[, c("Phenotype", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf",
                         "pval", "samplesize")]

fwrite(outcome, paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/05_two_sample_mr/external_analysis_for_paper_revision/CAD/outcome.CAD.SAS.txt"),  sep = " ", quote = F)

# read the outcome summary data
outcome <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/QC/CAD.QC.gz"), sep = "\t")
outcome <- as.data.frame(outcome)
outcome <- outcome[, c("SNP", "A1", "A2", "SE", "P", "BETA", "EAF")]
colnames(outcome) <- c("SNP", "effect_allele", "other_allele", 
                          "se", "pval", "beta", "eaf")
outcome$Phenotype <- rep("CAD", dim(outcome)[1])
outcome <- outcome[, c("Phenotype", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf",
                         "pval")]

fwrite(outcome, paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/05_two_sample_mr/external_analysis_for_paper_revision/CAD/outcome.CAD.txt"),  sep = " ", quote = F)