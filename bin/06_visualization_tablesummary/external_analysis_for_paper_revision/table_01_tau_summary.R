#' This code is to plot the histogram of the tau for different traits
#' 
#' @author Shi Yujia
#' @date 2023.02.22

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(latex2exp)
    library(patchwork)
})

traits <- c("IGF1", "CRP")
diseases <- c("CAD")
traits_name <- c("IGF-1", "C-Reactive Protein")
models <- c("model1b", "model3")

matrix.ivcf <- matrix(nrow = 4, ncol = 14) # previous 16 since we have IST.
matrix.cf <- matrix(nrow = 4, ncol = 14)
matrix.driv <- matrix(nrow = 4, ncol = 14)
n <- 1

for (t in 1:length(traits)){
  for (d in diseases){
    for (md in models){

      # continuous W  
      iv.tau.continuousW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_", md, ".csv"))
      cf.tau.continuousW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/cf_full_continuousW_pvalue_", md, ".csv"))
      driv.tau.continuousW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/driv_full_continuousW_te_ul_", md, ".csv"))

      # cor.test(iv.tau.continuousW$tau, driv.tau.continuousW$point)
      # cor.test(iv.tau.continuousW$tau, cf.tau.continuousW$tau)

      # plot(iv.tau.continuousW$tau, cf.tau.continuousW$tau)

      # # binary W
      # iv.tau.binaryW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_", md, ".csv"))
      # cf.tau.binaryW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_", md, ".csv"))
      # driv.tau.binaryW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_", md, ".csv"))

      # cor.test(iv.tau.binaryW$tau, driv.tau.binaryW$point)
      # cor.test(iv.tau.binaryW$tau, cf.tau.binaryW$tau)

      # plot(iv.tau.binaryW$tau, cf.tau.binaryW$tau)

      # Save continuous W results
      matrix.ivcf[n, 1:6] <- summary(iv.tau.continuousW$tau) 
      nominator <- ifelse(traits[t] == "HDL", length(which(iv.tau.continuousW$tau < 0)), length(which(iv.tau.continuousW$tau > 0)))
      matrix.ivcf[n, 7] <- nominator/length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n, 8] <- length(which(iv.tau.continuousW$tau.p.adjust < 0.1))/length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n, 9] <- length(which(iv.tau.continuousW$tau.pval < 0.1))/length(iv.tau.continuousW$tau.pval)
      matrix.ivcf[n, 10] <- "Continuous W"
      matrix.ivcf[n, 11] <- traits_name[t]
      matrix.ivcf[n, 12] <- d
      matrix.ivcf[n, 13] <- length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n, 14] <- md

      matrix.cf[n, 1:6] <- summary(cf.tau.continuousW$tau) 
      nominator <- ifelse(traits[t] == "HDL", length(which(cf.tau.continuousW$tau < 0)), length(which(cf.tau.continuousW$tau > 0)))
      matrix.cf[n, 7] <- nominator/length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n, 8] <- length(which(cf.tau.continuousW$tau.p.adjust < 0.1))/length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n, 9] <- length(which(cf.tau.continuousW$tau.pval < 0.1))/length(cf.tau.continuousW$tau.pval)
      matrix.cf[n, 10] <- "Continuous W"
      matrix.cf[n, 11] <- traits_name[t]
      matrix.cf[n, 12] <- d
      matrix.cf[n, 13] <- length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n, 14] <- md

      matrix.driv[n, 1:6] <- summary(driv.tau.continuousW$point) 
      nominator <- ifelse(traits[t] == "HDL", length(which(driv.tau.continuousW$point < 0)), length(which(driv.tau.continuousW$point > 0)))
      matrix.driv[n, 7] <- nominator/length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n, 8] <- length(which(driv.tau.continuousW$p_value_corrected < 0.1))/length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n, 9] <- length(which(driv.tau.continuousW$p_value < 0.1))/length(driv.tau.continuousW$p_value)
      matrix.driv[n, 10] <- "Continuous W"
      matrix.driv[n, 11] <- traits_name[t]
      matrix.driv[n, 12] <- d
      matrix.driv[n, 13] <- length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n, 14] <- md

      n <- n+1

    }
  }
}

colnames(matrix.ivcf) <- c("Min", "1st Qu", "Median", "Mean", "3rd Qu", "Max", 
                           "Positive Tau Percentage", "Significant Tau Percentage (P.val.adjusted)", "Significant Tau Percentage (P.val)", 
                           "Treatment Type", "Treatment", "Disease", "Sample Size", "Model")
colnames(matrix.driv) <- c("Min", "1st Qu", "Median", "Mean", "3rd Qu", "Max", 
                           "Positive Tau Percentage", "Significant Tau Percentage (P.val.adjusted)", "Significant Tau Percentage (P.val)",
                           "Treatment Type", "Treatment", "Disease", "Sample Size", "Model")
colnames(matrix.cf) <- c("Min", "1st Qu", "Median", "Mean", "3rd Qu", "Max", 
                         "Positive Tau Percentage", "Significant Tau Percentage (P.val.adjusted)", "Significant Tau Percentage (P.val)",
                         "Treatment Type", "Treatment", "Disease", "Sample Size", "Model")

matrix.ivcf <- as.data.frame(matrix.ivcf)
matrix.driv <- as.data.frame(matrix.driv)
matrix.cf <- as.data.frame(matrix.cf)

matrix.ivcf <- matrix.ivcf[, c("Treatment Type", "Treatment", "Disease", "Model", "Sample Size",
                               "Positive Tau Percentage", "Significant Tau Percentage (P.val)", 
                               "Min", "1st Qu", "Median", "3rd Qu", "Max", "Mean", 
                               "Significant Tau Percentage (P.val.adjusted)")]
matrix.driv <- matrix.driv[, c("Treatment Type", "Treatment", "Disease", "Model", "Sample Size",
                               "Positive Tau Percentage", "Significant Tau Percentage (P.val)", 
                               "Min", "1st Qu", "Median", "3rd Qu", "Max", "Mean", 
                               "Significant Tau Percentage (P.val.adjusted)")]
matrix.cf <- matrix.cf[, c("Treatment Type", "Treatment", "Disease", "Model", "Sample Size",
                               "Positive Tau Percentage", "Significant Tau Percentage (P.val)", 
                               "Min", "1st Qu", "Median", "3rd Qu", "Max", "Mean", 
                               "Significant Tau Percentage (P.val.adjusted)")]

matrix.cf$Method <- "CF"
matrix.driv$Method <- "DRIV"
matrix.ivcf$Method <- "IVCF"

fwrite(matrix.ivcf, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/02_summary_statistics/01_taus/method_IVCF_taus.csv")
fwrite(matrix.cf, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/02_summary_statistics/01_taus/method_CF_taus.csv")
fwrite(matrix.driv, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/02_summary_statistics/01_taus/method_DRIV_taus.csv")

matrix_final <- rbind(matrix.cf, matrix.driv, matrix.ivcf)
fwrite(matrix_final, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/02_summary_statistics/01_taus/method_taus.csv")
