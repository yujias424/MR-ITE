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

traits <- c("LDL", "TC", "HDL", "TG")
diseases <- c("CAD")
traits_name <- c("LDL-C", "Total Cholesterol", "HDL-C", "Triglycerides")
models <- c("model1b", "model3")

matrix.ivcf <- matrix(nrow = 16, ncol = 14) # previous 16 since we have IST.
matrix.cf <- matrix(nrow = 16, ncol = 14)
matrix.driv <- matrix(nrow = 16, ncol = 14)
n <- 1

for (t in 1:length(traits)){
  for (d in diseases){
    for (md in models){

      # continuous W  
      iv.tau.continuousW <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_", md, ".csv"))
      cf.tau.continuousW <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/cf_full_continuousW_pvalue_", md, ".csv"))
      driv.tau.continuousW <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/driv_full_continuousW_te_ul_", md, ".csv"))

      # cor.test(iv.tau.continuousW$tau, driv.tau.continuousW$point)
      # cor.test(iv.tau.continuousW$tau, cf.tau.continuousW$tau)

      # plot(iv.tau.continuousW$tau, cf.tau.continuousW$tau)

      # binary W
      iv.tau.binaryW <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_", md, ".csv"))
      cf.tau.binaryW <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_", md, ".csv"))
      driv.tau.binaryW <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_", md, ".csv"))

      # cor.test(iv.tau.binaryW$tau, driv.tau.binaryW$point)
      # cor.test(iv.tau.binaryW$tau, cf.tau.binaryW$tau)

      # plot(iv.tau.binaryW$tau, cf.tau.binaryW$tau)

      # Save continuous W results
      matrix.ivcf[n*2-1, 1:6] <- summary(iv.tau.continuousW$tau) 
      nominator <- ifelse(traits[t] == "HDL", length(which(iv.tau.continuousW$tau < 0)), length(which(iv.tau.continuousW$tau > 0)))
      matrix.ivcf[n*2-1, 7] <- nominator/length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n*2-1, 8] <- length(which(iv.tau.continuousW$tau.p.adjust < 0.1))/length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n*2-1, 9] <- length(which(iv.tau.continuousW$tau.pval < 0.1))/length(iv.tau.continuousW$tau.pval)
      matrix.ivcf[n*2-1, 10] <- "Continuous W"
      matrix.ivcf[n*2-1, 11] <- traits_name[t]
      matrix.ivcf[n*2-1, 12] <- d
      matrix.ivcf[n*2-1, 13] <- length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n*2-1, 14] <- md

      matrix.cf[n*2-1, 1:6] <- summary(cf.tau.continuousW$tau) 
      nominator <- ifelse(traits[t] == "HDL", length(which(cf.tau.continuousW$tau < 0)), length(which(cf.tau.continuousW$tau > 0)))
      matrix.cf[n*2-1, 7] <- nominator/length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n*2-1, 8] <- length(which(cf.tau.continuousW$tau.p.adjust < 0.1))/length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n*2-1, 9] <- length(which(cf.tau.continuousW$tau.pval < 0.1))/length(cf.tau.continuousW$tau.pval)
      matrix.cf[n*2-1, 10] <- "Continuous W"
      matrix.cf[n*2-1, 11] <- traits_name[t]
      matrix.cf[n*2-1, 12] <- d
      matrix.cf[n*2-1, 13] <- length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n*2-1, 14] <- md

      matrix.driv[n*2-1, 1:6] <- summary(driv.tau.continuousW$point) 
      nominator <- ifelse(traits[t] == "HDL", length(which(driv.tau.continuousW$point < 0)), length(which(driv.tau.continuousW$point > 0)))
      matrix.driv[n*2-1, 7] <- nominator/length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n*2-1, 8] <- length(which(driv.tau.continuousW$p_value_corrected < 0.1))/length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n*2-1, 9] <- length(which(driv.tau.continuousW$p_value < 0.1))/length(driv.tau.continuousW$p_value)
      matrix.driv[n*2-1, 10] <- "Continuous W"
      matrix.driv[n*2-1, 11] <- traits_name[t]
      matrix.driv[n*2-1, 12] <- d
      matrix.driv[n*2-1, 13] <- length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n*2-1, 14] <- md

      # Save binary W results
      matrix.ivcf[n*2, 1:6] <- summary(iv.tau.binaryW$tau) 
      matrix.ivcf[n*2, 7] <- length(which(iv.tau.binaryW$tau < 0))/length(iv.tau.binaryW$tau.p.adjust)
      matrix.ivcf[n*2, 8] <- length(which(iv.tau.binaryW$tau.p.adjust < 0.1))/length(iv.tau.binaryW$tau.p.adjust)
      matrix.ivcf[n*2, 9] <- length(which(iv.tau.binaryW$tau.pval < 0.1))/length(iv.tau.binaryW$tau.pval)
      matrix.ivcf[n*2, 10] <- "Binary W"
      matrix.ivcf[n*2, 11] <- traits_name[t]
      matrix.ivcf[n*2, 12] <- d
      matrix.ivcf[n*2, 13] <- length(iv.tau.binaryW$tau.p.adjust)
      matrix.ivcf[n*2, 14] <- md

      matrix.cf[n*2, 1:6] <- summary(cf.tau.binaryW$tau) 
      matrix.cf[n*2, 7] <- length(which(cf.tau.binaryW$tau < 0))/length(cf.tau.binaryW$tau.p.adjust)
      matrix.cf[n*2, 8] <- length(which(cf.tau.binaryW$tau.p.adjust < 0.1))/length(cf.tau.binaryW$tau.p.adjust)
      matrix.cf[n*2, 9] <- length(which(cf.tau.binaryW$tau.pval < 0.1))/length(cf.tau.binaryW$tau.pval)
      matrix.cf[n*2, 10] <- "Binary W"
      matrix.cf[n*2, 11] <- traits_name[t]
      matrix.cf[n*2, 12] <- d
      matrix.cf[n*2, 13] <- length(cf.tau.binaryW$tau.p.adjust)
      matrix.cf[n*2, 14] <- md

      matrix.driv[n*2, 1:6] <- summary(driv.tau.binaryW$point) 
      matrix.driv[n*2, 7] <- length(which(driv.tau.binaryW$point < 0))/length(driv.tau.binaryW$p_value_corrected)
      matrix.driv[n*2, 8] <- length(which(driv.tau.binaryW$p_value_corrected < 0.1))/length(driv.tau.binaryW$p_value_corrected)
      matrix.driv[n*2, 9] <- length(which(driv.tau.binaryW$p_value < 0.1))/length(driv.tau.binaryW$p_value)
      matrix.driv[n*2, 10] <- "Binary W"
      matrix.driv[n*2, 11] <- traits_name[t]
      matrix.driv[n*2, 12] <- d
      matrix.driv[n*2, 13] <- length(driv.tau.binaryW$p_value_corrected)
      matrix.driv[n*2, 14] <- md

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

fwrite(matrix.ivcf, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/05_summary_statistics/02_taus/method_IVCF_taus.csv")
fwrite(matrix.cf, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/05_summary_statistics/02_taus/method_CF_taus.csv")
fwrite(matrix.driv, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/05_summary_statistics/02_taus/method_DRIV_taus.csv")

matrix_final <- rbind(matrix.cf, matrix.driv, matrix.ivcf)
fwrite(matrix_final, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/05_summary_statistics/02_taus/method_taus.csv")
