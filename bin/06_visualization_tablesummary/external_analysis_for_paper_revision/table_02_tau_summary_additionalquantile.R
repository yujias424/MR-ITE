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

matrix.ivcf <- matrix(nrow = 4, ncol = 26) # previous 16 since we have IST.
matrix.cf <- matrix(nrow = 4, ncol = 26)
matrix.driv <- matrix(nrow = 4, ncol = 26)
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
      matrix.ivcf[n, 1:9] <- as.numeric(quantile(iv.tau.continuousW$tau, c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1))) 
      matrix.ivcf[n, 10] <- as.numeric(matrix.ivcf[n, 6]) - as.numeric(matrix.ivcf[n, 4]) # IQR
      matrix.ivcf[n, 11] <- as.numeric(matrix.ivcf[n, 7]) - as.numeric(matrix.ivcf[n, 3]) # 90% - 10%
      matrix.ivcf[n, 12] <- as.numeric(matrix.ivcf[n, 8]) - as.numeric(matrix.ivcf[n, 2]) # 95% - 5%
      matrix.ivcf[n, 13] <- summary(iv.tau.continuousW$tau)[4] 
      nominator <- ifelse(traits[t] == "HDL", length(which(iv.tau.continuousW$tau < 0)), length(which(iv.tau.continuousW$tau > 0)))
      matrix.ivcf[n, 14] <- nominator/length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n, 15] <- length(which(iv.tau.continuousW$tau.p.adjust < 0.1))/length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n, 16] <- length(which(iv.tau.continuousW$tau.pval < 0.1))/length(iv.tau.continuousW$tau.pval)
      matrix.ivcf[n, 17] <- "Continuous W"
      matrix.ivcf[n, 18] <- traits_name[t]
      matrix.ivcf[n, 19] <- d
      matrix.ivcf[n, 20] <- length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n, 21] <- md
      matrix.ivcf[n, 22] <- length(which(iv.tau.continuousW$tau.p.adjust < 0.05))/length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n, 23] <- length(which(iv.tau.continuousW$tau.pval < 0.05))/length(iv.tau.continuousW$tau.pval)
      matrix.ivcf[n, 24] <- abs((as.numeric(matrix.ivcf[n, 8]) - as.numeric(matrix.ivcf[n, 2]))/as.numeric(matrix.ivcf[n, 2]))
      matrix.ivcf[n, 25] <- abs((as.numeric(matrix.ivcf[n, 7]) - as.numeric(matrix.ivcf[n, 3]))/as.numeric(matrix.ivcf[n, 3]))
      matrix.ivcf[n, 26] <- abs((as.numeric(matrix.ivcf[n, 6]) - as.numeric(matrix.ivcf[n, 4]))/as.numeric(matrix.ivcf[n, 4]))

      matrix.cf[n, 1:9] <- as.numeric(quantile(cf.tau.continuousW$tau, c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1))) 
      matrix.cf[n, 10] <- as.numeric(matrix.cf[n, 6]) - as.numeric(matrix.cf[n, 4]) # IQR
      matrix.cf[n, 11] <- as.numeric(matrix.cf[n, 7]) - as.numeric(matrix.cf[n, 3]) # 90% - 10%
      matrix.cf[n, 12] <- as.numeric(matrix.cf[n, 8]) - as.numeric(matrix.cf[n, 2]) # 95% - 5%
      matrix.cf[n, 13] <- summary(cf.tau.continuousW$tau)[4] 
      nominator <- ifelse(traits[t] == "HDL", length(which(cf.tau.continuousW$tau < 0)), length(which(cf.tau.continuousW$tau > 0)))
      matrix.cf[n, 14] <- nominator/length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n, 15] <- length(which(cf.tau.continuousW$tau.p.adjust < 0.1))/length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n, 16] <- length(which(cf.tau.continuousW$tau.pval < 0.1))/length(cf.tau.continuousW$tau.pval)
      matrix.cf[n, 17] <- "Continuous W"
      matrix.cf[n, 18] <- traits_name[t]
      matrix.cf[n, 19] <- d
      matrix.cf[n, 20] <- length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n, 21] <- md
      matrix.cf[n, 22] <- length(which(cf.tau.continuousW$tau.p.adjust < 0.05))/length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n, 23] <- length(which(cf.tau.continuousW$tau.pval < 0.05))/length(cf.tau.continuousW$tau.pval)
      matrix.cf[n, 24] <- abs((as.numeric(matrix.cf[n, 8]) - as.numeric(matrix.cf[n, 2]))/as.numeric(matrix.cf[n, 2]))
      matrix.cf[n, 25] <- abs((as.numeric(matrix.cf[n, 7]) - as.numeric(matrix.cf[n, 3]))/as.numeric(matrix.cf[n, 3]))
      matrix.cf[n, 26] <- abs((as.numeric(matrix.cf[n, 6]) - as.numeric(matrix.cf[n, 4]))/as.numeric(matrix.cf[n, 4]))

      matrix.driv[n, 1:9] <- as.numeric(quantile(driv.tau.continuousW$point, c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1)))
      matrix.driv[n, 10] <- as.numeric(matrix.driv[n, 6]) - as.numeric(matrix.driv[n, 4]) # IQR
      matrix.driv[n, 11] <- as.numeric(matrix.driv[n, 7]) - as.numeric(matrix.driv[n, 3]) # 90% - 10%
      matrix.driv[n, 12] <- as.numeric(matrix.driv[n, 8]) - as.numeric(matrix.driv[n, 2]) # 95% - 5%
      matrix.driv[n, 13] <- summary(driv.tau.continuousW$point)[4] 
      nominator <- ifelse(traits[t] == "HDL", length(which(driv.tau.continuousW$point < 0)), length(which(driv.tau.continuousW$point > 0)))
      matrix.driv[n, 14] <- nominator/length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n, 15] <- length(which(driv.tau.continuousW$p_value_corrected < 0.1))/length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n, 16] <- length(which(driv.tau.continuousW$p_value < 0.1))/length(driv.tau.continuousW$p_value)
      matrix.driv[n, 17] <- "Continuous W"
      matrix.driv[n, 18] <- traits_name[t]
      matrix.driv[n, 19] <- d
      matrix.driv[n, 20] <- length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n, 21] <- md
      matrix.driv[n, 22] <- length(which(driv.tau.continuousW$p_value_corrected < 0.05))/length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n, 23] <- length(which(driv.tau.continuousW$p_value < 0.05))/length(driv.tau.continuousW$p_value)
      matrix.driv[n, 24] <- abs((as.numeric(matrix.driv[n, 8]) - as.numeric(matrix.driv[n, 2]))/as.numeric(matrix.driv[n, 2]))
      matrix.driv[n, 25] <- abs((as.numeric(matrix.driv[n, 7]) - as.numeric(matrix.driv[n, 3]))/as.numeric(matrix.driv[n, 3]))
      matrix.driv[n, 26] <- abs((as.numeric(matrix.driv[n, 6]) - as.numeric(matrix.driv[n, 4]))/as.numeric(matrix.driv[n, 4]))

      n <- n+1

    }
  }
}

colnames(matrix.ivcf) <- c("Min", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "Max", "50% Coverage (25% - 75%)", "80% Coverage (10% - 90%)", "90% Coverage (5% - 95%)", "Mean",
                           "Positive Tau Percentage", "Significant Tau Percentage (P.val.adjusted < 0.1)", "Significant Tau Percentage (P.val < 0.1)",
                           "Treatment Type", "Treatment", "Disease", "Sample Size", "Model", 
                           "Significant Tau Percentage (P.val.adjusted < 0.05)", "Significant Tau Percentage (P.val < 0.05)", "absolute % difference (95th vs 5th percentile)", "absolute % difference (90th vs 10th percentile)", "absolute % difference (75th vs 25th percentile)")
colnames(matrix.driv) <- c("Min", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "Max", "50% Coverage (25% - 75%)", "80% Coverage (10% - 90%)", "90% Coverage (5% - 95%)", "Mean",
                           "Positive Tau Percentage", "Significant Tau Percentage (P.val.adjusted < 0.1)", "Significant Tau Percentage (P.val < 0.1)",
                           "Treatment Type", "Treatment", "Disease", "Sample Size", "Model",
                           "Significant Tau Percentage (P.val.adjusted < 0.05)", "Significant Tau Percentage (P.val < 0.05)", "absolute % difference (95th vs 5th percentile)", "absolute % difference (90th vs 10th percentile)", "absolute % difference (75th vs 25th percentile)")
colnames(matrix.cf) <- c("Min", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "Max", "50% Coverage (25% - 75%)", "80% Coverage (10% - 90%)", "90% Coverage (5% - 95%)", "Mean",
                         "Positive Tau Percentage", "Significant Tau Percentage (P.val.adjusted < 0.1)", "Significant Tau Percentage (P.val < 0.1)",
                         "Treatment Type", "Treatment", "Disease", "Sample Size", "Model",
                         "Significant Tau Percentage (P.val.adjusted < 0.05)", "Significant Tau Percentage (P.val < 0.05)", "absolute % difference (95th vs 5th percentile)", "absolute % difference (90th vs 10th percentile)", "absolute % difference (75th vs 25th percentile)")

matrix.ivcf <- as.data.frame(matrix.ivcf)
matrix.driv <- as.data.frame(matrix.driv)
matrix.cf <- as.data.frame(matrix.cf)

matrix.ivcf <- matrix.ivcf[, c("Treatment Type", "Treatment", "Disease", "Model", "Sample Size",
                               "Positive Tau Percentage", "Significant Tau Percentage (P.val < 0.1)", "Significant Tau Percentage (P.val < 0.05)", 
                               "Min", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "Max", "50% Coverage (25% - 75%)", "80% Coverage (10% - 90%)", "90% Coverage (5% - 95%)", "Mean",
                               "Significant Tau Percentage (P.val.adjusted < 0.1)", "Significant Tau Percentage (P.val.adjusted < 0.05)", "absolute % difference (95th vs 5th percentile)", "absolute % difference (90th vs 10th percentile)", "absolute % difference (75th vs 25th percentile)")]
matrix.driv <- matrix.driv[, c("Treatment Type", "Treatment", "Disease", "Model", "Sample Size",
                               "Positive Tau Percentage", "Significant Tau Percentage (P.val < 0.1)", "Significant Tau Percentage (P.val < 0.05)", 
                               "Min", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "Max", "50% Coverage (25% - 75%)", "80% Coverage (10% - 90%)", "90% Coverage (5% - 95%)", "Mean",
                               "Significant Tau Percentage (P.val.adjusted < 0.1)", "Significant Tau Percentage (P.val.adjusted < 0.05)", "absolute % difference (95th vs 5th percentile)", "absolute % difference (90th vs 10th percentile)", "absolute % difference (75th vs 25th percentile)")]
matrix.cf <- matrix.cf[, c("Treatment Type", "Treatment", "Disease", "Model", "Sample Size",
                            "Positive Tau Percentage", "Significant Tau Percentage (P.val < 0.1)", "Significant Tau Percentage (P.val < 0.05)", 
                            "Min", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "Max", "50% Coverage (25% - 75%)", "80% Coverage (10% - 90%)", "90% Coverage (5% - 95%)", "Mean",
                            "Significant Tau Percentage (P.val.adjusted < 0.1)", "Significant Tau Percentage (P.val.adjusted < 0.05)", "absolute % difference (95th vs 5th percentile)", "absolute % difference (90th vs 10th percentile)", "absolute % difference (75th vs 25th percentile)")]

matrix.cf$Method <- "CF"
matrix.driv$Method <- "DRIV"
matrix.ivcf$Method <- "IVCF"

fwrite(matrix.ivcf, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/02_summary_statistics/01_taus/method_IVCF_taus_additionquantile_2.csv")
fwrite(matrix.cf, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/02_summary_statistics/01_taus/method_CF_taus_additionquantile_2.csv")
fwrite(matrix.driv, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/02_summary_statistics/01_taus/method_DRIV_taus_additionquantile_2.csv")

matrix_final <- rbind(matrix.cf, matrix.driv, matrix.ivcf)
fwrite(matrix_final, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/02_summary_statistics/01_taus/method_taus_additionquantile_2.csv")
