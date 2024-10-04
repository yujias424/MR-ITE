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

matrix.ivcf <- matrix(nrow = 16, ncol = 25) # previous 16 since we have IST.
matrix.cf <- matrix(nrow = 16, ncol = 25)
matrix.driv <- matrix(nrow = 16, ncol = 25)
n <- 1

for (t in 1:length(traits)){
  for (d in diseases){
    for (md in models){

      # continuous W  
      iv.tau.continuousW <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_", md, ".csv"))
      cf.tau.continuousW <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/cf_full_continuousW_pvalue_", md, ".csv"))
      driv.tau.continuousW <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/driv_full_continuousW_te_ul_", md, ".csv"))

      # prob of CAD
      prob_cad <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date3_james2022.gz")
      prob_cad <- prob_cad[prob_cad$IID %in% driv.tau.continuousW$IID, ]
      prob_cad <- table(prob_cad$CAD)[2]/(table(prob_cad$CAD)[1]+table(prob_cad$CAD)[2])

      iv.tau.continuousW$tau <- (iv.tau.continuousW$tau*10)/prob_cad
      cf.tau.continuousW$tau <- (cf.tau.continuousW$tau*10)/prob_cad
      driv.tau.continuousW$point <- (driv.tau.continuousW$point*10)/prob_cad

      # binary W
      iv.tau.binaryW <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_", md, ".csv"))
      cf.tau.binaryW <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_", md, ".csv"))
      driv.tau.binaryW <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_", md, ".csv"))

      iv.tau.binaryW$tau <- iv.tau.binaryW$tau/prob_cad
      cf.tau.binaryW$tau <- cf.tau.binaryW$tau/prob_cad
      driv.tau.binaryW$point <- driv.tau.binaryW$point/prob_cad

      # Save continuous W results
      matrix.ivcf[n*2-1, 1:9] <- as.numeric(quantile(iv.tau.continuousW$tau, c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1))) 
      matrix.ivcf[n*2-1, 10] <- as.numeric(matrix.ivcf[n*2-1, 6]) - as.numeric(matrix.ivcf[n*2-1, 4]) # IQR
      matrix.ivcf[n*2-1, 11] <- as.numeric(matrix.ivcf[n*2-1, 7]) - as.numeric(matrix.ivcf[n*2-1, 3]) # 90% - 10%
      matrix.ivcf[n*2-1, 12] <- as.numeric(matrix.ivcf[n*2-1, 8]) - as.numeric(matrix.ivcf[n*2-1, 2]) # 95% - 5%
      matrix.ivcf[n*2-1, 13] <- summary(iv.tau.continuousW$tau)[4] 
      nominator <- ifelse(traits[t] == "HDL", length(which(iv.tau.continuousW$tau < 0)), length(which(iv.tau.continuousW$tau > 0)))
      matrix.ivcf[n*2-1, 14] <- nominator/length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n*2-1, 15] <- length(which(iv.tau.continuousW$tau.p.adjust < 0.1))/length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n*2-1, 16] <- length(which(iv.tau.continuousW$tau.pval < 0.1))/length(iv.tau.continuousW$tau.pval)
      matrix.ivcf[n*2-1, 17] <- length(which(iv.tau.continuousW$tau.p.adjust < 0.05))/length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n*2-1, 18] <- length(which(iv.tau.continuousW$tau.pval < 0.05))/length(iv.tau.continuousW$tau.pval)
      matrix.ivcf[n*2-1, 19] <- "Continuous W"
      matrix.ivcf[n*2-1, 20] <- traits_name[t]
      matrix.ivcf[n*2-1, 21] <- d
      matrix.ivcf[n*2-1, 22] <- length(iv.tau.continuousW$tau.p.adjust)
      matrix.ivcf[n*2-1, 23] <- md
      matrix.ivcf[n*2-1, 24] <- abs((as.numeric(matrix.ivcf[n*2-1, 8]) - as.numeric(matrix.ivcf[n*2-1, 2]))/as.numeric(matrix.ivcf[n*2-1, 2]))
      matrix.ivcf[n*2-1, 25] <- abs((as.numeric(matrix.ivcf[n*2-1, 7]) - as.numeric(matrix.ivcf[n*2-1, 3]))/as.numeric(matrix.ivcf[n*2-1, 3]))

      matrix.cf[n*2-1, 1:9] <- as.numeric(quantile(cf.tau.continuousW$tau, c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1))) 
      matrix.cf[n*2-1, 10] <- as.numeric(matrix.cf[n*2-1, 6]) - as.numeric(matrix.cf[n*2-1, 4]) # IQR
      matrix.cf[n*2-1, 11] <- as.numeric(matrix.cf[n*2-1, 7]) - as.numeric(matrix.cf[n*2-1, 3]) # 90% - 10%
      matrix.cf[n*2-1, 12] <- as.numeric(matrix.cf[n*2-1, 8]) - as.numeric(matrix.cf[n*2-1, 2]) # 95% - 5%
      matrix.cf[n*2-1, 13] <- summary(cf.tau.continuousW$tau)[4] 
      nominator <- ifelse(traits[t] == "HDL", length(which(cf.tau.continuousW$tau < 0)), length(which(cf.tau.continuousW$tau > 0)))
      matrix.cf[n*2-1, 14] <- nominator/length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n*2-1, 15] <- length(which(cf.tau.continuousW$tau.p.adjust < 0.1))/length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n*2-1, 16] <- length(which(cf.tau.continuousW$tau.pval < 0.1))/length(cf.tau.continuousW$tau.pval)
      matrix.cf[n*2-1, 17] <- length(which(cf.tau.continuousW$tau.p.adjust < 0.05))/length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n*2-1, 18] <- length(which(cf.tau.continuousW$tau.pval < 0.05))/length(cf.tau.continuousW$tau.pval)
      matrix.cf[n*2-1, 19] <- "Continuous W"
      matrix.cf[n*2-1, 20] <- traits_name[t]
      matrix.cf[n*2-1, 21] <- d
      matrix.cf[n*2-1, 22] <- length(cf.tau.continuousW$tau.p.adjust)
      matrix.cf[n*2-1, 23] <- md
      matrix.cf[n*2-1, 24] <- abs((as.numeric(matrix.cf[n*2-1, 8]) - as.numeric(matrix.cf[n*2-1, 2]))/as.numeric(matrix.cf[n*2-1, 2]))
      matrix.cf[n*2-1, 25] <- abs((as.numeric(matrix.cf[n*2-1, 7]) - as.numeric(matrix.cf[n*2-1, 3]))/as.numeric(matrix.cf[n*2-1, 3]))

      matrix.driv[n*2-1, 1:9] <- as.numeric(quantile(driv.tau.continuousW$point, c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1)))
      matrix.driv[n*2-1, 10] <- as.numeric(matrix.driv[n*2-1, 6]) - as.numeric(matrix.driv[n*2-1, 4]) # IQR
      matrix.driv[n*2-1, 11] <- as.numeric(matrix.driv[n*2-1, 7]) - as.numeric(matrix.driv[n*2-1, 3]) # 90% - 10%
      matrix.driv[n*2-1, 12] <- as.numeric(matrix.driv[n*2-1, 8]) - as.numeric(matrix.driv[n*2-1, 2]) # 95% - 5%
      matrix.driv[n*2-1, 13] <- summary(driv.tau.continuousW$point)[4] 
      nominator <- ifelse(traits[t] == "HDL", length(which(driv.tau.continuousW$point < 0)), length(which(driv.tau.continuousW$point > 0)))
      matrix.driv[n*2-1, 14] <- nominator/length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n*2-1, 15] <- length(which(driv.tau.continuousW$p_value_corrected < 0.1))/length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n*2-1, 16] <- length(which(driv.tau.continuousW$p_value < 0.1))/length(driv.tau.continuousW$p_value)
      matrix.driv[n*2-1, 17] <- length(which(driv.tau.continuousW$p_value_corrected < 0.05))/length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n*2-1, 18] <- length(which(driv.tau.continuousW$p_value < 0.05))/length(driv.tau.continuousW$p_value)
      matrix.driv[n*2-1, 19] <- "Continuous W"
      matrix.driv[n*2-1, 20] <- traits_name[t]
      matrix.driv[n*2-1, 21] <- d
      matrix.driv[n*2-1, 22] <- length(driv.tau.continuousW$p_value_corrected)
      matrix.driv[n*2-1, 23] <- md
      matrix.driv[n*2-1, 24] <- abs((as.numeric(matrix.driv[n*2-1, 8]) - as.numeric(matrix.driv[n*2-1, 2]))/as.numeric(matrix.driv[n*2-1, 2]))
      matrix.driv[n*2-1, 25] <- abs((as.numeric(matrix.driv[n*2-1, 7]) - as.numeric(matrix.driv[n*2-1, 3]))/as.numeric(matrix.driv[n*2-1, 3]))

      # Save binary W results
      matrix.ivcf[n*2, 1:9] <- as.numeric(quantile(iv.tau.binaryW$tau, c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1))) 
      matrix.ivcf[n*2, 10] <- as.numeric(matrix.ivcf[n*2, 6]) - as.numeric(matrix.ivcf[n*2, 4]) # IQR
      matrix.ivcf[n*2, 11] <- as.numeric(matrix.ivcf[n*2, 7]) - as.numeric(matrix.ivcf[n*2, 3]) # 90% - 10%
      matrix.ivcf[n*2, 12] <- as.numeric(matrix.ivcf[n*2, 8]) - as.numeric(matrix.ivcf[n*2, 2]) # 95% - 5%
      matrix.ivcf[n*2, 13] <- summary(iv.tau.binaryW$tau)[4] 
      matrix.ivcf[n*2, 14] <- length(which(iv.tau.binaryW$tau < 0))/length(iv.tau.binaryW$tau.p.adjust)
      matrix.ivcf[n*2, 15] <- length(which(iv.tau.binaryW$tau.p.adjust < 0.1))/length(iv.tau.binaryW$tau.p.adjust)
      matrix.ivcf[n*2, 16] <- length(which(iv.tau.binaryW$tau.pval < 0.1))/length(iv.tau.binaryW$tau.pval)
      matrix.ivcf[n*2, 17] <- length(which(iv.tau.binaryW$tau.p.adjust < 0.05))/length(iv.tau.binaryW$tau.p.adjust)
      matrix.ivcf[n*2, 18] <- length(which(iv.tau.binaryW$tau.pval < 0.05))/length(iv.tau.binaryW$tau.pval)
      matrix.ivcf[n*2, 19] <- "Binary W"
      matrix.ivcf[n*2, 20] <- traits_name[t]
      matrix.ivcf[n*2, 21] <- d
      matrix.ivcf[n*2, 22] <- length(iv.tau.binaryW$tau.p.adjust)
      matrix.ivcf[n*2, 23] <- md
      matrix.ivcf[n*2, 24] <- abs((as.numeric(matrix.ivcf[n*2, 8]) - as.numeric(matrix.ivcf[n*2, 2]))/as.numeric(matrix.ivcf[n*2, 2]))
      matrix.ivcf[n*2, 25] <- abs((as.numeric(matrix.ivcf[n*2, 7]) - as.numeric(matrix.ivcf[n*2, 3]))/as.numeric(matrix.ivcf[n*2, 3]))

      matrix.cf[n*2, 1:9] <- as.numeric(quantile(cf.tau.binaryW$tau, c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1))) 
      matrix.cf[n*2, 10] <- as.numeric(matrix.cf[n*2, 6]) - as.numeric(matrix.cf[n*2, 4]) # IQR
      matrix.cf[n*2, 11] <- as.numeric(matrix.cf[n*2, 7]) - as.numeric(matrix.cf[n*2, 3]) # 90% - 10%
      matrix.cf[n*2, 12] <- as.numeric(matrix.cf[n*2, 8]) - as.numeric(matrix.cf[n*2, 2]) # 95% - 5%
      matrix.cf[n*2, 13] <- summary(cf.tau.binaryW$tau)[4] 
      matrix.cf[n*2, 14] <- length(which(cf.tau.binaryW$tau < 0))/length(cf.tau.binaryW$tau.p.adjust)
      matrix.cf[n*2, 15] <- length(which(cf.tau.binaryW$tau.p.adjust < 0.1))/length(cf.tau.binaryW$tau.p.adjust)
      matrix.cf[n*2, 16] <- length(which(cf.tau.binaryW$tau.pval < 0.1))/length(cf.tau.binaryW$tau.pval)
      matrix.cf[n*2, 17] <- length(which(cf.tau.binaryW$tau.p.adjust < 0.05))/length(cf.tau.binaryW$tau.p.adjust)
      matrix.cf[n*2, 18] <- length(which(cf.tau.binaryW$tau.pval < 0.05))/length(cf.tau.binaryW$tau.pval)
      matrix.cf[n*2, 19] <- "Binary W"
      matrix.cf[n*2, 20] <- traits_name[t]
      matrix.cf[n*2, 21] <- d
      matrix.cf[n*2, 22] <- length(cf.tau.binaryW$tau.p.adjust)
      matrix.cf[n*2, 23] <- md
      matrix.cf[n*2, 24] <- abs((as.numeric(matrix.cf[n*2, 8]) - as.numeric(matrix.cf[n*2, 2]))/as.numeric(matrix.cf[n*2, 2]))
      matrix.cf[n*2, 25] <- abs((as.numeric(matrix.cf[n*2, 7]) - as.numeric(matrix.cf[n*2, 3]))/as.numeric(matrix.cf[n*2, 3]))

      matrix.driv[n*2, 1:9] <- as.numeric(quantile(driv.tau.binaryW$point, c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 1)))
      matrix.driv[n*2, 10] <- as.numeric(matrix.driv[n*2, 6]) - as.numeric(matrix.driv[n*2, 4]) # IQR
      matrix.driv[n*2, 11] <- as.numeric(matrix.driv[n*2, 7]) - as.numeric(matrix.driv[n*2, 3]) # 90% - 10%
      matrix.driv[n*2, 12] <- as.numeric(matrix.driv[n*2, 8]) - as.numeric(matrix.driv[n*2, 2]) # 95% - 5%
      matrix.driv[n*2, 13] <- summary(driv.tau.binaryW$point)[4]
      matrix.driv[n*2, 14] <- length(which(driv.tau.binaryW$point < 0))/length(driv.tau.binaryW$p_value_corrected)
      matrix.driv[n*2, 15] <- length(which(driv.tau.binaryW$p_value_corrected < 0.1))/length(driv.tau.binaryW$p_value_corrected)
      matrix.driv[n*2, 16] <- length(which(driv.tau.binaryW$p_value < 0.1))/length(driv.tau.binaryW$p_value)
      matrix.driv[n*2, 17] <- length(which(driv.tau.binaryW$p_value_corrected < 0.05))/length(driv.tau.binaryW$p_value_corrected)
      matrix.driv[n*2, 18] <- length(which(driv.tau.binaryW$p_value < 0.05))/length(driv.tau.binaryW$p_value)
      matrix.driv[n*2, 19] <- "Binary W"
      matrix.driv[n*2, 20] <- traits_name[t]
      matrix.driv[n*2, 21] <- d
      matrix.driv[n*2, 22] <- length(driv.tau.binaryW$p_value_corrected)
      matrix.driv[n*2, 23] <- md
      matrix.driv[n*2, 24] <- abs((as.numeric(matrix.driv[n*2, 8]) - as.numeric(matrix.driv[n*2, 2]))/as.numeric(matrix.driv[n*2, 2]))
      matrix.driv[n*2, 25] <- abs((as.numeric(matrix.driv[n*2, 7]) - as.numeric(matrix.driv[n*2, 3]))/as.numeric(matrix.driv[n*2, 3]))

      n <- n+1

    }
  }
}

colnames(matrix.ivcf) <- c("Min", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "Max", "50% Coverage (25% - 75%)", "80% Coverage (10% - 90%)", "90% Coverage (5% - 95%)", "Mean",
                           "Positive Tau Percentage", "Significant Tau Percentage (P.val.adjusted < 0.1)", "Significant Tau Percentage (P.val < 0.1)", "Significant Tau Percentage (P.val.adjusted < 0.05)", "Significant Tau Percentage (P.val < 0.05)", 
                           "Treatment Type", "Treatment", "Disease", "Sample Size", "Model", "absolute % difference (95th vs 5th percentile)", "absolute % difference (90th vs 10th percentile)")
colnames(matrix.driv) <- c("Min", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "Max", "50% Coverage (25% - 75%)", "80% Coverage (10% - 90%)", "90% Coverage (5% - 95%)", "Mean",
                           "Positive Tau Percentage", "Significant Tau Percentage (P.val.adjusted < 0.1)", "Significant Tau Percentage (P.val < 0.1)", "Significant Tau Percentage (P.val.adjusted < 0.05)", "Significant Tau Percentage (P.val < 0.05)", 
                           "Treatment Type", "Treatment", "Disease", "Sample Size", "Model", "absolute % difference (95th vs 5th percentile)", "absolute % difference (90th vs 10th percentile)")
colnames(matrix.cf) <- c("Min", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "Max", "50% Coverage (25% - 75%)", "80% Coverage (10% - 90%)", "90% Coverage (5% - 95%)", "Mean",
                         "Positive Tau Percentage", "Significant Tau Percentage (P.val.adjusted < 0.1)", "Significant Tau Percentage (P.val < 0.1)", "Significant Tau Percentage (P.val.adjusted < 0.05)", "Significant Tau Percentage (P.val < 0.05)", 
                         "Treatment Type", "Treatment", "Disease", "Sample Size", "Model", "absolute % difference (95th vs 5th percentile)", "absolute % difference (90th vs 10th percentile)")

matrix.ivcf <- as.data.frame(matrix.ivcf)
matrix.driv <- as.data.frame(matrix.driv)
matrix.cf <- as.data.frame(matrix.cf)

matrix.ivcf <- matrix.ivcf[, c("Treatment Type", "Treatment", "Disease", "Model", "Sample Size",
                               "Positive Tau Percentage", "Significant Tau Percentage (P.val < 0.1)", "Significant Tau Percentage (P.val < 0.05)", 
                               "Min", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "Max", "50% Coverage (25% - 75%)", "80% Coverage (10% - 90%)", "90% Coverage (5% - 95%)", "Mean",
                               "Significant Tau Percentage (P.val.adjusted < 0.1)", "Significant Tau Percentage (P.val.adjusted < 0.05)", "absolute % difference (95th vs 5th percentile)", "absolute % difference (90th vs 10th percentile)")]
matrix.driv <- matrix.driv[, c("Treatment Type", "Treatment", "Disease", "Model", "Sample Size",
                               "Positive Tau Percentage", "Significant Tau Percentage (P.val < 0.1)", "Significant Tau Percentage (P.val < 0.05)", 
                               "Min", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "Max", "50% Coverage (25% - 75%)", "80% Coverage (10% - 90%)", "90% Coverage (5% - 95%)", "Mean",
                               "Significant Tau Percentage (P.val.adjusted < 0.1)", "Significant Tau Percentage (P.val.adjusted < 0.05)", "absolute % difference (95th vs 5th percentile)", "absolute % difference (90th vs 10th percentile)")]
matrix.cf <- matrix.cf[, c("Treatment Type", "Treatment", "Disease", "Model", "Sample Size",
                            "Positive Tau Percentage", "Significant Tau Percentage (P.val < 0.1)", "Significant Tau Percentage (P.val < 0.05)", 
                            "Min", "5%", "10%", "25%", "50%", "75%", "90%", "95%", "Max", "50% Coverage (25% - 75%)", "80% Coverage (10% - 90%)", "90% Coverage (5% - 95%)", "Mean",
                            "Significant Tau Percentage (P.val.adjusted < 0.1)", "Significant Tau Percentage (P.val.adjusted < 0.05)", "absolute % difference (95th vs 5th percentile)", "absolute % difference (90th vs 10th percentile)")]

matrix.cf$Method <- "CF"
matrix.driv$Method <- "DRIV"
matrix.ivcf$Method <- "IVCF"

fwrite(matrix.ivcf, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/05_summary_statistics/02_taus/method_IVCF_propchange_additionquantile.csv")
fwrite(matrix.cf, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/05_summary_statistics/02_taus/method_CF_propchange_additionquantile.csv")
fwrite(matrix.driv, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/05_summary_statistics/02_taus/method_DRIV_propchange_additionquantile.csv")

matrix_final <- rbind(matrix.cf, matrix.driv, matrix.ivcf)
fwrite(matrix_final, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/05_summary_statistics/02_taus/method_propchange_additionquantile.csv")
