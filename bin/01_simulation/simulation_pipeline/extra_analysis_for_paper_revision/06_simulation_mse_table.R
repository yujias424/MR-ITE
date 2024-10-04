#' This code is to visualize the mse results from our latest simulation on mse using a table.
#' 
#' @author Yujia Shi
#' @date 2023.08.25

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(rstatix)
  library(ggsignif)
  library(patchwork)
  library(data.table)
  library(Metrics)
  library(reshape2)
})

pleiotropy_scene <- c(1,2,3)
# valid_snps <- c(20,40,60)
valid_snps <- c(100,200,300)
tau_scene <- seq(1,8)
sim_times <- seq(1,30)
top_porp <- c(1)
pleiotropy_scene <- c(1,2,3)

first_ps <- T

for (tp in top_porp){

    message(paste0("Current running proportion cutoff at ", tp*100, "%."))

    for (ps in pleiotropy_scene){

        pleotropic.type  <- ps

        for (vs in valid_snps){
            
            invalid.snps.count <- vs

            for (ts in tau_scene){

                ate_cf_mse_c_top <- c()
                late_all_mse_c_top <- c()
                late_few_mse_c_top <- c()
                ite_cf_mse_c_top <- c()
                ite_ivcf_all_mse_c_top <- c()
                ite_ivcf_few_mse_c_top <- c()

                for (st in sim_times){

                    tau_estimate <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffSNPs_diffN/01_simulation/tau/estimate/", 
                                                    ps, "_", vs, "_", ts, "_", st, "_estimate.csv.gz"))
                    tau_std <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffSNPs_diffN/01_simulation/tau/std/", 
                                                    ps, "_", vs, "_", ts, "_", st, "_std.csv.gz"))

                    ate_cf_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$ATE_CF, sd = 0), nrow=50000))
                    late_all_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$LATE_all, sd = 0), nrow=50000))
                    late_few_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$LATE_few, sd = 0), nrow=50000))
                    ite_cf_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$ITE_CF, sd = 0), nrow=50000))
                    ite_ivcf_all_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$ITE_IVCF_all, sd = 0), nrow=50000))
                    ite_ivcf_few_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$ITE_IVCF_few, sd = 0), nrow=50000))
                    true_tau_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$true_tau, sd = 0), nrow=50000))

                    # sort by ite cf
                    ite_cf_mat$dummy <- NA
                    ite_ivcf_all_mat$dummy <- NA
                    ite_ivcf_few_mat$dummy <- NA
                    ite_cf_mat <- ite_cf_mat[rownames(true_tau_mat), ]
                    ite_ivcf_all_mat <- ite_ivcf_all_mat[rownames(true_tau_mat), ]
                    ite_ivcf_few_mat <- ite_ivcf_few_mat[rownames(true_tau_mat), ]

                    selected_p_c <- 50000 * tp

                    # calculate the mse (top 10%)
                    ate_cf_mse_top <- mse(ate_cf_mat[1:selected_p_c,1], true_tau_mat[1:selected_p_c,1])
                    late_all_mse_top <- mse(late_all_mat[1:selected_p_c,1], true_tau_mat[1:selected_p_c,1])
                    late_few_mse_top <- mse(late_few_mat[1:selected_p_c,1], true_tau_mat[1:selected_p_c,1])
                    ite_cf_mse_top <- mse(ite_cf_mat[1:selected_p_c,1], true_tau_mat[1:selected_p_c,1])
                    ite_ivcf_all_mse_top <- mse(ite_ivcf_all_mat[1:selected_p_c,1], true_tau_mat[1:selected_p_c,1])
                    ite_ivcf_few_mse_top <- mse(ite_ivcf_few_mat[1:selected_p_c,1], true_tau_mat[1:selected_p_c,1])

                    ate_cf_mse_c_top <- c(ate_cf_mse_c_top, abs(ate_cf_mse_top))
                    late_all_mse_c_top <- c(late_all_mse_c_top, abs(late_all_mse_top))
                    late_few_mse_c_top <- c(late_few_mse_c_top, abs(late_few_mse_top))
                    ite_cf_mse_c_top <- c(ite_cf_mse_c_top, abs(ite_cf_mse_top))
                    ite_ivcf_all_mse_c_top <- c(ite_ivcf_all_mse_c_top, abs(ite_ivcf_all_mse_top))
                    ite_ivcf_few_mse_c_top <- c(ite_ivcf_few_mse_c_top, abs(ite_ivcf_few_mse_top))

                }

                print(summary(true_tau_mat$V1))

                mse_df <- data.frame("ITE" = ite_cf_mse_c_top,
                                      "MR-ITE (keep pleiotropy)" = ite_ivcf_all_mse_c_top,
                                      "MR-ITE (remove pleiotropy)" = ite_ivcf_few_mse_c_top,
                                      "CF-ATE" = ate_cf_mse_c_top,
                                      "MR-ATE (keep pleiotropy)" = late_all_mse_c_top,
                                      "MR-ATE (remove pleiotropy)" = late_few_mse_c_top)
                
                mse_df$id <- rownames(mse_df)

                melt_mse_df <- melt(mse_df)
                melt_mse_df$variable <- c(rep("ITE", 30),
                                           rep("MR-ITE (keep pleiotropy)", 30),
                                           rep("MR-ITE (remove pleiotropy)", 30),
                                           rep("ATE", 30),
                                           rep("MR-ATE (keep pleiotropy)", 30),
                                           rep("MR-ATE (remove pleiotropy)", 30))
                
                melt_mse_df$tauscene <- ts
                melt_mse_df$pleiotropyscene <- ps

                if (ts == 1){
                    print(ts)
                    mse_df_res <- melt_mse_df
                } else {
                    print(ts)
                    mse_df_res <- rbind(mse_df_res, melt_mse_df)
                }
                
            }

            colnames(mse_df_res)[2:3] <- c("Methods", "MSE")

        }

        if (first_ps){
            mse_table_final <- mse_df_res
            first_ps <- F
        } else {
            mse_table_final <- rbind(mse_table_final, mse_df_res)
        }

    }
    
}

mse_table <- mse_table_final %>%
  group_by(Methods, tauscene, pleiotropyscene) %>%
  summarise_at(vars(MSE), list(name = mean))
mse_table <- as.data.frame(mse_table)
colnames(mse_table) <- c("Methods", "Tau Scenario", "Pleiotropy Scenario", "Mean Squared Error")

dir.create("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffSNPs_diffN/01_simulation/table/", recursive = T, showWarnings = F)
fwrite(mse_table, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffSNPs_diffN/01_simulation/table/mse_table.csv")
