#' This code is to visualize the MSE results from our latest simulation on bias using a table.
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
valid_snps <- c(20,40,60)
# valid_snps <- c(100,200,300)
tau_scene <- seq(1,8)
sim_times <- seq(1,30)
top_porp <- c(0.1, 0.5)
pleiotropy_scene <- c(1,2,3)

for (tp in top_porp){

    message(paste0("Current running proportion cutoff at ", tp*100, "%."))
    first_ps <- T

    for (ps in pleiotropy_scene){

        pleotropic.type  <- ps

        for (vs in valid_snps){
            
            invalid.snps.count <- vs

            for (ts in tau_scene){

                ate_cf_bias_c_tail <- c()
                late_all_bias_c_tail <- c()
                late_few_bias_c_tail <- c()
                ite_cf_bias_c_tail <- c()
                ite_ivcf_all_bias_c_tail <- c()
                ite_ivcf_few_bias_c_tail <- c()

                for (st in sim_times){

                    # tau_estimate <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffSNPs_diffN/01_simulation/tau/estimate/", 
                    #                                 ps, "_", vs, "_", ts, "_", st, "_estimate.csv.gz"))
                    # tau_std <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffSNPs_diffN/01_simulation/tau/std/", 
                    #                                 ps, "_", vs, "_", ts, "_", st, "_std.csv.gz"))

                    tau_estimate <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffN/01_simulation/tau/estimate/", 
                                                    ps, "_", vs, "_", ts, "_", st, "_estimate.csv.gz"))
                    tau_std <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffN/01_simulation/tau/std/", 
                                                    ps, "_", vs, "_", ts, "_", st, "_std.csv.gz"))

                    ate_cf_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$ATE_CF, sd = 0), nrow=50000))
                    late_all_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$LATE_all, sd = 0), nrow=50000))
                    late_few_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$LATE_few, sd = 0), nrow=50000))
                    ite_cf_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$ITE_CF, sd = 0), nrow=50000))
                    ite_ivcf_all_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$ITE_IVCF_all, sd = 0), nrow=50000))
                    ite_ivcf_few_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$ITE_IVCF_few, sd = 0), nrow=50000))
                    true_tau_mat <- as.data.frame(matrix(rnorm(50000 * 1, mean = tau_estimate$true_tau, sd = 0), nrow=50000))

                    true_tau_mat$dummy <- NA
                    true_tau_mat <- true_tau_mat[order(true_tau_mat$V1, decreasing = T), ]
                    
                    # sort by ite cf
                    ite_cf_mat$dummy <- NA
                    ite_ivcf_all_mat$dummy <- NA
                    ite_ivcf_few_mat$dummy <- NA
                    ite_cf_mat <- ite_cf_mat[rownames(true_tau_mat), ]
                    ite_ivcf_all_mat <- ite_ivcf_all_mat[rownames(true_tau_mat), ]
                    ite_ivcf_few_mat <- ite_ivcf_few_mat[rownames(true_tau_mat), ]

                    selected_p_c <- 50000 * tp

                    # calculate the bias (tail 10%)
                    ate_cf_bias_tail <- bias(ate_cf_mat[(50000-selected_p_c):50000,1], true_tau_mat[(50000-selected_p_c):50000,1])
                    late_all_bias_tail <- bias(late_all_mat[(50000-selected_p_c):50000,1], true_tau_mat[(50000-selected_p_c):50000,1])
                    late_few_bias_tail <- bias(late_few_mat[(50000-selected_p_c):50000,1], true_tau_mat[(50000-selected_p_c):50000,1])
                    ite_cf_bias_tail <- bias(ite_cf_mat[(50000-selected_p_c):50000,1], true_tau_mat[(50000-selected_p_c):50000,1])
                    ite_ivcf_all_bias_tail <- bias(ite_ivcf_all_mat[(50000-selected_p_c):50000,1], true_tau_mat[(50000-selected_p_c):50000,1])
                    ite_ivcf_few_bias_tail <- bias(ite_ivcf_few_mat[(50000-selected_p_c):50000,1], true_tau_mat[(50000-selected_p_c):50000,1])

                    ate_cf_bias_c_tail <- c(ate_cf_bias_c_tail, abs(ate_cf_bias_tail))
                    late_all_bias_c_tail <- c(late_all_bias_c_tail, abs(late_all_bias_tail))
                    late_few_bias_c_tail <- c(late_few_bias_c_tail, abs(late_few_bias_tail))
                    ite_cf_bias_c_tail <- c(ite_cf_bias_c_tail, abs(ite_cf_bias_tail))
                    ite_ivcf_all_bias_c_tail <- c(ite_ivcf_all_bias_c_tail, abs(ite_ivcf_all_bias_tail))
                    ite_ivcf_few_bias_c_tail <- c(ite_ivcf_few_bias_c_tail, abs(ite_ivcf_few_bias_tail))

                    # break
                }

                print(summary(true_tau_mat$V1))

                bias_df <- data.frame("ITE" = ite_cf_bias_c_tail,
                                      "MR-ITE (keep pleiotropy)" = ite_ivcf_all_bias_c_tail,
                                      "MR-ITE (remove pleiotropy)" = ite_ivcf_few_bias_c_tail,
                                      "ATE" = ate_cf_bias_c_tail,
                                      "MR-ATE (keep pleiotropy)" = late_all_bias_c_tail,
                                      "MR-ATE (remove pleiotropy)" = late_few_bias_c_tail)
                
                bias_df$id <- rownames(bias_df)

                melt_bias_df <- melt(bias_df)
                melt_bias_df$variable <- c(rep("ITE", 30),
                                           rep("MR-ITE (keep pleiotropy)", 30),
                                           rep("MR-ITE (remove pleiotropy)", 30),
                                           rep("ATE", 30),
                                           rep("MR-ATE (keep pleiotropy)", 30),
                                           rep("MR-ATE (remove pleiotropy)", 30))
                
                melt_bias_df$tauscene <- ts
                melt_bias_df$pleiotropyscene <- ps

                if (ts == 1){
                    print(ts)
                    bias_df_res <- melt_bias_df
                } else {
                    print(ts)
                    bias_df_res <- rbind(bias_df_res, melt_bias_df)
                }
                
            }

            colnames(bias_df_res)[2:3] <- c("Methods", "Bias")

        }

        if (first_ps){
            bias_table_final <- bias_df_res
            first_ps <- F
        } else {
            bias_table_final <- rbind(bias_table_final, bias_df_res)
        }
    }

    bias_table <- bias_table_final %>%
        group_by(Methods, tauscene, pleiotropyscene) %>%
        summarise_at(vars(Bias), list(name = mean))
    bias_table <- as.data.frame(bias_table)
    colnames(bias_table) <- c("Methods", "Tau Scenario", "Pleiotropy Scenario", "Bias")

    # dir.create("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffSNPs_diffN/01_simulation/table/", recursive = T, showWarnings = F)
    # fwrite(bias_table, paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffSNPs_diffN/01_simulation/table/bias_table_", tp, ".csv"))

    dir.create("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffN/01_simulation/table/", recursive = T, showWarnings = F)
    fwrite(bias_table, paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffN/01_simulation/table/bias_table_", tp, ".csv"))

}


