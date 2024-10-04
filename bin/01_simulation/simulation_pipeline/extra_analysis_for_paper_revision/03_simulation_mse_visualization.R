#' This code is to visualize the mse results from our latest simulation on mse.
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
valid_snps <- c(100,200,300)
tau_scene <- seq(1,8)
sim_times <- seq(1,30)
top_porp <- c(1)
pleiotropy_scene <- c(1,2,3)

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

                    tau_estimate <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision_diffSNPs_diffN/01_simulation/tau/estimate/", 
                                                    ps, "_", vs, "_", ts, "_", st, "_estimate.csv.gz"))
                    tau_std <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision_diffSNPs_diffN/01_simulation/tau/std/", 
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

                    # break
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

                if (ts == 1){
                    print(ts)
                    mse_df_res <- melt_mse_df
                } else {
                    print(ts)
                    mse_df_res <- rbind(mse_df_res, melt_mse_df)
                }
                
            }

            colnames(mse_df_res)[2:3] <- c("Methods", "MSE")
            bxp <- ggboxplot(
                mse_df_res, x = "tauscene", y = "MSE", fill = "Methods",
                palette = "npg"
            )

            stat.test.ite.greater <- mse_df_res %>%
                group_by(tauscene) %>%
                t_test(MSE ~ Methods, paired = T, alternative = "greater",
                        comparisons = list(c("MR-ITE (keep pleiotropy)", "MR-ITE (remove pleiotropy)")))

            stat.test.ate <- mse_df_res %>%
                group_by(tauscene) %>%
                t_test(MSE ~ Methods, paired = T, alternative = "greater", # two.sided, greater, less
                        comparisons = list(c("MR-ITE (remove pleiotropy)", "MR-ATE (remove pleiotropy)"), c("MR-ITE (keep pleiotropy)", "MR-ATE (keep pleiotropy)")))
                                        
            stat.test <- bind_rows(stat.test.ite.greater, stat.test.ate)

            # Box plots with p-values
            stat.test <- stat.test %>%
                arrange(stat.test) %>% 
                add_xy_position(x = "tauscene", dodge = 0.8, step.increase = 0.05) # 0.12

            for (l in 1:(dim(stat.test)[1]/3)){
                stat.test[l*3 - 2, 12] <- stat.test[l*3, 12]  # 0.05
                stat.test[l*3, 12] <- NA
            }

            for (l in 1:dim(stat.test)[1]){
                if (is.na(stat.test[l, 12])){
                    stat.test[l, 12] <- stat.test[l-1, 12] + 0.03 # 0.05
                }
            }

            head(as.data.frame(stat.test))

            bxp_final <- bxp + 
                        stat_pvalue_manual(
                            stat.test, label = "p.adj.signif", tip.length = 0.02, size = 10) + 
                        xlab("Tau Scenarios") + ylab("MSE") + ggtitle(paste0("Pleiotropy Scene: ", pleotropic.type, "; Invalid Snp Counts (Total SNPs: 500): ", invalid.snps.count)) +
                        theme(axis.text = element_text(size=35), 
                                axis.title = element_text(size=35),
                                plot.title = element_text(hjust = 0.5, size=35),
                                legend.text = element_text(size=35),
                                legend.title = element_text(size=35),
                                legend.position="right")

            assign(paste0("plot_", pleotropic.type, "_", invalid.snps.count), bxp_final)

        }
        
    }

    p.final <- ((plot_1_100 | plot_1_200 | plot_1_300) /
                (plot_2_100 | plot_2_200 | plot_2_300) /
                (plot_3_100 | plot_3_200 | plot_3_300)) +
                plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 35))

    ggsave(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/plot/extra_analysis_for_paper_revision/mse/all_mse_", tp * 100,".png"), p.final, dpi=300, 
            width = 100, height = 60, units = "cm", scale = 1.6, limitsize = F)
    
}

