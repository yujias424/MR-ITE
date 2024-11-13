#' This code is to visualize the MSE results from our latest simulation on bias.
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
tau_scene <- seq(1,8)
sim_times <- seq(1,30)

# tail_porp <- c(0.05, 0.1, 0.25, 0.5, 0.75, 1)
tail_porp <- c(0.1)

for (tp in tail_porp){

    message(paste0("Current running proportion cutoff at ", tp*100, "%."))

    for (ps in pleiotropy_scene){

        pleotropic.type  <- ps

        for (vs in valid_snps){
            
            invalid.snps.count <- vs

            for (ts in tau_scene){

                ite_cf_bias_c_tail <- c()
                ite_ivcf_few_bias_c_tail <- c()
                ite_driv_few_bias_c_tail <- c()

                for (st in sim_times){

                    tau_estimate <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/01_simulation/tau/estimate/", 
                                                    ps, "_", vs, "_", ts, "_", st, "_estimate.csv.gz"))
                    tau_estimate_driv <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/01_simulation/tau/estimate/", 
                                                ps, "_", vs, "_", ts, "_", st, "_estimate_driv.csv.gz"))
                    tau_std <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/01_simulation/tau/std/", 
                                                    ps, "_", vs, "_", ts, "_", st, "_std.csv.gz"))

                    ite_cf_mat <- as.data.frame(matrix(rnorm(10000 * 1, mean = tau_estimate$ITE_CF, sd = 0), nrow=10000))
                    ite_ivcf_few_mat <- as.data.frame(matrix(rnorm(10000 * 1, mean = tau_estimate$ITE_IVCF_few, sd = 0), nrow=10000))
                    ite_driv_few_mat <- as.data.frame(matrix(rnorm(10000 * 1, mean = tau_estimate_driv$ITE_DRIV_few, sd = 0), nrow=10000))
                    true_tau_mat <- as.data.frame(matrix(rnorm(10000 * 1, mean = tau_estimate$true_tau, sd = 0), nrow=10000))

                    true_tau_mat$dummy <- NA
                    true_tau_mat <- true_tau_mat[order(true_tau_mat$V1, decreasing = T), ]

                    # sort by ite cf
                    ite_cf_mat$dummy <- NA
                    ite_ivcf_few_mat$dummy <- NA
                    ite_driv_few_mat$dummy <- NA
                    ite_cf_mat <- ite_cf_mat[rownames(true_tau_mat), ]
                    ite_ivcf_few_mat <- ite_ivcf_few_mat[rownames(true_tau_mat), ]
                    ite_driv_few_mat <- ite_driv_few_mat[rownames(true_tau_mat), ]

                    selected_p_c <- 10000 * tp

                    # calculate the bias (tail 10%)
                    ite_cf_bias_tail <- bias(ite_cf_mat[(10000-selected_p_c):10000,1], true_tau_mat[(10000-selected_p_c):10000,1])
                    ite_ivcf_few_bias_tail <- bias(ite_ivcf_few_mat[(10000-selected_p_c):10000,1], true_tau_mat[(10000-selected_p_c):10000,1])
                    ite_driv_few_bias_tail <- bias(ite_driv_few_mat[(10000-selected_p_c):10000,1], true_tau_mat[(10000-selected_p_c):10000,1])

                    ite_cf_bias_c_tail <- c(ite_cf_bias_c_tail, abs(ite_cf_bias_tail))
                    ite_ivcf_few_bias_c_tail <- c(ite_ivcf_few_bias_c_tail, abs(ite_ivcf_few_bias_tail))
                    ite_driv_few_bias_c_tail <- c(ite_driv_few_bias_c_tail, abs(ite_driv_few_bias_tail))

                }

                print(summary(true_tau_mat$V1))

                bias_df <- data.frame("ITE" = ite_cf_bias_c_tail,
                                      "MR-ITE-IVCF (remove pleiotropy)" = ite_ivcf_few_bias_c_tail,
                                      "MR-ITE-DRIV (remove pleiotropy)" = ite_driv_few_bias_c_tail)
                
                bias_df$id <- rownames(bias_df)

                melt_bias_df <- melt(bias_df)
                melt_bias_df$variable <- c(rep("ITE", 30),
                                           rep("MR-ITE-IVCF (remove pleiotropy)", 30),
                                           rep("MR-ITE-DRIV (remove pleiotropy)", 30))
                
                melt_bias_df$tauscene <- ts

                if (ts == 1){
                    print(ts)
                    bias_df_res <- melt_bias_df
                } else {
                    print(ts)
                    bias_df_res <- rbind(bias_df_res, melt_bias_df)
                }
                
            }

            colnames(bias_df_res)[2:3] <- c("Methods", "Bias")
            bxp <- ggboxplot(
                bias_df_res, x = "tauscene", y = "Bias", fill = "Methods",
                palette = "npg"
            )

            stat.test.ite.greater <- bias_df_res %>%
                group_by(tauscene) %>%
                t_test(Bias ~ Methods, paired = T, alternative = "two.sided",
                        comparisons = list(c("ITE", "MR-ITE-DRIV (remove pleiotropy)"), c("MR-ITE-DRIV (remove pleiotropy)", "MR-ITE-IVCF (remove pleiotropy)")))
                                        
            stat.test <- bind_rows(stat.test.ite.greater)

            # Box plots with p-values
            stat.test <- stat.test %>%
                arrange(stat.test) %>% 
                add_xy_position(x = "tauscene", dodge = 0.8, step.increase = 0.05) # 0.12

            for (l in 1:dim(stat.test)[1]){
                if (is.na(stat.test[l, 12])){
                    stat.test[l, 12] <- stat.test[l-1, 12] + 0.05 # 0.05
                }
            }

            head(as.data.frame(stat.test))

            bxp_final <- bxp + 
                        stat_pvalue_manual(
                            stat.test, label = "p.adj.signif", tip.length = 0.02, size = 10) + 
                        xlab("Tau Scenarios") + ylab("Bias") + ggtitle(paste0("Pleiotropy Scene: ", pleotropic.type, "; Invalid Snp Counts: ", invalid.snps.count)) +
                        theme(axis.text = element_text(size=35), 
                                axis.title = element_text(size=35),
                                plot.title = element_text(hjust = 0.5, size=35),
                                legend.text = element_text(size=35),
                                legend.title=element_text(size=35),
                                legend.position="right")

            assign(paste0("plot_", pleotropic.type, "_", invalid.snps.count), bxp_final)

        }

    }

    p.final <- ((plot_1_20 | plot_1_40 | plot_1_60) /
                (plot_2_20 | plot_2_40 | plot_2_60) /
                (plot_3_20 | plot_3_40 | plot_3_60)) +
                plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 45))

    ggsave(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/plot/bias/all_bias_", tp * 100,"_ITEbenchmark.png"), p.final, dpi=300, 
            width = 100, height = 60, units = "cm", scale = 1.6, limitsize = F)
    
}

