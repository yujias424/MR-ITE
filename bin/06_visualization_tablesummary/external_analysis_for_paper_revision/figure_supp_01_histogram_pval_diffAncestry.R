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

diseases <- c("CAD") 
traits <- c("LDL", "TC")
traits_name <- c("LDL-C", "Total Cholesterol")
models <- c("model3")

for (ml in models){
  for (t in 1:length(traits)){
    for (d in diseases){
      
      ml <- models[1]
      t <- 1
      d <- diseases[1]

      # continuous W
      iv.tau.continuousW.white <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_", ml, ".csv"))
      iv.tau.continuousW.black <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_", ml, "_black.csv"))
      iv.tau.continuousW.asian <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_", ml, "_southasian.csv"))

      # # binary W
      # iv.tau.binaryW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_", ml, ".csv"))
      # cf.tau.binaryW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_", ml, ".csv"))
      # driv.tau.binaryW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_", ml, ".csv"))

      # making plot (continuous W)
      iv.tau.white <- iv.tau.continuousW.white$tau
      iv.tau.sig.white <- iv.tau.continuousW.white$tau[iv.tau.continuousW.white$tau.pval<0.1]
      iv.tau.sig.05.white <- iv.tau.continuousW.white$tau[iv.tau.continuousW.white$tau.pval<0.05]
      if (length(iv.tau.sig.white) == 0){
        iv.tau.sig.white <- c(mean(iv.tau.white))
      } else {
        iv.tau.sig.white <- iv.tau.sig.white
      }

      iv.tau.black <- iv.tau.continuousW.black$tau
      iv.tau.sig.black <- iv.tau.continuousW.black$tau[iv.tau.continuousW.black$tau.pval<0.1]
      iv.tau.sig.05.black <- iv.tau.continuousW.black$tau[iv.tau.continuousW.black$tau.pval<0.05]
      if (length(iv.tau.sig.black) == 0){
        iv.tau.sig.black <- c(mean(iv.tau.black))
      } else {
        iv.tau.sig.black <- iv.tau.sig.black
      }

      iv.tau.asian <- iv.tau.continuousW.asian$tau
      iv.tau.sig.asian <- iv.tau.continuousW.asian$tau[iv.tau.continuousW.asian$tau.pval<0.1]
      iv.tau.sig.05.asian <- iv.tau.continuousW.asian$tau[iv.tau.continuousW.asian$tau.pval<0.05]
      if (length(iv.tau.sig.asian) == 0){
        iv.tau.sig.asian <- c(mean(iv.tau.asian))
      } else {
        iv.tau.sig.asian <- iv.tau.sig.asian
      }

      # taus.iv.plot.continuousW <- c(iv.tau.white, iv.tau.black, iv.tau.asian, 
      #                               iv.tau.sig.white, iv.tau.sig.black, iv.tau.sig.asian,
      #                               iv.tau.sig.05.white, iv.tau.sig.05.black, iv.tau.sig.05.asian)
      # var.iv.plot.continuousW <- c(rep("European", length(iv.tau.white)), rep("African", length(iv.tau.black)), rep("South-Asian", length(iv.tau.asian)), 
      #                              rep("European", length(iv.tau.sig.white)), rep("African", length(iv.tau.sig.black)), rep("South-Asian", length(iv.tau.sig.asian)),
      #                              rep("European", length(iv.tau.sig.05.white)), rep("African", length(iv.tau.sig.05.black)), rep("South-Asian", length(iv.tau.sig.05.asian))) 
      # sig.iv.plot.continuousW <- c(rep("Origin", length(iv.tau.white)), rep("Origin", length(iv.tau.black)), rep("Origin", length(iv.tau.asian)),
      #                              rep("P.val < 0.1", length(iv.tau.sig.white)), rep("P.val < 0.1", length(iv.tau.sig.black)), rep("P.val < 0.1", length(iv.tau.sig.asian)), 
      #                              rep("P.val < 0.05", length(iv.tau.sig.05.white)), rep("P.val < 0.05", length(iv.tau.sig.05.black)), rep("P.val < 0.05", length(iv.tau.sig.05.asian))) 

      # plot.taus.iv.continuousW <- data.frame("tau" = taus.iv.plot.continuousW, "variable" = var.iv.plot.continuousW, "significance" = sig.iv.plot.continuousW)
      # plot.taus.iv.continuousW$significance <- factor(plot.taus.iv.continuousW$significance, levels = c("Origin", "P.val < 0.1", "P.val < 0.05"))
      # plot.taus.iv.continuousW$variable <- factor(plot.taus.iv.continuousW$variable, levels = c("European", "African", "South-Asian"))

      # p.iv <- ggplot(plot.taus.iv.continuousW, aes(x=tau, fill=significance)) +
      #             geom_histogram(color="black",  bins = 50, position = 'identity', alpha=0.6, key_glyph = "polygon") + 
      #             scale_fill_manual(name=TeX("$\\tau$ type"), values=c("Origin"="lightblue", "P.val < 0.1"="orange", "P.val < 0.05"="red"), labels=c("Origin"=TeX("$\\tau$"), "P.val < 0.1"=TeX("$\\tau$ with P.val < 0.1"), "P.val < 0.05"=TeX("$\\tau$ with P.val < 0.05"))) +
      #             facet_grid(rows = vars(variable), scales = "free_y", labeller = label_parsed) + ggtitle(traits_name[t]) + 
      #             theme_classic() + 
      #             theme(text = element_text(size = 35), plot.title = element_text(hjust = 0.5), 
      #                   axis.text = element_text(size = 35), axis.line=element_line(color='black'),
      #                   legend.text = element_text(size=35),
      #                   legend.title=element_text(size=35),
      #                   panel.grid.major = element_blank(),
      #                   panel.grid.minor = element_blank(),
      #                   strip.background = element_blank(),
      #                   legend.position="right", legend.text.align = 0) + 
      #             geom_vline(xintercept=0, color="red", linetype="dashed") + 
      #             xlab("Individual Treatment Effect") + ylab("Count")

      taus.iv.plot.continuousW <- c(iv.tau.black, iv.tau.asian, 
                                    iv.tau.sig.black, iv.tau.sig.asian,
                                    iv.tau.sig.05.black, iv.tau.sig.05.asian)
      var.iv.plot.continuousW <- c( rep("African", length(iv.tau.black)), rep("South-Asian", length(iv.tau.asian)), 
                                   rep("African", length(iv.tau.sig.black)), rep("South-Asian", length(iv.tau.sig.asian)),
                                   rep("African", length(iv.tau.sig.05.black)), rep("South-Asian", length(iv.tau.sig.05.asian))) 
      sig.iv.plot.continuousW <- c(rep("Origin", length(iv.tau.black)), rep("Origin", length(iv.tau.asian)),
                                   rep("P.val < 0.1", length(iv.tau.sig.black)), rep("P.val < 0.1", length(iv.tau.sig.asian)), 
                                  rep("P.val < 0.05", length(iv.tau.sig.05.black)), rep("P.val < 0.05", length(iv.tau.sig.05.asian))) 

      plot.taus.iv.continuousW <- data.frame("tau" = taus.iv.plot.continuousW, "variable" = var.iv.plot.continuousW, "significance" = sig.iv.plot.continuousW)
      plot.taus.iv.continuousW$significance <- factor(plot.taus.iv.continuousW$significance, levels = c("Origin", "P.val < 0.1", "P.val < 0.05"))
      plot.taus.iv.continuousW$variable <- factor(plot.taus.iv.continuousW$variable, levels = c("European", "African", "South-Asian"))

      p.iv <- ggplot(plot.taus.iv.continuousW, aes(x=tau, fill=significance)) +
                  geom_histogram(color="black",  bins = 50, position = 'identity', alpha=0.6, key_glyph = "polygon") + 
                  scale_fill_manual(name=TeX("$\\tau$ type"), values=c("Origin"="lightblue", "P.val < 0.1"="orange", "P.val < 0.05"="red"), labels=c("Origin"=TeX("$\\tau$"), "P.val < 0.1"=TeX("$\\tau$ with P.val < 0.1"), "P.val < 0.05"=TeX("$\\tau$ with P.val < 0.05"))) +
                  facet_grid(rows = vars(variable), scales = "free_y", labeller = label_parsed) + ggtitle(traits_name[t]) + 
                  theme_classic() + 
                  theme(text = element_text(size = 35), plot.title = element_text(hjust = 0.5), 
                        axis.text = element_text(size = 35), axis.line=element_line(color='black'),
                        legend.text = element_text(size=35),
                        legend.title=element_text(size=35),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        strip.background = element_blank(),
                        legend.position="right", legend.text.align = 0) + 
                  geom_vline(xintercept=0, color="red", linetype="dashed") + 
                  xlab("Individual Treatment Effect") + ylab("Count")

      # p.iv
      p.iv.nolegend <- p.iv + guides(fill = "none")
      assign(paste0("p.continuousW.", traits[t], ".", d), p.iv.nolegend)

    }
  }
}
p.continuousW.LDL.CAD

histogram_p <- (p1_title) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 35), legend.position = "right")
ggsave(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/extra_analysis_for_paper_revision/figure_supp_fig1_diffAnestry_model3.png"), p.continuousW.LDL.CAD, dpi=300, 
        width = 10, height = 8, units = "cm",scale = 3, limitsize = FALSE)
# break
