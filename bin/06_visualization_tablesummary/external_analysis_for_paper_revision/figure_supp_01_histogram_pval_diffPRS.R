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
traits <- c("LDL", "TC", "IGF1")
traits_name <- c("LDL-C", "Total Cholesterol", "IGF1", "C-reative protein")
models <- c("model3")

for (ml in models){
  for (t in 1:length(traits)){
    for (d in diseases){

      # continuous W
      iv.tau.continuousW.PRSice2 <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_", ml, "_white.csv"))
      iv.tau.continuousW.LDPred2 <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_", ml, "_LDPred2.csv"))

      cor.test(iv.tau.continuousW.PRSice2$tau, iv.tau.continuousW.LDPred2$tau)

      # making plot (continuous W)
      iv.tau.PRSice2 <- iv.tau.continuousW.PRSice2$tau
      iv.tau.sig.PRSice2 <- iv.tau.continuousW.PRSice2$tau[iv.tau.continuousW.PRSice2$tau.pval<0.1]
      iv.tau.sig.05.PRSice2 <- iv.tau.continuousW.PRSice2$tau[iv.tau.continuousW.PRSice2$tau.pval<0.05]
      if (length(iv.tau.sig.PRSice2) == 0){
        iv.tau.sig.PRSice2 <- c(mean(iv.tau.PRSice2))
      } else {
        iv.tau.sig.PRSice2 <- iv.tau.sig.PRSice2
      }

      iv.tau.LDPred2 <- iv.tau.continuousW.LDPred2$tau
      iv.tau.sig.LDPred2 <- iv.tau.continuousW.LDPred2$tau[iv.tau.continuousW.LDPred2$tau.pval<0.1]
      iv.tau.sig.05.LDPred2 <- iv.tau.continuousW.LDPred2$tau[iv.tau.continuousW.LDPred2$tau.pval<0.05]
      if (length(iv.tau.sig.LDPred2) == 0){
        iv.tau.sig.LDPred2 <- c(mean(iv.tau.LDPred2))
      } else {
        iv.tau.sig.LDPred2 <- iv.tau.sig.LDPred2
      }


      taus.iv.plot.continuousW <- c(iv.tau.PRSice2, iv.tau.LDPred2, 
                                    iv.tau.sig.PRSice2, iv.tau.sig.LDPred2, 
                                    iv.tau.sig.05.PRSice2, iv.tau.sig.05.LDPred2)
      var.iv.plot.continuousW <- c(rep("PRSice2", length(iv.tau.PRSice2)), rep("LDPred2", length(iv.tau.LDPred2)), 
                                   rep("PRSice2", length(iv.tau.sig.PRSice2)), rep("LDPred2", length(iv.tau.sig.LDPred2)), 
                                   rep("PRSice2", length(iv.tau.sig.05.PRSice2)), rep("LDPred2", length(iv.tau.sig.05.LDPred2))) 
      sig.iv.plot.continuousW <- c(rep("Origin", length(iv.tau.PRSice2)), rep("Origin", length(iv.tau.LDPred2)), 
                                   rep("P.val < 0.1", length(iv.tau.sig.PRSice2)), rep("P.val < 0.1", length(iv.tau.sig.LDPred2)), 
                                   rep("P.val < 0.05", length(iv.tau.sig.05.PRSice2)), rep("P.val < 0.05", length(iv.tau.sig.05.LDPred2))) 

      plot.taus.iv.continuousW <- data.frame("tau" = taus.iv.plot.continuousW, "variable" = var.iv.plot.continuousW, "significance" = sig.iv.plot.continuousW)
      plot.taus.iv.continuousW$significance <- factor(plot.taus.iv.continuousW$significance, levels = c("Origin", "P.val < 0.1", "P.val < 0.05"))
      plot.taus.iv.continuousW$variable <- factor(plot.taus.iv.continuousW$variable, levels = c("PRSice2", "LDPred2"))

      p.iv <- ggplot(plot.taus.iv.continuousW, aes(x=tau, fill=significance)) +
                  geom_histogram(color="black",  bins = 50, position = 'identity', alpha=0.6, key_glyph = "polygon") + 
                  scale_fill_manual(name=TeX("$\\tau$ type"), values=c("Origin"="lightblue", "P.val < 0.1"="orange", "P.val < 0.05"="red"), labels=c("Origin"=TeX("$\\tau$"), "P.val < 0.1"=TeX("$\\tau$ with P.val < 0.1"), "P.val < 0.05"=TeX("$\\tau$ with P.val < 0.05"))) +
                  facet_grid(rows = vars(variable), scales = "free", labeller = label_parsed) + ggtitle(traits_name[t]) + 
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

      # taus.iv.plot.continuousW <- c(iv.tau.PRSice2, iv.tau.LDPred2)
      # var.iv.plot.continuousW <- c(rep("PRSice2", length(iv.tau.PRSice2)), rep("LDPred2", length(iv.tau.LDPred2))) 
      # sig.iv.plot.continuousW <- c(rep("Origin", length(iv.tau.PRSice2)), rep("Origin", length(iv.tau.LDPred2))) 

      # plot.taus.iv.continuousW <- data.frame("tau" = taus.iv.plot.continuousW, "variable" = var.iv.plot.continuousW, "significance" = sig.iv.plot.continuousW)
      # plot.taus.iv.continuousW$significance <- factor(plot.taus.iv.continuousW$significance, levels = c("Origin"))
      # plot.taus.iv.continuousW$variable <- factor(plot.taus.iv.continuousW$variable, levels = c("PRSice2", "LDPred2"))

      # p.iv <- ggplot(plot.taus.iv.continuousW, aes(x=tau, fill=significance)) +
      #             geom_histogram(color="black",  bins = 50, position = 'identity', alpha=0.6, key_glyph = "polygon") + 
      #             scale_fill_manual(name=TeX("$\\tau$ type"), values=c("Origin"="lightblue"), labels=c("Origin"=TeX("$\\tau$"))) +
      #             facet_grid(rows = vars(variable), scales = "free", labeller = label_parsed) + ggtitle(traits_name[t]) + 
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

      # p.iv
      p.iv.nolegend <- p.iv + guides(fill = "none")
      assign(paste0("p.continuousW.", traits[t], ".", d), p.iv.nolegend)

    }
  }
}

cor.test(iv.tau.PRSice2, iv.tau.LDPred2)
p1 <- (p.continuousW.LDL.CAD | p.continuousW.TC.CAD | p.continuousW.IGF1.CAD)
p1_title <- wrap_elements(p1 + plot_annotation("Coronary Artery Disease\n Continuous Treatment", theme=theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 35)))) 

histogram_p <- (p1_title) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 45), legend.position = "right")
# ggsave(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/extra_analysis_for_paper_revision/figure_supp_fig1_diffPRS_", ml, ".png"), histogram_p, dpi=300, 
#         width = 30, height = 20, units = "cm",scale = 3, limitsize = FALSE)
ggsave(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/extra_analysis_for_paper_revision/figure_supp_fig1_diffPRS_", ml, "_withPvalue.png"), histogram_p, dpi=300, 
        width = 30, height = 20, units = "cm",scale = 3, limitsize = FALSE)
# break
