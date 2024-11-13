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
# traits <- c("LDL", "HDL", "TC", "TG")
# traits_name <- c("LDL-C", "HDL-C", "Total Cholesterol", "Triglycerides")
traits <- c("LDL","TC")
traits_name <- c("LDL-C", "Total Cholesterol")
# models <- c("model1", "model2", "model3")
# models <- c("model1b", "model2b", "model3")
models <- c("model1b", "model3")

for (ml in models){
  for (t in 1:length(traits)){
    for (d in diseases){

      # continuous W
      iv.tau.continuousW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_", ml, ".csv"))
      cf.tau.continuousW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/cf_full_continuousW_pvalue_", ml, ".csv"))
      driv.tau.continuousW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/driv_full_continuousW_te_ul_", ml, ".csv"))

      # binary W
      iv.tau.binaryW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/ivgrf_full_binaryW_continuousZ_pvalue_", ml, ".csv"))
      cf.tau.binaryW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/cf_full_binaryW_continuousZ_pvalue_", ml, ".csv"))
      driv.tau.binaryW <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", d, "/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_", ml, ".csv"))

      # making plot (continuous W)
      iv.tau <- iv.tau.continuousW$tau
      iv.tau.sig <- iv.tau.continuousW$tau[iv.tau.continuousW$tau.pval<0.1]
      iv.tau.sig.05 <- iv.tau.continuousW$tau[iv.tau.continuousW$tau.pval<0.05]
      if (length(iv.tau.sig) == 0){
        iv.tau.sig <- c(mean(iv.tau))
      } else {
        iv.tau.sig <- iv.tau.sig
      }

      driv.tau <- driv.tau.continuousW$point
      driv.tau.sig <- driv.tau.continuousW$point[driv.tau.continuousW$p_value<0.1]
      driv.tau.sig.05 <- driv.tau.continuousW$point[driv.tau.continuousW$p_value<0.05]
      if (length(driv.tau.sig) == 0){
        driv.tau.sig <- c(mean(driv.tau))
      } else {
        driv.tau.sig <- driv.tau.sig
      }

      cf.tau <- cf.tau.continuousW$tau
      cf.tau.sig <- cf.tau.continuousW$tau[cf.tau.continuousW$tau.pval<0.1]
      cf.tau.sig.05 <- cf.tau.continuousW$tau[cf.tau.continuousW$tau.pval<0.05]
      if (length(cf.tau.sig) == 0){
        cf.tau.sig <- c(mean(cf.tau))
      } else {
        cf.tau.sig <- cf.tau.sig
      }

      taus.iv.plot.continuousW <- c(iv.tau, driv.tau, iv.tau.sig, driv.tau.sig, iv.tau.sig.05, driv.tau.sig.05) 
      taus.cf.plot.continuousW <- c(cf.tau, cf.tau.sig, cf.tau.sig.05) 
      var.iv.plot.continuousW <- c(rep("IV-CF", length(iv.tau)), rep("DRIV", length(driv.tau)), rep("IV-CF", length(iv.tau.sig)), rep("DRIV", length(driv.tau.sig)), rep("IV-CF", length(iv.tau.sig.05)), rep("DRIV", length(driv.tau.sig.05))) 
      var.cf.plot.continuousW <- c(rep("CF", length(cf.tau)), rep("CF", length(cf.tau.sig)), rep("CF", length(cf.tau.sig.05))) 
      sig.iv.plot.continuousW <- c(rep("Origin", length(iv.tau)), rep("Origin", length(driv.tau)), rep("P.val < 0.1", length(iv.tau.sig)), rep("P.val < 0.1", length(driv.tau.sig)), rep("P.val < 0.05", length(iv.tau.sig.05)), rep("P.val < 0.05", length(driv.tau.sig.05))) 
      sig.cf.plot.continuousW <- c(rep("Origin", length(cf.tau)), rep("P.val < 0.1", length(cf.tau.sig)), rep("P.val < 0.05", length(cf.tau.sig.05))) 

      plot.taus.iv.continuousW <- data.frame("tau" = taus.iv.plot.continuousW, "variable" = var.iv.plot.continuousW, "significance" = sig.iv.plot.continuousW)
      plot.taus.iv.continuousW$significance <- factor(plot.taus.iv.continuousW$significance, levels = c("Origin", "P.val < 0.1", "P.val < 0.05"))
      plot.taus.iv.continuousW$variable <- factor(plot.taus.iv.continuousW$variable, levels = c("IV-CF", "DRIV"))
      plot.taus.cf.continuousW <- data.frame("tau" = taus.cf.plot.continuousW, "variable" = var.cf.plot.continuousW, "significance" = sig.cf.plot.continuousW)
      plot.taus.cf.continuousW$significance <- factor(plot.taus.cf.continuousW$significance, levels = c("Origin", "P.val < 0.1", "P.val < 0.05"))
      plot.taus.cf.continuousW$variable <- factor(plot.taus.cf.continuousW$variable, levels = c("CF"))

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
        xlab("Individual Treatment Effect") + ylab("Count") + # TeX("$\\tau$ (CLATE)")
        xlim(min(plot.taus.cf.continuousW$tau, plot.taus.iv.continuousW$tau), max(plot.taus.cf.continuousW$tau, plot.taus.iv.continuousW$tau))
      
      p.cf <- ggplot(plot.taus.cf.continuousW, aes(x=tau, fill=significance)) +
        geom_histogram(color="black", bins = 50, position = 'identity', alpha=0.6) + 
        scale_fill_manual(name=TeX("$\\tau$ type"), values=c("Origin"="lightblue", "P.val < 0.1"="orange", "P.val < 0.05"="red"), labels=c("Origin"=TeX("$\\tau$"), "P.val < 0.1"=TeX("$\\tau$ with P.val < 0.1"), "P.val < 0.05"=TeX("$\\tau$ with P.val < 0.05"))) +
        facet_grid(rows = vars(variable), scales = "free", labeller = label_parsed) +
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
        xlab("Individual Treatment Effect") + ylab("Count") + # xlim(layer_scales(p.iv)$x$range$range[1], layer_scales(p.iv)$x$range$range[2]) # TeX("$\\tau$ (CATE)")
        xlim(min(plot.taus.cf.continuousW$tau, plot.taus.iv.continuousW$tau), max(plot.taus.cf.continuousW$tau, plot.taus.iv.continuousW$tau))

      p.iv.nolegend <- p.iv + guides(fill = "none")
      assign(paste0("p.continuousW.", traits[t], ".", d), (p.iv.nolegend/p.cf) + plot_layout(heights = c(8,4), guides = "collect") & theme(legend.position = "right"))

      # making plot (binary W)
      iv.tau <- iv.tau.binaryW$tau
      iv.tau.sig <- iv.tau.binaryW$tau[iv.tau.binaryW$tau.pval<0.1]
      iv.tau.sig.05 <- iv.tau.binaryW$tau[iv.tau.binaryW$tau.pval<0.05]
      if (length(iv.tau.sig) == 0){
        iv.tau.sig <- c(mean(iv.tau))
      } else {
        iv.tau.sig <- iv.tau.sig
      }

      driv.tau <- driv.tau.binaryW$point
      driv.tau.sig <- driv.tau.binaryW$point[driv.tau.binaryW$p_value<0.1]
      driv.tau.sig.05 <- driv.tau.binaryW$point[driv.tau.binaryW$p_value<0.05]
      if (length(driv.tau.sig) == 0){
        driv.tau.sig <- c(mean(driv.tau))
      } else {
        driv.tau.sig <- driv.tau.sig
      }

      cf.tau <- cf.tau.binaryW$tau
      cf.tau.sig <- cf.tau.binaryW$tau[cf.tau.binaryW$tau.pval<0.1]
      cf.tau.sig.05 <- cf.tau.binaryW$tau[cf.tau.binaryW$tau.pval<0.05]
      if (length(cf.tau.sig) == 0){
        cf.tau.sig <- c(mean(cf.tau))
      } else {
        cf.tau.sig <- cf.tau.sig
      }

      taus.iv.plot.binaryW <- c(iv.tau, driv.tau, iv.tau.sig, driv.tau.sig, iv.tau.sig.05, driv.tau.sig.05) 
      taus.cf.plot.binaryW <- c(cf.tau, cf.tau.sig, cf.tau.sig.05) 
      var.iv.plot.binaryW <- c(rep("IV-CF", length(iv.tau)), rep("DRIV", length(driv.tau)), rep("IV-CF", length(iv.tau.sig)), rep("DRIV", length(driv.tau.sig)), rep("IV-CF", length(iv.tau.sig.05)), rep("DRIV", length(driv.tau.sig.05))) 
      var.cf.plot.binaryW <- c(rep("CF", length(cf.tau)), rep("CF", length(cf.tau.sig)), rep("CF", length(cf.tau.sig.05))) 
      sig.iv.plot.binaryW <- c(rep("Origin", length(iv.tau)), rep("Origin", length(driv.tau)), rep("P.val < 0.1", length(iv.tau.sig)), rep("P.val < 0.1", length(driv.tau.sig)), rep("P.val < 0.05", length(iv.tau.sig.05)), rep("P.val < 0.05", length(driv.tau.sig.05))) 
      sig.cf.plot.binaryW <- c(rep("Origin", length(cf.tau)), rep("P.val < 0.1", length(cf.tau.sig)), rep("P.val < 0.05", length(cf.tau.sig.05))) 

      plot.taus.iv.binaryW <- data.frame("tau" = taus.iv.plot.binaryW, "variable" = var.iv.plot.binaryW, "significance" = sig.iv.plot.binaryW)
      plot.taus.iv.binaryW$significance <- factor(plot.taus.iv.binaryW$significance, levels = c("Origin", "P.val < 0.1", "P.val < 0.05"))
      plot.taus.iv.binaryW$variable <- factor(plot.taus.iv.binaryW$variable, levels = c("IV-CF", "DRIV"))
      plot.taus.cf.binaryW <- data.frame("tau" = taus.cf.plot.binaryW, "variable" = var.cf.plot.binaryW, "significance" = sig.cf.plot.binaryW)
      plot.taus.cf.binaryW$significance <- factor(plot.taus.cf.binaryW$significance, levels = c("Origin", "P.val < 0.1", "P.val < 0.05"))
      plot.taus.cf.binaryW$variable <- factor(plot.taus.cf.binaryW$variable, levels = c("CF"))

      p.iv <- ggplot(plot.taus.iv.binaryW, aes(x=tau, fill=significance)) +
        geom_histogram(color="black",  bins = 50, position = 'identity', alpha=0.6) + 
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
        xlab("Individual Treatment Effect") + ylab("Count") + # TeX("$\\tau$ (CLATE)")
        xlim(min(plot.taus.cf.binaryW$tau, plot.taus.iv.binaryW$tau), max(plot.taus.cf.binaryW$tau, plot.taus.iv.binaryW$tau))

      p.cf <- ggplot(plot.taus.cf.binaryW, aes(x=tau, fill=significance)) +
        geom_histogram(color="black", bins = 50, position = 'identity', alpha=0.6) + 
        scale_fill_manual(name=TeX("$\\tau$ type"), values=c("Origin"="lightblue", "P.val < 0.1"="orange", "P.val < 0.05"="red"), labels=c("Origin"=TeX("$\\tau$"), "P.val < 0.1"=TeX("$\\tau$ with P.val < 0.1"), "P.val < 0.05"=TeX("$\\tau$ with P.val < 0.05"))) +
        facet_grid(rows = vars(variable), scales = "free", labeller = label_parsed) +
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
        xlab("Individual Treatment Effect") + ylab("Count") + # xlim(layer_scales(p.iv)$x$range$range[1], layer_scales(p.iv)$x$range$range[2]) # TeX("$\\tau$ (CATE)")
        xlim(min(plot.taus.cf.binaryW$tau, plot.taus.iv.binaryW$tau), max(plot.taus.cf.binaryW$tau, plot.taus.iv.binaryW$tau))

      p.iv.nolegend <- p.iv + guides(fill = "none")
      assign(paste0("p.binaryW.", traits[t], ".", d), (p.iv.nolegend/p.cf) + plot_layout(heights = c(8,4), guides = "collect") & theme(legend.position = "right")) 

    }
  }

  p1 <- (p.continuousW.LDL.CAD | p.continuousW.TC.CAD) 
  p2 <- (p.binaryW.LDL.CAD | p.binaryW.TC.CAD)
  p1_title <- wrap_elements(p1 + plot_annotation("Coronary Artery Disease\n Continuous Treatment", theme=theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 35)))) 
  p2_title <- wrap_elements(p2 + plot_annotation("Binary Treatment", theme=theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 35))))

  histogram_p <- (p1_title / p2_title) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 45), legend.position = "right")
  ggsave(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/figure1_normal_pvalue_", ml, "_james.png"), histogram_p, dpi=300, 
          width = 50, height = 40, units = "cm",scale = 3, limitsize = FALSE)
  # break
}
