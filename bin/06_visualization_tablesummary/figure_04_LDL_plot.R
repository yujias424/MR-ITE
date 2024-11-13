#' This code is to plot the figure for LDL results, incorporating both binary and continuous traits.
#' 
#' @author Shi Yujia
#' @date 2024.08.22

suppressPackageStartupMessages({
    library(patchwork)
    library(ggplot2)
    library(stringr)
})

# Load beeswarm plot
load("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/beeswarm.RData")

# load LDL continuous W plots
load("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/figure3_variable_importance_LDL_CAD_continuousW_boxplot.RData")

for (i in ls()){

    if (str_detect(i, "p.shap.feature")){
        assign(paste0("continousW.LDL.", i), get(i))
    }

    if (str_detect(i, "p.tau.feature")){
        assign(paste0("continousW.LDL.", i), get(i))
    }

}

# load LDL binary W plots
load("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/figure3_variable_importance_LDL_CAD_binaryW_boxplot.RData")

for (i in ls()){

    if (startsWith(i, "p.shap.feature")){
        assign(paste0("binaryW.LDL.", i), get(i))
    }

    if (startsWith(i, "p.tau.feature")){
        assign(paste0("binaryW.LDL.", i), get(i))
    }

}

design <- "ABCDEF
           AGHIJK
           LMNOPQ
           LRSTUV"

final_plot <- ggobject.LDL.CAD.continuousW + 
              continousW.LDL.p.shap.feature.1 + continousW.LDL.p.shap.feature.2 + continousW.LDL.p.shap.feature.3 + continousW.LDL.p.shap.feature.4 + continousW.LDL.p.shap.feature.5 + 
              continousW.LDL.p.tau.feature.1 + continousW.LDL.p.tau.feature.2 + continousW.LDL.p.tau.feature.3 + continousW.LDL.p.tau.feature.4 + continousW.LDL.p.tau.feature.5 + 
              ggobject.LDL.CAD.binaryW + 
              binaryW.LDL.p.shap.feature.1 + binaryW.LDL.p.shap.feature.2 + binaryW.LDL.p.shap.feature.3 + binaryW.LDL.p.shap.feature.4 + binaryW.LDL.p.shap.feature.5 + 
              binaryW.LDL.p.tau.feature.1 + binaryW.LDL.p.tau.feature.2 + binaryW.LDL.p.tau.feature.3 + binaryW.LDL.p.tau.feature.4 + binaryW.LDL.p.tau.feature.5 + 
              plot_layout(design = design) +
              plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 45)) 
              
ggsave("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/figure5_LDL.png", final_plot, dpi=300, 
        width = 17, height = 15, units = "cm", scale = 10, limitsize = FALSE)

design <- "ABCDE
           FGHIJ
           KLMNO
           PQRST"

final_plot <- continousW.LDL.p.shap.feature.6 + continousW.LDL.p.shap.feature.7 + continousW.LDL.p.shap.feature.8 + continousW.LDL.p.shap.feature.9 + continousW.LDL.p.shap.feature.10 + 
              continousW.LDL.p.tau.feature.6 + continousW.LDL.p.tau.feature.7 + continousW.LDL.p.tau.feature.8 + continousW.LDL.p.tau.feature.9 + continousW.LDL.p.tau.feature.10 +
              binaryW.LDL.p.shap.feature.6 + binaryW.LDL.p.shap.feature.7 + binaryW.LDL.p.shap.feature.8 + binaryW.LDL.p.shap.feature.9 + binaryW.LDL.p.shap.feature.10 + 
              binaryW.LDL.p.tau.feature.6 + binaryW.LDL.p.tau.feature.7 + binaryW.LDL.p.tau.feature.8 + binaryW.LDL.p.tau.feature.9 + binaryW.LDL.p.tau.feature.10 + 
              plot_layout(design = design) +
              plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 45)) 
ggsave("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/figure5_LDL_supp.png", final_plot, dpi=300, 
        width = 17, height = 15, units = "cm", scale = 10, limitsize = FALSE)
