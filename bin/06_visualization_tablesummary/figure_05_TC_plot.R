#' This code is to plot the figure for TC results, incorporating both binary and continuous traits.
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

# load TC continuous W plots
load("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/figure3_variable_importance_TC_CAD_continuousW_boxplot.RData")

for (i in ls()){

    if (str_detect(i, "p.shap.feature")){
        assign(paste0("continousW.TC.", i), get(i))
    }

    if (str_detect(i, "p.tau.feature")){
        assign(paste0("continousW.TC.", i), get(i))
    }

}

# load TC binary W plots
load("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/figure3_variable_importance_TC_CAD_binaryW_boxplot.RData")

for (i in ls()){

    if (startsWith(i, "p.shap.feature")){
        assign(paste0("binaryW.TC.", i), get(i))
    }

    if (startsWith(i, "p.tau.feature")){
        assign(paste0("binaryW.TC.", i), get(i))
    }

}

# design <- "ABCDEF
#            AGHIJK
#            LMNOPQ
#            LRSTUV"

# final_plot <- ggobject.TC.CAD.continuousW + 
#               continousW.TC.p.shap.feature.1 + continousW.TC.p.shap.feature.2 + continousW.TC.p.shap.feature.3 + continousW.TC.p.shap.feature.4 + continousW.TC.p.shap.feature.5 + 
#               continousW.TC.p.tau.feature.1 + continousW.TC.p.tau.feature.2 + continousW.TC.p.tau.feature.3 + continousW.TC.p.tau.feature.4 + continousW.TC.p.tau.feature.5 + 
#               ggobject.TC.CAD.binaryW + 
#               binaryW.TC.p.shap.feature.1 + binaryW.TC.p.shap.feature.2 + binaryW.TC.p.shap.feature.3 + binaryW.TC.p.shap.feature.4 + binaryW.TC.p.shap.feature.5 + 
#               binaryW.TC.p.tau.feature.1 + binaryW.TC.p.tau.feature.2 + binaryW.TC.p.tau.feature.3 + binaryW.TC.p.tau.feature.4 + binaryW.TC.p.tau.feature.5 + 
#               plot_layout(design = design) +
#               plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 30)) 
              
# ggsave("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/figure5_TC.png", final_plot, dpi=300, 
#         width = 17, height = 15, units = "cm", scale = 10, limitsize = FALSE)

design <- "ABCDE
           FGHIJ
           KLMNO
           PQRST"

final_plot <- continousW.TC.p.shap.feature.6 + continousW.TC.p.shap.feature.7 + continousW.TC.p.shap.feature.8 + continousW.TC.p.shap.feature.9 + continousW.TC.p.shap.feature.10 + 
              continousW.TC.p.tau.feature.6 + continousW.TC.p.tau.feature.7 + continousW.TC.p.tau.feature.8 + continousW.TC.p.tau.feature.9 + continousW.TC.p.tau.feature.10 +
              binaryW.TC.p.shap.feature.6 + binaryW.TC.p.shap.feature.7 + binaryW.TC.p.shap.feature.8 + binaryW.TC.p.shap.feature.9 + binaryW.TC.p.shap.feature.10 + 
              binaryW.TC.p.tau.feature.6 + binaryW.TC.p.tau.feature.7 + binaryW.TC.p.tau.feature.8 + binaryW.TC.p.tau.feature.9 + binaryW.TC.p.tau.feature.10 + 
              plot_layout(design = design) +
              plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 30)) 
ggsave("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/figure5_TC_supp.png", final_plot, dpi=300, 
        width = 17, height = 15, units = "cm", scale = 10, limitsize = FALSE)