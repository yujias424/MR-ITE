#' This code is to plot the figure for IGF1 results, incorporating both binary and continuous traits.
#' 
#' @author Shi Yujia
#' @date 2024.08.22

suppressPackageStartupMessages({
    library(patchwork)
    library(ggplot2)
    library(stringr)
})

# Load beeswarm plot
load("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/extra_analysis_for_paper_revision/beeswarm.RData")

# load IGF1 continuous W plots
load("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/extra_analysis_for_paper_revision/figure3_variable_importance_IGF1_CAD_continuousW_boxplot.RData")

for (i in ls()){

    if (str_detect(i, "p.shap.feature")){
        assign(paste0("continousW.IGF1.", i), get(i))
    }

    if (str_detect(i, "p.tau.feature")){
        assign(paste0("continousW.IGF1.", i), get(i))
    }

}

design <- "ABC
           AGH
           DEF
           IJK"

final_plot <- ggobject.IGF1.CAD.continuousW + 
              continousW.IGF1.p.shap.feature.1 + continousW.IGF1.p.shap.feature.2 + continousW.IGF1.p.shap.feature.3 + continousW.IGF1.p.shap.feature.4 + continousW.IGF1.p.shap.feature.5 + 
              continousW.IGF1.p.tau.feature.1 + continousW.IGF1.p.tau.feature.2 + continousW.IGF1.p.tau.feature.3 + continousW.IGF1.p.tau.feature.4 + continousW.IGF1.p.tau.feature.5 + 
              plot_layout(design = design) +
              plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 30)) 
              
ggsave("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/extra_analysis_for_paper_revision/figure5_IGF1.png", final_plot, dpi=300, 
        width = 17, height = 15, units = "cm", scale = 10, limitsize = FALSE)

# body fat; BMI; SBP; Testosterone