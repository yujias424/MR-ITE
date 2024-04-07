#' This is the code to merge the variable importance plot from python.
#' 
#' @author Yujia Shi
#' @date 2022.11.02

suppressPackageStartupMessages({
  library(gridExtra)
  library(ggplot2)
  library(grid)
  library(gridGraphics)
  library(ggplotify)
  library(cowplot)
})

# Without TG
file_name <- c("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/LDL/heart_disease/variable_importance/binaryWZ/driv/beeswarm_shap.png",
               "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/LDL/ischemic_stroke/variable_importance/binaryWZ/driv/beeswarm_shap.png",
               "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/HDL/heart_disease/variable_importance/binaryWZ/driv_trainmodel/beeswarm_shap.png",
               "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/HDL/ischemic_stroke/variable_importance/binaryWZ/driv_trainmodel/beeswarm_shap.png",
               "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TC/heart_disease/variable_importance/binaryWZ/driv/beeswarm_shap.png",
               "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TC/ischemic_stroke/variable_importance/binaryWZ/driv/beeswarm_shap.png")

taus_file <- grep("tau_*", file_name, ignore.case = TRUE, value = TRUE)
interaction_file <-  grep("_interact.png", file_name, ignore.case = TRUE, value = TRUE)
origin_file <- setdiff(setdiff(file_name, interaction_file), taus_file)

rl = lapply(sprintf(c(origin_file, taus_file)), png::readPNG)
gl = lapply(rl, grid::rasterGrob)
gp <- lapply(gl, as.ggplot)

merged_plot <- plot_grid(gp[[1]], gp[[2]], gp[[3]], gp[[4]], gp[[5]],
                         gp[[6]],
                         ncol = 2, labels = "AUTO", label_size = 20, vjust = 1, label_fontface = "plain")
merged_plot

ggsave(plot = merged_plot, filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/variable_importance/beeswarm_overall_withoutTG.png",
       dpi = 300, height = 18, width = 20, bg = 'white')

# With TG
file_name <- c("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/LDL/heart_disease/variable_importance/binaryWZ/driv/beeswarm_shap_withoutPC.png",
               "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/LDL/ischemic_stroke/variable_importance/binaryWZ/driv/beeswarm_shap_withoutPC.png",
               "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/HDL/heart_disease/variable_importance/binaryWZ/driv_trainmodel/beeswarm_shap_withoutPC.png",
               "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/HDL/ischemic_stroke/variable_importance/binaryWZ/driv_trainmodel/beeswarm_shap_withoutPC.png",
               "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TC/heart_disease/variable_importance/binaryWZ/driv/beeswarm_shap_withoutPC.png",
               "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TC/ischemic_stroke/variable_importance/binaryWZ/driv/beeswarm_shap_withoutPC.png",
               "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TG/heart_disease/variable_importance/binaryWZ/driv_trainmodel/beeswarm_shap_withoutPC.png",
               "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TG/ischemic_stroke/variable_importance/binaryWZ/driv/beeswarm_shap_withoutPC.png")

taus_file <- grep("tau_*", file_name, ignore.case = TRUE, value = TRUE)
interaction_file <-  grep("_interact.png", file_name, ignore.case = TRUE, value = TRUE)
origin_file <- setdiff(setdiff(file_name, interaction_file), taus_file)

rl = lapply(sprintf(c(origin_file, taus_file)), png::readPNG)
gl = lapply(rl, grid::rasterGrob)
gp <- lapply(gl, as.ggplot)

gp[[1]] <- gp[[1]] + ggtitle("Coronary Artery Disease\nLDL-C") + theme(plot.title = element_text(hjust = 0.5, size = 20))
gp[[2]] <- gp[[2]] + ggtitle("Ischemic Stroke\nLDL-C") + theme(plot.title = element_text(hjust = 0.5, size = 20))
gp[[3]] <- gp[[3]] + ggtitle("HDL-C") + theme(plot.title = element_text(hjust = 0.5, size = 20))
gp[[4]] <- gp[[4]] + ggtitle("HDL-C") + theme(plot.title = element_text(hjust = 0.5, size = 20))
gp[[5]] <- gp[[5]] + ggtitle("Total Cholesterol") + theme(plot.title = element_text(hjust = 0.5, size = 20))
gp[[6]] <- gp[[6]] + ggtitle("Total Cholesterol") + theme(plot.title = element_text(hjust = 0.5, size = 20))
gp[[7]] <- gp[[7]] + ggtitle("Triglycerides") + theme(plot.title = element_text(hjust = 0.5, size = 20))
gp[[8]] <- gp[[8]] + ggtitle("Triglycerides") + theme(plot.title = element_text(hjust = 0.5, size = 20))

merged_plot <- plot_grid(gp[[1]], gp[[2]], gp[[3]], gp[[4]], gp[[5]],
                         gp[[6]], gp[[7]], gp[[8]],
                         ncol = 2, labels = "AUTO", label_size = 20, vjust = 1, label_fontface = "plain")
# merged_plot

ggsave(plot = merged_plot, filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/variable_importance/beeswarm_overall_withTG.png",
       dpi = 300, height = 25, width = 25, bg = 'white')
