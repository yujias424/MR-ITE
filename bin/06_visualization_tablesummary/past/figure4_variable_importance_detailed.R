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

#' =========== LDL =========== 

#' ===========
#' CAD
#' ===========

file_name <- list.files("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/LDL/heart_disease/variable_importance/binaryWZ/driv/variable/", include.dirs = T)

for (i in 1:length(file_name)){
  file_name[i] <- paste0("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/LDL/heart_disease/variable_importance/binaryWZ/driv/variable/", file_name[i])
}

taus_file <- grep("tau_*", file_name, ignore.case = TRUE, value = TRUE)
interaction_file <-  grep("_interact.png", file_name, ignore.case = TRUE, value = TRUE)
origin_file <- setdiff(setdiff(file_name, interaction_file), taus_file)

rl = lapply(sprintf(c(origin_file, taus_file)), png::readPNG)
gl = lapply(rl, grid::rasterGrob)
gp <- lapply(gl, as.ggplot)

merged_plot <- plot_grid(gp[[1]], gp[[2]], gp[[3]], gp[[4]], gp[[5]],
                         gp[[6]], gp[[7]], gp[[8]], gp[[9]], gp[[10]],
                         gp[[11]], gp[[12]], gp[[13]], gp[[14]], gp[[15]],
                         gp[[16]], gp[[17]], gp[[18]], gp[[19]], gp[[20]],
                         ncol = 5, labels = "AUTO", label_size = 20, vjust = 1, label_fontface = "plain")

ggsave(plot = merged_plot, filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/variable_importance/LDL_CAD_binary.png",
       dpi = 300, height = 20, width = 20, bg = 'white')

#' ===========
#' IS
#' ===========

file_name <- list.files("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/LDL/ischemic_stroke/variable_importance/binaryWZ/driv/variable/", include.dirs = T)

for (i in 1:length(file_name)){
  file_name[i] <- paste0("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/LDL/ischemic_stroke/variable_importance/binaryWZ/driv/variable/", file_name[i])
}

taus_file <- grep("tau_*", file_name, ignore.case = TRUE, value = TRUE)
interaction_file <-  grep("_interact.png", file_name, ignore.case = TRUE, value = TRUE)
origin_file <- setdiff(setdiff(file_name, interaction_file), taus_file)

rl = lapply(sprintf(c(origin_file, taus_file)), png::readPNG)
gl = lapply(rl, grid::rasterGrob)
gp <- lapply(gl, as.ggplot)

merged_plot <- plot_grid(gp[[1]], gp[[2]], gp[[3]], gp[[4]], gp[[5]],
                         gp[[6]], gp[[7]], gp[[8]], gp[[9]], gp[[10]],
                         gp[[11]], gp[[12]], gp[[13]], gp[[14]], gp[[15]],
                         gp[[16]], gp[[17]], gp[[18]], gp[[19]], gp[[20]],
                         ncol = 5, labels = "AUTO", label_size = 20, vjust = 1, label_fontface = "plain")

ggsave(plot = merged_plot, filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/variable_importance/LDL_IS_binary.png",
       dpi = 300, height = 20, width = 20, bg = 'white')

#' =========== HDL =========== 

#' ===========
#' CAD
#' ===========

file_name <- list.files("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/HDL/heart_disease/variable_importance/binaryWZ/driv_trainmodel/variable/", include.dirs = T)

for (i in 1:length(file_name)){
  file_name[i] <- paste0("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/HDL/heart_disease/variable_importance/binaryWZ/driv_trainmodel/variable/", file_name[i])
}

taus_file <- grep("tau_*", file_name, ignore.case = TRUE, value = TRUE)
interaction_file <-  grep("_interact.png", file_name, ignore.case = TRUE, value = TRUE)
origin_file <- setdiff(setdiff(file_name, interaction_file), taus_file)

rl = lapply(sprintf(c(origin_file, taus_file)), png::readPNG)
gl = lapply(rl, grid::rasterGrob)
gp <- lapply(gl, as.ggplot)

merged_plot <- plot_grid(gp[[1]], gp[[2]], gp[[3]], gp[[4]], gp[[5]],
                         gp[[6]], gp[[7]], gp[[8]], gp[[9]], gp[[10]],
                         gp[[11]], gp[[12]], gp[[13]], gp[[14]], gp[[15]],
                         gp[[16]], gp[[17]], gp[[18]], gp[[19]], gp[[20]],
                         ncol = 5, labels = "AUTO", label_size = 20, vjust = 1, label_fontface = "plain")

ggsave(plot = merged_plot, filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/variable_importance/HDL_CAD_binary.png",
       dpi = 300, height = 20, width = 20, bg = 'white')

#' ===========
#' IS
#' ===========

file_name <- list.files("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/HDL/ischemic_stroke/variable_importance/binaryWZ/driv_trainmodel/variable/", include.dirs = T)

for (i in 1:length(file_name)){
  file_name[i] <- paste0("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/HDL/ischemic_stroke/variable_importance/binaryWZ/driv_trainmodel/variable/", file_name[i])
}

taus_file <- grep("tau_*", file_name, ignore.case = TRUE, value = TRUE)
interaction_file <-  grep("_interact.png", file_name, ignore.case = TRUE, value = TRUE)
origin_file <- setdiff(setdiff(file_name, interaction_file), taus_file)

rl = lapply(sprintf(c(origin_file, taus_file)), png::readPNG)
gl = lapply(rl, grid::rasterGrob)
gp <- lapply(gl, as.ggplot)

merged_plot <- plot_grid(gp[[1]], gp[[2]], gp[[3]], gp[[4]], gp[[5]],
                         gp[[6]], gp[[7]], gp[[8]], gp[[9]], gp[[10]],
                         gp[[11]], gp[[12]], gp[[13]], gp[[14]], gp[[15]],
                         gp[[16]], gp[[17]], gp[[18]], gp[[19]], gp[[20]],
                         ncol = 5, labels = "AUTO", label_size = 20, vjust = 1, label_fontface = "plain")

ggsave(plot = merged_plot, filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/variable_importance/HDL_IS_binary.png",
       dpi = 300, height = 20, width = 20, bg = 'white')

#' =========== TC =========== 

#' ===========
#' CAD
#' ===========

file_name <- list.files("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TC/heart_disease/variable_importance/binaryWZ/driv/variable/", include.dirs = T)

for (i in 1:length(file_name)){
  file_name[i] <- paste0("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TC/heart_disease/variable_importance/binaryWZ/driv/variable/", file_name[i])
}

taus_file <- grep("tau_*", file_name, ignore.case = TRUE, value = TRUE)
interaction_file <-  grep("_interact.png", file_name, ignore.case = TRUE, value = TRUE)
origin_file <- setdiff(setdiff(file_name, interaction_file), taus_file)

rl = lapply(sprintf(c(origin_file, taus_file)), png::readPNG)
gl = lapply(rl, grid::rasterGrob)
gp <- lapply(gl, as.ggplot)

merged_plot <- plot_grid(gp[[1]], gp[[2]], gp[[3]], gp[[4]], gp[[5]],
                         gp[[6]], gp[[7]], gp[[8]], gp[[9]], gp[[10]],
                         gp[[11]], gp[[12]], gp[[13]], gp[[14]], gp[[15]],
                         gp[[16]], gp[[17]], gp[[18]], gp[[19]], gp[[20]],
                         ncol = 5, labels = "AUTO", label_size = 20, vjust = 1, label_fontface = "plain")

ggsave(plot = merged_plot, filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/variable_importance/TC_CAD_binary.png",
       dpi = 300, height = 20, width = 20, bg = 'white')

#' ===========
#' IS
#' ===========

file_name <- list.files("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TC/ischemic_stroke/variable_importance/binaryWZ/driv/variable/", include.dirs = T)

for (i in 1:length(file_name)){
  file_name[i] <- paste0("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TC/ischemic_stroke/variable_importance/binaryWZ/driv/variable/", file_name[i])
}

taus_file <- grep("tau_*", file_name, ignore.case = TRUE, value = TRUE)
interaction_file <-  grep("_interact.png", file_name, ignore.case = TRUE, value = TRUE)
origin_file <- setdiff(setdiff(file_name, interaction_file), taus_file)

rl = lapply(sprintf(c(origin_file, taus_file)), png::readPNG)
gl = lapply(rl, grid::rasterGrob)
gp <- lapply(gl, as.ggplot)

merged_plot <- plot_grid(gp[[1]], gp[[2]], gp[[3]], gp[[4]], gp[[5]],
                         gp[[6]], gp[[7]], gp[[8]], gp[[9]], gp[[10]],
                         gp[[11]], gp[[12]], gp[[13]], gp[[14]], gp[[15]],
                         gp[[16]], gp[[17]], gp[[18]], gp[[19]], gp[[20]],
                         ncol = 5, labels = "AUTO", label_size = 20, vjust = 1, label_fontface = "plain")

ggsave(plot = merged_plot, filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/variable_importance/TC_IS_binary.png",
       dpi = 300, height = 20, width = 20, bg = 'white')

#' =========== TG =========== 

#' ===========
#' CAD
#' ===========

file_name <- list.files("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TG/heart_disease/variable_importance/binaryWZ/driv_trainmodel/variable/", include.dirs = T)

for (i in 1:length(file_name)){
  file_name[i] <- paste0("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TG/heart_disease/variable_importance/binaryWZ/driv_trainmodel/variable/", file_name[i])
}

taus_file <- grep("tau_*", file_name, ignore.case = TRUE, value = TRUE)
interaction_file <-  grep("_interact.png", file_name, ignore.case = TRUE, value = TRUE)
origin_file <- setdiff(setdiff(file_name, interaction_file), taus_file)

rl = lapply(sprintf(c(origin_file, taus_file)), png::readPNG)
gl = lapply(rl, grid::rasterGrob)
gp <- lapply(gl, as.ggplot)

merged_plot <- plot_grid(gp[[1]], gp[[2]], gp[[3]], gp[[4]], gp[[5]],
                         gp[[6]], gp[[7]], gp[[8]], gp[[9]], gp[[10]],
                         gp[[11]], gp[[12]], gp[[13]], gp[[14]], gp[[15]],
                         gp[[16]], gp[[17]], gp[[18]], gp[[19]], gp[[20]],
                         ncol = 5, labels = "AUTO", label_size = 20, vjust = 1, label_fontface = "plain")

ggsave(plot = merged_plot, filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/variable_importance/TG_CAD_binary.png",
       dpi = 300, height = 20, width = 20, bg = 'white')

#' ===========
#' IS
#' ===========

file_name <- list.files("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TG/ischemic_stroke/variable_importance/binaryWZ/driv/variable/", include.dirs = T)

for (i in 1:length(file_name)){
  file_name[i] <- paste0("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TG/ischemic_stroke/variable_importance/binaryWZ/driv/variable/", file_name[i])
}

taus_file <- grep("tau_*", file_name, ignore.case = TRUE, value = TRUE)
interaction_file <-  grep("_interact.png", file_name, ignore.case = TRUE, value = TRUE)
origin_file <- setdiff(setdiff(file_name, interaction_file), taus_file)

rl = lapply(sprintf(c(origin_file, taus_file)), png::readPNG)
gl = lapply(rl, grid::rasterGrob)
gp <- lapply(gl, as.ggplot)

merged_plot <- plot_grid(gp[[1]], gp[[2]], gp[[3]], gp[[4]], gp[[5]],
                         gp[[6]], gp[[7]], gp[[8]], gp[[9]], gp[[10]],
                         gp[[11]], gp[[12]], gp[[13]], gp[[14]], gp[[15]],
                         gp[[16]], gp[[17]], gp[[18]], gp[[19]], gp[[20]],
                         ncol = 5, labels = "AUTO", label_size = 20, vjust = 1, label_fontface = "plain")

ggsave(plot = merged_plot, filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/variable_importance/TG_IS_binary.png",
       dpi = 300, height = 20, width = 20, bg = 'white')
