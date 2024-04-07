#' This is the code to merge the simultation plot from python.
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

#' =========== merge simulation plot =========== 

file_name <- list.files("~/Project/2021-11-10-individual_MR/res/plot/simulation/", include.dirs = T)
file_name <- file_name[1:9]

for (i in 1:length(file_name)){
  file_name[i] <- paste0("~/Project/2021-11-10-individual_MR/res/plot/simulation/", file_name[i])
}

taus_file <- grep("tau_*", file_name, ignore.case = TRUE, value = TRUE)
interaction_file <-  grep("_interact.png", file_name, ignore.case = TRUE, value = TRUE)
origin_file <- setdiff(setdiff(file_name, interaction_file), taus_file)

rl = lapply(sprintf(c(origin_file, taus_file)), png::readPNG)


rl = lapply(sprintf(c(file_name)), png::readPNG)
gl = lapply(rl, grid::rasterGrob)
gp <- lapply(gl, as.ggplot)

merged_plot <- plot_grid(gp[[1]], gp[[2]], gp[[3]], gp[[4]], gp[[5]],
                         gp[[6]], gp[[7]], gp[[8]], gp[[9]],
                         ncol = 3, labels = "AUTO", label_size = 20, vjust = 1, label_fontface = "plain")

merged_plot

ggsave(plot = merged_plot, filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/figure1_simulation.png",
       dpi = 300, height = 20, width = 20, bg = 'white')
