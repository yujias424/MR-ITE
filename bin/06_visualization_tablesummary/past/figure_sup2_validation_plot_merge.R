#' This code is to merge the validation plot and create supplement figure 2.
#' 
#' @author Shi Yujia
#' @date 2022.12.10

suppressPackageStartupMessages({
  library(gridExtra)
  library(ggplot2)
  library(grid)
  library(gridGraphics)
  library(ggplotify)
  library(cowplot)
  library(rsvg) # this one is not installed on the server, thus the code needs to be ran locally.
})

p1 <- png::readPNG("../Downloads/binaryLDL_CAD_boxplot.png")
p1_gl <- as.ggplot(grid::rasterGrob(p1))
p1_gl <- p1_gl # + ggtitle("LDL-CAD model") + theme(plot.title = element_text(hjust = 0.5, size = 25))

rsvg_png(
  "../Downloads/binaryLDL_CAD_besttree.svg", "../Downloads/binaryLDL_CAD_besttree.png",
  width = 1800, height = 1800,
)

p2 <- png::readPNG("../Downloads/binaryLDL_CAD_besttree.png")
p2_gl <- as.ggplot(grid::rasterGrob(p2))
p2_gl <- p2_gl + ggtitle("LDL-CAD model") + theme(plot.title = element_text(hjust = 0.5, size = 25))

p3 <- png::readPNG("../Downloads/binaryHDL_IS_boxplot.png")
p3_gl <- as.ggplot(grid::rasterGrob(p3))
p3_gl <- p3_gl # + ggtitle("HDL-IS model") + theme(plot.title = element_text(hjust = 0.5, size = 25))

rsvg_png(
  "../Downloads/binaryHDL_IS_besttree.svg", "../Downloads/binaryHDL_IS_besttree.png",
  width = 1800, height = 1800,
)

p4 <- png::readPNG("../Downloads/binaryHDL_IS_besttree.png")
p4_gl <- as.ggplot(grid::rasterGrob(p4))
p4_gl <- p4_gl + ggtitle("HDL-IS model") + theme(plot.title = element_text(hjust = 0.5, size = 25))

merged_plot <- plot_grid(p2_gl, p1_gl, p4_gl, p3_gl, 
                         nrow = 2, ncol = 2,
                         labels = "AUTO", label_size = 20, 
                         vjust = 1, label_fontface = "plain")
merged_plot

ggsave(plot = merged_plot, filename = "../Downloads/supp_figure4.png",
       dpi = 300, height = 25, width = 25, bg = 'white')
