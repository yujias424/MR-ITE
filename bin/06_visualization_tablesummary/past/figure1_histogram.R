#' This code is to merge the tau summary plot into a single plot.
#' 
#' @author Yujia Shi
#' @date 2022.10.02

# CAD
load("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/HDL/heart_disease/taus_summary/HDL_taus_plot.RData")
load("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/LDL/heart_disease/taus_summary/LDL_taus_plot.RData")
load("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TC/heart_disease/taus_summary/TC_taus_plot.RData")
load("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TG/heart_disease/taus_summary/TG_taus_plot.RData")

tmp_merge <- (p.continuousZW.LDL.cad | p.continuousZW.HDL.cad | p.continuousZW.TC.cad | p.continuousZW.TG.cad) /
              (p.binaryZW.LDL.cad | p.binaryZW.HDL.cad | p.binaryZW.TC.cad | p.binaryZW.TG.cad) + plot_annotation(tag_levels = list(c("A", "", "", "", "", "", "", "", "B", "", "", "", "", "", "", "")))
tmp_merge

ggsave(filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/heart_disease/CAD_taus_summary_ver5.png",
       plot = tmp_merge,
       dpi = "retina", width = 32, height = 18)

# IS
load("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/HDL/ischemic_stroke/taus_summary/HDL_taus_plot.RData")
load("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/LDL/ischemic_stroke/taus_summary/LDL_taus_plot.RData")
load("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TC/ischemic_stroke/taus_summary/TC_taus_plot.RData")
load("~/Project/2022-09-01-individual_MR/res/plot/Manuscript/TG/ischemic_stroke/taus_summary/TG_taus_plot.RData")

tmp_merge <- (p.continuousZW.LDL.is | p.continuousZW.HDL.is | p.continuousZW.TC.is | p.continuousZW.TG.is) /
  (p.binaryZW.LDL.is | p.binaryZW.HDL.is | p.binaryZW.TC.is | p.binaryZW.TG.is) + plot_annotation(tag_levels = list(c("A", "", "", "", "", "", "", "", "B", "", "", "", "", "", "", "")))
tmp_merge

ggsave(filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/ischemic_stroke/IS_taus_summary_ver4.png",
       plot = tmp_merge,
       dpi = "retina", width = 32, height = 18)

tmp_merge <- (p.continuousZW.LDL.cad | p.continuousZW.HDL.cad | p.continuousZW.TC.cad | p.continuousZW.TG.cad) /
  (p.binaryZW.LDL.cad | p.binaryZW.HDL.cad | p.binaryZW.TC.cad | p.binaryZW.TG.cad) / 
  (p.continuousZW.LDL.is | p.continuousZW.HDL.is | p.continuousZW.TC.is | p.continuousZW.TG.is) /
  (p.binaryZW.LDL.is | p.binaryZW.HDL.is | p.binaryZW.TC.is | p.binaryZW.TG.is) + plot_annotation(tag_levels = list(c("A", "", "", "", "", "", "", "", "B", "", "", "", "", "", "", "", "C", "", "", "", "", "", "", "", "D", "", "", "", "", "", "", "")))
tmp_merge

library(ggpubr)
figure.sub1 <- ggarrange(p.continuousZW.LDL.cad, p.continuousZW.HDL.cad, p.continuousZW.TC.cad, p.continuousZW.TG.cad,
                         ncol = 4, nrow = 1)
figure.sub1 <- annotate_figure(figure.sub1, top = text_grob("Coronary Artery Disease\n Continuous Treatment", 
                                      color = "black",  size = 25))

figure.sub2 <- ggarrange(p.binaryZW.LDL.cad, p.binaryZW.HDL.cad, p.binaryZW.TC.cad, p.binaryZW.TG.cad,
                         ncol = 4, nrow = 1)
figure.sub2 <- annotate_figure(figure.sub2, top = text_grob("Binary Treatment", 
                                                            color = "black",  size = 25))

figure.sub3 <- ggarrange(p.continuousZW.LDL.is, p.continuousZW.HDL.is, p.continuousZW.TC.is, p.continuousZW.TG.is,
                         ncol = 4, nrow = 1)
figure.sub3 <- annotate_figure(figure.sub3, top = text_grob("Ischemic Stroke\n Continuous Treatment", 
                                                            color = "black",  size = 25))

figure.sub4 <- ggarrange(p.binaryZW.LDL.is, p.binaryZW.HDL.is, p.binaryZW.TC.is, p.binaryZW.TG.is,
                         ncol = 4, nrow = 1)
figure.sub4 <- annotate_figure(figure.sub4, top = text_grob("Binary Treatment", 
                                                            color = "black",  size = 25))

tmp_merge <- figure.sub1 /
             figure.sub2 / 
             figure.sub3 /
             figure.sub4 + plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) & theme(plot.tag = element_text(size = 25, face = 'plain'))
tmp_merge

ggsave(filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/withTG/Figure2_withTG.png",
       plot = tmp_merge,
       dpi = "retina", width = 32, height = 40)

# ==========================================================================
library(ggpubr)
figure.sub1 <- ggarrange(p.continuousZW.LDL.cad, p.continuousZW.HDL.cad, p.continuousZW.TC.cad, 
                         ncol = 3, nrow = 1)
figure.sub1 <- annotate_figure(figure.sub1, top = text_grob("Coronary Artery Disease\n Continuous Treatment", 
                                                            color = "black",  size = 25))

figure.sub2 <- ggarrange(p.binaryZW.LDL.cad, p.binaryZW.HDL.cad, p.binaryZW.TC.cad, 
                         ncol = 3, nrow = 1)
figure.sub2 <- annotate_figure(figure.sub2, top = text_grob("Binary Treatment", 
                                                            color = "black",  size = 25))

figure.sub3 <- ggarrange(p.continuousZW.LDL.is, p.continuousZW.HDL.is, p.continuousZW.TC.is, 
                         ncol = 3, nrow = 1)
figure.sub3 <- annotate_figure(figure.sub3, top = text_grob("Ischemic Stroke\n Continuous Treatment", 
                                                            color = "black",  size = 25))

figure.sub4 <- ggarrange(p.binaryZW.LDL.is, p.binaryZW.HDL.is, p.binaryZW.TC.is,
                         ncol = 3, nrow = 1)
figure.sub4 <- annotate_figure(figure.sub4, top = text_grob("Binary Treatment", 
                                                            color = "black",  size = 25))

tmp_merge <- figure.sub1 /
  figure.sub2 / 
  figure.sub3 /
  figure.sub4 + plot_annotation(tag_levels = list(c("A", "B", "C", "D"))) & theme(plot.tag = element_text(size = 25, face = 'plain'))
tmp_merge

ggsave(filename = "~/Project/2022-09-01-individual_MR/res/plot/Manuscript/Figure/Figure2_withoutTG.png",
       plot = tmp_merge,
       dpi = "retina", width = 32, height = 40)
