#' This code is to plot the beeswarm for top SHAP variables.
#' 
#' @author Shi Yujia
#' @date 2023.02.15

suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(patchwork)
    library(ggbeeswarm)
    library(tidyr)
    library(ggeasy)
})

traits <- c("IGF1") 
diseases <- c("CAD")
diseases_name <- c("Coronary Artery Disease") 
traits_type <- c("continuous")
traits_name <- c("Insulin-like growth factor 1") 

for (t in 1:length(traits)){
    for (d in 1:length(diseases)){
        for (wtype in 1:1){
            
            # t <- 1; d <- 1; wtype <- 1

            shap_values <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", diseases[d], "/03_variable_importance/driv_", traits_type[wtype], "W_shap_model3.csv.gz"))
            feature_values <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", diseases[d], "/03_variable_importance/driv_", traits_type[wtype], "W_value_model3.csv.gz"))

            # get the order of top important variables.
            shaporder <- names(sort(sapply(abs(shap_values), mean), decreasing = T))

            # remove the PC and WHR
            shaporder <- shaporder[!startsWith(shaporder, "PC")]
            shaporder <- shaporder[!startsWith(shaporder, "Genotype")]
            print(shaporder)
            # # remove the PC variables from the top order
            # shaporder.noPC.top10 <- shaporder[!startsWith(shaporder, "PC")][1:10]

            # select top 10
            shaporder.noPC.top10 <- shaporder[1:10]

            # standardize the feature value for each feature.
            feature_values_scaled <- scale(feature_values)
            feature_values <- data.table(feature_values_scaled)

            # select top 10 variables
            shap_values.top10 <- shap_values[, ..shaporder.noPC.top10]
            feature_values.top10 <- feature_values[, ..shaporder.noPC.top10]

            shap_values.gathered <- gather(shap_values.top10, variable_name, shap_value)
            feature_values.gathered <- gather(feature_values.top10, variable_name, feature_value)

            newdf <- cbind(shap_values.gathered, feature_values.gathered)
            newdf <- newdf[, c(1,2,4)]
            newdf$variable_name <- factor(newdf$variable_name, levels = rev(shaporder.noPC.top10))

            # we need to sampling the data otherwise it is impossible to make the final plot.
            selected_index <- caret::createDataPartition(newdf$variable_name, p = 0.01, list = FALSE)
            newdf_selected <- newdf[selected_index,]

            assign(paste0("ggobject.", traits[t], ".", diseases[d], ".", traits_type[wtype], "W"), newdf_selected %>% ggplot(aes(shap_value, variable_name, color = feature_value)) +
                                                                        geom_beeswarm(aes(color=feature_value), dodge.width = .1, cex=0.1) +
                                                                        scale_color_gradient2(
                                                                            low = "blue",
                                                                            mid = "white",
                                                                            high = "red",
                                                                            midpoint = 0, 
                                                                            limits=c(-5,5),
                                                                            labels=c("Low","High"),
                                                                            breaks=c(-5,5)
                                                                        ) +
                                                                        geom_vline(xintercept=0, color="black", linetype="dashed") +
                                                                        theme_classic() +
                                                                        xlab("SHAP value (impact on model output)") + 
                                                                        ggtitle(paste0(traits_name[t], "\n", diseases_name[d], "\n", R.utils::capitalize(traits_type[wtype]), " Treatment")) + 
                                                                        theme(
                                                                            plot.title = element_text(hjust = 0.5, size = 35),
                                                                            axis.title.x = element_text(size = 35),
                                                                            axis.title.y = element_blank(),
                                                                            axis.text = element_text(size = 35),
                                                                            legend.title = element_text(angle = -90, size = 30),
                                                                            legend.title.align = 0.5,
                                                                            legend.text = element_text(size = 30)
                                                                        ) + 
                                                                        guides(colour = guide_colourbar(title = "Feature Value", 
                                                                                                        title.position = "right",
                                                                                                        label = T)) + 
                                                                        easy_remove_y_axis(what = c("line"), teach = FALSE)
            )

        }
        # break
    }
    # break
}

save(ggobject.IGF1.CAD.continuousW,
     file = "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/extra_analysis_for_paper_revision/beeswarm.RData")

# p.final <- (ggobject.LDL.CAD.continuousW | ggobject.TC.CAD.continuousW) / (ggobject.LDL.CAD.binaryW | ggobject.TC.CAD.binaryW) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 35))
# p.final <- (ggobject.LDL.CAD.binaryW | ggobject.TC.CAD.binaryW) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 35))
# p.final <- (ggobject.TC.CAD.continuousW | ggobject.TC.CAD.binaryW) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 35))

# ggobject.IGF1.CAD.continuousW

# ggsave('/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/extra_analysis_for_paper_revision/figure2_variable_importance_beeswarm_IGF1.png', ggobject.IGF1.CAD.continuousW, dpi=300, 
#         width = 30, height = 20, units = "cm", scale = 2)
