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

traits <- c("LDL", "TC") 
diseases <- c("CAD")
diseases_name <- c("Coronary Artery Disease") 
traits_type <- c("binary", "continuous")
traits_name <- c("LDL-C", "Total Cholesterol") 

for (t in 1:length(traits)){
    for (d in 1:length(diseases)){
        for (wtype in 1:2){
            
            # t <- 1; d <- 1; wtype <- 1

            shap_values <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", diseases[d], "/03_variable_importance/driv_", traits_type[wtype], "W_shap_model3.csv.gz"))
            feature_values <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", diseases[d], "/03_variable_importance/driv_", traits_type[wtype], "W_value_model3.csv.gz"))

            if (traits[t] == "LDL"){
                shap_values  <- shap_values  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0",   
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0", 
                         "BMI"="21001-0.0", "Body Fat Percentage"="23099-0.0", "Weight"="21002-0.0", # "Waist-hip-ratio"="whr", 
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "Apolipoprotein A"="30630-0.0", "HDL-C"="30760-0.0", "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", "Calcium"="30680-0.0", 
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0", "Total protein"="30860-0.0", 
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")
                
                feature_values  <- feature_values  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0",   
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0",
                         "BMI"="21001-0.0", "Body Fat Percentage"="23099-0.0", "Weight"="21002-0.0",  # "Waist-hip-ratio"="whr", 
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "Apolipoprotein A"="30630-0.0", "HDL-C"="30760-0.0", "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", "Calcium"="30680-0.0",
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0", "Total protein"="30860-0.0", 
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")
            } else {
                shap_values <- shap_values  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0",   
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0",
                         "BMI"="21001-0.0", "Body Fat Percentage"="23099-0.0", "Weight"="21002-0.0", # "Waist-hip-ratio"="whr",  
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", "Calcium"="30680-0.0",
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0", "Total protein"="30860-0.0", 
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")

                feature_values <- feature_values  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0",   
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0",
                         "BMI"="21001-0.0", "Body Fat Percentage"="23099-0.0", "Weight"="21002-0.0", # "Waist-hip-ratio"="whr", 
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", "Calcium"="30680-0.0",
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0", "Total protein"="30860-0.0", 
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")
            }

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

            # replace glucode with HbA1c
            if (t == 1 && wtype == 1){
                shaporder.noPC.top10 <- shaporder.noPC.top10[!shaporder.noPC.top10 %in% c("Glucose", "Creatinine", "Urea")]
                shaporder.noPC.top10 <- c(shaporder.noPC.top10, shaporder[12:14])
                shaporder.noPC.top10[10] <- "HbA1c"
                print(shaporder.noPC.top10)
            } else if (t == 1 && wtype == 2){
                shaporder.noPC.top10 <- shaporder.noPC.top10[!shaporder.noPC.top10 %in% c("HDL-C")]
                shaporder.noPC.top10 <- c(shaporder.noPC.top10, shaporder[13])
                print(shaporder.noPC.top10)
            }

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

p.final <- (ggobject.LDL.CAD.continuousW | ggobject.TC.CAD.continuousW) / (ggobject.LDL.CAD.binaryW | ggobject.TC.CAD.binaryW) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 35))
# p.final <- (ggobject.LDL.CAD.binaryW | ggobject.TC.CAD.binaryW) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 35))
# p.final <- (ggobject.TC.CAD.continuousW | ggobject.TC.CAD.binaryW) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 35))

ggsave('/home/yujia/Project/2023-07-20-individual_MR/res/03_plot/figure2_variable_importance_beeswarm.png', p.final, dpi=300, 
        width = 60, height = 40, units = "cm", scale = 2)