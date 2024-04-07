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
    library(ggExtra)
    library(latex2exp)
    library(ggpointdensity)
    library(viridis)
})

traits <- c("LDL" , "TC") 
diseases <- c("CAD")
diseases_name <- c("Coronary Artery Disease") 
traits_type <- c("binary", "continuous")
traits_name <- c("LDL-C", "Total Cholesterol")
set.seed(1)

for (t in 1:length(traits)){
    for (d in 1:length(diseases)){
        for (wtype in 1:2){
            
            # t <- 1; d <- 1; wtype <- 1
            if (wtype == 1){
                shap_values <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", diseases[d], "/03_variable_importance/driv_", traits_type[wtype], "W_shap_model3.csv.gz"))
                feature_values <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", diseases[d], "/03_variable_importance/driv_", traits_type[wtype], "W_value_model3.csv.gz"))
                taus <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", diseases[d], "/01_individual_treatment_effect/driv_full_", traits_type[wtype], "W_continuousZ_te_ul_model3.csv"))
                # taus <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", diseases[d], "/01_individual_treatment_effect/ivgrf_full_", traits_type[wtype], "W_continuousZ_pvalue_model3.csv"))
            
            } else {
                shap_values <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", diseases[d], "/03_variable_importance/driv_", traits_type[wtype], "W_shap_model3.csv.gz"))
                feature_values <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", diseases[d], "/03_variable_importance/driv_", traits_type[wtype], "W_value_model3.csv.gz"))
                taus <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", diseases[d], "/01_individual_treatment_effect/driv_full_", traits_type[wtype], "W_te_ul_model3.csv"))
                # taus <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/", diseases[d], "/01_individual_treatment_effect/ivgrf_full_", traits_type[wtype], "W_pvalue_model3.csv"))
            }
            
            if (traits[t] == "LDL"){
                shap_values  <- shap_values  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0",   
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0",
                         "Body Fat Percentage"="23099-0.0", "BMI"="21001-0.0", "Weight"="21002-0.0",  # "Waist-hip-ratio"="whr", 
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "Apolipoprotein A"="30630-0.0", "HDL-C"="30760-0.0", "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Calcium"="30680-0.0","Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", 
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Total protein"="30860-0.0", "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0",
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")
                
                feature_values  <- feature_values  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0",   
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0", 
                         "Body Fat Percentage"="23099-0.0", "BMI"="21001-0.0", "Weight"="21002-0.0", # "Waist-hip-ratio"="whr", 
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "Apolipoprotein A"="30630-0.0", "HDL-C"="30760-0.0", "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Calcium"="30680-0.0","Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", 
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Total protein"="30860-0.0", "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0",
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")
            } else {
                shap_values <- shap_values  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0",   
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0",
                         "Body Fat Percentage"="23099-0.0", "BMI"="21001-0.0", "Weight"="21002-0.0", # "Waist-hip-ratio"="whr", 
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Calcium"="30680-0.0","Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", 
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Total protein"="30860-0.0", "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0",
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")

                feature_values <- feature_values  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0",   
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0",
                         "Body Fat Percentage"="23099-0.0", "BMI"="21001-0.0", "Weight"="21002-0.0", # "Waist-hip-ratio"="whr", 
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Calcium"="30680-0.0","Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", 
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Total protein"="30860-0.0", "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0",
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")
            }

            # get the order of top important variables. (based on DRIV)
            shaporder <- names(sort(sapply(abs(shap_values), mean), decreasing = T))
            
            # # can also plot top important variable based on ivgrf
            # if (wtype == 1){
            #     shaporder <- read.csv(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/CAD/03_variable_importance/ivgrf_binaryW_continuousZ_variable_importance_model3.csv"))
            #     shaporder <- shaporder[order(shaporder$Importance, decreasing = T), ]$Variable
            # } else {
            #     shaporder <- read.csv(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", traits[t], "/CAD/03_variable_importance/ivgrf_continuousW_variable_importance_model3.csv"))
            #     shaporder <- shaporder[order(shaporder$Importance, decreasing = T), ]$Variable
            # }

            # # remove the PC variables from the top order
            # shaporder.noPC.top10 <- shaporder[!startsWith(shaporder, "PC")][1:10]

            # remove the PC and WHR
            shaporder <- shaporder[!startsWith(shaporder, "PC")]
            shaporder <- shaporder[!startsWith(shaporder, "Genotype")]

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

            # select top 10 variables
            shap_values.top10 <- shap_values[, ..shaporder.noPC.top10]
            feature_values.top10 <- feature_values[, ..shaporder.noPC.top10]

            shap_values.top10 <- as.data.frame(shap_values.top10)
            feature_values.top10 <- as.data.frame(feature_values.top10)

            for (i in 1:length(shaporder.noPC.top10)){
                
                p1.dat <- data.frame("shap_value" = shap_values.top10[, c(shaporder.noPC.top10[i])], "feature_value" = feature_values.top10[, c(shaporder.noPC.top10[i])])
                random_sample <- sort(sample(dim(p1.dat)[1], 5000))

                p1 <- ggplot(p1.dat, aes(y=shap_value, x=feature_value)) +
                        geom_point(size = 2, color = "#0091ff") +
                        geom_smooth(data = p1.dat[random_sample, ], se = F, linetype = "dashed", color = "red", size = 2) +
                        theme_classic() +
                        xlab(shaporder.noPC.top10[i]) + ylab(paste0("SHAP value for\n", shaporder.noPC.top10[i])) +
                        theme(
                            axis.text = element_text(size = 35),
                            axis.title.x = element_text(size = 35),
                            axis.title.y = element_text(size = 35),
                            aspect.ratio=1
                        ) 

                # p2 <- ggMarginal(p1, margins = "x", type = "histogram", size = 10, col = "black", fill = "gray", bins = 100)
                # assign(paste0("p.shap.feature.", i), patchwork::wrap_elements(p2))
                # rm(p2)

                assign(paste0("p.shap.feature.", i), p1)
                
                # # p3.dat <- data.frame("tau_value" = taus$point, "feature_value" = feature_values.top10[, c(shaporder.noPC.top10[i])]) # driv
                # p3.dat <- data.frame("tau_value" = taus$tau, "feature_value" = feature_values.top10[, c(shaporder.noPC.top10[i])]) # ivgrf
                # p3 <- ggplot(p3.dat, aes(y=tau_value, x=feature_value)) +
                #         geom_pointdensity(adjust = .001) +
                #         scale_color_viridis() +
                #         geom_smooth(data = p3.dat[random_sample, ], se = F, linetype = "dashed", color = "red", size = 2) +
                #         theme_classic() +
                #         xlab(shaporder.noPC.top10[i]) + ylab(TeX("Individual Treatment Effect")) +
                #         theme(
                #             axis.text = element_text(size = 35),
                #             axis.title.x = element_text(size = 35),
                #             axis.title.y = element_text(size = 35),
                #             legend.position = "none",
                #             aspect.ratio=1
                #         ) 

                p3.dat <- data.frame("tau_value" = taus$point, "feature_value" = feature_values.top10[, c(shaporder.noPC.top10[i])])
                p4.dat <- p3.dat[order(p3.dat$feature_value, decreasing = FALSE), ]
                # we split the data into 10 quantiles
                unordered <- rep(
                    seq_len(10),
                    each = dim(p4.dat)[1] / 10
                )
                unordered <- c(unordered, rep(10, 4))
                p4.dat$group <- unordered
                p4.dat$group <- as.factor(p4.dat$group)

                random_sample_p4 <- caret::createDataPartition(p4.dat$group, p = 0.01)[[1]]
                
                p4 <- ggplot(p4.dat, aes(y=tau_value, x=group)) +
                        geom_boxplot() +
                        geom_smooth(aes(group=1), se = F, linetype = "solid", color = "red", size = 2) +
                        theme_classic() +
                        xlab(shaporder.noPC.top10[i]) + ylab(TeX("Individual Treatment Effect")) +
                        theme(
                            axis.text = element_text(size = 35),
                            axis.title.x = element_text(size = 35),
                            axis.title.y = element_text(size = 35),
                            legend.position = "none",
                            aspect.ratio=1
                        )

                assign(paste0("p.tau.feature.", i), p4)
                # break
            }

            # option 1:
            # p.final.1 <- (p.shap.feature.1 | p.shap.feature.2 | p.shap.feature.3 | p.shap.feature.4 | p.shap.feature.5)
            # p.final.2 <- (p.shap.feature.6 | p.shap.feature.7 | p.shap.feature.8 | p.shap.feature.9 | p.shap.feature.10)
            # p.final.3 <- (p.tau.feature.1 | p.tau.feature.2 | p.tau.feature.3 | p.tau.feature.4 | p.tau.feature.5)
            # p.final.4 <- (p.tau.feature.6 | p.tau.feature.7 | p.tau.feature.8 | p.tau.feature.9 | p.tau.feature.10)

            # p.final <- p.final.1 / p.final.2 / p.final.3 / p.final.4 
            
            # option 2:
            # p.final <- ((p.shap.feature.1 | p.shap.feature.2 | p.shap.feature.3 | p.shap.feature.4 | p.shap.feature.5) / 
            #             (p.shap.feature.6 | p.shap.feature.7 | p.shap.feature.8 | p.shap.feature.9 | p.shap.feature.10) / 
            #             (p.tau.feature.1 | p.tau.feature.2 | p.tau.feature.3 | p.tau.feature.4 | p.tau.feature.5) / 
            #             (p.tau.feature.6 | p.tau.feature.7 | p.tau.feature.8 | p.tau.feature.9 | p.tau.feature.10)) + 
            #             plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 25))
            
            # option 3:
            p.final <- (p.shap.feature.1 / p.tau.feature.1 / p.shap.feature.6 / p.tau.feature.6 |
                        p.shap.feature.2 / p.tau.feature.2 / p.shap.feature.7 / p.tau.feature.7 |
                        p.shap.feature.3 / p.tau.feature.3 / p.shap.feature.8 / p.tau.feature.8 |
                        p.shap.feature.4 / p.tau.feature.4 / p.shap.feature.9 / p.tau.feature.9 |
                        p.shap.feature.5 / p.tau.feature.5 / p.shap.feature.10 / p.tau.feature.10) + 
                        plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 35))
                
            ggsave(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/03_plot/figure3_variable_importance_", traits[t], "_", diseases[d], "_", traits_type[wtype], "W_boxplot.png"), p.final, dpi=300, 
                        width = 70, height = 60, units = "cm", scale = 2, limitsize = FALSE)
            message("\n\n=========================\n\n")
            # break
            
        }
        # break
    }
    # break
}
