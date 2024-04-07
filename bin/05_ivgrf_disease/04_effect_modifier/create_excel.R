#' This file is to summarize the effect modifier results from multiple dataset.
#' 
#' @author Shi Yujia
#' @date 2022.12.21

suppressPackageStartupMessages({
    library(xlsx)
})

diseases <- c("CAD")
lipids_traits <- c("LDL", "TC")
trait_type <- c("continuousW", "binaryW")
first_index <- TRUE

for (lt in lipids_traits){
    for (d in diseases){
        for (tt in trait_type){
        
            tmp_df <- read.csv(paste0("/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/04_effct_modifiers/segmented_LS/", lt, "_", d, "_", tt, ".csv"))
        
            if (first_index){
                write.xlsx(tmp_df, file = paste0("/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/04_effct_modifiers/segmented_LS/effect_modifier_segmented_LS_", tt, ".xlsx"), sheetName = paste0(lt, "_", d), 
                            col.names = TRUE, row.names = F, append = FALSE)
                first_index <- F
            } else {
                write.xlsx(tmp_df, file = paste0("/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/04_effct_modifiers/segmented_LS/effect_modifier_segmented_LS_", tt, ".xlsx"), sheetName = paste0(lt, "_", d), 
                            col.names = TRUE, row.names = F, append = TRUE)
            }

        }
    }
}

first_index <- TRUE

for (lt in lipids_traits){
    for (d in diseases){
        for (tt in trait_type){

            tmp_df <- read.csv(paste0("/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/04_effct_modifiers/OLS/", lt, "_", d, "_", tt, ".csv"))
        
            if (first_index){
                write.xlsx(tmp_df, file = paste0("/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/04_effct_modifiers/OLS/effect_modifier_OLS_", tt, ".xlsx"), sheetName = paste0(lt, "_", d), 
                            col.names = TRUE, row.names = F, append = FALSE)
                first_index <- F
            } else {
                write.xlsx(tmp_df, file = paste0("/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/04_effct_modifiers/OLS/effect_modifier_OLS_", tt, ".xlsx"), sheetName = paste0(lt, "_", d), 
                            col.names = TRUE, row.names = F, append = TRUE)
            }
            
        }
    }
}
