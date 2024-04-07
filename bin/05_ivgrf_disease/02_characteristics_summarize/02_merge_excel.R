#' This code is to merge the csv file into a single excel file.
#' 
#' @author Shi Yujia
#' @date 2022.12.12

suppressPackageStartupMessages({
  library(xlsx)
})

traits <- c("LDL", "TC", "TG") # "HDL", 
diseases <- c("CAD")

# continous
first_index <- T
for (i in traits){
  for (j in diseases){
    
    tmp_df <- read.csv(paste0("/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/05_summary_statistics/", i, "_", j, "_continuousW.csv"))
    
    if (first_index){
      write.xlsx(tmp_df, file = "/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/05_summary_statistics/continuous_covariates_summary_statistics.xlsx", sheetName = paste0(i, "_", j), 
                 col.names = TRUE, row.names = F, append = FALSE)
      first_index <- F
    } else {
      write.xlsx(tmp_df, file = "/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/05_summary_statistics/continuous_covariates_summary_statistics.xlsx", sheetName = paste0(i, "_", j), 
                 col.names = TRUE, row.names = F, append = TRUE)
    }
    
    # break
     
  }
}

# binary
first_index <- T
for (i in traits){
  for (j in diseases){
    
    tmp_df <- read.csv(paste0("/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/05_summary_statistics/", i, "_", j, "_binaryW.csv"))
    
    if (first_index){
      write.xlsx(tmp_df, file = "/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/05_summary_statistics/binary_covariates_summary_statistics.xlsx", sheetName = paste0(i, "_", j), 
                 col.names = TRUE, row.names = F, append = FALSE)
      first_index <- F
    } else {
      write.xlsx(tmp_df, file = "/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/05_summary_statistics/binary_covariates_summary_statistics.xlsx", sheetName = paste0(i, "_", j), 
                 col.names = TRUE, row.names = F, append = TRUE)
    }
    
    # break
    
  }
}

# first_index <- T
# for (i in traits){
#   for (j in diseases){
    
#     tmp_df <- read.csv(paste0("/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/05_summary_statistics/", i, "_", j, "_continuousW.csv"))
    
#     if (first_index){
#       final_df <- tmp_df
#       first_index <- F
#     } else {
#       final_df <- dplyr::left_join(final_df, tmp_df)
#     }
     
#   }
# }


