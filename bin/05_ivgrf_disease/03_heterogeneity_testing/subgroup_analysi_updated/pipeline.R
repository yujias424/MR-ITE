filename <- c("/home/yujia/Project/2023-07-20-individual_MR/bin/05_ivgrf_disease/03_heterogeneity_testing/subgroup_analysi_updated/validation_subgroup_analysis_binaryLDL_CAD.R",
              "/home/yujia/Project/2023-07-20-individual_MR/bin/05_ivgrf_disease/03_heterogeneity_testing/subgroup_analysi_updated/validation_subgroup_analysis_binaryTC_CAD.R",
              "/home/yujia/Project/2023-07-20-individual_MR/bin/05_ivgrf_disease/03_heterogeneity_testing/subgroup_analysi_updated/validation_subgroup_analysis_continuousLDL_CAD.R",
              "/home/yujia/Project/2023-07-20-individual_MR/bin/05_ivgrf_disease/03_heterogeneity_testing/subgroup_analysi_updated/validation_subgroup_analysis_continuousTC_CAD.R")

for (i in filename){

    system(paste0("Rscript ", i))

    message(" ")
    print(paste0("Finihsed ", i))
    message("\n =============================================== \n")
    
}