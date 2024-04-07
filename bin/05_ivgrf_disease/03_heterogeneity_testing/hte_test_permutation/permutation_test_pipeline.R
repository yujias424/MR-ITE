#' This code is to testing the presence of heterogeneity.
#' 
#' @author: Yujia Shi
#' @date: 2022.09.01

traits <- c("LDL", "TC") # "LDL", "HDL", "TC", "TG"

for (i in traits){

    message(paste0("Currently running ", i, "."))
    commands <- paste0("Rscript /home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/01_ITE_analysis/", i, "/01_ivgrf_ITE_HTE_testing/02_ivgrf_test_heterogeneity_", i, ".R")
    system(command = commands, wait = TRUE)
    message("==================================================\n")

}