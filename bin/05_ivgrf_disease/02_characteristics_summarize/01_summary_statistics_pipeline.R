#' This code is to testing the run the pipeline to get the summary statistics of all covariates.
#' 
#' @author: Yujia Shi
#' @date: 2022.09.01

traits <- c("LDL", "HDL", "TC", "TG") 

for (i in traits){

    message(paste0("Currently running ", i, "."))
    commands <- paste0("Rscript /home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/02_characteristics_summarize/", i, "/ivgrf_CAD_prs_iv_treatment.R")
    system(command = commands, wait = TRUE)
    # commands <- paste0("Rscript /home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/02_characteristics_summarize/", i, "/ivgrf_ischemic_stroke_prs_iv_treatment.R")
    # system(command = commands, wait = TRUE)
    message("==================================================\n")

}