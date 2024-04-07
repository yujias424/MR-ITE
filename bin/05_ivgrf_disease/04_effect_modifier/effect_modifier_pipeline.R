#' This code is to run the pipeline to test for effect modifier of top 10 important covariates.
#' 
#' @author Shi Yujia
#' @date 2022.12.31

lipid_traits <- c("LDL", "TC")
diseases <- c("CAD")

for (lt in lipid_traits){
    for (d in diseases){
        
        message(paste0("Currently running disease ", d, " with lipid trait ", lt))

        commands <- paste0("Rscript ~/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/04_effect_modifier/",
                            lt, "/ivgrf_", d, "_em.R")
        system(command = commands, wait = TRUE)
        
        message("Finished one.\n===========================\n")

    }
}