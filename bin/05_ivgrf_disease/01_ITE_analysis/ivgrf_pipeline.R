traits <- c("LDL", "TC", "HDL", "TG")

for (t in traits){

    message(paste0("\nCurrent running trait ", t))
    runcommand <- paste0("Rscript /home/yujia/Project/2023-07-20-individual_MR/bin/05_ivgrf_disease/01_ITE_analysis/", t, "/01_ivgrf_ITE_HTE_testing/CAD/01_ITE_CAD_", t, ".R")
    system(runcommand, wait = T)
    message(paste0("Finished running trait ", t, "\n"))

}