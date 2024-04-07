traits <- c("LDL", "TC", "HDL", "TG")

for (t in traits){

    command_arg <- paste0("Rscript /home/yujia/Project/2023-07-20-individual_MR/bin/02_data_preprocessing/04_pheno_covar_processing/CAD/01_ITE_CAD_", t, ".R")
    system(command_arg, wait = T)
    print(paste0("Finished running trait ", t, "."))

}