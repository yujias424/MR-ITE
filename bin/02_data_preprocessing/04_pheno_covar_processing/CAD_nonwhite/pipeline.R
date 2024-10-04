traits <- c("LDL", "TC", "HDL", "TG", "IGF1", "CRP")

for (t in traits){

    command_arg <- paste0("Rscript /mnt/md0/yujia/project/2023-07-20-individual_MR/bin/02_data_preprocessing/04_pheno_covar_processing/CAD_nonwhite/01_ITE_CAD_", t, ".R")
    system(command_arg, wait = T)
    print(paste0("Finished running trait ", t, "."))

}