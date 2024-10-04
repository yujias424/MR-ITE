traits <- c("LDL", "HDL", "TC", "TG", "CRP", "IGF1")

for (i in traits){

    message(paste0("Currently running the trait ", i, "."))
    commands <- paste0("Rscript /mnt/md0/yujia/project/2023-07-20-individual_MR/bin/05_ivgrf_disease/01_ITE_analysis/extra_analysis_for_paper_revision/PRSice2_PRS/", i, "/01_ivgrf_ITE_HTE_testing/CAD/01_ITE_CAD_", i, ".R")
    system(commands)

}