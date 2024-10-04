command_1 <- "Rscript /mnt/md0/yujia/project/2023-07-20-individual_MR/bin/05_ivgrf_disease/01_ITE_analysis/IGF1/01_ivgrf_ITE_HTE_testing/CAD/01_ITE_CAD_IGF1.R"
command_2 <- "python -u /mnt/md0/yujia/project/2023-07-20-individual_MR/bin/05_ivgrf_disease/01_ITE_analysis/IGF1/01_ivgrf_ITE_HTE_testing/CAD/01_ITE_CAD_IGF1_driv.py"
# command_3 <- "Rscript /mnt/md0/yujia/project/2023-07-20-individual_MR/bin/05_ivgrf_disease/01_ITE_analysis/IGF1/01_ivgrf_ITE_HTE_testing/CAD/02_ITE_CAD_IGF1_testHTE.R"

# commands <- c(command_1, command_2, command_3)
commands <- c(command_1, command_2)

for (cmd in commands){

    message(cmd)
    system(cmd, wait = T)
    message(paste0("Finished running command ", cmd, "\n =============================================================== \n"))

}