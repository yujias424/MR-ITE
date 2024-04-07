traits <- c("LDL", "HDL", "TC", "TG")

for (t in traits){

    if (t == "HDL"){
        commands <- paste0("Rscript /home/yujia/Project/2023-07-20-individual_MR/bin/03_pleiotropy_removal/HDL/01_harmonize_HDL_GLGC.R")
        system(commands, wait = T)

        commands <- paste0("Rscript /home/yujia/Project/2023-07-20-individual_MR/bin/03_pleiotropy_removal/HDL/02_conmix_pleiotropy_removal_GLGC.R")
        system(commands, wait = T)

        print("==================\nFinished HDL.\n==================")

    } else {
        commands <- paste0("Rscript /home/yujia/Project/2023-07-20-individual_MR/bin/03_pleiotropy_removal/", t, "/01_harmonize_", t, ".R")
        system(commands, wait = T)

        commands <- paste0("Rscript /home/yujia/Project/2023-07-20-individual_MR/bin/03_pleiotropy_removal/", t, "/02_conmix_pleiotropy_removal.R")
        system(commands, wait = T)

        print(paste0("==================\nFinished ", t, ".\n=================="))
    }
}