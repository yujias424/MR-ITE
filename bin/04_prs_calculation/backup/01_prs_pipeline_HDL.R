#' This code is to run the PRS calculation automatically for all seven different disease we are interested in.
#' 
#' @author Shi Yujia
#' @date 2023.01.19

####### lipid traits #######
traits <- c("HDL")
base_data <- c("GLGC")
diseases <- c("CAD")

p.value.cutoff <- c(1e-8, 1e-7, 1e-6, 1e-5)
barlevels <- c("1e-8,5e-8,7e-8,9e-8,1e-7,5e-7,7e-7,9e-7,1e-6,5e-6,7e-6,9e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5 ")

for (t in traits){
    for (i in diseases){
        for (bd in base_data){

            message(paste0("Currently running ", t, " on ", i, "."))

            for (p in 1:length(p.value.cutoff)){

                message(paste0("Running P-value cut-off ", p.value.cutoff[p]))
                dir.create(paste0("/home/yujia/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/", t, "/", i, "/", p.value.cutoff[p]), recursive = T, showWarnings = F)
                setwd(paste0("/home/yujia/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/", t, "/", i, "/", p.value.cutoff[p]))

                commands <- paste0("/home/yujia/software/PRSice-2/PRSice_linux ",
                                    "--a1 A1 ",
                                    "--a2 A2 ",
                                    "--bp BP ",
                                    "--chr CHR ", "--bar-levels ", barlevels[1], "--fastscore ", 
                                    "--base /home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/", t, "/QC/", t, ".", bd, ".QC.gz ", 
                                    "--beta ",  
                                    "--binary-target F ", 
                                    "--clump-kb 250kb ",
                                    "--clump-p 1.000000 ",
                                    "--clump-r2 0.100000 ", # follow PRSice tutorial
                                    "--cov /home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.prs ",
                                    "--cov-factor Genotype_Batch,Gender ", # Assessment_Centre,
                                    "--num-auto 22 ",
                                    "--out ", t, "_prs_", bd, " ",
                                    "--pheno /home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/", t, "/", i, "/ukbb.phenotype.", t, ".mgdL ", 
                                    "--print-snp ",
                                    "--pvalue P ",
                                    "--seed 2389825711 ",
                                    "--snp SNP ",
                                    "--stat BETA ",
                                    "--target /mnt/data/yujia/GWAS_ukbb_maf0.001_info0.8/bed/allchrs/ukb_allchrs ", # /mnt/data/yujia/GWAS_ukbb_maf0.001_info0.3/chr#_QC ",
                                    "--thread 15 ",
                                    "--score con-std ",
                                    "--extract /home/yujia/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/", t, "/", i, "/selected_snps_", p.value.cutoff[p], "_", bd, ".csv")
                system(command = commands, wait = TRUE)

                message("==================================================\n")

            }
        }
    }
}
