#' This code is to run the PRS calculation automatically for all seven different disease we are interested in.
#' 
#' @author Shi Yujia
#' @date 2023.01.19

####### lipid traits #######
# traits <- c("LDL", "TC", "TG")
traits <- c("TG")
diseases <- c("CAD")

p.value.cutoff <- c(1e-8) # , 1e-7, 1e-6, 1e-5
barlevels <- c("1e-8,5e-8,7e-8,9e-8,1e-7,5e-7,7e-7,9e-7,1e-6,5e-6,7e-6,9e-6,1e-5,5e-5,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2,0.1,0.5 ")
barlevels <- c("1e-10,5e-10,1e-9,5e-9,1e-8 ")
set_arr <- c("set1", "set2", "set3")

for (t in traits){
    for (sa in set_arr){
        for (i in diseases){

            message(paste0("Currently running ", t, " on ", i, "."))

            for (p in 1:length(p.value.cutoff)){

                message(paste0("Running P-value cut-off ", p.value.cutoff[p]))
                dir.create(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/TG_GLGC_2021/", i, "/", sa, "/", p.value.cutoff[p]), recursive = T, showWarnings = F)
                setwd(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/TG_GLGC_2021/", i, "/", sa, "/", p.value.cutoff[p]))

                if (sa == "set1"){
                    commands <- paste0("/home/yujia/Software/PRSice/PRSice_linux ",
                                        "--a1 A1 ",
                                        "--a2 A2 ",
                                        "--bp BP ",
                                        "--chr CHR ", # "--bar-levels ", barlevels[1], "--fastscore ", 
                                        "--lower 1e-12 ",
                                        "--base /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/", t, ".QC.gz ", # .QC.wthPos.
                                        "--beta ",  
                                        "--binary-target F ", 
                                        "--clump-kb 250kb ",
                                        "--clump-p 1.000000 ",
                                        "--clump-r2 0.100000 ", # follow PRSice tutorial
                                        "--cov /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/", t, "/", i, "/ukbb.covariate.prs.", t, ".", sa, " ", 
                                        "--cov-factor Genotype_Batch,Gender ", # Assessment_Centre,
                                        "--num-auto 22 ",
                                        "--out ", t, "_prs ",
                                        "--pheno /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/", t, "/", i, "/ukbb.phenotype.", t, ".mgdL ", 
                                        "--print-snp ",
                                        "--pvalue P ",
                                        "--seed 2389825711 ",
                                        "--snp SNP ",
                                        "--stat BETA ",
                                        "--target /mnt/md0/public_data/UKBB_data/genotype/GWAS_ukbb_maf0.001_info0.8/bed/allchrs/ukb_allchrs ", # /mnt/data/yujia/GWAS_ukbb_maf0.001_info0.3/chr#_QC ",
                                        "--thread 30 ",
                                        "--score con-std ", 
                                        "--extract /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/TG_GLGC_2021/", i, "/selected_snps_", p.value.cutoff[p],".csv")
                } else if (sa == "set2") {
                    commands <- paste0("/home/yujia/Software/PRSice/PRSice_linux ",
                                        "--a1 A1 ",
                                        "--a2 A2 ",
                                        "--bp BP ",
                                        "--chr CHR ", # "--bar-levels ", barlevels[1], "--fastscore ", 
                                        "--lower 1e-12 ",
                                        "--base /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/", t, ".QC.gz ", 
                                        "--beta ",  
                                        "--binary-target F ", 
                                        "--clump-kb 250kb ",
                                        "--clump-p 1.000000 ",
                                        "--clump-r2 0.100000 ", # follow PRSice tutorial
                                        "--cov /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/", t, "/", i, "/ukbb.covariate.prs.", t, ".", sa, " ", 
                                        "--cov-factor Genotype_Batch,Gender,Non_smoker,Previous_smoker,Current_smoker,",
                                        "Non_alcohol_drinker,Previous_alcohol_drinker,Current_alcohol_drinker,",
                                        "Cholesterol_lowering_medication,Blood_pressure_medication,Insulin,", # Assessment_Centre, # No_medication,
                                        "htn,t2dm,heart_failure,hemorrhage_stroke,ischemic_stroke ",
                                        "--num-auto 22 ",
                                        "--out ", t, "_prs ",
                                        "--pheno /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/", t, "/", i, "/ukbb.phenotype.", t, ".mgdL ", 
                                        "--print-snp ",
                                        "--pvalue P ",
                                        "--seed 2389825711 ",
                                        "--snp SNP ",
                                        "--stat BETA ",
                                        "--target /mnt/md0/public_data/UKBB_data/genotype/GWAS_ukbb_maf0.001_info0.8/bed/allchrs/ukb_allchrs ", # /mnt/data/yujia/GWAS_ukbb_maf0.001_info0.3/chr#_QC ",
                                        "--thread 30 ",
                                        "--score con-std ", 
                                        "--extract /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/TG_GLGC_2021/", i, "/selected_snps_", p.value.cutoff[p],".csv")
                } else {
                    commands <- paste0("/home/yujia/Software/PRSice/PRSice_linux ",
                                        "--a1 A1 ",
                                        "--a2 A2 ",
                                        "--bp BP ",
                                        "--chr CHR ", # "--bar-levels ", barlevels[1], "--fastscore ", 
                                        "--lower 1e-12 ",
                                        "--base /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/", t, ".QC.gz ", 
                                        "--beta ",  
                                        "--binary-target F ", 
                                        "--clump-kb 250kb ",
                                        "--clump-p 1.000000 ",
                                        "--clump-r2 0.100000 ", # follow PRSice tutorial
                                        "--cov /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/", t, "/", i, "/ukbb.covariate.prs.", t, ".", sa, " ", 
                                        "--cov-factor Genotype_Batch,Gender,Non_smoker,Previous_smoker,Current_smoker,",
                                        "Non_alcohol_drinker,Previous_alcohol_drinker,Current_alcohol_drinker,",
                                        "Cholesterol_lowering_medication,Blood_pressure_medication,Insulin,", # Assessment_Centre, # No_medication,
                                        "htn,t2dm,heart_failure,hemorrhage_stroke,ischemic_stroke ",
                                        "--num-auto 22 ",
                                        "--out ", t, "_prs ",
                                        "--pheno /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/", t, "/", i, "/ukbb.phenotype.", t, ".mgdL ", 
                                        "--print-snp ",
                                        "--pvalue P ",
                                        "--seed 2389825711 ",
                                        "--snp SNP ",
                                        "--stat BETA ",
                                        "--target /mnt/md0/public_data/UKBB_data/genotype/GWAS_ukbb_maf0.001_info0.8/bed/allchrs/ukb_allchrs ", # /mnt/data/yujia/GWAS_ukbb_maf0.001_info0.3/chr#_QC ",
                                        "--thread 30 ",
                                        "--score con-std ", 
                                        "--extract /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/TG_GLGC_2021/", i, "/selected_snps_", p.value.cutoff[p],".csv")
                }
                system(command = commands, wait = TRUE)

                message("==================================================\n")

            }
        }
    }
}
