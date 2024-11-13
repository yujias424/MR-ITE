#' This code is to perform pleiotropy removal using package Mendelian Randomization with Conmix approach.
#' An aim of this code is to identify those potential invalid instruments.
#' 
#' @author Yujia Shi
#' @date 2022.09.01

suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
  library(dplyr)
})

diseases <- c("CAD")
p.value.cutoff <- c(1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5)

for (ds in diseases){
  for (p in p.value.cutoff){
    
    # merge file is obtained from harmonize_IGF1.R file.
    ssdat <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/05_two_sample_mr/", ds, "/IGF1/merged_", ds, "_", p, ".txt.gz"))
    ssdat <- as.data.frame(ssdat)

    MRInputObject <- mr_input(snps = ssdat$SNP,
                              bx = ssdat$beta.exposure,
                              bxse = ssdat$se.exposure,
                              by = ssdat$beta.outcome,
                              byse = ssdat$se.outcome,
                              effect_allele = ssdat$effect_allele.exposure,
                              other_allele = ssdat$other_allele.exposure
                              )

    MRAllObject_conmix <- mr_conmix(MRInputObject)
    message(paste0("Currently running the disease ", ds, " on p value cutoff ", p, "."))
    message(length(MRAllObject_conmix@ValidSNPs))
    message("\n==============================\n")

    ssdat_filter <- ssdat[ssdat$SNP %in% MRAllObject_conmix@ValidSNPs, ]

    selected_snps <- data.frame("SNP" = MRAllObject_conmix@ValidSNPs)
    write.csv(selected_snps, file = paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/IGF1/", ds, "/selected_snps_", p, ".csv"), 
              row.names = F, quote = F)
    # break
  }

  message("\n============================================================\n")

}

message("Finished the pipeline.")