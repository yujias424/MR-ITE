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

diseases <- c("CAD", "IST")
p.value.cutoff <- c(1e-8, 1e-7, 1e-6, 1e-5)

for (ds in diseases){
  for (p in p.value.cutoff){

    # merge file is obtained from harmonize_HDL.R file.
    ssdat <- fread(paste0("~/Project/2022-09-01-individual_MR/dat/05_two_sample_mr/HDL/merged_", ds, "_", p, "_WFGS.txt.gz"))
    ssdat <- as.data.frame(ssdat)

    MRInputObject <- mr_input(snps = ssdat$SNP,
                              bx = ssdat$beta.exposure,
                              bxse = ssdat$se.exposure,
                              by = ssdat$beta.outcome,
                              byse = ssdat$se.outcome)

    MRAllObject_conmix <- mr_conmix(MRInputObject)
    message(paste0("Currently running the disease ", ds, " on p value cutoff ", p, "."))
    message(length(MRAllObject_conmix@ValidSNPs))
    message("\n==============================\n")

    ssdat_filter <- ssdat[ssdat$SNP %in% MRAllObject_conmix@ValidSNPs, ]

    selected_snps <- data.frame("SNP" = MRAllObject_conmix@ValidSNPs)
    write.csv(selected_snps, file = paste0("~/Project/2022-09-01-individual_MR/dat/06_PRS_calculation/HDL/", ds, "/selected_snps_", p, "_WFGS.csv"), 
              row.names = F, quote = F)

  }

  message("\n============================================================\n")

}

message("Finished the pipeline.")