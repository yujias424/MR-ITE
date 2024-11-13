#' This code is to using twosampleMR to harmonize the data (including outcome and exposure) for further pleiotropy removal
#' 
#' @author Yujia Shi
#' @date 2022.09.01

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

trait <- "HDL"

# read the exposure summary data (done)
# exposure <- fread(paste0("/home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/", trait, "/QC/", trait, ".GLGC.QC.gz"), sep = "\t")
# exposure <- as.data.frame(exposure)
# exposure <- exposure[, c("SNP", "A1", "A2", "N", "SE", "P", "BETA", "EAF")]
# colnames(exposure) <- c("SNP", "effect_allele", "other_allele", "samplesize", 
#                           "se", "pval", "beta", "eaf")
# exposure$Phenotype <- rep("HDL", dim(exposure)[1])
# exposure$units <- rep("mmol/L", dim(exposure)[1])
# exposure <- exposure[, c("Phenotype", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf",
#                          "pval", "units", "samplesize")]

# write.table(exposure, paste0("~/Project/2023-07-20-individual_MR/dat/05_two_sample_mr/exposure.GLGC.", trait, ".txt"), row.names = F, sep = " ", quote = F)

p.value.cutoff <- c(1e-8, 1e-7, 1e-6, 1e-5) 

trait <- "HDL"
hdl_exp_dat.origin <- read_exposure_data(paste0("~/Project/2023-07-20-individual_MR/dat/05_two_sample_mr/exposure.GLGC.", trait, ".txt"))

for (p in p.value.cutoff){

    message(paste0("Currently we run on pvalue cutoff ", p, "."))

    # Set a cutoff for selecting significant P-value snps in base data.
    hdl_exp_dat <- hdl_exp_dat.origin[which(hdl_exp_dat.origin$pval.exposure < p), ] 

    # We can also try not setting any filtering here.
    hdl_exp_dat <- hdl_exp_dat

    row.names(hdl_exp_dat) <- NULL
    print(dim(hdl_exp_dat))

    #' ==========================
    #' using IEU OpenGWAS dataset
    #' ==========================
    # # Heart Disease (CAD)
    # chd_out_dat <- extract_outcome_data(
    #     snps = hdl_exp_dat$SNP,
    #     outcomes = 'ebi-a-GCST003116' # ieu-a-7: from CARDIoGRAMplusC4D, "A comprehensive 1000 Genomes–based genome-wide association meta-analysis of coronary artery disease", no UKBB sample included, mixed population; ebi-a-GCST003116: from CARDIoGRAMplusC4D, European.
    # )

    # harmonise_dat <- harmonise_data(
    #     exposure_dat = hdl_exp_dat, 
    #     outcome_dat = chd_out_dat
    # )

    # harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
    #                                           "eaf.exposure", 
    #                                           "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
    #                                           "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")] 
    # print(dim(harmonise_dat_merged)) 
    # message("\n\n")

    # fwrite(harmonise_dat_merged, paste0("~/Project/2023-07-20-individual_MR/dat/05_two_sample_mr/CAD/", trait, "/merged_CAD_", p, "_GLGC.txt.gz"), row.names = F, sep = " ", quote = F)
    
    # Distolic Blood Pressure
    dbp_out_dat <- extract_outcome_data(
        snps = hdl_exp_dat$SNP,
        outcomes = 'ieu-b-39' # ieu-b-39 International Consortium of Blood Pressure
    )

    harmonise_dat <- harmonise_data(
        exposure_dat = hdl_exp_dat, 
        outcome_dat = dbp_out_dat
    )

    harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
                                              "eaf.exposure", 
                                              "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
                                              "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")] 
    print(dim(harmonise_dat_merged)) 
    message("\n\n")

    fwrite(harmonise_dat_merged, paste0("~/Project/2023-07-20-individual_MR/dat/05_two_sample_mr/DBP/", trait, "/merged_DBP_", p, "_GLGC.txt.gz"), row.names = F, sep = " ", quote = F)

    # Systolic Blood Pressure
    sbp_out_dat <- extract_outcome_data(
        snps = hdl_exp_dat$SNP,
        outcomes = 'ieu-b-38' # ieu-b-38 International Consortium of Blood Pressure
    )

    harmonise_dat <- harmonise_data(
        exposure_dat = hdl_exp_dat, 
        outcome_dat = sbp_out_dat
    )

    harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
                                              "eaf.exposure", 
                                              "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
                                              "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")] 
    print(dim(harmonise_dat_merged)) 
    message("\n\n")

    fwrite(harmonise_dat_merged, paste0("~/Project/2023-07-20-individual_MR/dat/05_two_sample_mr/SBP/", trait, "/merged_SBP_", p, "_GLGC.txt.gz"), row.names = F, sep = " ", quote = F)

    message(paste0("Finished running on pvalue cutoff ", p, ".\n\n=================================\n\n"))

}
