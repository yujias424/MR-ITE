#' This code is to using twosampleMR to harmonize the data (including outcome and exposure) for further pleiotropy removal
#' 
#' @author Yujia Shi
#' @date 2022.09.01

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

trait <- "CRP"

# read the exposure summary data
exposure <- fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/", trait, "/QC/", trait, ".QC.gz"), sep = "\t")
exposure <- as.data.frame(exposure)
exposure <- exposure[, c("SNP", "A1", "A2", "SE", "P", "BETA")]
colnames(exposure) <- c("SNP", "effect_allele", "other_allele",
                          "se", "pval", "beta")
exposure$Phenotype <- rep("CRP", dim(exposure)[1])
exposure$units <- rep("mmol/L", dim(exposure)[1])
exposure <- exposure[, c("Phenotype", "SNP", "beta", "se", "effect_allele", "other_allele", 
                         "pval", "units")]

write.table(exposure, paste0("~/Project/2023-07-20-individual_MR/dat/05_two_sample_mr/exposure.", trait, ".txt"), row.names = F, sep = " ", quote = F)

p.value.cutoff <- c(1e-8, 1e-7, 1e-6, 1e-5)

trait <- "CRP"
CRP_exp_dat.origin <- read_exposure_data(paste0("~/Project/2023-07-20-individual_MR/dat/05_two_sample_mr/exposure.", trait, ".txt"))

for (p in p.value.cutoff){

    message(paste0("Currently we run on pvalue cutoff ", p, "."))

    # Set a cutoff for selecting significant P-value snps in base data.
    CRP_exp_dat <- CRP_exp_dat.origin[which(CRP_exp_dat.origin$pval.exposure < p), ]

    # We can also try not setting any filtering here.
    CRP_exp_dat <- CRP_exp_dat

    row.names(CRP_exp_dat) <- NULL
    print(dim(CRP_exp_dat))

    #' ==========================
    #' using IEU OpenGWAS dataset
    #' ==========================
    # Heart Disease (CAD)
    chd_out_dat <- extract_outcome_data(
        snps = CRP_exp_dat$SNP,
        outcomes = 'ebi-a-GCST003116' # ieu-a-7: from CARDIoGRAMplusC4D, "A comprehensive 1000 Genomes–based genome-wide association meta-analysis of coronary artery disease", no UKBB sample included, mixed population; ebi-a-GCST003116: from CARDIoGRAMplusC4D, European.
    )

    harmonise_dat <- harmonise_data(
        exposure_dat = CRP_exp_dat, 
        outcome_dat = chd_out_dat
    )

    harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
                                              "beta.exposure", "se.exposure", "pval.exposure",
                                              "beta.outcome", "se.outcome", "pval.outcome")]
    print(dim(harmonise_dat_merged)) 
    head(harmonise_dat_merged)
    message("\n\n")

    fwrite(harmonise_dat_merged, paste0("~/Project/2023-07-20-individual_MR/dat/05_two_sample_mr/CAD/", trait, "/merged_CAD_", p, ".txt.gz"), row.names = F, sep = " ", quote = F)

    # # Distolic Blood Pressure
    # dbp_out_dat <- extract_outcome_data(
    #     snps = CRP_exp_dat$SNP,
    #     outcomes = 'ieu-b-39' # ieu-b-39 International Consortium of Blood Pressure
    # )

    # harmonise_dat <- harmonise_data(
    #     exposure_dat = CRP_exp_dat, 
    #     outcome_dat = dbp_out_dat
    # )

    # harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
    #                                           "eaf.exposure", 
    #                                           "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
    #                                           "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")] 
    # print(dim(harmonise_dat_merged)) 
    # message("\n\n")

    # fwrite(harmonise_dat_merged, paste0("~/Project/2023-07-20-individual_MR/dat/05_two_sample_mr/DBP/", trait, "/merged_DBP_", p, ".txt.gz"), row.names = F, sep = " ", quote = F)

    # # Systolic Blood Pressure
    # sbp_out_dat <- extract_outcome_data(
    #     snps = CRP_exp_dat$SNP,
    #     outcomes = 'ieu-b-38' # ieu-b-38 International Consortium of Blood Pressure https://www.ahajournals.org/doi/full/10.1161/HYPERTENSIONAHA.120.16138
    # )

    # harmonise_dat <- harmonise_data(
    #     exposure_dat = CRP_exp_dat, 
    #     outcome_dat = sbp_out_dat
    # )

    # harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
    #                                           "eaf.exposure", 
    #                                           "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
    #                                           "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")] 
    # print(dim(harmonise_dat_merged)) 
    # message("\n\n")

    # fwrite(harmonise_dat_merged, paste0("~/Project/2023-07-20-individual_MR/dat/05_two_sample_mr/SBP/", trait, "/merged_SBP_", p, ".txt.gz"), row.names = F, sep = " ", quote = F)

    message(paste0("Finished running on pvalue cutoff ", p, ".\n\n=================================\n\n"))

}
