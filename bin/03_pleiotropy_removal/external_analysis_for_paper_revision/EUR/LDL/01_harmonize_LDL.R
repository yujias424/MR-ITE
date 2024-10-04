#' This code is to using twosampleMR to harmonize the data (including outcome and exposure) for further pleiotropy removal
#' 
#' @author Yujia Shi
#' @date 2022.09.01

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

trait <- "LDL"

# # read the exposure summary data
# exposure <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/", trait, "/QC/", trait, ".QC.gz"), sep = "\t")
# exposure <- as.data.frame(exposure)
# exposure <- exposure[, c("SNP", "A1", "A2", "N", "SE", "P", "BETA", "EAF")]
# colnames(exposure) <- c("SNP", "effect_allele", "other_allele", "samplesize", 
#                           "se", "pval", "beta", "eaf")
# exposure$Phenotype <- rep("LDL", dim(exposure)[1])
# exposure$units <- rep("mmol/L", dim(exposure)[1])
# exposure <- exposure[, c("Phenotype", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf",
#                          "pval", "units", "samplesize")]

# write.table(exposure, paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/05_two_sample_mr/exposure.", trait, ".txt"), row.names = F, sep = " ", quote = F)

p.value.cutoff <- c(0.2)

trait <- "LDL"
ldl_exp_dat.origin <- read_exposure_data(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/05_two_sample_mr/exposure.", trait, ".txt"))
outcome_dat.origin <- read_outcome_data(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/05_two_sample_mr/external_analysis_for_paper_revision/CAD/outcome.CAD.txt"))

for (p in p.value.cutoff){

    message(paste0("Currently we run on pvalue cutoff ", p, "."))

    # Set a cutoff for selecting significant P-value snps in base data.
    ldl_exp_dat <- ldl_exp_dat.origin[which(ldl_exp_dat.origin$pval.exposure < p), ]

    # We can also try not setting any filtering here.
    ldl_exp_dat <- ldl_exp_dat

    row.names(ldl_exp_dat) <- NULL
    print(dim(ldl_exp_dat))

    #' ==========================
    #' using IEU OpenGWAS dataset
    #' ==========================
    # # Heart Disease (CAD)
    # options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
    # chd_out_dat <- extract_outcome_data(
    #     snps = ldl_exp_dat$SNP,
    #     outcomes = 'ebi-a-GCST003116' # ieu-a-7: from CARDIoGRAMplusC4D, "A comprehensive 1000 Genomes–based genome-wide association meta-analysis of coronary artery disease", no UKBB sample included, mixed population; ebi-a-GCST003116: from CARDIoGRAMplusC4D, European.
    # )

    chd_out_dat <- outcome_dat.origin[outcome_dat.origin$SNP %in% ldl_exp_dat$SNP, ]

    harmonise_dat <- harmonise_data(
        exposure_dat = ldl_exp_dat, 
        outcome_dat = chd_out_dat
    )

    harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
                                              "eaf.exposure", 
                                              "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
                                              "beta.outcome", "se.outcome", "pval.outcome")] # , "samplesize.outcome"
    print(dim(harmonise_dat_merged)) 
    message("\n\n")

    fwrite(harmonise_dat_merged, paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/05_two_sample_mr/external_analysis_for_paper_revision/CAD/", trait, "/merged_CAD_", p, ".txt.gz"), row.names = F, sep = " ", quote = F)

    message(paste0("Finished running on pvalue cutoff ", p, ".\n\n=================================\n\n"))

}
