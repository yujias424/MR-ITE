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

# read the exposure summary data
exposure <- fread(paste0("/home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/", trait, "/QC/", trait, ".WFGS.QC.gz"), sep = "\t")
exposure <- as.data.frame(exposure)
exposure <- exposure[, c("SNP", "A1", "A2", "N", "SE", "P", "BETA")]
colnames(exposure) <- c("SNP", "effect_allele", "other_allele", "samplesize", 
                          "se", "pval", "beta")
exposure$Phenotype <- rep("HDL", dim(exposure)[1])
exposure$units <- rep("mg/dL", dim(exposure)[1])
exposure <- exposure[, c("Phenotype", "SNP", "beta", "se", "effect_allele", "other_allele",
                         "pval", "units", "samplesize")]

write.table(exposure, paste0("~/Project/2022-09-01-individual_MR/dat/05_two_sample_mr/", trait, "/exposure.WFGS.txt"), row.names = F, sep = " ", quote = F)

p.value.cutoff <- c(1e-8, 1e-7, 1e-6, 1e-5) 

trait <- "HDL"
hdl_exp_dat.origin <- read_exposure_data(paste0("~/Project/2022-09-01-individual_MR/dat/05_two_sample_mr/", trait, "/exposure.WFGS.txt"))

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
    # Heart Disease (CAD)
    chd_out_dat <- extract_outcome_data(
        snps = hdl_exp_dat$SNP,
        outcomes = 'ebi-a-GCST003116' # ieu-a-7: from CARDIoGRAMplusC4D, "A comprehensive 1000 Genomes–based genome-wide association meta-analysis of coronary artery disease", no UKBB sample included, mixed population; ebi-a-GCST003116: from CARDIoGRAMplusC4D, European.
    )

    harmonise_dat <- harmonise_data(
        exposure_dat = hdl_exp_dat, 
        outcome_dat = chd_out_dat
    )

    harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                                              "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
                                              "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")] 
    print(dim(harmonise_dat_merged)) 
    message("\n\n")

    fwrite(harmonise_dat_merged, paste0("~/Project/2022-09-01-individual_MR/dat/05_two_sample_mr/", trait, "/merged_CAD_", p, "_WFGS.txt.gz"), row.names = F, sep = " ", quote = F)

    # Ischemic stroke
    is_out_dat <- extract_outcome_data(
        snps = hdl_exp_dat$SNP,
        outcomes = 'ebi-a-GCST006908' # ebi-a-GCST006908: European; ukb-d-I9_STR_EXH ieu-a-1108
    )

    harmonise_dat <- harmonise_data(
        exposure_dat = hdl_exp_dat, 
        outcome_dat = is_out_dat
    )

    harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                                              "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
                                              "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")]
    print(dim(harmonise_dat_merged))

    fwrite(harmonise_dat_merged, paste0("~/Project/2022-09-01-individual_MR/dat/05_two_sample_mr/", trait, "/merged_IST_", p, "_WFGS.txt.gz"), row.names = F, sep = " ", quote = F)
    message(paste0("Finished running on pvalue cutoff ", p, ".\n\n=================================\n\n"))

}

#' ====================================================================================================================================================================================

# #' =======================================================================
# #' CAD using https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.117.312086
# #' =======================================================================

# chd_out_dat <- read_outcome_data(snps = hdl_exp_dat$SNP, filename = "~/Project/2022-09-01-individual_MR/dat/analysis/pheno_data/GWAS_CAD/CAD_UKBIOBANK.gz", sep = " ",
#                                  snp_col = "oldID",
#                                  beta_col = "beta",
#                                  se_col = "se",
#                                  effect_allele_col = "a1",
#                                  other_allele_col = "a2",
#                                  eaf_col = "af",
#                                  pval_col = "pval",
#                                  samplesize_col = "N",
#                                  chr_col = "chr",
#                                  pos_col = "bp")
# chd_out_dat$outcome <- rep("CAD", dim(chd_out_dat)[1])
# harmonise_dat <- harmonise_data(
#     exposure_dat = hdl_exp_dat, 
#     outcome_dat = chd_out_dat
# )

# harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
#                                           "eaf.exposure", 
#                                           "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
#                                           "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")]
# dim(harmonise_dat_merged) # 3128 SNPs

# write.table(harmonise_dat_merged, paste0("~/Project/2022-09-01-individual_MR/dat/analysis/two_sample_mr/", trait, "/merged_heart_disease.txt"), row.names = F, sep = " ", quote = F)

# #' =======================================================================
# #' Any Stroke using https://www.megastroke.org/download.html
# #' =======================================================================

# as_out_dat <- read_outcome_data(snps = hdl_exp_dat$SNP, filename = "~/Project/2022-09-01-individual_MR/dat/analysis/pheno_data/GWAS_IS/MEGASTROKE.1.AS.TRANS.out", sep = " ",
#                                  snp_col = "MarkerName",
#                                  beta_col = "Effect",
#                                  se_col = "StdErr",
#                                  effect_allele_col = "Allele1",
#                                  other_allele_col = "Allele2",
#                                  eaf_col = "Freq1",
#                                  pval_col = "P-value")
# as_out_dat$outcome <- rep("AS", dim(as_out_dat)[1])
# harmonise_dat <- harmonise_data(
#     exposure_dat = hdl_exp_dat, 
#     outcome_dat = as_out_dat
# )

# harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
#                                           "eaf.exposure", 
#                                           "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
#                                           "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")]
# dim(harmonise_dat_merged) # 3119 SNPs

# write.table(harmonise_dat_merged, paste0("~/Project/2022-09-01-individual_MR/dat/analysis/two_sample_mr/", trait, "/merged_stroke.txt"), row.names = F, sep = " ", quote = F)

# #' =======================================================================
# #' Ischemic Stroke using https://www.megastroke.org/download.html
# #' =======================================================================

# is_out_dat <- read_outcome_data(snps = hdl_exp_dat$SNP, filename = "~/Project/2022-09-01-individual_MR/dat/analysis/pheno_data/GWAS_IS/MEGASTROKE.2.AIS.TRANS.out", sep = " ",
#                                  snp_col = "MarkerName",
#                                  beta_col = "Effect",
#                                  se_col = "StdErr",
#                                  effect_allele_col = "Allele1",
#                                  other_allele_col = "Allele2",
#                                  eaf_col = "Freq1",
#                                  pval_col = "P-value")
# is_out_dat$outcome <- rep("IS", dim(is_out_dat)[1])
# harmonise_dat <- harmonise_data(
#     exposure_dat = hdl_exp_dat, 
#     outcome_dat = is_out_dat
# )

# harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
#                                           "eaf.exposure", 
#                                           "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
#                                           "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")]
# dim(harmonise_dat_merged) # 3117 SNPs

# write.table(harmonise_dat_merged, paste0("~/Project/2022-09-01-individual_MR/dat/analysis/two_sample_mr/", trait, "/merged_ischemic_stroke.txt"), row.names = F, sep = " ", quote = F)

#' ====================================================================================================================================================================================

# # Stroke (temporarily deprecated)
# chd_out_dat <- extract_outcome_data(
#     snps = hdl_exp_dat$SNP,
#     outcomes = 'ebi-a-GCST005838'
# )

# harmonise_dat <- harmonise_data(
#     exposure_dat = hdl_exp_dat, 
#     outcome_dat = chd_out_dat
# )

# harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
#                                           "eaf.exposure", 
#                                           "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure",
#                                           "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome")]

# write.table(harmonise_dat_merged, paste0("~/Project/2022-09-01-individual_MR/dat/analysis/two_sample_mr/", trait, "/merged_stroke.txt"), row.names = F, sep = " ", quote = F)