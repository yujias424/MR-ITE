#' This code is to using twosampleMR to harmonize the data (including outcome and exposure) for further pleiotropy removal
#' 
#' @author Yujia Shi
#' @date 2022.09.01

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

trait <- "IGF1"

# # read the exposure summary data
# exposure <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/", trait, "/QC/", trait, ".QC.gz"), sep = "\t")
# exposure <- as.data.frame(exposure)
# exposure <- exposure[, c("SNP", "A1", "A2", "SE", "P", "BETA", "EAF")]
# colnames(exposure) <- c("SNP", "effect_allele", "other_allele",
#                           "se", "pval", "beta", "eaf")
# exposure$Phenotype <- rep("IGF1", dim(exposure)[1])
# exposure$units <- rep("nmol/L", dim(exposure)[1])
# exposure <- exposure[, c("Phenotype", "SNP", "beta", "se", "effect_allele", "other_allele", 
#                          "pval", "eaf", "units")]

# write.table(exposure, paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/05_two_sample_mr/exposure.", trait, ".txt"), row.names = F, sep = " ", quote = F)

trait <- "IGF1"
IGF1_exp_dat.origin <- read_exposure_data(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/05_two_sample_mr/exposure.", trait, ".txt"))

# Set a cutoff for selecting significant P-value snps in base data.
# IGF1_exp_dat <- IGF1_exp_dat.origin[which(IGF1_exp_dat.origin$pval.exposure < 0.1), ]
IGF1_exp_dat <- IGF1_exp_dat.origin

# We can also try not setting any filtering here.
IGF1_exp_dat <- IGF1_exp_dat

row.names(IGF1_exp_dat) <- NULL
print(dim(IGF1_exp_dat))

#' ==========================
#' using IEU OpenGWAS dataset
#' ==========================
# Heart Disease (CAD)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
chd_out_dat <- extract_outcome_data(
    snps = IGF1_exp_dat$SNP,
    splitsize = 10000,
    outcomes = 'ebi-a-GCST003116' # ieu-a-7: from CARDIoGRAMplusC4D, "A comprehensive 1000 Genomes–based genome-wide association meta-analysis of coronary artery disease", no UKBB sample included, mixed population; ebi-a-GCST003116: from CARDIoGRAMplusC4D, European.
)

harmonise_dat <- harmonise_data(
    exposure_dat = IGF1_exp_dat, 
    outcome_dat = chd_out_dat
)

harmonise_dat_merged <- harmonise_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure", 
                                            "beta.exposure", "se.exposure", "pval.exposure",
                                            "beta.outcome", "se.outcome", "pval.outcome")]
print(dim(harmonise_dat_merged)) 
head(harmonise_dat_merged)
message("\n\n")

# fwrite(harmonise_dat_merged, paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/05_two_sample_mr/external_analysis_for_paper_revision/CAD/", trait, "/merged_CAD_0.1.txt.gz"), row.names = F, sep = " ", quote = F)
fwrite(harmonise_dat_merged, paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/05_two_sample_mr/external_analysis_for_paper_revision/CAD/", trait, "/merged_CAD_1.txt.gz"), row.names = F, sep = " ", quote = F)

message(paste0("Finished running on pvalue cutoff ", p, ".\n\n=================================\n\n"))

