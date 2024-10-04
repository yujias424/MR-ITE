# https://github.com/privefl/bigstatsr/issues/90

suppressPackageStartupMessages({
    library(data.table)
    library(bigsnpr)
})

# ALready Done: Generate the bed file
message(paste0("Currently Running chr", chr))
snp_readBed2(paste0("/mnt/md0/public_data/UKBB_data/genotype/GWAS_ukbb_maf0.001_info0.8/bed/allchrs/ukb_allchrs.bed"),
            backingfile = paste0("/mnt/md1/share/ukbb_geno_rds_allchr/ukb_allchrs"), ncores = 40)