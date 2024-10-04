#' This code is to map the SNP location to the file. mamba activate bsgenome
#' 
#' @author Shi Yujia
#' @date 2024.08.04

suppressPackageStartupMessages({
    library(data.table)
    library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
})

snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37

IGF1 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1.QC.gz")

my_rsids <- IGF1$SNP[1:10]
my_snps <- snpsById(snps, my_rsids)

my_rsids <- IGF1$SNP
my_snps <- snpsById(snps, my_rsids, ifnotfound="drop")
my_snps_df <- as.data.frame(my_snps)
my_snps_df <- my_snps_df[, c("seqnames", "pos", "RefSNP_id")]
colnames(my_snps_df) <- c("CHR", "BP", "SNP")

IGF1_merge <- merge(IGF1[, c("SNP", "A1", "A2", "N", "SE", "P", "BETA", "MAF", "EAF")], my_snps_df, by = "SNP")
dim(IGF1_merge)
dim(my_snps_df)

IGF1_merge <- IGF1_merge[, c("CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "BETA", "MAF", "EAF")]
IGF1_merge$CHR <- as.integer(IGF1_merge$CHR)

fwrite(IGF1_merge, "/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1.QC.withPos.gz", sep = "\t")
