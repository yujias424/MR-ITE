#' This code is to map the SNP location to the file.
#' 
#' @author Shi Yujia
#' @date 2024.08.04

suppressPackageStartupMessages({
    library(data.table)
    library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
})

snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37

CystatinC <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC.QC.gz")

my_rsids <- CystatinC$SNP[1:10]
my_snps <- snpsById(snps, my_rsids)

my_rsids <- CystatinC$SNP
my_snps <- snpsById(snps, my_rsids, ifnotfound="drop")
my_snps_df <- as.data.frame(my_snps)
my_snps_df <- my_snps_df[, c("seqnames", "pos", "RefSNP_id")]
colnames(my_snps_df) <- c("CHR", "BP", "SNP")

CystatinC_merge <- merge(CystatinC, my_snps_df, by = "SNP")
dim(CystatinC_merge)
dim(my_snps_df)

CystatinC_merge <- CystatinC_merge[, c("CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "BETA", "MAF", "N")]
CystatinC_merge$CHR <- as.integer(CystatinC_merge$CHR)

fwrite(CystatinC_merge, "/home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC.QC.withPos.gz", sep = "\t")

# aaa <- data.table::fread("/home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC.QC.withPos.gz")
