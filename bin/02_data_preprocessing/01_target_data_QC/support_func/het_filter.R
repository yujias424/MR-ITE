#' This is the code to perform QC on heterozygosity rates
#' individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean will be removed.
#' This code will be used in 01_plink_qc.py
#'
#' @author: Choi Shingwan
#' @date: 2020

suppressPackageStartupMessages({
  library(tidyverse)
})

setwd("/home/yujia/Project/2022-09-01-individual_MR/dat/01_target_data/GWAS/GWAS_dat/")

for (i in 1:23) {
  if (!file.exists(paste0("./chr", i, "_QC_valid.sample"))) {
    if (file.exists(paste0("./chr", i, "_QC.het"))) {
      dat <- read.table(paste0("./chr", i, "_QC.het"), header = F)
      m <- mean(dat$V6) # Calculate the mean
      s <- sd(dat$V6) # Calculate the SD
      valid <- subset(dat, V6 <= m + 3 * s & V6 >= m - 3 * s) # Get any samples with F coefficient within 3 SD of the population mean

      valid_write <- valid[, c(1, 2)]
      colnames(valid_write) <- c("ID_1", "ID_2")

      tmp.sample <- read.table(paste0("./chr", i, "_QC.sample"), header = T)
      valid_write <- rbind(data.frame("ID_1" = c(0), "ID_2" = c(0)), valid_write)
      row.names(valid_write) <- NULL
      tmp.sample.valid <- tmp.sample[tmp.sample$ID_1 %in% valid_write$ID_1, ]
      tmp.sample.invalid <- tmp.sample[!(tmp.sample$ID_1 %in% valid_write$ID_1), ]
      tmp.sample.invalid <- rbind(data.frame("ID_1" = c(0), "ID_2" = c(0), "missing" = 0, "sex" = "D"), tmp.sample.invalid)

      print(paste0("Raw patients num is ", dim(dat)[1], "."))
      print(paste0("HR QC patients num is ", (dim(tmp.sample.valid)[1]-1), "."))
      print(paste0("Removed patients num is ", dim(dat)[1] - (dim(tmp.sample.valid)[1]-1), "."))
      print(paste0("Removed patients num is ", (dim(tmp.sample.invalid)[1]-1), "."))

      write.table(tmp.sample.valid, paste0("./chr", i, "_QC_valid.sample"), quote = F, row.names = F) # print FID and IID for valid samples
      write.table(tmp.sample.invalid, paste0("./chr", i, "_QC_invalid.sample"), quote = F, row.names = F) # print FID and IID for valid samples
    }
  }
}