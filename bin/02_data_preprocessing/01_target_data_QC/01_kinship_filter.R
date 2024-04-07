#' This code is to filter the patients with high kinship (relatedness).
#' https://wiki.arrayserver.com/wiki/index.php/KinshipCutoff
#' A default kinship value of 0.0442 is set to identify pairs of subjects with 3rd-degree or closer relationships.
#' an estimated kinship coefficient range >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] corresponds to duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree relationships respectively
#' According to Choi's tutorial, we use the cutoff 0.0884 
#' 
#'
#' @author Shi Yujia
#' @date 2022.09.01

suppressPackageStartupMessages({
    library(tidyverse)
})

ukb.dat <- read.csv("/home/yujia/Project/2022-09-01-individual_MR/dat/01_target_data/relatedness_kinship/ukb_rel_a28732_s488191.dat", header = T, sep = " ")
head(ukb.dat)
dim(ukb.dat)

ukb.dat.invalid <- ukb.dat[ukb.dat$Kinship > 0.0884 | ukb.dat$Kinship<0, ]
row.names(ukb.dat.invalid) <- NULL
colnames(ukb.dat.invalid)[1:2] <- c("ID_1", "ID_2")
head(ukb.dat.invalid)
dim(ukb.dat.invalid)

tmp.sample <- read.table("/home/yujia/Project/2022-09-01-individual_MR/dat/01_target_data/sample/ukb_imp_chr1_v3.sample", header = T)

ukb.dat.invalid <- rbind(data.frame("ID_1" = c(0), "ID_2" = c(0)), ukb.dat.invalid[, c("ID_1", "ID_2")])
row.names(ukb.dat.invalid) <- NULL

tmp.sample.invalid <- tmp.sample[tmp.sample$ID_2 %in% ukb.dat.invalid$ID_2, ]
tmp.sample.valid <- tmp.sample[!(tmp.sample$ID_2 %in% ukb.dat.invalid$ID_2), ]
tmp.sample.valid <- rbind(data.frame("ID_1" = c(0), "ID_2" = c(0), "missing" = 0, "sex" = "D"), tmp.sample.valid)
row.names(tmp.sample.valid) <- NULL
row.names(tmp.sample.invalid) <- NULL
head(tmp.sample.valid)
head(tmp.sample.invalid)
dim(tmp.sample.valid)
dim(tmp.sample.invalid)

write.table(tmp.sample.invalid, file = "/home/yujia/Project/2022-09-01-individual_MR/dat/01_target_data/relatedness_kinship/relatedness_invalid.sample", quote = F, row.names = F)
