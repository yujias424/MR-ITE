#' This code is to test the correlation between IV-CF and DRIV.
#' 
#' @author Shi Yujia
#' @date 2024.09.19

suppressPackageStartupMessages({
    library(data.table)
})

traits <- c("LDL", "HDL", "TC", "TG", "IGF1", "CRP")

res <- matrix(nrow = 12, ncol = 6)

index <- 0

for (i in traits){

    ivcf_continuous_model1 <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", i, "/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model1b.csv"))
    driv_continuous_model1 <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", i, "/CAD/01_individual_treatment_effect/driv_full_continuousW_te_ul_model1b.csv"))

    ivcf_continuous_model3 <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", i, "/CAD/01_individual_treatment_effect/ivgrf_full_continuousW_pvalue_model3.csv"))
    driv_continuous_model3 <- fread(paste0("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/", i, "/CAD/01_individual_treatment_effect/driv_full_continuousW_te_ul_model3.csv"))

    res.model1.pearson <- cor.test(ivcf_continuous_model1$tau, driv_continuous_model1$point, method = "pearson")
    res.model3.pearson <- cor.test(ivcf_continuous_model3$tau, driv_continuous_model3$point, method = "pearson")

    res.model1.spearman <- cor.test(ivcf_continuous_model1$tau, driv_continuous_model1$point, method = "spearman")
    res.model3.spearman <- cor.test(ivcf_continuous_model3$tau, driv_continuous_model3$point, method = "spearman")

    res[(index*2+1), 1] <- i
    res[(index*2+2), 1] <- i
    res[(index*2+1), 2] <- "Model 1"
    res[(index*2+2), 2] <- "Model 2"
    res[(index*2+1), 3] <- res.model1.pearson$estimate
    res[(index*2+2), 3] <- res.model3.pearson$estimate
    res[(index*2+1), 4] <- res.model1.pearson$`p.value`
    res[(index*2+2), 4] <- res.model3.pearson$`p.value`
    res[(index*2+1), 5] <- res.model1.spearman$estimate
    res[(index*2+2), 5] <- res.model3.spearman$estimate
    res[(index*2+1), 6] <- res.model1.spearman$`p.value`
    res[(index*2+2), 6] <- res.model3.spearman$`p.value`

    index <- index + 1
}

res <- as.data.frame(res)
colnames(res) <- c("Trait", "Covarite Set", "Pearson Correlation", "P-value", "Spearman Correlation", "P-value")
fwrite(res, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/PearsonSpearmanCorr.csv")
