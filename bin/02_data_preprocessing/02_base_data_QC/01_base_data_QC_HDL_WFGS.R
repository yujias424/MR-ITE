#' This code is to download the WFGS HDL-C dataset. https://gwas.mrcieu.ac.uk/datasets/ieu-b-4844/
#' 
#' @author Yujia Shi
#' @date 2022.09.01

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

hdl_exp_dat <- extract_instruments(
    outcomes = 'ieu-b-4844', # ieu-b-4844: from Within family GWAS consortium, European. 
    p1 = 0.1,
    clump = F
)

fwrite(hdl_exp_dat, file = "/home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/HDL/download/HDL_WFGS.txt.gz", sep="\t", row.names = F)
