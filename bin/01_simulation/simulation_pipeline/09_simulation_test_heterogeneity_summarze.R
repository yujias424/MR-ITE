#' This code is to summarize the results from test heterogeneity.
#' 
#' @author Shi Yujia
#' @date 2022.09.01

suppressPackageStartupMessages({
    library(data.table)
})

ps_scenario <- c(2,3,4)
tau_scenario <- c(1,2,3,4,5,6,7,8)

res_mat <- matrix(0, ncol = 6, nrow = 8)

for (ps in ps_scenario){
    for (ts in tau_scenario){

        tmp_dat <- as.data.frame(fread(paste0("/home/yujia/Project/2023-07-20-individual_MR/res/01_simulation/framework_sims_HTE/", ps, "_40_", ts, "_res.csv")))

        var_test <- sum(tmp_dat$`perm-var`<0.05)/dim(tmp_dat)[1]
        taurisk_test <- sum(tmp_dat$`perm-tau-risk`<0.05)/dim(tmp_dat)[1]
        
        res_mat[ts, ps-1] <- var_test
        res_mat[ts, ps-1+3] <- taurisk_test

    }
}

