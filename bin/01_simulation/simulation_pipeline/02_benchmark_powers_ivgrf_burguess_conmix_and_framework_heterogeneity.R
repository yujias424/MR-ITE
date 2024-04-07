#' This code is to evaluate the performance of our proposed framework in heterogeneity testing.
#' 
#' @author Yujia Shi
#' @date 2023.08.25

suppressPackageStartupMessages({
  
  source("~/Project/2023-07-20-individual_MR/bin/01_simulation/support_func_data_generation/simulation_powers.R")
  source("~/Project/2023-07-20-individual_MR/bin/01_simulation/support_func_test_heterogeneity/test_heterogeneity_update.R")

  library(grf)
  library(DescTools)
  library(data.table)
  library(Metrics)
  library(MendelianRandomization)
  library(glmnet)
  library(caret)
  library(foreach)
  library(doParallel)
  library(R.utils)

})

# parameter settings
features.p <- c(50, 50, 40, 40, 30, 30, 20, 20, 50, 50, 40, 40, 30, 30, 20, 20)
mu.mu <- c(function.8, function.5, function.4, function.7, function.3, function.1, function.2, function.6, function.8, function.5, function.4, function.7, function.3, function.1, function.2, function.6)
tau.tau <- c(function.1, function.2, function.3, function.4, function.5, function.6, function.7, function.8, function.1, function.2, function.3, function.4, function.5, function.6, function.7, function.8)

n.experiments <- 50
seed <- 1
set.seed(seed)

# pleiotropy_scen <- c(2, 3, 4) # must be 2,3,4 sinece ps=1 means xi=0, a=0
pleiotropy_scen <- c(2, 3, 4)
# invalid_snps_n_set <- c(20, 40, 60)
invalid_snps_n_set <- c(40)

# test #
for (ps in pleiotropy_scen){
  for (isns in invalid_snps_n_set){
    
    # ps <- pleiotropy_scen[2]

    message(paste0("Running pleiotropy scenario ", ps))
    message(paste0("With invalid SNPs ", isns))
   
    for (i in 1:8){
      
      # i <- 7
      message(paste0("Running scenario ", i))
      res.mat <- matrix(0, nrow = n.experiments, ncol = 2)
      
      ne <- 1
      while (ne <= n.experiments){
        
        message(paste0(n.experiments-ne, " exps left."))
        skip_to_next <- FALSE

        dat1 <- scenarios.burgess(n = 40000, p = features.p[i], 
                                  mu = mu.mu[i][[1]], tau = tau.tau[i][[1]], 
                                  invalid_snps_n = isns, snps_n = 100, pleiotropy_scen = ps)

        td.1 <- cbind(as.data.frame(dat1$Z[1:10000,]), dat1$T[1:10000])
        colnames(td.1)[101] <- "X"
        zx <- lm(X ~ ., data = td.1)
        out.zx <- summary(zx)

        td.2 <- cbind(as.data.frame(dat1$Z[10001:20000,]), dat1$Y[10001:20000])
        colnames(td.2)[101] <- "Y"
        zy <- lm(Y ~ ., data = td.2)
        out.zy <- summary(zy)

        message("Running Conmix.")
        MRInputObject_1 <- mr_input(snps = colnames(dat1$Z),
                                    bx = out.zx$coefficients[, 1],
                                    bxse = out.zx$coefficients[, 2],
                                    by =  out.zy$coefficients[, 1],
                                    byse = out.zy$coefficients[, 2])
        # MRAllObject_conmix <- mr_conmix(MRInputObject_1)

        tryCatch({
          MRAllObject_conmix <- withTimeout({
            mr_conmix(MRInputObject_1)
          }, timeout = 240, onTimeout = "warning")
        }, TimeoutException = function(e){
          skip_to_next <- TRUE 
        })

        if (is.null(MRAllObject_conmix)){
          print("Timeout in running conmix.")
          next 
        } 
        
        message("Finished running Conmix.")
        message("                        ")

        regression_forest_fewsnps <- regression_forest(X = dat1$Z[20001:30000, colnames(dat1$Z)[colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs]], Y = dat1$T[20001:30000], 
                                                        num.threads = 40, 
                                                        sample.fraction = 0.1, num.trees = 3000, min.node.size = 10)

        # denominator.fewsnps <- rowSums(as.data.frame(dat5$Z[30001:40000, colnames(dat5$Z)[colnames(dat5$Z) %in% MRAllObject_conmix@ValidSNPs]]))

        Z_fewsnps <- predict(regression_forest_fewsnps, newdata = dat1$Z[30001:40000, colnames(dat1$Z)[colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs]], num.threads = 40)$predictions

        # # avg
        # prs.fewsnps <- Z_fewsnps/denominator.fewsnps
        # prs.allsnps <- Z_allsnps/denominator.allsnps

        # # sum
        # prs.fewsnps <- Z_fewsnps
        # prs.allsnps <- Z_allsnps

        # std
        prs.fewsnps <- (Z_fewsnps-mean(Z_fewsnps))/sd(Z_fewsnps)

        iv.tau.forest.fewsnps <- instrumental_forest(X = dat1$`X`[30001:40000, ], Y = dat1$`Y`[30001:40000], 
                                                    W = dat1$`T`[30001:40000], Z = prs.fewsnps, 
                                                    sample.fraction = 0.1, num.threads = 40, 
                                                    num.trees = 3000, min.node.size = 30)
        iv.pred.fewsnps <- predict(iv.tau.forest.fewsnps)$predictions

        tau.var <- var(iv.pred.fewsnps)
        Y.hat <- iv.tau.forest.fewsnps[["Y.hat"]]
        W.hat <- iv.tau.forest.fewsnps[["W.hat"]]
        Z.hat <- iv.tau.forest.fewsnps[["Z.hat"]]

        fixed.YW.tau.risk <- assess.explained.tau.fixed.YW.risk(fitted.obj = iv.tau.forest.fewsnps, 
                                                                Y = dat1$Y[30001:40000], 
                                                                Y.hat = Y.hat, 
                                                                W = dat1$`T`[30001:40000], 
                                                                W.hat = W.hat,
                                                                Z = prs.fewsnps,
                                                                Z.hat = Z.hat)

        perm.pvals <- adaptive.permutate.covariates.testing(X = dat1$X[30001:40000, ],
                                                            Y = dat1$Y[30001:40000], 
                                                            Y.hat = Y.hat,
                                                            W = dat1$`T`[30001:40000],
                                                            W.hat = W.hat,
                                                            Z = prs.fewsnps,
                                                            Z.hat = Z.hat,
                                                            fixed.YW.tau.risk = fixed.YW.tau.risk,
                                                            tau.var = tau.var,
                                                            min.node.size = 30, 
                                                            sample.fraction = 0.1,
                                                            num.trees = 1000,
                                                            num.strap = 50)

        iv.pred.fewsnps <- predict(iv.tau.forest.fewsnps)$predictions

        perm_var_p <- perm.pvals[1]
        perm_risk_p <- perm.pvals[2]
        res.p <- c(perm_var_p, perm_risk_p)

        message(paste0("pval of variance by permutation covariates:", perm.pvals[1]))
        message(paste0("pval of tau.risk (fixed YW) by permutation covariates:", perm.pvals[2]))  
        message(paste0("\n------------------------------\n"))  

        res.mat[ne, ] <- res.p
        ne <- ne + 1
        # break
      }
      # break
      colnames(res.mat) <- c("perm-var", "perm-tau-risk")
      res.mat <- as.data.frame(res.mat)
      fwrite(res.mat, paste0("~/Project/2023-07-20-individual_MR/res/01_simulation/framework_sims_HTE/", ps, "_", isns, "_", i, "_res.csv"))

    }
    # break
  }
  # break
}
