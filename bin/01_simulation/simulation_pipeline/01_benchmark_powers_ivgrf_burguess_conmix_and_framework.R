#' This code is to benchmark the performance between our proposed framework and the traditional causal forest using simulation dataset,
#' we aim to check the correlation between true tau and the predicted tau.
#' 
#' @author Yujia Shi
#' @date 2023.08.25

suppressPackageStartupMessages({
  
  source("/home/yujia/Project/2023-07-20-individual_MR/01_simulation/support_func_data_generation/simulation_powers.R")
  source("/home/yujia/Project/2023-07-20-individual_MR/01_simulation/support_func_test_heterogeneity/test_heterogeneity.R")

  library(grf)
  library(DescTools)
  library(data.table)
  library(Metrics)
  library(MendelianRandomization)
  library(glmnet)
  library(caret)
  library(ivreg)
  library(R.utils)
  library(reticulate)

})

use_condaenv("mr")
np <- import("numpy")
econml <- import("econml")

# parameter settings
features.p <- c(50, 50, 40, 40, 30, 30, 20, 20, 50, 50, 40, 40, 30, 30, 20, 20)
mu.mu <- c(function.8, function.5, function.4, function.7, function.3, function.1, function.2, function.6, function.8, function.5, function.4, function.7, function.3, function.1, function.2, function.6)
tau.tau <- c(function.1, function.2, function.3, function.4, function.5, function.6, function.7, function.8, function.1, function.2, function.3, function.4, function.5, function.6, function.7, function.8)

n.experiments <- 30
seed <- 1
set.seed(seed)

pleiotropy_scen <- c(2, 3, 4) # must be 2,3,4 sinece ps=1 means xi=0, a=0
invalid_snps_n_set <- c(20, 40, 60)

for (ps in pleiotropy_scen){
  for (isns in invalid_snps_n_set){
    
    scenario.index <- ps-1
    message(paste0("Running pleiotropy scenario ", scenario.index))
    message(paste0("With invalid SNPs ", isns))
    message("              ")
  
    for (i in 1:8){ 

      message(paste0("Running scenario ", i))
      res.mat <- matrix(0, nrow = n.experiments, ncol = 11)
      
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
        
        # calculate the conmix accuracy.
        predicted_value <- rep(0, 100)
        predicted_value[which(colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs)] <- 1
        predicted_value[which(!colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs)] <- 0
        predicted_value <- factor(predicted_value)
        expected_value <- factor(c(rep(0, isns), rep(1, 100-isns)))
        confusion_mat <- confusionMatrix(data = predicted_value, reference = expected_value)

        # IVW keep pleiotropy 
        MRInputObject_ivw_keeppleiotropy <- mr_input(snps = colnames(dat1$Z),
                                                      bx = out.zx$coefficients[, 1],
                                                      bxse = out.zx$coefficients[, 2],
                                                      by =  out.zy$coefficients[, 1],
                                                      byse = out.zy$coefficients[, 2])
        ivw_keeppleiotropy <- mr_ivw(MRInputObject_ivw_keeppleiotropy)

        # IVW remove pleiotropy 
        MRInputObject_ivw_removepleiotropy <- mr_input(snps = colnames(dat1$Z)[which(colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs)],
                                                        bx = out.zx$coefficients[which(colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs), 1],
                                                        bxse = out.zx$coefficients[which(colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs), 2],
                                                        by =  out.zy$coefficients[which(colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs), 1],
                                                        byse = out.zy$coefficients[which(colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs), 2])
        ivw_removepleiotropy <- mr_ivw(MRInputObject_ivw_removepleiotropy)

        # Regular Approach (Causal Forest)
        # Causal Forest 
        causal.tau.forest <- causal_forest(X = dat1$`X`[30001:40000, ], Y = dat1$`Y`[30001:40000], 
                                          W = dat1$`T`[30001:40000], 
                                          sample.fraction = 0.1, num.trees = 3000, min.node.size = 5,
                                          num.threads = 40)
        cf.pred <- predict(causal.tau.forest, estimate.variance = T, num.threads = 40)
        cf.pred.estimate <- cf.pred$predictions
        cf.pred.std <- sqrt(cf.pred$variance.estimates)

        # propensity score
        regression_forest_propensity <- regression_forest(X = dat1$`X`[30001:40000, ], Y = dat1$T[30001:40000], 
                                                          num.threads = 40, 
                                                          sample.fraction = 0.1, num.trees = 3000, min.node.size = 5)
        ps.pred <- predict(regression_forest_propensity, estimate.variance = T, num.threads = 40)
        ps.pred.estimate <- ps.pred$predictions
        ps.pred.std <- sqrt(ps.pred$variance.estimates)

        # ATE
        ate.cf <- average_treatment_effect(causal.tau.forest)
        ate.cf.estimate <- ate.cf[1]
        ate.cf.std <- ate.cf[2]*100 # SD = SE*sqrt(n)       

        # Calculate the PRS (Regression Forest)
        regression_forest_allsnps <- regression_forest(X = dat1$Z[20001:30000, ], Y = dat1$T[20001:30000], 
                                                        num.threads = 40, 
                                                        sample.fraction = 0.1, num.trees = 3000, min.node.size = 5)
        regression_forest_fewsnps <- regression_forest(X = dat1$Z[20001:30000, colnames(dat1$Z)[colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs]], Y = dat1$T[20001:30000], 
                                                        num.threads = 40, 
                                                        sample.fraction = 0.1, num.trees = 3000, min.node.size = 5)

        Z_allsnps <- predict(regression_forest_allsnps, newdata = dat1$`Z`[30001:40000, ], num.threads = 40)$predictions
        Z_fewsnps <- predict(regression_forest_fewsnps, newdata = dat1$Z[30001:40000, colnames(dat1$Z)[colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs]], num.threads = 40)$predictions

        # Calculate the PRS score.
        # std
        prs.fewsnps <- (Z_fewsnps-mean(Z_fewsnps))/sd(Z_fewsnps)
        prs.allsnps <- (Z_allsnps-mean(Z_allsnps))/sd(Z_allsnps)

        # MR-ITE (IVgrf approach)
        iv.tau.forest.allsnps <- instrumental_forest(X = dat1$`X`[30001:40000, ], Y = dat1$`Y`[30001:40000], 
                                                    W = dat1$`T`[30001:40000], Z = prs.allsnps, 
                                                    sample.fraction = 0.1, num.threads = 40, 
                                                    num.trees = 3000, min.node.size = 5)

        iv.tau.forest.fewsnps <- instrumental_forest(X = dat1$`X`[30001:40000, ], Y = dat1$`Y`[30001:40000], 
                                                    W = dat1$`T`[30001:40000], Z = prs.fewsnps, 
                                                    sample.fraction = 0.1, num.threads = 40, 
                                                    num.trees = 3000, min.node.size = 5)

        iv.tau.forest.prs <- instrumental_forest(X = dat1$`X`[30001:40000, ], Y = dat1$`Y`[30001:40000], 
                                                W = dat1$`T`[30001:40000], Z = dat1$`prs`[30001:40000],
                                                sample.fraction = 0.1, num.threads = 40, 
                                                num.trees = 3000, min.node.size = 5)

        iv.pred.allsnps <- predict(iv.tau.forest.allsnps, estimate.variance = T, num.threads = 40)
        iv.pred.fewsnps <- predict(iv.tau.forest.fewsnps, estimate.variance = T, num.threads = 40)
        iv.pred.prs <- predict(iv.tau.forest.prs, estimate.variance = T, num.threads = 40)

        iv.pred.allsnps.estimate <- iv.pred.allsnps$predictions
        iv.pred.fewsnps.estimate <- iv.pred.fewsnps$predictions
        iv.pred.prs.estimate <- iv.pred.prs$predictions

        iv.pred.allsnps.std <- sqrt(iv.pred.allsnps$variance.estimates)
        iv.pred.fewsnps.std <- sqrt(iv.pred.fewsnps$variance.estimates)
        iv.pred.prs.std <- sqrt(iv.pred.prs$variance.estimates)
        
        # MR-ITE (DRIV approach)
        X_py <- as.matrix(dat1$`X`[30001:40000, ])
        rownames(X_py) <- NULL
        colnames(X_py) <- NULL
        Y_py <- dat1$`Y`[30001:40000]
        W_py <- dat1$`T`[30001:40000]
        Z_py <- dat1$`prs`[30001:40000]
        Z_py_few <- prs.fewsnps
        Z_py_all <- prs.allsnps

        est_driv_continuousW <- econml$iv$dr$ForestDRIV(projection=F, discrete_treatment=F, discrete_instrument=F,
                                                        n_estimators=as.integer(3000), min_samples_leaf=as.integer(5), max_samples=0.1, cov_clip=1,
                                                        random_state=as.integer(309), n_jobs=as.integer(40))
        est_driv_continuousW$fit(Y_py, W_py, Z=Z_py_all, X=X_py)
        point_driv_prsall <- est_driv_continuousW$effect_inference(X_py)
        point_driv_prsall_pred <- point_driv_prsall$point_estimate
        point_driv_prsall_std <- sqrt(point_driv_prsall$var)

        est_driv_continuousW$fit(Y_py, W_py, Z=Z_py_few, X=X_py)
        point_driv_prsfew <- est_driv_continuousW$effect_inference(X_py)
        point_driv_prsfew_pred <- point_driv_prsfew$point_estimate
        point_driv_prsfew_std <- sqrt(point_driv_prsfew$var)

        est_driv_continuousW$fit(Y_py, W_py, Z=Z_py, X=X_py)
        point_driv_prstrue <- est_driv_continuousW$effect_inference(X_py)
        point_driv_prstrue_pred <- point_driv_prstrue$point_estimate
        point_driv_prstrue_std <- sqrt(point_driv_prstrue$var)

        taus.mat.estimate.driv <- data.frame("ITE_DRIV_all" = point_driv_prsall_pred, 
                                              "ITE_DRIV_few" = point_driv_prsfew_pred, 
                                              "ITE_DRIV_trueZ" = point_driv_prstrue_pred, 
                                              "true_tau" = dat1$tau[30001:40000])
        taus.mat.std.driv <- data.frame("ITE_DRIV_all" = point_driv_prsall_std, 
                                        "ITE_DRIV_few" = point_driv_prsfew_std, 
                                        "ITE_DRIV_trueZ" = point_driv_prstrue_std)

        fwrite(taus.mat.estimate.driv, paste0("/home/yujia/Project/2023-07-20-individual_MR/res/01_simulation/tau/estimate/", scenario.index, "_", isns, "_", i, "_", ne,"_estimate_driv.csv.gz"))
        fwrite(taus.mat.std.driv, paste0("/home/yujia/Project/2023-07-20-individual_MR/res/01_simulation/tau/std/", scenario.index, "_", isns, "_", i, "_", ne,"_std_driv.csv.gz")) 

        # LATE
        tmp.ivreg.allsnps <- dat1$`X`[30001:40000, ]
        tmp.ivreg.allsnps$Y <- dat1$`Y`[30001:40000]
        tmp.ivreg.allsnps$W <- dat1$`T`[30001:40000]
        tmp.ivreg.allsnps$Z <- prs.allsnps
        ivreg.res.as <- ivreg(Y ~ . - W - Z| W | Z, data = tmp.ivreg.allsnps)
        late.allsnps <- summary(ivreg.res.as)$coefficients
        late.allsnps.estimate <- late.allsnps[2,1]
        late.allsnps.std <- late.allsnps[2,2] * 100 # SD = SE * sqrt(n)

        tmp.ivreg.fewsnps <- dat1$`X`[30001:40000, ]
        tmp.ivreg.fewsnps$Y <- dat1$`Y`[30001:40000]
        tmp.ivreg.fewsnps$W <- dat1$`T`[30001:40000]
        tmp.ivreg.fewsnps$Z <- prs.fewsnps
        ivreg.res.fs <- ivreg(Y ~ . - W - Z| W | Z, data = tmp.ivreg.fewsnps)
        late.fewsnps <- summary(ivreg.res.fs)$coefficients
        late.fewsnps.estimate <- late.fewsnps[2,1]
        late.fewsnps.std <- late.fewsnps[2,2] * 100 # SD = SE * sqrt(n)

        tmp.ivreg.trueZ <- dat1$`X`[30001:40000, ]
        tmp.ivreg.trueZ$Y <- dat1$`Y`[30001:40000]
        tmp.ivreg.trueZ$W <- dat1$`T`[30001:40000]
        tmp.ivreg.trueZ$Z <- dat1$`prs`[30001:40000]
        ivreg.res.fs <- ivreg(Y ~ . - W - Z| W | Z, data = tmp.ivreg.trueZ)
        late.trueZ <- summary(ivreg.res.fs)$coefficients
        late.trueZ.estimate <- late.trueZ[2,1]
        late.trueZ.std <- late.trueZ[2,2] * 100 # SD = SE * sqrt(n)

        message("MSE: ")
        message(paste0("ATE: ", mse(rep(ate.cf.estimate, length(cf.pred.estimate)), dat1$tau[30001:40000])))
        message(paste0("LATE IVReg (keep pleiotropy): ", mse(rep(late.allsnps.estimate, length(cf.pred.estimate)), dat1$tau[30001:40000])))
        message(paste0("LATE IVReg (remove pleiotropy): ", mse(rep(late.fewsnps.estimate, length(cf.pred.estimate)), dat1$tau[30001:40000])))
        message(paste0("LATE IVReg (true Z): ", mse(rep(late.trueZ.estimate, length(cf.pred.estimate)), dat1$tau[30001:40000])))
        message(paste0("LATE IVW (keep pleiotropy): ", mse(rep(ivw_keeppleiotropy@Estimate, length(cf.pred.estimate)), dat1$tau[30001:40000])))
        message(paste0("LATE IVW (remove pleiotropy): ", mse(rep(ivw_removepleiotropy@Estimate, length(cf.pred.estimate)), dat1$tau[30001:40000])))
        message(paste0("ITE: ", mse(cf.pred.estimate, dat1$tau[30001:40000])))
        message(paste0("ITE-IVCF (keep pleiotropy): ", mse(iv.pred.allsnps.estimate, dat1$tau[30001:40000])))
        message(paste0("ITE-IVCF (remove pleiotropy): ", mse(iv.pred.fewsnps.estimate, dat1$tau[30001:40000])))
        message(paste0("ITE-IVCF (true Z): ", mse(iv.pred.prs.estimate, dat1$tau[30001:40000])))
        message(paste0("ITE-DRIV (keep pleiotropy): ", mse(point_driv_prsall_pred, dat1$tau[30001:40000])))
        message(paste0("ITE-DRIV (remove pleiotropy): ", mse(point_driv_prsfew_pred, dat1$tau[30001:40000])))
        message(paste0("ITE-DRIV (true Z): ", mse(point_driv_prstrue_pred, dat1$tau[30001:40000])))
        message("\nCorrelation Test:   ")
        message(paste0("ITE: ", cor.test(dat1$tau[30001:40000], cf.pred.estimate)$estimate))
        message(paste0("ITE-IVCF (keep pleiotropy): ", cor.test(dat1$tau[30001:40000], iv.pred.allsnps.estimate)$estimate))
        message(paste0("ITE-IVCF (remove pleiotropy): ", cor.test(dat1$tau[30001:40000], iv.pred.fewsnps.estimate)$estimate))
        message(paste0("ITE-IVCF (true Z): ", cor.test(dat1$tau[30001:40000], iv.pred.prs.estimate)$estimate))
        message(paste0("ITE-DRIV (keep pleiotropy): ", cor.test(dat1$tau[30001:40000], point_driv_prsall_pred)$estimate))
        message(paste0("ITE-DRIV (remove pleiotropy): ", cor.test(dat1$tau[30001:40000], point_driv_prsfew_pred)$estimate))
        message(paste0("ITE-DRIV (true Z): ", cor.test(dat1$tau[30001:40000], point_driv_prstrue_pred)$estimate))
        message("              ")
        message(paste0("The conmix accuracy is ", confusion_mat$overall[1]))
        message("              ")

        res.mat[ne, 1] <- mse(cf.pred.estimate, dat1$tau[30001:40000])
        res.mat[ne, 2] <- mse(iv.pred.allsnps.estimate, dat1$tau[30001:40000])
        res.mat[ne, 3] <- mse(iv.pred.fewsnps.estimate, dat1$tau[30001:40000])
        res.mat[ne, 4] <- mse(point_driv_prsall_pred, dat1$tau[30001:40000])
        res.mat[ne, 5] <- mse(point_driv_prsfew_pred, dat1$tau[30001:40000])
        res.mat[ne, 6] <- mse(rep(ate.cf.estimate, length(cf.pred.estimate)), dat1$tau[30001:40000])
        res.mat[ne, 7] <- mse(rep(late.allsnps.estimate, length(cf.pred.estimate)), dat1$tau[30001:40000])
        res.mat[ne, 8] <- mse(rep(late.fewsnps.estimate, length(cf.pred.estimate)), dat1$tau[30001:40000])
        res.mat[ne, 9] <- mse(rep(ivw_keeppleiotropy@Estimate, length(cf.pred.estimate)), dat1$tau[30001:40000])
        res.mat[ne, 10] <- mse(rep(ivw_removepleiotropy@Estimate, length(cf.pred.estimate)), dat1$tau[30001:40000])
        res.mat[ne, 11] <- confusion_mat$overall[1]

        taus.mat.estimate <- data.frame("ATE_CF" = rep(ate.cf.estimate, length(cf.pred.estimate)), 
                                        "LATE_all" = rep(late.allsnps.estimate, length(cf.pred.estimate)), "LATE_few" = rep(late.fewsnps.estimate, length(cf.pred.estimate)), "LATE_trueZ" = rep(late.trueZ.estimate, length(cf.pred.estimate)),
                                        "ITE_CF" = cf.pred.estimate, 
                                        "ITE_IVCF_all" = iv.pred.allsnps.estimate, "ITE_IVCF_few" = iv.pred.fewsnps.estimate, "ITE_IVCF_trueZ" = iv.pred.prs.estimate, 
                                        "true_tau" = dat1$tau[30001:40000],
                                        "true_propensity_score" = dat1$lt[30001:40000],
                                        "propensity" = ps.pred.estimate)
        taus.mat.std <- data.frame("ATE_CF" = rep(ate.cf.std, length(cf.pred.estimate)), 
                                  "LATE_all" = rep(late.allsnps.std, length(cf.pred.estimate)), "LATE_few" = rep(late.fewsnps.std, length(cf.pred.estimate)), "LATE_trueZ" = rep(late.trueZ.std, length(cf.pred.estimate)),
                                  "ITE_CF" = cf.pred.std, 
                                  "ITE_IVCF_all" = iv.pred.allsnps.std, "ITE_IVCF_few" = iv.pred.fewsnps.std, "ITE_IVCF_trueZ" = iv.pred.prs.std,
                                  "propensity" = ps.pred.std)

        fwrite(taus.mat.estimate, paste0("/home/yujia/Project/2023-07-20-individual_MR/res/01_simulation/tau/estimate/", scenario.index, "_", isns, "_", i, "_", ne,"_estimate.csv.gz"))
        fwrite(taus.mat.std, paste0("/home/yujia/Project/2023-07-20-individual_MR/res/01_simulation/tau/std/", scenario.index, "_", isns, "_", i, "_", ne,"_std.csv.gz")) 

        ne <- ne + 1
        rm(MRAllObject_conmix) # since we use MRAllObject_conmix == NULL as the criterion to jump

      }

      res.mat <- as.data.frame(res.mat)
      colnames(res.mat) <- c("CFgrf_MSE", 
                              "IVgrf_All_SNPs_MSE", "IVgrf_Valid_SNPs_MSE",
                              "DRIV_All_SNPs_MSE", "DRIV_Valid_SNPs_MSE",
                              "ATE_MSE", 
                              "LATE_IVreg_All_SNPs_MSE", "LATE_IVreg_Valid_SNPs_MSE",
                              "LATE_IVW_All_SNPs_MSE", "LATE_IVW_Valid_SNPs_MSE", "Conmix_aacuracy")

      fwrite(res.mat, paste0("/home/yujia/Project/2023-07-20-individual_MR/res/01_simulation/mse/", scenario.index, "_", isns, "_", i, "_res.csv"))

    }
  }
}

