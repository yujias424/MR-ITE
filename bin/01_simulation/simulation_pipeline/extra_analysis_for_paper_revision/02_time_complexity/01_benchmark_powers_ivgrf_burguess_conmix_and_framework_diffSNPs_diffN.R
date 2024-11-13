#' This code is to perform the time-complexity analysis.
#' 
#' @author Yujia Shi
#' @date 2024.08.25

suppressPackageStartupMessages({
  
  source("/mnt/md0/yujia/project/2023-07-20-individual_MR/bin/01_simulation/support_func_data_generation/simulation_powers_diffSNP.R")
  source("/mnt/md0/yujia/project/2023-07-20-individual_MR/bin/01_simulation/support_func_test_heterogeneity/test_heterogeneity.R")

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

  library(memtime)

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
n_snps <- 500
invalid_snps_n_set <- c(100, 200, 300)

time_elapse <- c()
memory_resource <- c()
invalid_snps_array <- c()
pleiotropy_scen_array <- c()
tau_scenario_array <- c()

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
      tmp_time <- memtime({
        while (ne <= n.experiments){
          message(paste0(n.experiments-ne, " exps left."))
          skip_to_next <- FALSE

          dat1 <- scenarios.burgess(n = 200000, p = features.p[i], 
                                    mu = mu.mu[i][[1]], tau = tau.tau[i][[1]], 
                                    invalid_snps_n = isns, snps_n = n_snps, pleiotropy_scen = ps)

          td.1 <- cbind(as.data.frame(dat1$Z[1:50000,]), dat1$T[1:50000])
          colnames(td.1)[n_snps+1] <- "X"
          zx <- lm(X ~ ., data = td.1)
          out.zx <- summary(zx)

          td.2 <- cbind(as.data.frame(dat1$Z[50001:100000,]), dat1$Y[50001:100000])
          colnames(td.2)[n_snps+1] <- "Y"
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
            }, timeout = 60, onTimeout = "warning")
          }, TimeoutException = function(e){
            skip_to_next <- TRUE 
          }, error = function(e) {
            skip_to_next <- TRUE 
          })
          
          if (is.null(MRAllObject_conmix)){
            print("Timeout in running conmix.")
            next 
          } 
          
          message("Finished running Conmix.")
          message("                        ")
          
          # calculate the conmix accuracy.
          predicted_value <- rep(0, n_snps)
          predicted_value[which(colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs)] <- 1
          predicted_value[which(!colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs)] <- 0
          predicted_value <- factor(predicted_value)
          expected_value <- factor(c(rep(0, isns), rep(1, n_snps-isns)))
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
          causal.tau.forest <- causal_forest(X = dat1$`X`[150001:200000, ], Y = dat1$`Y`[150001:200000], 
                                            W = dat1$`T`[150001:200000], 
                                            sample.fraction = 0.1, num.trees = 3000, min.node.size = 25,
                                            num.threads = 40)
          cf.pred <- predict(causal.tau.forest, estimate.variance = T, num.threads = 40)
          cf.pred.estimate <- cf.pred$predictions
          cf.pred.std <- sqrt(cf.pred$variance.estimates)

          # propensity score
          regression_forest_propensity <- regression_forest(X = dat1$`X`[150001:200000, ], Y = dat1$T[150001:200000], 
                                                            num.threads = 40, 
                                                            sample.fraction = 0.1, num.trees = 3000, min.node.size = 25)
          ps.pred <- predict(regression_forest_propensity, estimate.variance = T, num.threads = 40)
          ps.pred.estimate <- ps.pred$predictions
          ps.pred.std <- sqrt(ps.pred$variance.estimates)

          # ATE
          ate.cf <- average_treatment_effect(causal.tau.forest)
          ate.cf.estimate <- ate.cf[1]
          ate.cf.std <- ate.cf[2]*100 # SD = SE*sqrt(n)       

          # Calculate the PRS (Regression Forest)
          regression_forest_allsnps <- regression_forest(X = dat1$Z[100001:150000, ], Y = dat1$T[100001:150000], 
                                                          num.threads = 40, 
                                                          sample.fraction = 0.1, num.trees = 3000, min.node.size = 25)
          regression_forest_fewsnps <- regression_forest(X = dat1$Z[100001:150000, colnames(dat1$Z)[colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs]], Y = dat1$T[100001:150000], 
                                                          num.threads = 40, 
                                                          sample.fraction = 0.1, num.trees = 3000, min.node.size = 25)
          
          Z_allsnps <- predict(regression_forest_allsnps, newdata = dat1$`Z`[150001:200000, ], num.threads = 40)$predictions
          Z_fewsnps <- predict(regression_forest_fewsnps, newdata = dat1$Z[150001:200000, colnames(dat1$Z)[colnames(dat1$Z) %in% MRAllObject_conmix@ValidSNPs]], num.threads = 40)$predictions

          # Calculate the PRS score.
          # std
          prs.fewsnps <- (Z_fewsnps-mean(Z_fewsnps))/sd(Z_fewsnps)
          prs.allsnps <- (Z_allsnps-mean(Z_allsnps))/sd(Z_allsnps)

          # MR-ITE (IVgrf approach)
          iv.tau.forest.allsnps <- instrumental_forest(X = dat1$`X`[150001:200000, ], Y = dat1$`Y`[150001:200000], 
                                                      W = dat1$`T`[150001:200000], Z = prs.allsnps, 
                                                      sample.fraction = 0.1, num.threads = 40, 
                                                      num.trees = 3000, min.node.size = 25)

          iv.tau.forest.fewsnps <- instrumental_forest(X = dat1$`X`[150001:200000, ], Y = dat1$`Y`[150001:200000], 
                                                      W = dat1$`T`[150001:200000], Z = prs.fewsnps, 
                                                      sample.fraction = 0.1, num.threads = 40, 
                                                      num.trees = 3000, min.node.size = 25)

          iv.tau.forest.prs <- instrumental_forest(X = dat1$`X`[150001:200000, ], Y = dat1$`Y`[150001:200000], 
                                                  W = dat1$`T`[150001:200000], Z = dat1$`prs`[150001:200000],
                                                  sample.fraction = 0.1, num.threads = 40, 
                                                  num.trees = 3000, min.node.size = 25)

          iv.pred.allsnps <- predict(iv.tau.forest.allsnps, estimate.variance = T, num.threads = 40)
          iv.pred.fewsnps <- predict(iv.tau.forest.fewsnps, estimate.variance = T, num.threads = 40)
          iv.pred.prs <- predict(iv.tau.forest.prs, estimate.variance = T, num.threads = 40)

          iv.pred.allsnps.estimate <- iv.pred.allsnps$predictions
          iv.pred.fewsnps.estimate <- iv.pred.fewsnps$predictions
          iv.pred.prs.estimate <- iv.pred.prs$predictions

          iv.pred.allsnps.std <- sqrt(iv.pred.allsnps$variance.estimates)
          iv.pred.fewsnps.std <- sqrt(iv.pred.fewsnps$variance.estimates)
          iv.pred.prs.std <- sqrt(iv.pred.prs$variance.estimates)

          # LATE
          tmp.ivreg.allsnps <- dat1$`X`[150001:200000, ]
          tmp.ivreg.allsnps$Y <- dat1$`Y`[150001:200000]
          tmp.ivreg.allsnps$W <- dat1$`T`[150001:200000]
          tmp.ivreg.allsnps$Z <- prs.allsnps
          ivreg.res.as <- ivreg(Y ~ . - W - Z| W | Z, data = tmp.ivreg.allsnps)
          late.allsnps <- summary(ivreg.res.as)$coefficients
          late.allsnps.estimate <- late.allsnps[2,1]
          late.allsnps.std <- late.allsnps[2,2] * 100 # SD = SE * sqrt(n)

          tmp.ivreg.fewsnps <- dat1$`X`[150001:200000, ]
          tmp.ivreg.fewsnps$Y <- dat1$`Y`[150001:200000]
          tmp.ivreg.fewsnps$W <- dat1$`T`[150001:200000]
          tmp.ivreg.fewsnps$Z <- prs.fewsnps
          ivreg.res.fs <- ivreg(Y ~ . - W - Z| W | Z, data = tmp.ivreg.fewsnps)
          late.fewsnps <- summary(ivreg.res.fs)$coefficients
          late.fewsnps.estimate <- late.fewsnps[2,1]
          late.fewsnps.std <- late.fewsnps[2,2] * 100 # SD = SE * sqrt(n)

          tmp.ivreg.trueZ <- dat1$`X`[150001:200000, ]
          tmp.ivreg.trueZ$Y <- dat1$`Y`[150001:200000]
          tmp.ivreg.trueZ$W <- dat1$`T`[150001:200000]
          tmp.ivreg.trueZ$Z <- dat1$`prs`[150001:200000]
          ivreg.res.fs <- ivreg(Y ~ . - W - Z| W | Z, data = tmp.ivreg.trueZ)
          late.trueZ <- summary(ivreg.res.fs)$coefficients
          late.trueZ.estimate <- late.trueZ[2,1]
          late.trueZ.std <- late.trueZ[2,2] * 100 # SD = SE * sqrt(n)

          message("MSE: ")
          message(paste0("ATE: ", mse(rep(ate.cf.estimate, length(cf.pred.estimate)), dat1$tau[150001:200000])))
          message(paste0("LATE IVReg (keep pleiotropy): ", mse(rep(late.allsnps.estimate, length(cf.pred.estimate)), dat1$tau[150001:200000])))
          message(paste0("LATE IVReg (remove pleiotropy): ", mse(rep(late.fewsnps.estimate, length(cf.pred.estimate)), dat1$tau[150001:200000])))
          message(paste0("LATE IVReg (true Z): ", mse(rep(late.trueZ.estimate, length(cf.pred.estimate)), dat1$tau[150001:200000])))
          message(paste0("LATE IVW (keep pleiotropy): ", mse(rep(ivw_keeppleiotropy@Estimate, length(cf.pred.estimate)), dat1$tau[150001:200000])))
          message(paste0("LATE IVW (remove pleiotropy): ", mse(rep(ivw_removepleiotropy@Estimate, length(cf.pred.estimate)), dat1$tau[150001:200000])))
          message(paste0("ITE: ", mse(cf.pred.estimate, dat1$tau[150001:200000])))
          message(paste0("ITE-IVCF (keep pleiotropy): ", mse(iv.pred.allsnps.estimate, dat1$tau[150001:200000])))
          message(paste0("ITE-IVCF (remove pleiotropy): ", mse(iv.pred.fewsnps.estimate, dat1$tau[150001:200000])))
          message(paste0("ITE-IVCF (true Z): ", mse(iv.pred.prs.estimate, dat1$tau[150001:200000])))
          message("\nCorrelation Test:   ")
          message(paste0("ITE: ", cor.test(dat1$tau[150001:200000], cf.pred.estimate)$estimate))
          message(paste0("ITE-IVCF (keep pleiotropy): ", cor.test(dat1$tau[150001:200000], iv.pred.allsnps.estimate)$estimate))
          message(paste0("ITE-IVCF (remove pleiotropy): ", cor.test(dat1$tau[150001:200000], iv.pred.fewsnps.estimate)$estimate))
          message(paste0("ITE-IVCF (true Z): ", cor.test(dat1$tau[150001:200000], iv.pred.prs.estimate)$estimate))
          message("              ")
          message(paste0("The conmix accuracy is ", confusion_mat$overall[1]))
          message("              ")

          res.mat[ne, 1] <- mse(cf.pred.estimate, dat1$tau[150001:200000])
          res.mat[ne, 2] <- mse(iv.pred.allsnps.estimate, dat1$tau[150001:200000])
          res.mat[ne, 3] <- mse(iv.pred.fewsnps.estimate, dat1$tau[150001:200000])
          res.mat[ne, 6] <- mse(rep(ate.cf.estimate, length(cf.pred.estimate)), dat1$tau[150001:200000])
          res.mat[ne, 7] <- mse(rep(late.allsnps.estimate, length(cf.pred.estimate)), dat1$tau[150001:200000])
          res.mat[ne, 8] <- mse(rep(late.fewsnps.estimate, length(cf.pred.estimate)), dat1$tau[150001:200000])
          res.mat[ne, 9] <- mse(rep(ivw_keeppleiotropy@Estimate, length(cf.pred.estimate)), dat1$tau[150001:200000])
          res.mat[ne, 10] <- mse(rep(ivw_removepleiotropy@Estimate, length(cf.pred.estimate)), dat1$tau[150001:200000])
          res.mat[ne, 11] <- confusion_mat$overall[1]

          taus.mat.estimate <- data.frame("ATE_CF" = rep(ate.cf.estimate, length(cf.pred.estimate)), 
                                          "LATE_all" = rep(late.allsnps.estimate, length(cf.pred.estimate)), "LATE_few" = rep(late.fewsnps.estimate, length(cf.pred.estimate)), "LATE_trueZ" = rep(late.trueZ.estimate, length(cf.pred.estimate)),
                                          "ITE_CF" = cf.pred.estimate, 
                                          "ITE_IVCF_all" = iv.pred.allsnps.estimate, "ITE_IVCF_few" = iv.pred.fewsnps.estimate, "ITE_IVCF_trueZ" = iv.pred.prs.estimate, 
                                          "true_tau" = dat1$tau[150001:200000],
                                          "true_propensity_score" = dat1$lt[150001:200000],
                                          "propensity" = ps.pred.estimate)
          taus.mat.std <- data.frame("ATE_CF" = rep(ate.cf.std, length(cf.pred.estimate)), 
                                    "LATE_all" = rep(late.allsnps.std, length(cf.pred.estimate)), "LATE_few" = rep(late.fewsnps.std, length(cf.pred.estimate)), "LATE_trueZ" = rep(late.trueZ.std, length(cf.pred.estimate)),
                                    "ITE_CF" = cf.pred.std, 
                                    "ITE_IVCF_all" = iv.pred.allsnps.std, "ITE_IVCF_few" = iv.pred.fewsnps.std, "ITE_IVCF_trueZ" = iv.pred.prs.std,
                                    "propensity" = ps.pred.std)

          ne <- ne + 1
          rm(MRAllObject_conmix) # since we use MRAllObject_conmix == NULL as the criterion to jump
          break
        }
      })


      time_elapse <- c(time_elapse, tmp_time$time[3])
      memory_resource <- c(memory_resource, tmp_time$memory[2,2])
      invalid_snps_array <- c(invalid_snps_array, isns)
      pleiotropy_scen_array <- c(pleiotropy_scen_array, ps)
      tau_scenario_array <- c(tau_scenario_array, i)
    
    }
  }
}

res.mat <- data.frame("Pleiotropy_Scenario" = pleiotropy_scen_array, "Invalid_SNPs" = invalid_snps_array, "Tau_Scenario" = tau_scenario_array, "Time_Elapse_Seconds" = time_elapse, "Memory_Requirement" = memory_resource)

fwrite(res.mat, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/extra_analysis_for_paper_revision/01_simulation/extra_analysis_for_paper_revision_diffSNPs_diffN/01_simulation/table/time_memory_complexity.csv")