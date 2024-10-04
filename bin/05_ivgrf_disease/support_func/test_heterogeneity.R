#' This code is to save the testing method that kai developed before,
#' our goal is to test whether these methods perform well in instrumental setting.
#' 
#' @author: Yujia Shi
#' @date: 2021.12.01

#' ===================
#' split half testing
#' ===================

split_half_testing <- function(covariates, Y, Z,
                               treatment,
                               binary = F,
                               is_save = T,
                               save_split = T,
                               file_prefix = NULL,
                               col_names = NULL,
                               no_repeats = 10,
                               seed = NULL, # seed will be set within run.hte call
                               n_core = 30,
                               u = 2, 
                               num.trees = 50000,
                               min_node_size = 10,
                               ci.group.size = 5,
                               partial_simes_thre = 0.05,
                               random_rep_seed = TRUE) {
  
  #' Following Kai's idea, we keep using split half testing to assess the model stability, the idea is that if significant ite is detected in 
  #' one model, it should also be detected in another model.
  #' This code is to perform split half testing
  #' 
  #' @param covariates: the X covariates matrix
  #' @param Y: outcome
  #' @param Z: instrument
  #' @param treatment: treatment assignment
  #' @param binary: binary treatment or not
  #' @param is_save: whether to save the results
  #' @param col_names: c("simes.pval", "partial.simes.pval", "pearson.estimate", "pearson.pvalue", "kendall.estimate", "kendall.pvalue", "spearman.estimate", "spearman.pvalue", "fisher.pval", "t.test.a.pval", "t.test.b.pval")
  #' @param no_repeats: how many times you want to perform the test.
  #' @param seed: seed
  #' @param n_core: core number you want to use.
  #' @param random_rep_seed: whether change seed in each repeat.
  #' @param num.trees: number of trees you want to use in the ivgrf
  #' @param min_node_size: min.node.size setting for ivgrf
  #' @param ci.group.size: ci.group.size setting for causal forest 
  #' @param u: partial simes test u setting. can be either a number or a vector of u.
  #' 
  #' @return the aggregate results
  
  cf.estimator <- ivgrf
  
  # note: if is_save is TRUE, then file_prefix and col_names must be provided.
  # when is_save turns on, results for every repeat will be saved. 
  # if save_split is true, then observation level result will be stored.
  # Here file names is a file prefix for saving
  if ((is_save | save_split) && is.null(file_prefix)) stop("if is_save or save_split is T, then file_prefix must be provided.")
  
  no.obs <- dim(covariates)[1]
   
  # For loop
  correlation_matrix <- NULL
  for (i in seq(no_repeats)) {
    observation_result.a <- matrix(0, nrow = no.obs, ncol = 4)
    observation_result.b <- matrix(0, nrow = no.obs, ncol = 4)
    
    if (binary) {
      treat <- seq(no.obs)[treatment == 1]
      control <- seq(no.obs)[treatment == 0]
      sampled_treat <- floor(length(treat)/2)
      trainId <- c(sample(treat, sampled_treat), sample(control, floor(no.obs/2) - sampled_treat))
    }else{
      trainId <- sample(1: no.obs, floor(no.obs/2), replace = FALSE)
    }
    
    no.obs.train <- length(trainId)
    no.obs.test <- no.obs - no.obs.train
    
    tau.forest.train <- cf.estimator(X = covariates[trainId, ], Y = Y[trainId], W = treatment[trainId], Z = Z[trainId], num.trees = num.trees, 
                                      min.node.size = min_node_size, ci.group.size = ci.group.size, seed = seed) # we would prefer a larger num.trees here, to guarantee a better CI as we don't want CI to be too conservative.
    tau.pred.train.a <- predict(tau.forest.train, estimate.variance = T, num.threads = n_core)
    tau.pred.test.a <- predict(tau.forest.train, newdata = covariates[-trainId, ], estimate.variance = T, num.threads = n_core)
    
    tau.forest.test <- cf.estimator(X = covariates[-trainId, ], Y = Y[-trainId], W = treatment[-trainId], Z = Z[-trainId], num.trees = num.trees, 
                                      min.node.size = min_node_size, ci.group.size = ci.group.size, seed = seed) 
    tau.pred.train.b <- predict(tau.forest.test, estimate.variance = T, num.threads = n_core)
    tau.pred.test.b <- predict(tau.forest.test, newdata = covariates[trainId, ], estimate.variance = T, num.threads = n_core)
    
    # compute z-score, pvalues, and ajusted.pvalues
    tau.a.train.stats <- compute_stats(tau.pred.train.a)
    tau.a.stats <- compute_stats(tau.pred.test.a)
    tau.b.train.stats <- compute_stats(tau.pred.train.b)
    tau.b.stats <- compute_stats(tau.pred.test.b)
    
    if (save_split) {
      observation_result.a[trainId, ] <- as.matrix(tau.a.train.stats)
      observation_result.a[-trainId, ] <- as.matrix(tau.a.stats)
      write.csv(observation_result.a, file = paste0(file_prefix, "_observation_", i, "_result_a.csv"))
      observation_result.b[-trainId, ] <- as.matrix(tau.b.train.stats)
      observation_result.b[trainId, ] <- as.matrix(tau.b.stats)
      write.csv(observation_result.b, file = paste0(file_prefix, "_observation_", i, "_result_b.csv"))
      
      # Extract varimp from each of the split half forest (Jun 13, 2020 @alex)
      varImp <- variable_importance(tau.forest.train, max.depth = 4)
      train.varimp <- data.frame(variable = colnames(covariates), varImp)
      write.csv(train.varimp, file = paste0(file_prefix, "_observation_", i, "_varimp_train.csv"), row.names = F, quote = F)
      
      varImp <- variable_importance(tau.forest.test, max.depth = 4)
      test.varimp <- data.frame(variable = colnames(covariates), varImp)
      write.csv(test.varimp, file = paste0(file_prefix, "_observation_", i, "_varimp_test.csv"), row.names = F, quote = F)
      #message(paste0("Varimp for observation ", i, " saved."))
    }
    
    # use simes test to combine the p-value, which aims to test whether at least one patient's ITE is significantly different from 0.
    simes_pval.a <- simes.test(tau.a.stats[, 3])
    simes_pval.b <- simes.test(tau.b.stats[, 3])
    
    partial_simes_pval.a <- simes.partial(floor(no.obs.test * partial_simes_thre), tau.a.stats[, 3])
    partial_simes_pval.b <- simes.partial(floor(no.obs.train * partial_simes_thre), tau.b.stats[, 3])
    
    # check the correlation between two predictions from two datasets, we use this to assess the model stability.
    test_res.a <- correlation_test(tau.a.train.stats[, 1], tau.b.stats[, 1], methods = c("pearson", "kendall", "spearman"))
    test_res.b <- correlation_test(tau.a.stats[, 1], tau.b.train.stats[, 1], methods = c("pearson", "kendall", "spearman"))
    
    # fisher.pval.a <- fisher.exact.test(tau.a.train.stats[, 3], tau.b.stats[, 3])
    # fisher.pval.b <- fisher.exact.test(tau.a.stats[, 3], tau.b.train.stats[, 3])
    
    # t.test.pval.a <- quantile.t.test(tau.a.train.stats[, 1], tau.b.stats[, 1])
    # t.test.pval.b <- quantile.t.test(tau.a.stats[, 1], tau.b.train.stats[, 1]) 
    
    # correlation_rslt <- rbind(c(simes_pval.a, partial_simes_pval.a, test_res.a, fisher.pval.a, t.test.pval.a), 
    #                           c(simes_pval.b, partial_simes_pval.b, test_res.b, fisher.pval.b, t.test.pval.b))
    
    correlation_rslt <- rbind(c(simes_pval.a, partial_simes_pval.a, test_res.a), 
                              c(simes_pval.b, partial_simes_pval.b, test_res.b))
    
    correlation_matrix <- rbind(correlation_matrix, correlation_rslt)
  }
  
  if (is_save) {
    colnames(correlation_matrix) <- col_names 
    write.csv(correlation_matrix, file = paste0(file_prefix, '_split_half.csv'), row.names = F, quote = F)
  }

  # print(correlation_matrix)
  
  # change to partial_simes_pval 20190929
  if (length(u) == 1){

    aggregated_rslt <- sapply(seq(dim(correlation_matrix)[2]), function(i) {

      removed_na <- na.omit(correlation_matrix[, i])
      # simes_pval <- ifelse(length(removed_na) > 0, simes.test(removed_na), NA)
      if (!(i %in% c(3, 5, 7))) {
        pval <- ifelse(length(removed_na) > 0, simes.partial(u, removed_na), NA)
      } else {
        pval <- ifelse(length(removed_na) > 0, extract_binom_pval(removed_na), NA) 
      }
      return(pval)
    })

  } else {
    
    col.count <- 2 + length(u) * 6
    aggregated_rslt <- matrix(, nrow = 1, ncol = col.count)

    estimate.col <- c()
    for (i in 1:length(u)){
      estimate.col <- c(estimate.col, 3 + 6 * (i-1))
      estimate.col <- c(estimate.col, 5 + 6 * (i-1))
      estimate.col <- c(estimate.col, 7 + 6 * (i-1))
    }

    for (i in 1:2){
      removed_na <- na.omit(correlation_matrix[, i])
      pval <- ifelse(length(removed_na) > 0, simes.partial(2, removed_na), NA)
      aggregated_rslt[1, i] <- pval
    }

    for (u.select in 1:length(u)){
      
      for (j in 3:8){

        i <- j + 6 * (u.select - 1)
        removed_na <- na.omit(correlation_matrix[, j])
        if (!(i %in% estimate.col)) {
          pval <- ifelse(length(removed_na) > 0, simes.partial(u[u.select], removed_na), NA)
          aggregated_rslt[1, i] <- pval
        } else {
          pval <- ifelse(length(removed_na) > 0, extract_binom_pval(removed_na), NA) 
          aggregated_rslt[1, i] <- pval
        }
      }

    }

  }
  
  return(aggregated_rslt)
  # return(correlation_matrix)
  if (random_rep_seed) seed <- seed + 1 # whether different seeds are used for each replicate
  
}

#' ===================
#' support functions
#' ===================

compute_stats <- function(tau.pred) {
  
  #' This function is to compute zval, pval and ajusted.p using the tau prediction
  #' 
  #' @param tau.pred: the prediction return by grf::predict
  #' 
  #' @return: a dataframe containing four columns, including original tau, tau z value, tau p-value, tau adjusted p-value
  
  tau.zval <- tau.pred$predictions/sqrt(tau.pred$variance.estimates) 
  tau.pval <- pnorm(abs(tau.zval), lower.tail = FALSE) * 2
  tau.p.adjust <- p.adjust(tau.pval, method = "BH")
  stats <- cbind(tau.pred$predictions, tau.zval, tau.pval, tau.p.adjust)
  colnames(stats) <- c("tau", "tau.zval", "tau.pval", "tau.p.adjust")
  return(as.data.frame(stats))

}

correlation_test <- function(x, y, methods, alt = "greater") {
  
  #' This function is to test correlation between two data vector with same length
  #' 
  #' @param x: input data x
  #' @param y: input data y
  #' @param methods: the correlation testing method you want to use (options: "pearson", "kendall", "spearman")
  #' @param alt: alternative, default is "greater"
  #' 
  #' @return: a dataframe containing correlation test results, column is method, and rows are correlation and p-value
  
  fitted_res <- sapply(methods, function(mthd) {
    fitted <- cor.test(x, y, method = mthd, alternative = alt, use = "na.or.complete")
    c(fitted$estimate, fitted$p.value)
  })
  fitted_res
}

ivgrf <- function(X, Y, W, Z, Y.hat = NULL, W.hat = NULL, Z.hat = NULL, 
                  num.trees = 5000, min.node.size = 10, sample.fraction = 0.2, ci.group.size = 5, seed = NULL, num.threads = 10) {
  
  #' This function is for a easy way to apply grf::instrumental_forest inline
  #' should noticed that for permutation test, Y.hat and W.hat should be provided in order to fix Y and W.
  #' 
  #' @param X: covariates matrix X
  #' @param Y: outcome
  #' @param W: treatment assignment
  #' @param Z: instrument (since we are using instrumental forest set in grf, Z should be in 1 dimension)
  #' @param Y.hat: check grf doc
  #' @param W.hat: check grf doc
  #' @param Z.hat: check grf doc
  #' @param num.trees: default 5000, should increase in need.
  #' @param min.node.size: default 10, can increase in need.
  #' @param sample.fraction: default 0.2, can increase in need.
  #' @param ci.group.size: default 5, can increase in need.
  #' @param num.threads: default 10, can increase in need.
  #' 
  #' @return: the instrumental forest
  
  if ((is.null(W.hat) + is.null(W.hat)) == 1) stop("Y.hat and W.hat should be provided at the same time!")
  
  tau.forest <- grf::instrumental_forest(X = X,
                                         Y = Y,
                                         W = W,
                                         Z = Z,
                                         Y.hat = Y.hat,
                                         W.hat = W.hat,
                                         Z.hat = Z.hat,
                                         ci.group.size = ci.group.size,
                                         sample.fraction = sample.fraction,  # When confidence intervals are enabled, the sampling fraction must be less than 0.5.
                                         num.trees = num.trees,
                                         min.node.size = min.node.size,   # previous setting is 10
                                         num.threads = num.threads)
  return(tau.forest)
}

early.stop <- function(pvals) {
  
  #' function to indicate whether we should stop the iteration earilier, based on the given pvals.
  #' 
  #' @param pvals: a pval vector, we will stop eariler if two pvals are all smaller than 0.05.
  #' 
  #' @return: True if we want early stop, False if not.

  is_stop <- (sum(pvals > 0.05) == length(pvals))
  is_stop
}

get.lower.BinomCI <- function(perm.val, observed.val, min.perm, conf.level, BinomCI_method) {
  
  #' function to compute confidence intervals for binomial proportions following the most popular methods.
  #' 
  #' @param perm.val: variance of the tau estimated from permuted model
  #' @param observed.val: variance of the tau estimated from default model
  #' @param min.perm: number of permutation
  #' @param conf.level: confidence interval you wants
  #' @param BinomCI_method: which method to use, default is wilson.
  #' 
  #' @return: lower CI of the binomial dis
  
  CI <- BinomCI(x = sum(perm.val > observed.val),
                n = min.perm,
                conf.level = conf.level,
                method = BinomCI_method)
  CI[2] 
}

extract_binom_pval <- function(est) {
  
  #' function to perform the binomal test for the coefficient estimation.
  #' 
  #' @param est: estimation
  #' 
  #' @return: binom test p val.

  binom_obj <- binom.test(x = sum(est > 0),
                          n = 20,
                          p = 0.5,
                          alternative = "greater",
                          conf.level = 0.95)
  return(binom_obj$p.value)
}

simes.test <- function(x, returnstat = FALSE) {
  
  #' function to perform simes test
  #' 
  #' @param x: p value vector
  #' @param returnstat: whether we return the test statistics.
  #' 
  #' @return: vector containing test statistics and combined p-value / combined p-value only.
  
  r <- rank(x,  ties.method = "random")
  t <- min(length(x) * x / r)
  if (returnstat) c(t, t) else t
}

simes.partial <- function(u, pvec) {  # equation 2 in the above reference
  
  #' function to perform simes partial test
  #' 
  #' @param  u: samples size you believe to be significant
  #' @param pvec: pvalue vector
  #' 
  #' @return: simes partial test results
  
  n <- length(pvec)
  pvec <- pvec[order(pvec)]  # sort pval in ascending order
  
  corr.p <- numeric(n - u + 1)
  
  for (i in 1: (n - u + 1)) {
    corr.p[i] <- (n - u + 1)/i  * pvec[u - 1 + i]
  }
  
  partial.conj.p <- min(corr.p)
  return(partial.conj.p)
}

assess.explained.tau.fixed.YW.risk <- function(fitted.obj, newdata = NULL, Y, Y.hat, W, W.hat, Z, Z.hat, constant.tau = NULL, num.threads = 10) {
  
  #' function to calculate tau-risk R
  #' 
  #' @param fitted.obj: the forest 
  #' @param Y: outcome
  #' @param Y.hat: estimated outcome
  #' @param W: treatment assignment
  #' @param W.hat: estimated W, E[W|X]
  #' @param Z: instrument
  #' @param Z.hat estimated Z, E[Z|X]
  #' @param constant.tau: a constant tau to estimate the tau
  #' 
  #' @return: risk can be explained via the formula
  
  tau.pred <- predict(fitted.obj, newdata = newdata, estimate.variance = TRUE, num.threads = num.threads)
  center.Y <- Y - Y.hat
  center.Z <- Z - Z.hat
  
  if (is.null(constant.tau)) {
    mean.tau <- (W - W.hat) * mean(tau.pred$predictions)
    tau.2 <- (W - W.hat) * tau.pred$predictions
  } else {
    mean.tau <- (W - W.hat) * mean(constant.tau)
    tau.2 <- (W - W.hat) * tau.pred$predictions
  }
  
  risk.total <- (center.Z * (center.Y - mean.tau))^2 
  risk.unexplained <- (center.Z * (center.Y - tau.2))^2

  risk.explained <- mean(risk.total - risk.unexplained)
  
  return(risk.explained)
}

paralleled.perm.cf <- function(X, Y, W,
                               Y.hat, W.hat,
                               Z, Z.hat,
                               num.trees,
                               num.perm,
                               min.node.size,
                               sample.fraction = 0.5,
                               ci.group.size = 5,
                               num.threads = 10, 
                               seed = NULL,
                               cluster_id = NULL,
                               fitted_object = None) {
  
  #' function to perform permutation variance test and permutation risk test
  #' 
  #' @param X: the X covariates
  #' @param Y: outcome
  #' @param W: treatment assignment W
  #' @param Y.hat: estimated Y, E[Y|X]
  #' @param W.hat: estimated W, E[W|X]
  #' @param Z: instrument
  #' @param Z.hat: estimated instrument Z, , E[Z|X]
  #' @param num.trees: tree number
  #' @param num.perm: number of permutations
  #' @param min.node.size: min.node.size, default is 5
  #' @param sample.fraction: sample fraction used to estimate each tree. default: 0.5
  #' @param ci.group.size: ci.group.size setting for causal forest
  #' @param seed: random seed
  #' @param is_save: whether to save the results of each permutation
  #' @param cluster_id: whether we need to do a stratified permutation to keep cluster heterogeneity present.
  #' @param reestimate: whether to re-estimate the Y.hat, Z.hat and W.hat using the permuted covariates.
  #' @param updatedRloss: whether to use updated R-loss based on neuralps paper.
  #' 
  #' @return: results containing permutation variance p-value and permutation risk variance.  
  
  cf.estimator <- ivgrf 
  
  perm.risk.mat <- foreach(i = seq(num.perm), .combine = "rbind", .options.multicore = list(set.seed = T)) %do%  {
    
    # if cluster_id is provided, we do a stratified permutation to keep cluter heterogeneity present
    if (is.null(cluster_id)) {
      samp <- sample(dim(X)[1], dim(X)[1], replace = F)
    } else {
      samp <- stratified.permutation(cluster_id)
    }
    
    sampled.X <- X[samp, ]

    # changed this line 
    perm.tau.forest <- cf.estimator(X = sampled.X,
                                    Y = Y, 
                                    Y.hat = Y.hat,
                                    W = W,
                                    W.hat = W.hat,
                                    Z = Z,
                                    Z.hat = Z.hat,
                                    seed = seed,
                                    ci.group.size = ci.group.size,
                                    sample.fraction = sample.fraction,
                                    num.trees = num.trees,
                                    min.node.size = min.node.size,
                                    num.threads = num.threads)
    perm.var <- var(perm.tau.forest$predictions)
    
    assess.explained.tau.fixed.YW.risk.inline <- assess.explained.tau.fixed.YW.risk
    perm.risk <- assess.explained.tau.fixed.YW.risk.inline(fitted.obj = perm.tau.forest, 
                                                           Y = Y, 
                                                           Y.hat = Y.hat, 
                                                           W = W, 
                                                           W.hat = W.hat, 
                                                           Z = Z,
                                                           Z.hat = Z.hat,
                                                           num.threads = num.threads)

    c(perm.var, perm.risk)
  }
  return(perm.risk.mat)
}

permutate.covariates.testing <- function(X, 
                                         Y,
                                         Y.hat,
                                         W,
                                         W.hat,
                                         Z,
                                         Z.hat,
                                         fixed.YW.tau.risk,
                                         tau.var,
                                         cluster_id = NULL,
                                         num.threads = 10, 
                                         num.trees = 1000,
                                         num.perm = 50,
                                         min.node.size = 20,
                                         ci.group.size = 5,
                                         seed = NULL) {
  
  #' function is to permute X with Y, W, Y.hat, and W.hat, to investigate whether covariates contribute to HTE
  #' Note: if is_save is T, then file_prefix must be provided. 
  #' file_prefix is a prefix for saving, since several results will be saved. 
  #' 
  #' @param X: the X covariates
  #' @param Y: outcome
  #' @param Y.hat: estimated Y, E[Y|X]
  #' @param W: treatment assignment
  #' @param W.hat: estimated W, E[W|X]
  #' @param Z: instrument
  #' @param Z.hat: estimated Z, E[Z|X]
  #' @param fixed.YW.tau.risk: tau risk calculated in the default model
  #' @param tau.var: tau variance calculated in the default model
  #' @param num_trees: tree number 
  #' @param ci.group.size: ci.group.size setting for causal forest
  #' @param seed: random seed
  #' @param cluster_id: whether we need to do a stratified permutation to keep cluster heterogeneity present.
  #' 
  #' @return: tau.risk for each permutation
  
  perm.risk.mat <- paralleled.perm.cf(X = X,
                                      Y = Y, 
                                      W = W,
                                      Z = Z,
                                      Y.hat = Y.hat, 
                                      W.hat = W.hat,
                                      Z.hat = Z.hat,
                                      num.trees = num.trees,
                                      num.perm = num.perm,
                                      min.node.size = min.node.size,
                                      ci.group.size = ci.group.size,
                                      seed = seed,
                                      cluster_id = NULL,
                                      num.threads = num.threads)
  
  # calculate the P-value, based on formula: Pr(Null) = #(Var_perm(tau) >= Var_observed(tau))/N
  permute.var.pval <- mean(perm.risk.mat[, 1] > tau.var) 
  fixed.YW.permutation.pval <- mean(perm.risk.mat[, 2] > fixed.YW.tau.risk)
  
  return(c(permute.var.pval, fixed.YW.permutation.pval))
}

adaptive.permutate.covariates.testing <- function(X, 
                                                  Y, Y.hat,
                                                  W, W.hat, 
                                                  Z, Z.hat,
                                                  fixed.YW.tau.risk,
                                                  tau.var,
                                                  cluster_id = NULL,
                                                  num.trees = 1000,
                                                  num.strap = 500,
                                                  min.perm = 10,
                                                  conf.level = 0.99,
                                                  min.node.size = 20,
                                                  ci.group.size = 5,
                                                  sample.fraction = 0.5,
                                                  BinomCI_method = "wilson",
                                                  num.threads = num.threads,
                                                  seed = NULL) {
  
  #' function to permute X to investigate whether X contribute to HTE
  #' 
  #' @param X: covariates
  #' @param Y: outcome
  #' @param Y.hat: estimated Y, E[Y|X]
  #' @param W: treatment assignment
  #' @param W.hat: estimated W, E[W|X]
  #' @param W.hat.XZ: estimated W, E[W|X,Z] default NA
  #' @param Z: instrument
  #' @param Z.hat: estimated Z, E[Z|X]
  #' @param fixed.YW.tau.risk: tau risk calculated in the default model
  #' @param tau.var: tau variance calculated in the default model
  #' @param cluster_id：whether we need to do a stratified permutation to keep cluster heterogeneity present, if so, provide the cluster id here.
  #' @param num.trees: tree number 
  #' @param num.strap: how many times of permutation we want to perform
  #' @param min.perm: minimum number of permutation we need to perform
  #' @param conf.level: confidence interval we need
  #' @param BinomCI_method: method to calculate the binomial CI, default is wilson
  #' @param is_save: whether to save the results of each permutation
  #' @param file_prefix: the file prefix to store the results of each permutation
  #' @param seed: random seed
  #' @param min.node.size: min node size in ivgrf
  #' @param ci.group.size: ci.group.size setting for causal forest
  #' @param sample.fraction: sample fraction used to estimate each tree. default: 0.5
  #' @param reestimate: whether to re-estimate the Y.hat, Z.hat and W.hat using the permuted covariates.
  #' @param updatedRloss: whether to use updated R-loss based on neuralps paper.
  #' 
  #' @return: tau.risk for each permutation

  old.perm.risk.mat <- NULL; stopping <- FALSE; perm_num <- min.perm; idx <- 0
  
  while (!stopping) {
    
    if (perm_num >= num.strap) break
    
    perm.risk.mat <- paralleled.perm.cf(X = X, 
                                        Y = Y,
                                        W = W,
                                        Z = Z,
                                        Y.hat = Y.hat,
                                        W.hat = W.hat,
                                        Z.hat = Z.hat,
                                        num.trees = num.trees,
                                        num.perm = min.perm,
                                        min.node.size = min.node.size,
                                        ci.group.size = ci.group.size,
                                        sample.fraction = sample.fraction,
                                        num.threads = num.threads,
                                        seed = NULL,
                                        cluster_id = NULL)
    # print(perm.risk.mat)
    if (is.null(old.perm.risk.mat)) {
      old.perm.risk.mat <- perm.risk.mat
    } else {
      old.perm.risk.mat <- rbind(old.perm.risk.mat, perm.risk.mat)
    }
    
    perm_num <- nrow(old.perm.risk.mat)
    
    lower.var.pval <- get.lower.BinomCI(old.perm.risk.mat[, 1], tau.var, perm_num, conf.level, BinomCI_method)
    lower.fixed.YW.pval <- get.lower.BinomCI(old.perm.risk.mat[, 2], fixed.YW.tau.risk, perm_num, conf.level, BinomCI_method)
    
    idx <- idx + 1
    stopping <- ifelse(early.stop(c(lower.var.pval, lower.fixed.YW.pval)), TRUE, FALSE)
  } # end of while
  
  permute.var.pval <- mean(old.perm.risk.mat[, 1] > tau.var) 
  fixed.YW.permute.pval <- mean(old.perm.risk.mat[, 2] > fixed.YW.tau.risk)
  
  return(c(permute.var.pval, fixed.YW.permute.pval))
}

test_heterogeneity <- function(X, Y, W, Z, 
                               Y.hat = NULL, W.hat = NULL, Z.hat = NULL, 
                               num.trees = 10000, seed = 309,
                               sample.fraction = 0.15, min.node.size = 100, ci.group.size = 10, min.perm = 30, num.threads = 10,
                               num.strap = 100, constant.tau = NULL){

  #' function to permute X to investigate whether X contribute to HTE
  #' 
  #' @param X: covariates
  #' @param Y: outcome
  #' @param W: treatment assignment
  #' @param Z: instrument
  #' @param num.trees: tree number 
  #' @param num.strap: how many times of permutation we want to perform
  #' @param min.perm: minimum number of permutation we need to perform
  #' @param min.node.size: min node size in ivgrf
  #' @param ci.group.size: ci.group.size setting for causal forest
  #' @param sample.fraction: sample fraction used to estimate each tree. default: 0.5
  #' 
  #' @return: p-value for permutation-variance test and permutation-risk test (perm-var, perm-risk)

  iv.forest.continuous <- instrumental_forest(X, Y, W, Z, 
                                              Y.hat = Y.hat,
                                              W.hat = W.hat,
                                              Z.hat = Z.hat,
                                              num.threads = num.threads, num.trees = num.trees, 
                                              sample.fraction = sample.fraction, min.node.size = min.node.size, ci.group.size = ci.group.size, seed = seed)
  iv.pred.continuous <- predict(iv.forest.continuous, estimate.variance = TRUE)
  iv.pred.stats.continuous <- compute_stats(iv.pred.continuous)
  summary(iv.pred.stats.continuous)
  
  fixed.YW.tau.risk.continuous <- assess.explained.tau.fixed.YW.risk(fitted.obj = iv.forest.continuous, 
                                                                     Y = Y, 
                                                                     Y.hat = iv.forest.continuous[["Y.hat"]], 
                                                                     W = W,  
                                                                     W.hat = iv.forest.continuous[["W.hat"]],
                                                                     Z = Z, 
                                                                     Z.hat = iv.forest.continuous[["Z.hat"]],
                                                                     constant.tau = constant.tau,
                                                                     num.threads = num.threads)
  
  print(fixed.YW.tau.risk.continuous)
  perm.pvals.continuous <- adaptive.permutate.covariates.testing(X = X,
                                                                 Y = Y,  
                                                                 Y.hat = iv.forest.continuous[["Y.hat"]],
                                                                 W = W, 
                                                                 W.hat = iv.forest.continuous[["W.hat"]],
                                                                 Z = Z, 
                                                                 Z.hat = iv.forest.continuous[["Z.hat"]],
                                                                 fixed.YW.tau.risk = fixed.YW.tau.risk.continuous,
                                                                 tau.var = var(iv.pred.continuous$predictions),
                                                                 min.node.size = min.node.size,
                                                                 sample.fraction = sample.fraction,
                                                                 ci.group.size = ci.group.size,
                                                                 num.trees = num.trees, 
                                                                 num.strap = num.strap,
                                                                 min.perm = min.perm,
                                                                 num.threads = num.threads)
  
  perm.pvals.continuous
}

#' ##################################
#' Test functions
#' ##################################

#' In this part, we want to write a new function which used the update R-loss to test for the heterogeneity.
#' Mainly based on the paper "Machine Learning Estimation of Heterogeneous Treatment Effects with Instruments"
#' 
#' We use the updated R-loss DMLIV. The main difference between this one and the one used by GRF is that the use of instrument.
#' 
#' DMLIV Partially orthogonal, convex loss:
#' L^1(\theat; q_0, h_0, p_0) := E[(Y - E[Y|X] - \theta(X)(E[T|Z,X] - E[T|X]))^2]
#' 
#' GRF:
#' (Z_i - \bar{Z}_P)((Y_i - \bar{Y}_P) - (W_i - \bar{W}_P) \hat{\tau}_P)
#' 
#' We can observe that in DMLIV, the instrument is used to estimate E[T|Z,X], while in GRF is used at the beginning.
#' An obviously issue is that how can we get the E[T|Z,X]?
#' Solution 1: We just use another algorithm to estimate it.

#' ##################################
#' functions currently not being used
#' comment will be ignored
#' ##################################

fisher.exact.test <- function(x, y, sig = 0.1) {
  crosstabular <- table(as.factor(x > sig), as.factor(y > sig))
  unname(crosstabular)
  if (dim(crosstabular)[1] == 2 && dim(crosstabular)[2] == 2) {
    fitted <- fisher.test(crosstabular, alternative = "greater")
    return(fitted$p.value)
  } else {
    return(NA)
  }
}

median.t.test <- function(x, y) {
  x.median <- median(x)
  y.median <- median(y)
  
  if (sum(x <= y.median) < 10 | sum(x > y.median) < 10) {
    x.p.value <- NA
  }else{
    x.fit <- t.test(x[x <= y.median], x[x > y.median], alternative = "less")
    x.p.value <- x.fit$p.value
  }
  
  if (sum(y <= x.median) < 10 | sum(y > x.median) < 10) {
    y.p.value = NA
  } else {
    y.fit <- t.test(y[y <= x.median], y[y > x.median], alternative = "less")
    y.p.value <- y.fit$p.value
  }
  return(c(x.p.value, y.p.value))
}

quantile.t.test <- function(x, y, quantile = 0.5) {
  x.quantile <- quantile(x, quantile)
  y.quantile <- quantile(y, quantile)
  
  if (sum(x <= y.quantile) < 5 | sum(x > y.quantile) < 5) {
    x.p.value <- NA
  } else {
    x.fit <- t.test(x[x <= y.quantile], x[x > y.quantile], alternative = "less", na.rm = T)
    x.p.value <- x.fit$p.value
  }
  
  if (sum(y <= x.quantile) < 5 | sum(y > x.quantile) < 5) {
    y.p.value <- NA
  } else {
    y.fit <- t.test(y[y <= x.quantile], y[y > x.quantile], alternative = "less", na.rm = T)
    y.p.value <- y.fit$p.value
  }
  return(c(x.p.value, y.p.value))
}

permutated.pval <- function(Y1, Y2, no.simulations = 10000) {
  # test fitting with permutating one outcome variable, and a pval will be returned
  # set.seed(0)
  risk <- mean((Y1 - Y2)^2)
  
  permutated.risk <- sapply(seq(no.simulations), function(x) {
    mean((Y1[sample(length(Y1))] - Y2)^2)
  })
  return(mean(permutated.risk < risk))
}