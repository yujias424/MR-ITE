#' This code is to generate the simulation data in instrumental setting.
#' The main idea comes from the power's paper "Some methods for heterogeneous treatment effect
#' estimation in high dimensions"
#' Also, Kai's thesis will be used as reference
#' 
#' @author: Yujia Shi
#' @date: 2022.01.07

#' ============================================================
#' 8 functions to generate the mean effect and treatment effect
#' ============================================================

function.1 <- function(x){
    res <- rep(0, dim(x)[1])
    res
}

function.2 <- function(x){
    res <- 5 * ifelse(x[, 5]>1, 1, 0) - 5 # 5 * ifelse(x[, 5]>1, 1, 0) - 5
    res
}

function.3 <- function(x){
    res <- 2 * x[, 1] - 4
    res
}

function.4 <- function(x){
    
    res <- x[, 2] * x[, 4] * x[, 6] + 2 * x[, 2] * x[, 4] * (1 - x[, 6]) + 
            3 * x[, 2] * (1 - x[, 4]) * x[, 6] +
            4 * x[, 2] * (1 - x[, 4]) * (1 - x[, 6]) + 
            5 * (1 - x[, 2]) * x[, 4] * x[, 6] + 
            6 * (1 - x[, 2]) * x[, 4] * (1 - x[, 6]) +
            7 * (1 - x[, 2]) * (1 - x[, 4]) * x[, 6] + 
            8 * (1 - x[, 2]) * (1 - x[, 4]) * (1 - x[, 6])
    res
}

function.5 <- function(x){
    res <- x[, 1] + x[, 3] + x[, 5] + x[, 7] + x[, 8] + x[, 9] - 2
    res
}

function.6 <- function(x){
    res <- 4 * ifelse(x[, 1] > 1, 1, 0) * ifelse(x[, 3] > 0, 1, 0) +
            4 * ifelse(x[, 5] > 1, 1, 0) * ifelse(x[, 7] > 0, 1, 0) +
            2 * x[, 8] * x[, 9]
    res
}

function.7 <- function(x){
    res <- x[, 1] ** 2 + x[, 2] + x[, 3] ** 2 + 
            x[, 4] + x[, 5] ** 2 + x[, 6] + x[, 7] ** 2 + x[, 8] + x[, 9] ** 2 - 11
    res <- 0.5 * res
    res 
}

function.8 <- function(x){
    res <- 1/sqrt(2) * (function.4(x) + function.5(x))
}

#' ====================
#' scenarios function
#' ====================

scenarios <- function(n, p, mu, tau, sigma, tprop = NULL, corr_Z_T = 0.6){
    
    #' this function is to create different scenarios for testing the
    #' performance of our test in instrumental way.
    #' 
    #' @param n number of observations
    #' @param p number of features
    #' @param mu mean effect function
    #' @param tau treatmnet effect function
    #' @param sigma outcome standard deviation (as outcome ~ normal)
    #' @param corr_Z_T correlation between instrument and treatment assignment
    #' 
    #' @return a list containing X, Z, Y, T and true tau.

    # simulate the x
    dat.x <- as.data.frame(matrix(nrow = n, ncol = p))
    
    odd.cols <- colnames(dat.x[, c(TRUE, FALSE)])
    even.cols <- colnames(dat.x[, c(FALSE, TRUE)])

    # simulate the odd columns
    for (i in odd.cols){
        dat.x[, i] <- rnorm(n = n)
    }

    # simulate the even columns
    for (i in even.cols){
        dat.x[, i] <- rbinom(n = n, size = 1, prob = 0.5)
    }

    # simulate the instrument
    instrument.z <- rnorm(n = n)

    # simulate the noise
    xi <- rnorm(n) # xi
    zeta <- rnorm(n) # zeta

    # simulate the treatment with instrument (this is already an endogenity case, as the treatment assignment is not randomized)
    latent.index <- corr_Z_T * instrument.z + 0.5 * zeta + sqrt(2/pi - 0.5**2 - corr_Z_T**2) * rnorm(n)
    if (is.null(tprop)){
        latent.index <- latent.index
    } else {
        sclaesize <- quantile(latent.index, tprop)
        latent.index <- latent.index + (-sclaesize)
    }

    treatment.t <-  0.5 * (sign(latent.index) + 1)

    # simulate the outcome
    mean.y <- mu(dat.x) + (treatment.t - 1/2) * tau(dat.x)
    outcome.y  <- rnorm(n, mean = mean.y, sigma) + xi + zeta

    # list("X" = dat.x, "Z" = instrument.z, "Y" = outcome.y, "T" = treatment.t - 0.5, "tau" = (treatment.t - 1/2) * tau(dat.x))
    list("X" = dat.x, "Z" = instrument.z, "Y" = outcome.y, "T" = treatment.t, "tau" = tau(dat.x), "latentT" = latent.index)
}

scenarios.burgess <- function(n, p, mu, tau, snps_n = 100, invalid_snps_n = 20, pleiotropy_scen = 1){
  
  #' this function is to create different scenarios for testing the
  #' performance of our test in instrumental way.
  #' 
  #' @param n number of observations
  #' @param p number of features
  #' @param mu mean effect function
  #' @param tau treatmnet effect function
  #' @param corr_Z_T correlation between instrument and treatment assignment
  #' 
  #' @return a list containing X, Z, Y, T and true tau.
  
  # simulate the x
  dat.x <- as.data.frame(matrix(nrow = n, ncol = p))
  
  odd.cols <- colnames(dat.x[, c(TRUE, FALSE)])
  even.cols <- colnames(dat.x[, c(FALSE, TRUE)])
  
  # simulate the odd columns
  for (i in odd.cols){
    dat.x[, i] <- rnorm(n = n)
  }
  
  # simulate the even columns
  for (i in even.cols){
    dat.x[, i] <- rbinom(n = n, size = 1, prob = 0.5)
    # dat.x[, i] <- rnorm(n = n)
  }
  
  # simulate the genetic variant
  gv.df <- matrix(nrow = n, ncol = snps_n)
  for (i in 1:snps_n){
    gv.df[, i] <- rbinom(n, 2, 0.3)
  }
  
  # simulate the noise
  epsilon_U <- rnorm(n)
  epsilon_T <- rnorm(n)
  epsilon_Y <- rnorm(n)
  
  # simulate the gamma
  # uniform distribution
  gamma_vec <- runif(snps_n, min = 0.03, max = 0.1)
  intercept <- (0.1+0.03)/2 * 0.6 * snps_n
  
  # normal distribution
  # gamma_vec <- rnorm(snps_n, 0.065, 0.02)
  # intercept <- 0.065 * 0.6 * snps_n
  
  # simulate the xi and alpha
  if (pleiotropy_scen == 1){
    xi_vec <- rep(0, snps_n)
    a_vec <- rep(0, snps_n)
  } else if (pleiotropy_scen == 2) {
    xi_vec <- rep(0, snps_n)
    a_vec <- c(runif(invalid_snps_n, min = -0.1, max = 0.1), rep(0, snps_n-invalid_snps_n))
  } else if (pleiotropy_scen == 3) {
    xi_vec <- rep(0, snps_n)
    a_vec <- c(runif(invalid_snps_n, min = 0, max = 0.1), rep(0, snps_n-invalid_snps_n))
  } else if (pleiotropy_scen == 4) {
    xi_vec <- c(runif(invalid_snps_n, min = -0.1, max = 0.1), rep(0, snps_n-invalid_snps_n))
    a_vec <- rep(0, snps_n)
  } else {
    message("Invalid pleiotropy scene.")
  }
  
  # simulate the U
  U.xi <- rowSums(t(t(gv.df) * xi_vec)) + epsilon_U
  Y.alpha <- rowSums(t(t(gv.df) * a_vec)) + epsilon_Y
  
  # 50% quantile
  # latent.index <- rowSums(t(t(gv.df) * gamma_vec)) + epsilon_T + U.xi
  # fifty.quantile <- quantile(latent.index, 0.5)
  # treatment.t <-  ifelse(latent.index > fifty.quantile, 1, 0)
  
  # zero mean
  # latent.index <- rowSums(t(t(gv.df) * gamma_vec)) + epsilon_T + U.xi
  # treatment.t <-  0.5 * (sign(latent.index - intercept) + 1)

  # continuous treatment t
  latent.index <- rowSums(t(t(gv.df) * gamma_vec)) + epsilon_T + U.xi
  treatment.t <- latent.index
  
  # simulate the outcome
  outcome.y <- treatment.t * (tau(dat.x)/10) + U.xi + Y.alpha # + mu(dat.x)
  # outcome.y <- treatment.t * tau(dat.x) + U.xi + Y.alpha
  list("X" = dat.x, "Z" = as.data.frame(gv.df), "Y" = outcome.y, "T" = treatment.t, 
       "tau" = (tau(dat.x)/10), "prs" = rowSums(t(t(gv.df) * gamma_vec)), "lt" = latent.index)

  # outcome.y <- treatment.t * 0.1 + U.xi + Y.alpha
  # list("X" = dat.x, "Z" = as.data.frame(gv.df), "Y" = outcome.y, "T" = treatment.t, "tau" = rep(0.1, n), "prs" = rowSums(t(t(gv.df) * gamma_vec)))
}
