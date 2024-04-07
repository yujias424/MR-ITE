#' This code is to perform propensity score matching in selecting PSM control for further analysis.
#' 
#' @author: Yujia Shi
#' @date: 2021.12.10

suppressPackageStartupMessages({
  library(MatchIt)
  library(tidyverse)
  library(data.table)
})
  
psm <- function(X, Y, method = "nearest", distance = "glm", link = "logit", ratio = 1) {
  
  #' This function is to perform PSM to get a better 1:1 ratio in final sample for further analysis.
  #' @param X
  #' @param Y: the id is already homonized, the data is converted, should noticed that id needs to be the same in all files
  #' 
  #' @return: matched id.
  
  Y <- as.data.frame(Y)
  X$treatment_w <- Y[, 2]
  X_nomiss <- X %>%  # MatchIt does not allow missing values
    na.omit()
  
  X_nomiss <- as.data.frame(X_nomiss)
  
  mod_match <- MatchIt::matchit(treatment_w ~ . - IID - FID,
                                method = method, distance = distance, data = X_nomiss, link = link, ratio = ratio) # IID and FID are first two column of X data.frame
  dta_m <- match.data(mod_match)
   
  dta_m$IID
   
}
