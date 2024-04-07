#' This code is to perform the split-half permutation validation test to check the presence of heterogeneity
#' 
#' @author Shi Yujia
#' @date 2022.09.01

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

files.list <- list.files("/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/06_validation_test_res/")
for (i in files.list) {
  
  file.list.name <- stringr::str_split(i, "_")[[1]]
  trait <- file.list.name[1]
  disease <- file.list.name[2]
  perm_type <- stringr::str_extract(file.list.name[3], ".*(?=\\.)")
  
  res <- read.csv(paste0("/home/yujia/Project/2022-09-01-individual_MR/res/02_ITE_analysis/06_validation_test_res/", i), row.names = 1)
  continuous_res <- simes.test(c(res$continuous_train, res$continuous_test))
  binary_res <- simes.test(c(res$binary_train, res$binary_test))
  
  message(paste0("Now we are analyzing results from ", trait, " on ", disease, " with perm type ", perm_type, "."))
  
  print(continuous_res)
  print(binary_res)
  
}
