#' This code is to perform extra validation on the heterogeneity testing results.
#' The idea comes from the discussion between Alex and Alejandro
#' details can be found in related email.
#' In breif, we have two approaches:
#' 1. Retrain a CF model and use the optimal tree as the partitioning model.
#' 2. Retrain a regular regression tree model where the target is the taus (i.e., Regression Forest and use the optimal tree as the partition model.)
#' 
#' Here we use option 2.
#' 
#' binary LDL -> CAD 
#'
#' @author Yujia Shi
#' @date 2022-12-10

suppressPackageStartupMessages({
  # library and functions for running analysis.
  library(grf)
  library(tidyverse)
  library(data.table)
  library(caret)
  library(ggplot2)
  library(xgboost)
  library(foreach)
  library(doParallel)
  library(DescTools)
  library(patchwork)
  library(latex2exp)
  library(AER)
  library(foreign)
  library(xlsx)
  library(ivreg)
  library(ggsignif) 
  library(evtree)
  
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/psm.R")
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/test_heterogeneity.R")

  # library and functions for testing.
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/03_heterogeneity_testing/subgroup_analysis/find_best_tree.R")
})

# LDL - CAD Binary trait

# ==================
# Data Preprocessing
# ==================
# X, Y, T, Z from their original source file
X <- fread("~/Project/2022-09-01-individual_MR/dat/04_pheno_covar_data/LDL/ukbb.covariate.LDL.complete", sep = "\t")
Z <- fread("~/Project/2022-09-01-individual_MR/dat/06_PRS_calculation/LDL/CAD/1e-08/LDL_prs.best", sep = " ") # score_constd
W <- fread("~/Project/2022-09-01-individual_MR/dat/04_pheno_covar_data/LDL/ukbb.phenotype.LDL.complete.mgdL", sep = "\t")
Y <- fread("~/Project/2022-09-01-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid", sep = "\t")

# homonize the ID included in the analysis
Z_id <- unique(Z$IID)
X_id <- unique(X$IID)
W_id <- unique(W$IID)
Y_id <- unique(Y$IID) 

selected_id <- Reduce(intersect, list(X_id, Y_id, W_id, Z_id))

X <- X[X$IID %in% selected_id, ]
Y <- Y[Y$IID %in% selected_id, ]
W <- W[W$IID %in% selected_id, ]
Z <- Z[Z$IID %in% selected_id, ]

# select covariates
selected.covariates <- c("FID", "IID", 
                         "22001-0.0", "4079-0.0", "4080-0.0", "189-0.0", "21022-0.0",                 
                         "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5",                 
                         "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10",
                         "whr", # "23099-0.0", "21001-0.0", 
                         "Blood_pressure_medication", "No_medication", "Insulin", # "Cholesterol_lowering_medication", 
                         "Non_alcohol_drinker", "Previous_alcohol_drinker", "Current_alcohol_drinker",
                         "Non_smoker", "Previous_smoker", "Current_smoker",
                         "30630-0.0", "30760-0.0", "30870-0.0", # lipid-related covariates
                         "30700-0.0", "30710-0.0", "30720-0.0", "30730-0.0", "30740-0.0", "30750-0.0", 
                         "30650-0.0", "30660-0.0", "30670-0.0", "30770-0.0",
                         "30810-0.0", "30830-0.0", "30850-0.0", "30860-0.0", "30880-0.0", "30890-0.0") 
X <- dplyr::select(X, all_of(selected.covariates))
X$`22001-0.0` <- as.integer(X$`22001-0.0`)
X$`21022-0.0` <- as.numeric(X$`21022-0.0`)
X$`4079-0.0` <- as.numeric(X$`4079-0.0`)
X$`4080-0.0` <- as.numeric(X$`4080-0.0`)

X <- cbind(X, Y[, c("copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "stroke", "ischemic_stroke")])
Y <- Y[, c("IID", "CAD")]

samples <- X$IID

W.mat <- W[W$IID %in% samples, c("30780-0.0")]
X.mat <- X[X$IID %in% samples, 3:dim(X)[2]]
Y.mat <- Y[Y$IID %in% samples, 2]
Z.mat <- Z[Z$IID %in% samples, 4]

X.mat <- X.mat %>% 
           rename("Gender"="22001-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0", 
                  "Age"="21022-0.0",             
                  "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                  "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10",
                  "Waist-hip-ratio"="whr", # "Body Fat Percentage"="23099-0.0", "BMI"="21001-0.0", 
                  "Blood pressure medication"="Blood_pressure_medication", "No medication"="No_medication", "Insulin"="Insulin", # "Cholesterol lowering medication"="Cholesterol_lowering_medication", 
                  "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                  "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                  "Apolipoprotein A"="30630-0.0", "HDL-C"="30760-0.0", "Triglycerides"="30870-0.0", # lipid-related covariates
                  "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", 
                  "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                  "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                  "Total protein"="30860-0.0", "Urate"="30880-0.0", "Vitamin D"="30890-0.0",
                  "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", 
                  "Venous thromboembolism history"="vte", "Chronic obstructive pulmonary disease history"="copd", "Atrial fibrillation history"="af", 
                  "Hemorrhage Stroke history"="hemorrhage_stroke", "Other stroke history"="stroke", "Ischemic Stroke history"="ischemic_stroke")

W.vector <- as.vector(W.mat$`30780-0.0`)
Y.vector <- as.vector(Y.mat$CAD)
Y.vector <- as.integer(Y.vector)
Z.vector <- as.vector(Z.mat$PRS)

W.vector.binary.cutoff <- as.integer(ifelse(W.vector <= 130, 1, 0)) # mimic usage of statin
W.vector.binary.cutoff <- as.integer(W.vector.binary.cutoff)

sample.id <- Y$IID

#' =======================
#' Start Analysis
#' =======================
# Run the validation test
set.seed(1)

# train a ivf model.
iv.forest.continuousZ.rawmodel <- instrumental_forest(X, Y.vector, W.vector.binary.cutoff, Z.vector, 
                                                      num.threads = 10, num.trees = 5000, 
                                                      sample.fraction = 0.15, min.node.size = 300, ci.group.size = 100) 
iv.pred.raw <- predict(iv.forest.continuousZ.rawmodel, estimate.variance = TRUE)
iv.pred.raw$IID <- sample.id

# split the data into training and testing set.
trainIndex <- caret::createDataPartition(iv.pred.raw$predictions, list = FALSE, p = .3)

iv.pred.train <- iv.pred.raw[trainIndex, c("IID", "predictions")]
iv.pred.test <- iv.pred.raw[-trainIndex, c("IID", "predictions")]

# train a evtree
X.train <- cbind(X[trainIndex, ], iv.pred.train)
X.train$IID <- NULL
colnames(X.train)[length(X.train)] <- "taus"

ev.tree <- evtree::evtree(taus ~ ., data = X.train, maxdepth=3, minbucket=10)
plot(ev.tree)

# test results on testing set
X.test <- cbind(X[-trainIndex, ], iv.pred.test)
X.test$IID <- NULL
colnames(X.test)[length(X.test)] <- "taus"

X.test <- X[-trainIndex, ]

evres.test <- predict(ev.tree, newdata = X.test)
evres.test <- data.frame("IID" = iv.pred.test$IID,
                         "subgroup" = predict(ev.tree, newdata = X.test))
subgroup.id <- sort(unique(evres.test$subgroup), decreasing = T)
evres.test$subgroup <- plyr::mapvalues(evres.test$subgroup, from = subgroup.id, to = seq(1,length(subgroup.id)))
evres.test$index <- as.numeric(row.names(evres.test))

Y.test <- Y[evres.test$index, c("heart_disease")]
W.test <- W.vector.binary.cutoff[evres.test$index]
X.test <- X[evres.test$index, ]
Z.test <- Z[evres.test$index]
subgroup.index.test <- evres.test$subgroup

dat.test <- cbind(Y.test, W.test, Z.test, X.test, subgroup.index.test)
colnames(dat.test)[1] <- "Y"
colnames(dat.test)[2] <- "W"
colnames(dat.test)[3] <- "Z"
colnames(dat.test)[dim(dat.test)[2]] <- "subgroup"

# ==============================

# Alejandro approach
results.list <- list()
for (i in 1:length(table(dat.test$subgroup))){
  for (j in 1:length(table(dat.test$subgroup))){
    
    if (i > j){
      group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
      group.index.b <- as.numeric(names(table(dat.test$subgroup))[j])
      
      ivreg.res.test.a <- ivreg(formula = Y ~ W + . - Z - subgroup  | Z + . - subgroup, 
                                data = dat.test[dat.test$subgroup == group.index.a,]) # https://www.econometrics-with-r.org/12-2-TGIVRM.html
      ivreg.test.summary.a <- summary(ivreg.res.test.a)
      ivreg.res.test.b <- ivreg(formula = Y ~ W + . - Z - subgroup  | Z + . - subgroup, 
                                data = dat.test[dat.test$subgroup == group.index.b,])
      ivreg.test.summary.b <- summary(ivreg.res.test.b)
      
      mean.a <- ivreg.test.summary.a$coefficients[2,1]
      std.error.a <- ivreg.test.summary.a$coefficients[2,2]  
      mean.b <- ivreg.test.summary.b$coefficients[2,1]
      std.error.b <- ivreg.test.summary.b$coefficients[2,2] 
      
      T_value.ab <- (mean.a-mean.b)/sqrt(std.error.a**2 + std.error.b**2)
      p_value.ab <- pt(T_value.ab, (table(dat.test$subgroup)[[i]]+table(dat.test$subgroup)[[i]]-2), lower.tail = F)
      
      T_value.ba <- (mean.b-mean.a)/sqrt(std.error.a**2 + std.error.b**2)
      p_value.ba <- pt(T_value.ba, (table(dat.test$subgroup)[[j]]+table(dat.test$subgroup)[[j]]-2), lower.tail = F)
      
      message("         ")
      print(paste0(group.index.a, " : ", group.index.b))
      print(paste0(mean.a, " : ", mean.b))
      print(p_value.ab)
      print(p_value.ba)
      
      if (min(p_value.ab, p_value.ba) < 0.05){
        # print(max(mean.a, mean.b) + 1.96*ifelse(mean.a > mean.b, std.error.a, std.error.b))
        results.list[[paste0(group.index.a, " : ", group.index.b)]] <- c(min(p_value.ab, p_value.ba), gtools::stars.pval(min(p_value.ab, p_value.ba)), 
                                                                         names(table(dat.test$subgroup))[i], names(table(dat.test$subgroup))[j], 
                                                                         max(mean.a, mean.b) + 1.96*ifelse(mean.a > mean.b, std.error.a, std.error.b))
      }
    }
    
  }
}

results.list <- results.list[order(sapply(results.list, function(x) x[5], simplify=TRUE), decreasing=TRUE)]

mean.cf <- c()
std.error.cf <- c()
class.cf <- c()
sample.size.cf <- c()

for (i in 1:length(table(dat.test$subgroup))){
  
  group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
  ivreg.res.test.a <- ivreg(formula = Y ~ W + . - Z - subgroup  | Z + . - subgroup, 
                            data = dat.test[dat.test$subgroup == group.index.a,])
  ivreg.test.summary.a <- summary(ivreg.res.test.a)
  
  mean.cf <- c(mean.cf, ivreg.test.summary.a$coefficients[2,1])
  std.error.cf <- c(std.error.cf, ivreg.test.summary.a$coefficients[2,2])
  sample.size.cf <- c(sample.size.cf, dim(dat.test[dat.test$subgroup == group.index.a,])[1])
  class.cf <- c(class.cf, names(table(dat.test$subgroup))[i])
}

sample.df <- data.frame(mean = mean.cf, std.error = std.error.cf, sample.size = sample.size.cf, class = class.cf)

comparsion_list <- list()
annotations_vector <- c()
y_position_vector <- c()
for (i in 1:length(results.list)){
  res.tmp <- results.list[[i]]
  annotations_vector <- c(annotations_vector, res.tmp[2])
  comparsion_list[[i]] <- c(res.tmp[3], res.tmp[4])
  y_position_vector <- c(y_position_vector, as.numeric(res.tmp[5]))
}

y_position_vector_tmp <- y_position_vector
for (i in 1:length(y_position_vector)){
  if (i == 1){
    y_position_vector_tmp[which(order(y_position_vector) == i)] <- max(y_position_vector_tmp) + 0.0002
  } else {
    y_position_vector_tmp[which(order(y_position_vector) == i)] <- max(y_position_vector_tmp) + 0.0004
  }
}

subgroup_matching_list <- list()
group_index <- 1
for (i in as.character(sort(as.numeric(names(table(sample.df$class)))))){
  subgroup_matching_list[[i]] <- letters[group_index]
  group_index <- group_index + 1
}

p.validation <- ggplot(sample.df, aes(factor(class), y=mean, ymin=mean - 1.96 * std.error, ymax = mean + 1.96 * std.error))+
  geom_pointrange() + geom_signif(comparisons = comparsion_list, annotations = annotations_vector, y_position = y_position_vector_tmp) +
  theme_classic() + xlab("Subgroup") + ylab("LATE") + ggtitle("Coronary Artery Disease") + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 20),
        axis.text = element_text(size = 25), axis.line=element_line(color='black')) + scale_x_discrete(labels = subgroup_matching_list)
p.validation

# ===========================================================================================================================
library("grid")
library("ggplotify")

plot(ev.tree)
ev.tree.plot.gg <- as.ggplot(plot(ev.tree))
ev.tree.plot + p.validation
p.validation

