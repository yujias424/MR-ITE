#' This code is to perform extra validation on the heterogeneity testing results.
#' The idea comes from the discussion between Alex and Alejandro,
#' details can be found in related email.
#' In breif, we have two approaches:
#' 1. Retrain a CF model and use the optimal tree as the partitioning model.
#' 2. Retrain a regular regression tree model where the target is the taus (i.e., Regression Forest and use the optimal tree as the partition model.)
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
  
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/psm.R")
  source("/home/yujia/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/support_func/test_heterogeneity.R")
  
  # library and functions for testing.
  source("~/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/03_heterogeneity_testing/subgroup_analysis/find_best_tree.R")
})

# LDL - CAD Binary trait
set.seed(1234)

# ==================
# Data Preprocessing
# ==================
X.set3 <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/LDL/CAD/ukbb.covariate.LDL.set3.gz", sep = "\t")
Z.set3 <- fread("~/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/LDL/CAD/set3/1e-08/LDL_prs.best", sep = " ") # score_constd
W <- fread("~/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/LDL/CAD/ukbb.phenotype.LDL.mgdL", sep = "\t")
Y.date3 <- fread("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date3_james2022.gz", sep = "\t")

# homonize the ID included in the analysis
selected_id_set3 <- Reduce(intersect, list(X.set3$IID, Z.set3$IID, W$IID, Y.date3$IID))

# model 3
X.model3 <- X.set3[X.set3$IID %in% selected_id_set3, ]
Y.model3 <- Y.date3[Y.date3$IID %in% selected_id_set3, ]
W.model3 <- W[W$IID %in% selected_id_set3, ]
Z.model3 <- Z.set3[Z.set3$IID %in% selected_id_set3, ]
X.model3 <- cbind(X.model3, Y.model3[, c("htn", "t2dm", "heart_failure", "hemorrhage_stroke", "ischemic_stroke")])
Y.model3 <- Y.model3[, c("IID", "CAD")]
message("Patient Info for model 3")
dim(X.model3); dim(Y.model3); dim(W.model3); dim(Z.model3);
print(table(Y.model3$CAD))

# generate mat file for model 3
W.model3.mat <- W.model3[, c("30780-0.0")]
X.model3.mat <- X.model3[, 3:dim(X.model3)[2]]
Y.model3.mat <- Y.model3[, c("CAD")]
Z.model3.mat <- Z.model3[, 4]
dim(W.model3.mat); dim(X.model3.mat); dim(Y.model3.mat); dim(Z.model3.mat)

# rename the mat colnames
X.model3.mat.PC  <- X.model3.mat  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0", 
                         "PC1"="22009-0.1", "PC2"="22009-0.2", "PC3"="22009-0.3", "PC4"="22009-0.4", "PC5"="22009-0.5",                 
                         "PC6"="22009-0.6", "PC7"="22009-0.7", "PC8"="22009-0.8", "PC9"="22009-0.9", "PC10"="22009-0.10", "Genotype Batch"="22000-0.0",
                         "BMI"="21001-0.0", "Body Fat Percentage"="23099-0.0", "Weight"="21002-0.0", # "Waist-hip-ratio"="whr", 
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "Apolipoprotein A"="30630-0.0", "HDL-C"="30760-0.0", "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", "Calcium"="30680-0.0", 
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0", "Total protein"="30860-0.0",
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")

X.model3.mat.noPC  <- X.model3.mat  %>% 
                  rename("Gender"="22001-0.0", "Age"="21022-0.0", "diastolic blood pressure"="4079-0.0", "systolic blood pressure"="4080-0.0", "Townsend deprivation index"="189-0.0", 
                         "BMI"="21001-0.0", "Body Fat Percentage"="23099-0.0", "Weight"="21002-0.0", # "Waist-hip-ratio"="whr", 
                         "Blood pressure medication"="Blood_pressure_medication", "Cholesterol lowering medication"="Cholesterol_lowering_medication", "Insulin"="Insulin", # "No medication"="No_medication", 
                         "Non-alcohol drinker"="Non_alcohol_drinker" , "Previous alcohol drinker"="Previous_alcohol_drinker", "Current alcohol drinker"="Current_alcohol_drinker",
                         "Non-smoker"="Non_smoker" , "Previous smoker"="Previous_smoker", "Current smoker"="Current_smoker",
                         "Apolipoprotein A"="30630-0.0", "HDL-C"="30760-0.0", "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", "Calcium"="30680-0.0", 
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0", "Total protein"="30860-0.0",
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")

X.model3.mat.noPC <- X.model3.mat.PC[, -c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Genotype Batch")]

# model 3
W.model3.vector <- as.vector(W.model3.mat$`30780-0.0`)
Y.model3.vector <- as.integer(as.vector(Y.model3.mat$CAD))
Z.model3.vector <- as.vector(Z.model3.mat$PRS)
Z.model3.vector.binary <- as.integer(ifelse(Z.model3.mat$PRS <= 0, 1, 0))
W.model3.vector.binary <- as.integer(ifelse(W.model3.vector <= 130, 1, 0)) # mimic usage of statin

#' =======================
#' Start Analysis
#' =======================
# Run the validation test

Y.rf.model <- regression_forest(X.model3.mat.PC, Y.model3.vector,
                                num.threads = 10, num.trees = 5000,
                                sample.fraction = 0.02, min.node.size = 200, ci.group.size = 100)
Y.hat <- predict(Y.rf.model)$predictions

W.rf.model <- regression_forest(X.model3.mat.PC, W.model3.vector.binary,
                                num.threads = 10, num.trees = 5000,
                                sample.fraction = 0.02, min.node.size = 200, ci.group.size = 100)
W.hat <- predict(W.rf.model)$predictions

Z.rf.model <- regression_forest(X.model3.mat.PC, Z.model3.vector,
                                num.threads = 10, num.trees = 5000,
                                sample.fraction = 0.02, min.node.size = 200, ci.group.size = 100)
Z.hat <- predict(Z.rf.model)$predictions

iv.forest.continuousZ.validatemodel <- instrumental_forest(X.model3.mat.noPC, Y.model3.vector, W.model3.vector.binary, Z.model3.vector, 
                                                           Y.hat = Y.hat,
                                                           W.hat = W.hat,
                                                           Z.hat = Z.hat,
                                                           num.threads = 10, num.trees = 5000,
                                                           sample.fraction = 0.02, min.node.size = 200, ci.group.size = 100)
# iv.pred.validate <- predict(iv.forest.continuousZ.validatemodel, estimate.variance = TRUE) # iv.forest.continuousZ.rawmodel
# iv.pred.stats.validate <- compute_stats(iv.pred.validate)

best.tree.index <- find_best_tree(iv.forest.continuousZ.validatemodel, type = "instrumental")
best.tree <- get_tree(iv.forest.continuousZ.validatemodel, best.tree.index$best_tree)
subgroup.index <- subgroup_patients(best.tree, X.model3.mat.noPC)

tree_plot <- plot(best.tree)
tree_plot

subgroup.index.train <- subgroup.index[best.tree$drawn_samples]
subgroup.index.test <- subgroup.index[-best.tree$drawn_samples]

Y.train <- Y.model3.vector[best.tree$drawn_samples]
W.train <- W.model3.vector.binary[best.tree$drawn_samples]
X.train <- X.model3.mat.noPC[best.tree$drawn_samples, ]
Z.train <- Z.model3.vector[best.tree$drawn_samples]

Y.test <- Y.model3.vector[-best.tree$drawn_samples]
W.test <- W.model3.vector.binary[-best.tree$drawn_samples]
X.test <- X.model3.mat.noPC[-best.tree$drawn_samples, ]
Z.test <- Z.model3.vector[-best.tree$drawn_samples]

dat.test <- cbind(Y.test, W.test, Z.test, X.test, subgroup.index.test)
dat.train <- cbind(Y.train, W.train, Z.train, X.train, subgroup.index.train)
colnames(dat.test)[1] <- "Y"
colnames(dat.test)[2] <- "W"
colnames(dat.test)[3] <- "Z"
colnames(dat.test)[dim(dat.test)[2]] <- "subgroup"
colnames(dat.train)[1] <- "Y"
colnames(dat.train)[2] <- "W"
colnames(dat.train)[3] <- "Z"
colnames(dat.train)[dim(dat.train)[2]] <- "subgroup"

table(dat.test$subgroup)
table(dat.train$subgroup)

# tau.test <- iv.pred.validate$predictions[-best.tree$drawn_samples]
# tau.train <- iv.pred.validate$predictions[best.tree$drawn_samples]
# tau.test.subgroup <- as.data.frame(cbind(tau.test, subgroup.index.test))
# tau.train.subgroup <- as.data.frame(cbind(tau.train, subgroup.index.train))

# Alejandro approach
results.list <- list()
results.list.individual <- list()
for (i in 1:length(table(dat.test$subgroup))){

  # for (j in 1:length(table(dat.test$subgroup))){
    
  #   if (i > j){
  #     group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
  #     group.index.b <- as.numeric(names(table(dat.test$subgroup))[j])
      
  #     ivreg.res.test.a <- ivreg(formula = Y ~ W + . - Z - subgroup | Z + . - subgroup,
  #                               data = dat.test[dat.test$subgroup == group.index.a,]) # https://www.econometrics-with-r.org/12-2-TGIVRM.html
  #     ivreg.test.summary.a <- summary(ivreg.res.test.a)
  #     ivreg.res.test.b <- ivreg(formula = Y ~ W + . - Z - subgroup | Z + . - subgroup, 
  #                               data = dat.test[dat.test$subgroup == group.index.b,])
  #     ivreg.test.summary.b <- summary(ivreg.res.test.b)
      
  #     mean.a <- ivreg.test.summary.a$coefficients[2,1]
  #     std.error.a <- ivreg.test.summary.a$coefficients[2,2]  
  #     mean.b <- ivreg.test.summary.b$coefficients[2,1]
  #     std.error.b <- ivreg.test.summary.b$coefficients[2,2] 
      
  #     T_value.ab <- (mean.a-mean.b)/sqrt(std.error.a**2 + std.error.b**2)
  #     p_value.ab <- pt(T_value.ab, (table(dat.test$subgroup)[[i]]+table(dat.test$subgroup)[[i]]-2), lower.tail = F)
      
  #     T_value.ba <- (mean.b-mean.a)/sqrt(std.error.a**2 + std.error.b**2)
  #     p_value.ba <- pt(T_value.ba, (table(dat.test$subgroup)[[j]]+table(dat.test$subgroup)[[j]]-2), lower.tail = F)
      
  #     message("         ")
  #     print(paste0(group.index.a, " : ", group.index.b))
  #     print(paste0(mean.a, " : ", mean.b))
  #     print(p_value.ab)
  #     print(p_value.ba)
      
  #     if (min(p_value.ab, p_value.ba) < 0.05){
  #       # print(max(mean.a, mean.b) + 1.96*ifelse(mean.a > mean.b, std.error.a, std.error.b))
  #       results.list[[paste0(group.index.a, " : ", group.index.b)]] <- c(min(p_value.ab, p_value.ba), gtools::stars.pval(min(p_value.ab, p_value.ba)), 
  #                                                                        names(table(dat.test$subgroup))[i], names(table(dat.test$subgroup))[j], 
  #                                                                        max(mean.a, mean.b) + 1.96*ifelse(mean.a > mean.b, std.error.a, std.error.b))
  #     }
  #   }
  # }

  group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
  # ivreg.res.test.a <- ivreg(formula = Y ~ W + . - Z - subgroup | Z + . - subgroup,
# data = dat.test[dat.test$subgroup == group.index.a,]) # https://www.econometrics-with-r.org/12-2-TGIVRM.html
  ivreg.res.test.a <- ivreg(formula = Y ~ . - W - Z - subgroup | W | Z ,
                                data = dat.test[dat.test$subgroup == group.index.a,]) # https://www.econometrics-with-r.org/12-2-TGIVRM.html

  
  ivreg.test.summary.a <- summary(ivreg.res.test.a)
  mean.a <- ivreg.test.summary.a$coefficients[2,1]
  std.error.a <- ivreg.test.summary.a$coefficients[2,2]
  results.list.individual[[i]] <- c(table(dat.test$subgroup)[i], mean.a, std.error.a)
    
}

print(results.list.individual)
results.mat <- matrix(nrow = length(results.list.individual), ncol = 4)

for (i in 1:length(results.list.individual)){
  results.mat[i, 1] <- i
  results.mat[i, 2] <- results.list.individual[[i]][1]
  results.mat[i, 3] <- results.list.individual[[i]][2]
  results.mat[i, 4] <- results.list.individual[[i]][3]
}

results.mat
# results.list <- results.list[order(sapply(results.list, function(x) x[5], simplify=TRUE), decreasing=TRUE)]

mean.cf <- c()
std.error.cf <- c()
class.cf <- c()
sample.size.cf <- c()

for (i in 1:length(table(dat.test$subgroup))){
  
  group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
  # ivreg.res.test.a <- ivreg(formula = Y ~ W + . - Z - subgroup  | Z + . - subgroup, 
  #                           data = dat.test[dat.test$subgroup == group.index.a,])
  ivreg.res.test.a <- ivreg(formula = Y ~ . - W - Z - subgroup | W | Z ,
                            data = dat.test[dat.test$subgroup == group.index.a,])
  ivreg.test.summary.a <- summary(ivreg.res.test.a)
  
  mean.cf <- c(mean.cf, ivreg.test.summary.a$coefficients[2,1])
  std.error.cf <- c(std.error.cf, ivreg.test.summary.a$coefficients[2,2])
  sample.size.cf <- c(sample.size.cf, dim(dat.test[dat.test$subgroup == group.index.a,])[1])
  class.cf <- c(class.cf, names(table(dat.test$subgroup))[i])
}

sample.df <- data.frame(mean = mean.cf, std.error = std.error.cf, sample.size = sample.size.cf, class = class.cf)
sample.df$std.dev <- sample.df$std.error * sqrt(sample.df$sample.size)
sample.df <- sample.df[order(sample.df$mean), ]
sample.df$class <- as.factor(c(1,2,3,4,5))
rownames(sample.df) <- NULL

# comparsion_list <- list()
# annotations_vector <- c()
# y_position_vector <- c()
# for (i in 1:length(results.list)){
#   res.tmp <- results.list[[i]]
#   annotations_vector <- c(annotations_vector, res.tmp[2])
#   comparsion_list[[i]] <- c(res.tmp[3], res.tmp[4])
#   y_position_vector <- c(y_position_vector, as.numeric(res.tmp[5]))
# }

# y_position_vector_tmp <- y_position_vector
# for (i in 1:length(y_position_vector)){
#   if (i == 1){
#     y_position_vector_tmp[which(order(y_position_vector) == i)] <- max(y_position_vector_tmp) + 0.001
#   } else {
#     y_position_vector_tmp[which(order(y_position_vector) == i)] <- max(y_position_vector_tmp) + 0.002
#   }
# }

subgroup_matching_list <- list()
group_index <- 1
for (i in as.character(sort(as.numeric(names(table(sample.df$class)))))){
  subgroup_matching_list[[i]] <- LETTERS[group_index]
  group_index <- group_index + 1
}

p.validation <- ggplot(sample.df, aes(factor(class), y=mean, ymin=mean - 1.96 * std.error, ymax = mean + 1.96 * std.error)) +
  geom_point(size = 5) + 
  geom_errorbar(width = 0.3, size = 1) + # geom_signif(comparisons = comparsion_list, annotations = annotations_vector, y_position = y_position_vector_tmp) +
  theme_classic() + xlab("Subgroup") + ylab("Local Average Treatment Effect") + ggtitle("LDL-C\nBinary Treatment") + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 20),
        axis.text = element_text(size = 25), axis.line=element_line(color='black')) + scale_x_discrete(labels = subgroup_matching_list)
p.validation

ggsave("/home/yujia/Project/2023-07-20-individual_MR/res/03_plot/supp_fig1_LDL_CAD_binary.png", plot = p.validation, dpi = 300,
       height = 10, width = 10, scale = 1)

# save.image("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/03_validation_update/binaryW_LDL_CAD.RData")
# ========================================================================================

# library("grid")
# library("ggplotify")

# # plotting
# tree_plot
# cat(DiagrammeRsvg::export_svg(tree_plot), file = '/home/yujia/Project/2023-07-20-individual_MR/res/03_plot/update_besttree/binaryLDL_CAD_besttree.svg')
# # ggsave(filename = "/home/yujia/Project/2023-07-20-individual_MR/res/03_plot/update_besttree/binaryLDL_CAD_boxplot.png", plot = p.validation, dpi = 300)
# sample.df

# ANOVA
# https://statpages.info/anova1sm.html ANOVA res: 0.0048

# ========================================================================================

# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# m0: the null value for the difference in means to be tested for. Default is 0. 
# equal.variance: whether or not to assume equal variance. Default is FALSE. 
t.test2 <- function(m1,m2,s1,s2,n1,n2, m0=0, equal.variance=FALSE){
    if( equal.variance==FALSE ) 
    {
        se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        # welch-satterthwaite df
        df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    } else
    {
        # pooled standard deviation, scaled by the sample sizes
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
        df <- n1+n2-2
    }      
    t <- (m1-m2-m0)/se 
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat) 
}

t.test2(sample.df$mean[1], sample.df$mean[5], sample.df$std.dev[1], sample.df$std.dev[5], sample.df$sample.size[1], sample.df$sample.size[5])
