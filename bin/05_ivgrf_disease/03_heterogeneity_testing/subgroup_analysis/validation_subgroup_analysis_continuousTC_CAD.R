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
  # library(xlsx)
  library(ivreg)
  library(ggsignif) 
  
  source("/mnt/md0/yujia/project/2023-07-20-individual_MR/bin/05_ivgrf_disease/support_func/psm.R")
  source("/mnt/md0/yujia/project/2023-07-20-individual_MR/bin/05_ivgrf_disease/support_func/test_heterogeneity.R")
  
  # library and functions for testing.
  source("/mnt/md0/yujia/project/2023-07-20-individual_MR/bin/05_ivgrf_disease/03_heterogeneity_testing/subgroup_analysis/find_best_tree.R")
})

# TC - CAD Binary trait

# ==================
# Data Preprocessing
# ==================
X.set3 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TC/CAD/ukbb.covariate.TC.set3.gz", sep = "\t")
Z.set3 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/TC/CAD/set3/1e-08/TC_prs.best", sep = " ") # score_constd
W <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TC/CAD/ukbb.phenotype.TC.mgdL", sep = "\t")
Y.date3 <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date3_james2022.gz", sep = "\t")

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
W.model3.mat <- W.model3[, c("30690-0.0")]
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
                         "Triglycerides"="30870-0.0", # lipid-related covariates
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
                         "Triglycerides"="30870-0.0", # lipid-related covariates
                         "Creatinine"="30700-0.0", "C-reactive protein"="30710-0.0", "Cystatin C"="30720-0.0", "Gamma glutamyltransferase"="30730-0.0", "Calcium"="30680-0.0",
                         "Glucose"="30740-0.0", "HbA1c"="30750-0.0", "Aspartate aminotransferase"="30650-0.0", "Direct bilirubin"="30660-0.0", 
                         "Urea"="30670-0.0", "IGF-1"="30770-0.0", "Phosphate"="30810-0.0", "SHBG"="30830-0.0", "Testosterone"="30850-0.0", 
                         "Urate"="30880-0.0", "Vitamin D"="30890-0.0", "Total bilirubin"="30840-0.0", "Total protein"="30860-0.0", 
                         "Type 2 diabetes history"="t2dm", "Hypertension history"="htn", "Heart failure history"="heart_failure", "Hemorrhage Stroke history"="hemorrhage_stroke", "Ischemic Stroke history"="ischemic_stroke")

X.model3.mat.noPC <- X.model3.mat.PC[, -c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Genotype Batch")]

# model 3
W.model3.vector <- as.vector(W.model3.mat$`30690-0.0`)
Y.model3.vector <- as.integer(as.vector(Y.model3.mat$CAD))
Z.model3.vector <- as.vector(Z.model3.mat$PRS)
Z.model3.vector.binary <- as.integer(ifelse(Z.model3.mat$PRS <= 0, 1, 0))
W.model3.vector.binary <- as.integer(ifelse(W.model3.vector <= 220, 1, 0)) # mimic usage of statin

#' =======================
#' Start Analysis
#' =======================
# Run the validation test

Y.rf.model <- regression_forest(X.model3.mat.PC, Y.model3.vector,
                                num.threads = 10, num.trees = 5000,
                                sample.fraction = 0.02, min.node.size = 200, ci.group.size = 100)
Y.hat <- predict(Y.rf.model)$predictions

W.rf.model <- regression_forest(X.model3.mat.PC, W.model3.vector,
                                num.threads = 10, num.trees = 5000,
                                sample.fraction = 0.02, min.node.size = 200, ci.group.size = 100)
W.hat <- predict(W.rf.model)$predictions

Z.rf.model <- regression_forest(X.model3.mat.PC, Z.model3.vector,
                                num.threads = 10, num.trees = 5000,
                                sample.fraction = 0.02, min.node.size = 200, ci.group.size = 100)
Z.hat <- predict(Z.rf.model)$predictions

iv.forest.continuousZ.validatemodel <- instrumental_forest(X.model3.mat.noPC, Y.model3.vector, W.model3.vector, Z.model3.vector, 
                                                           Y.hat = Y.hat,
                                                           W.hat = W.hat,
                                                           Z.hat = Z.hat,
                                                           num.threads = 10, num.trees = 5000, 
                                                           sample.fraction = 0.025, min.node.size = 200, ci.group.size = 100)
# iv.pred.validate <- predict(iv.forest.continuousZ.validatemodel, estimate.variance = TRUE) # iv.forest.continuousZ.rawmodel
# iv.pred.stats.validate <- compute_stats(iv.pred.validate)

best.tree.index <- find_best_tree(iv.forest.continuousZ.validatemodel, type = "instrumental")
best.tree <- get_tree(iv.forest.continuousZ.validatemodel, best.tree.index$best_tree)

tree_plot <- plot(best.tree)
tree_plot

subgroup.index <- subgroup_patients(best.tree, X.model3.mat.noPC)
subgroup.index.train <- subgroup.index[best.tree$drawn_samples]
subgroup.index.test <- subgroup.index[-best.tree$drawn_samples]

Y.train <- Y.model3.vector[best.tree$drawn_samples]
W.train <- W.model3.vector[best.tree$drawn_samples]
X.train <- X.model3.mat.noPC[best.tree$drawn_samples, ]
Z.train <- Z.model3.vector[best.tree$drawn_samples]

Y.test <- Y.model3.vector[-best.tree$drawn_samples]
W.test <- W.model3.vector[-best.tree$drawn_samples]
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
table(dat.train$subgroup)/2

# tau.test <- iv.pred.validate$predictions[-best.tree$drawn_samples]
# tau.train <- iv.pred.validate$predictions[best.tree$drawn_samples]
# tau.test.subgroup <- as.data.frame(cbind(tau.test, subgroup.index.test))
# tau.train.subgroup <- as.data.frame(cbind(tau.train, subgroup.index.train))

# Alejandro approach
results.list <- list()
results.list.individual <- list()
for (i in 1:length(table(dat.test$subgroup))){

  group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
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

mean.cf <- c()
std.error.cf <- c()
class.cf <- c()
sample.size.cf <- c()

for (i in 1:length(table(dat.test$subgroup))){
  
  group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
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
sample.df$class
sample.df$class <- as.factor(c(1,2,3,4,5,6))
rownames(sample.df) <- NULL

# 7 - 1 - A
# 11 - 2 - B
# 5 - 3 - C
# 10 - 4 - D
# 6 - 5 - E
# 8 - 6 - F

subgroup_matching_list <- list()
group_index <- 1
for (i in as.character(sort(as.numeric(names(table(sample.df$class)))))){
  subgroup_matching_list[[i]] <- LETTERS[group_index]
  group_index <- group_index + 1
}

p.validation <- ggplot(sample.df, aes(factor(class), y=mean, ymin=mean - 1.96 * std.error, ymax = mean + 1.96 * std.error)) +
  geom_point(size = 5) + 
  geom_errorbar(width = 0.3, size = 1) + # geom_signif(comparisons = comparsion_list, annotations = annotations_vector, y_position = y_position_vector_tmp) +
  theme_classic() + xlab("Subgroup") + ylab("Local Average Treatment Effect") + ggtitle("Total Cholesterol\nContinuous Treatment") + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 20),
        axis.text = element_text(size = 25), axis.line=element_line(color='black')) + scale_x_discrete(labels = subgroup_matching_list)
p.validation

ggsave("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/supp_fig1_TC_CAD_continuous.png", plot = p.validation, dpi = 300,
       height = 10, width = 10, scale = 1)

# save.image("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/03_validation/continuousW_TC_CAD.RData")
# ========================================================================================

# library("grid")
# library("ggplotify")

# # plotting
# tree_plot
# cat(DiagrammeRsvg::export_svg(tree_plot), file = '/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/update_besttree/continuousW_TC_CAD_besttree.svg')
# # tree_plot %>% export_svg %>% charToRaw %>% rsvg_pdf("'/mnt/md0/yujia/project/2023-07-20-individual_MR/res/03_plot/continuousW_TC_CAD_besttree.pdf")


