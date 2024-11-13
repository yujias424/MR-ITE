#' This code is to perform extra validation on the heterogeneity testing results.
#' The idea comes from the discussion between Alex and Alejandro,
#' details can be found in related email.
#' In breif, we have two approaches:
#' 1. Retrain a CF model and use the optimal tree as the partitioning model.
#' 2. Retrain a regular regression tree model where the target is the taus (i.e., Regression Forest and use the optimal tree as the partition model.)
#'
#' @author Yujia Shi
#' @date 2022-09-23

suppressPackageStartupMessages({
  
  library(grf)
  library(data.table)
  source("~/Project/2022-09-01-individual_MR/bin/05_ivgrf_disease/heterogeneity_testing/subgroup_analysis/find_best_tree.R")

})

# Old approach

#' =====================
#' Continuous W, Z
#' =====================
iv.forest.continuousZ.rawmodel <- instrumental_forest(X, Y.vector, W.vector, Z.vector, 
                                 num.threads = 30, num.trees = 10000, 
                                 sample.fraction = 0.15, min.node.size = 100, ci.group.size = 100) 
iv.pred.raw <- predict(iv.forest.continuousZ.rawmodel, estimate.variance = TRUE)
iv.pred.stats.raw <- compute_stats(iv.pred.raw)

iv.forest.continuousZ.validatemodel <- instrumental_forest(X, Y.vector, W.vector, Z.vector, 
                                             num.threads = 30, num.trees = 10000, 
                                             sample.fraction = 0.15, min.node.size = 100, ci.group.size = 100)
iv.pred.validate <- predict(iv.forest.continuousZ.rawmodel, estimate.variance = TRUE)
iv.pred.stats.validate <- compute_stats(iv.pred.validate)

best.tree.index <- find_best_tree(iv.forest.continuousZ.validatemodel, type = "instrumental")
best.tree <- get_tree(iv.forest.continuousZ.validatemodel, best.tree.index$best_tree)
subgroup.index <- subgroup_patients(best.tree, X)

subgroup.index.train <- subgroup.index[best.tree$drawn_samples]
subgroup.index.test <- subgroup.index[-best.tree$drawn_samples]

Y.train <- Y[best.tree$drawn_samples]
W.train <- W[best.tree$drawn_samples]
X.train <- X[best.tree$drawn_samples, ]
Z.train <- Z[best.tree$drawn_samples]

Y.test <- Y[-best.tree$drawn_samples]
W.test <- W[-best.tree$drawn_samples]
X.test <- X[-best.tree$drawn_samples, ]
Z.test <- Z[-best.tree$drawn_samples]

dat.test <- cbind(Y.test, W.test, Z.test, X.test, subgroup.index.test)
dat.train <- cbind(Y.train, W.train, Z.train, X.train, subgroup.index.train)
colnames(dat.test)[1] <- "Y"
colnames(dat.test)[2] <- "W"
colnames(dat.test)[3] <- "Z"
colnames(dat.test)[53] <- "subgroup"
colnames(dat.train)[1] <- "Y"
colnames(dat.train)[2] <- "W"
colnames(dat.train)[3] <- "Z"
colnames(dat.train)[53] <- "subgroup"

colnames(dat)
table(dat.test$subgroup)
table(dat.train$subgroup)

tau.test <- iv.pred.validate$predictions[-best.tree$drawn_samples]
tau.train <- iv.pred.validate$predictions[best.tree$drawn_samples]
tau.test.subgroup <- as.data.frame(cbind(tau.test, subgroup.index.test))
tau.train.subgroup <- as.data.frame(cbind(tau.train, subgroup.index.train))

table(dat.test$subgroup)

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
      
      if (min(p_value.ab, p_value.ba) < 0.1){
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

for (i in 1:length(table(dat.test$subgroup))){
  
  group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
  ivreg.res.test.a <- ivreg(formula = Y ~ W + . - Z - subgroup  | Z + . - subgroup, 
                            data = dat.test[dat.test$subgroup == group.index.a,])
  ivreg.test.summary.a <- summary(ivreg.res.test.a)
  
  mean.cf <- c(mean.cf, ivreg.test.summary.a$coefficients[2,1])
  std.error.cf <- c(std.error.cf, ivreg.test.summary.a$coefficients[2,2])
  class.cf <- c(class.cf, names(table(dat.test$subgroup))[i])
}

sample.df <- data.frame(mean = mean.cf, std.error = std.error.cf, class = class.cf)

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
    y_position_vector_tmp[which(order(y_position_vector) == i)] <- max(y_position_vector_tmp) + 0.005
  } else {
    y_position_vector_tmp[which(order(y_position_vector) == i)] <- max(y_position_vector_tmp) + 0.01
  }
}

subgroup_matching_list <- list()
group_index = 1
for (i in as.character(sort(as.numeric(names(table(sample.df$class)))))){
  subgroup_matching_list[[i]] <- letters[group_index]
  group_index = group_index + 1
}

p.validation <- ggplot(sample.df, aes(factor(class), y=mean, ymin=mean - 1.96 * std.error, ymax = mean + 1.96 * std.error))+
  geom_pointrange() + geom_signif(comparisons = comparsion_list, annotations = annotations_vector, y_position = y_position_vector_tmp) +
  theme_classic() + xlab("Subgroup") + ylab("LATE") + ggtitle("Coronary Artery Disease") + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 20),
        axis.text = element_text(size = 25), axis.line=element_line(color='black')) + scale_x_discrete(labels = subgroup_matching_list)
p.validation

#' =====================================================================================================================================

#' =====================
#' Binary W, Z
#' =====================
iv.forest.binaryWZ.rawmodel <- instrumental_forest(X, Y.vector, W.vector.binary.cutoff, Z.vector.binary, 
                                                      num.threads = 30, num.trees = 10000, 
                                                      sample.fraction = 0.15, min.node.size = 50, ci.group.size = 100) 
iv.pred.raw <- predict(iv.forest.binaryWZ.rawmodel, estimate.variance = TRUE)
iv.pred.stats.raw <- compute_stats(iv.pred.raw)

iv.forest.binaryWZ.validatemodel <- instrumental_forest(X, Y.vector, W.vector.binary.cutoff, Z.vector.binary, 
                                                           num.threads = 30, num.trees = 10000, 
                                                           sample.fraction = 0.15, min.node.size = 50, ci.group.size = 100)
iv.pred.validate <- predict(iv.forest.binaryWZ.validatemodel, estimate.variance = TRUE)
iv.pred.stats.validate <- compute_stats(iv.pred.validate)

best.tree.index <- find_best_tree(iv.forest.binaryWZ.validatemodel, type = "instrumental")
best.tree <- get_tree(iv.forest.binaryWZ.validatemodel, best.tree.index$best_tree)
subgroup.index <- subgroup_patients(best.tree, X)

subgroup.index.train <- subgroup.index[best.tree$drawn_samples]
subgroup.index.test <- subgroup.index[-best.tree$drawn_samples]

Y.train <- Y[best.tree$drawn_samples]
W.train <- W.vector.binary.cutoff[best.tree$drawn_samples]
X.train <- X[best.tree$drawn_samples, ]
Z.train <- Z.vector.binary[best.tree$drawn_samples]

Y.test <- Y[-best.tree$drawn_samples]
W.test <- W.vector.binary.cutoff[-best.tree$drawn_samples]
X.test <- X[-best.tree$drawn_samples, ]
Z.test <- Z.vector.binary[-best.tree$drawn_samples]

dat.test <- cbind(Y.test, W.test, Z.test, X.test, subgroup.index.test)
dat.train <- cbind(Y.train, W.train, Z.train, X.train, subgroup.index.train)
colnames(dat.test)[1] <- "Y"
colnames(dat.test)[2] <- "W"
colnames(dat.test)[3] <- "Z"
colnames(dat.test)[53] <- "subgroup"
colnames(dat.train)[1] <- "Y"
colnames(dat.train)[2] <- "W"
colnames(dat.train)[3] <- "Z"
colnames(dat.train)[53] <- "subgroup"

colnames(dat)
table(dat.test$subgroup)
table(dat.train$subgroup)

tau.test <- iv.pred.validate$predictions[-best.tree$drawn_samples]
tau.train <- iv.pred.validate$predictions[best.tree$drawn_samples]
tau.test.subgroup <- as.data.frame(cbind(tau.test, subgroup.index.test))
tau.train.subgroup <- as.data.frame(cbind(tau.train, subgroup.index.train))

table(dat.test$subgroup)

# Alejandro approach
results.list <- list()
for (i in 1:length(table(dat.test$subgroup))){
  for (j in 1:length(table(dat.test$subgroup))){
    
    if (i > j){
      
      group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
      group.index.b <- as.numeric(names(table(dat.test$subgroup))[j])
      
      Y.a <- dat.test[dat.test$subgroup == group.index.a,]$Y
      W.a <- dat.test[dat.test$subgroup == group.index.a,]$W
      Z.a <- dat.test[dat.test$subgroup == group.index.a,]$Z
      X.a <- dat.test[dat.test$subgroup == group.index.a, -c("Y", "W", "Z", "subgroup")]
      
      ivreg.res.test.a <- instrumental_forest(X.a, Y.a, W.a, Z.a, 
                                              num.threads = 30, num.trees = 10000, 
                                              sample.fraction = 0.1, min.node.size = 10, ci.group.size = 100)
      late.a <- average_treatment_effect(ivreg.res.test.a)
      
      Y.b <- dat.test[dat.test$subgroup == group.index.b,]$Y
      W.b <- dat.test[dat.test$subgroup == group.index.b,]$W
      Z.b <- dat.test[dat.test$subgroup == group.index.b,]$Z
      X.b <- dat.test[dat.test$subgroup == group.index.b, -c("Y", "W", "Z", "subgroup")]
      
      ivreg.res.test.b <- instrumental_forest(X.b, Y.b, W.b, Z.b, 
                                              num.threads = 30, num.trees = 10000, 
                                              sample.fraction = 0.1, min.node.size = 10, ci.group.size = 100)
      late.b <- average_treatment_effect(ivreg.res.test.b)
      
      T_value.ab <- (late.a[1]-late.b[1])/sqrt(late.a[2]**2 + late.b[2]**2)
      pt(T_value.ab, dim(X.a)[1] + dim(X.b)[2] - 2)
      
      mean.a <- late.a[1]
      std.error.a <- late.a[2] 
      mean.b <- late.b[1]
      std.error.b <- late.b[2]
      
      T_value.ab <- (mean.a-mean.b)/sqrt(std.error.a**2 + std.error.b**2)
      p_value.ab <- pt(T_value.ab, (table(dat.test$subgroup)[[i]]+table(dat.test$subgroup)[[i]]-2), lower.tail = F)
      
      T_value.ba <- (mean.b-mean.a)/sqrt(std.error.a**2 + std.error.b**2)
      p_value.ba <- pt(T_value.ba, (table(dat.test$subgroup)[[j]]+table(dat.test$subgroup)[[j]]-2), lower.tail = F)
      
      message("         ")
      print(paste0(group.index.a, " : ", group.index.b))
      print(paste0(mean.a, " : ", mean.b))
      print(p_value.ab)
      print(p_value.ba)
      
      if (min(p_value.ab, p_value.ba) < 0.1){
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

for (i in 1:length(table(dat.test$subgroup))){
  
  group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
  ivreg.res.test.a <- ivreg(formula = Y ~ W + . - Z - subgroup  | Z + . - subgroup, 
                            data = dat.test[dat.test$subgroup == group.index.a,])
  ivreg.test.summary.a <- summary(ivreg.res.test.a)
  
  mean.cf <- c(mean.cf, ivreg.test.summary.a$coefficients[2,1])
  std.error.cf <- c(std.error.cf, ivreg.test.summary.a$coefficients[2,2])
  class.cf <- c(class.cf, names(table(dat.test$subgroup))[i])
}

sample.df <- data.frame(mean = mean.cf, std.error = std.error.cf, class = class.cf)

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
    y_position_vector_tmp[which(order(y_position_vector) == i)] <- max(y_position_vector_tmp) + 0.005
  } else {
    y_position_vector_tmp[which(order(y_position_vector) == i)] <- max(y_position_vector_tmp) + 0.01
  }
}

subgroup_matching_list <- list()
group_index = 1
for (i in as.character(sort(as.numeric(names(table(sample.df$class)))))){
  subgroup_matching_list[[i]] <- letters[group_index]
  group_index = group_index + 1
}

p.validation <- ggplot(sample.df, aes(factor(class), y=mean, ymin=mean - 1.96 * std.error, ymax = mean + 1.96 * std.error))+
  geom_pointrange() + geom_signif(comparisons = comparsion_list, annotations = annotations_vector, y_position = y_position_vector_tmp) +
  theme_classic() + xlab("Subgroup") + ylab("LATE") + ggtitle("Coronary Artery Disease") + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 20),
        axis.text = element_text(size = 25), axis.line=element_line(color='black')) + scale_x_discrete(labels = subgroup_matching_list)
p.validation

#' ====================================================================================================================================
#' New approach
#' ====================================================================================================================================

#' =====================
#' Continuous W, Z
#' =====================
iv.forest.continuousZ.rawmodel <- instrumental_forest(X, Y.vector, W.vector, Z.vector, 
                                                      num.threads = 30, num.trees = 10000, 
                                                      sample.fraction = 0.15, min.node.size = 100, ci.group.size = 100) 
iv.pred.raw <- predict(iv.forest.continuousZ.rawmodel, estimate.variance = TRUE)
iv.pred.stats.raw <- compute_stats(iv.pred.raw)

iv.forest.continuousZ.validatemodel <- regression_forest(X, iv.pred.raw$predictions,  
                                                           num.threads = 30, num.trees = 10000, 
                                                           sample.fraction = 0.15, min.node.size = 450)
iv.pred.validate <- predict(iv.forest.continuousZ.validatemodel, estimate.variance = TRUE)
iv.pred.stats.validate <- compute_stats(iv.pred.validate)
get_tree(iv.forest.continuousZ.validatemodel, 1)

best.tree.index <- find_best_tree(iv.forest.continuousZ.validatemodel, type = "regression")
best.tree <- get_tree(iv.forest.continuousZ.validatemodel, best.tree.index$best_tree)
subgroup.index <- subgroup_patients(best.tree, X)

subgroup.index.train <- subgroup.index[best.tree$drawn_samples]
subgroup.index.test <- subgroup.index[-best.tree$drawn_samples]

Y.train <- Y.vector[best.tree$drawn_samples]
W.train <- W.vector[best.tree$drawn_samples]
X.train <- X[best.tree$drawn_samples, ]
Z.train <- Z.vector[best.tree$drawn_samples]

Y.test <- Y.vector[-best.tree$drawn_samples]
W.test <- W.vector[-best.tree$drawn_samples]
X.test <- X[-best.tree$drawn_samples, ]
Z.test <- Z.vector[-best.tree$drawn_samples]

dat.test <- cbind(Y.test, W.test, Z.test, X.test, subgroup.index.test)
dat.train <- cbind(Y.train, W.train, Z.train, X.train, subgroup.index.train)
colnames(dat.test)[1] <- "Y"
colnames(dat.test)[2] <- "W"
colnames(dat.test)[3] <- "Z"
colnames(dat.test)[53] <- "subgroup"
colnames(dat.train)[1] <- "Y"
colnames(dat.train)[2] <- "W"
colnames(dat.train)[3] <- "Z"
colnames(dat.train)[53] <- "subgroup"

colnames(dat)
table(dat.test$subgroup)
table(dat.train$subgroup)

tau.test <- iv.pred.validate$predictions[-best.tree$drawn_samples]
tau.train <- iv.pred.validate$predictions[best.tree$drawn_samples]
tau.test.subgroup <- as.data.frame(cbind(tau.test, subgroup.index.test))
tau.train.subgroup <- as.data.frame(cbind(tau.train, subgroup.index.train))

table(dat.test$subgroup)

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
      
      if (min(p_value.ab, p_value.ba) < 0.1){
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

for (i in 1:length(table(dat.test$subgroup))){
  
  group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
  ivreg.res.test.a <- ivreg(formula = Y ~ W + . - Z - subgroup  | Z + . - subgroup, 
                            data = dat.test[dat.test$subgroup == group.index.a,])
  ivreg.test.summary.a <- summary(ivreg.res.test.a)
  
  mean.cf <- c(mean.cf, ivreg.test.summary.a$coefficients[2,1])
  std.error.cf <- c(std.error.cf, ivreg.test.summary.a$coefficients[2,2])
  class.cf <- c(class.cf, names(table(dat.test$subgroup))[i])
  
}

sample.df <- data.frame(mean = mean.cf, std.error = std.error.cf, class = class.cf)

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
    y_position_vector_tmp[which(order(y_position_vector) == i)] <- max(y_position_vector_tmp) + 0.005
  } else {
    y_position_vector_tmp[which(order(y_position_vector) == i)] <- max(y_position_vector_tmp) + 0.01
  }
}

subgroup_matching_list <- list()
group_index = 1
for (i in as.character(sort(as.numeric(names(table(sample.df$class)))))){
  subgroup_matching_list[[i]] <- letters[group_index]
  group_index = group_index + 1
}

p.validation <- ggplot(sample.df, aes(factor(class), y=mean, ymin=mean - 1.96 * std.error, ymax = mean + 1.96 * std.error))+
  geom_pointrange() + geom_signif(comparisons = comparsion_list, annotations = annotations_vector, y_position = y_position_vector_tmp) +
  theme_classic() + xlab("Subgroup") + ylab("LATE") + ggtitle("Coronary Artery Disease") + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 20),
        axis.text = element_text(size = 25), axis.line=element_line(color='black')) + scale_x_discrete(labels = subgroup_matching_list)
p.validation

#' =====================
#' Binary W, Z
#' =====================
iv.forest.binaryWZ.rawmodel <- instrumental_forest(X, Y.vector, W.vector.binary.cutoff, Z.vector.binary, 
                                                   num.threads = 30, num.trees = 10000, 
                                                   sample.fraction = 0.15, min.node.size = 50, ci.group.size = 100) 
iv.pred.raw <- predict(iv.forest.binaryWZ.rawmodel, estimate.variance = TRUE)
iv.pred.stats.raw <- compute_stats(iv.pred.raw)

iv.forest.binaryWZ.validatemodel <- regression_forest(X, iv.pred.raw$predictions,  
                                                         num.threads = 30, num.trees = 10000, 
                                                         sample.fraction = 0.15, min.node.size = 450)
iv.pred.validate <- predict(iv.forest.binaryWZ.validatemodel, estimate.variance = TRUE)
iv.pred.stats.validate <- compute_stats(iv.pred.validate)
get_tree(iv.forest.continuousZ.validatemodel, 1)

best.tree.index <- find_best_tree(iv.forest.binaryWZ.validatemodel, type = "instrumental")
best.tree <- get_tree(iv.forest.binaryWZ.validatemodel, best.tree.index$best_tree)
subgroup.index <- subgroup_patients(best.tree, X)

subgroup.index.train <- subgroup.index[best.tree$drawn_samples]
subgroup.index.test <- subgroup.index[-best.tree$drawn_samples]

Y.train <- Y[best.tree$drawn_samples]
W.train <- W.vector.binary.cutoff[best.tree$drawn_samples]
X.train <- X[best.tree$drawn_samples, ]
Z.train <- Z.vector.binary[best.tree$drawn_samples]

Y.test <- Y[-best.tree$drawn_samples]
W.test <- W.vector.binary.cutoff[-best.tree$drawn_samples]
X.test <- X[-best.tree$drawn_samples, ]
Z.test <- Z.vector.binary[-best.tree$drawn_samples]

dat.test <- cbind(Y.test, W.test, Z.test, X.test, subgroup.index.test)
dat.train <- cbind(Y.train, W.train, Z.train, X.train, subgroup.index.train)
colnames(dat.test)[1] <- "Y"
colnames(dat.test)[2] <- "W"
colnames(dat.test)[3] <- "Z"
colnames(dat.test)[53] <- "subgroup"
colnames(dat.train)[1] <- "Y"
colnames(dat.train)[2] <- "W"
colnames(dat.train)[3] <- "Z"
colnames(dat.train)[53] <- "subgroup"

colnames(dat)
table(dat.test$subgroup)
table(dat.train$subgroup)

tau.test <- iv.pred.validate$predictions[-best.tree$drawn_samples]
tau.train <- iv.pred.validate$predictions[best.tree$drawn_samples]
tau.test.subgroup <- as.data.frame(cbind(tau.test, subgroup.index.test))
tau.train.subgroup <- as.data.frame(cbind(tau.train, subgroup.index.train))

table(dat.test$subgroup)

# Alejandro approach
results.list <- list()
for (i in 1:length(table(dat.test$subgroup))){
  for (j in 1:length(table(dat.test$subgroup))){
    
    if (i > j){
      
      group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
      group.index.b <- as.numeric(names(table(dat.test$subgroup))[j])
      
      Y.a <- dat.test[dat.test$subgroup == group.index.a,]$Y
      W.a <- dat.test[dat.test$subgroup == group.index.a,]$W
      Z.a <- dat.test[dat.test$subgroup == group.index.a,]$Z
      X.a <- dat.test[dat.test$subgroup == group.index.a, -c("Y", "W", "Z", "subgroup")]
      
      ivreg.res.test.a <- instrumental_forest(X.a, Y.a, W.a, Z.a, 
                                              num.threads = 30, num.trees = 10000, 
                                              sample.fraction = 0.1, min.node.size = 10, ci.group.size = 100)
      late.a <- average_treatment_effect(ivreg.res.test.a)
      
      Y.b <- dat.test[dat.test$subgroup == group.index.b,]$Y
      W.b <- dat.test[dat.test$subgroup == group.index.b,]$W
      Z.b <- dat.test[dat.test$subgroup == group.index.b,]$Z
      X.b <- dat.test[dat.test$subgroup == group.index.b, -c("Y", "W", "Z", "subgroup")]
      
      ivreg.res.test.b <- instrumental_forest(X.b, Y.b, W.b, Z.b, 
                                              num.threads = 30, num.trees = 10000, 
                                              sample.fraction = 0.1, min.node.size = 10, ci.group.size = 100)
      late.b <- average_treatment_effect(ivreg.res.test.b)
      
      T_value.ab <- (late.a[1]-late.b[1])/sqrt(late.a[2]**2 + late.b[2]**2)
      pt(T_value.ab, dim(X.a)[1] + dim(X.b)[2] - 2)
      
      mean.a <- late.a[1]
      std.error.a <- late.a[2] 
      mean.b <- late.b[1]
      std.error.b <- late.b[2]
      
      T_value.ab <- (mean.a-mean.b)/sqrt(std.error.a**2 + std.error.b**2)
      p_value.ab <- pt(T_value.ab, (table(dat.test$subgroup)[[i]]+table(dat.test$subgroup)[[i]]-2), lower.tail = F)
      
      T_value.ba <- (mean.b-mean.a)/sqrt(std.error.a**2 + std.error.b**2)
      p_value.ba <- pt(T_value.ba, (table(dat.test$subgroup)[[j]]+table(dat.test$subgroup)[[j]]-2), lower.tail = F)
      
      message("         ")
      print(paste0(group.index.a, " : ", group.index.b))
      print(paste0(mean.a, " : ", mean.b))
      print(p_value.ab)
      print(p_value.ba)
      
      if (min(p_value.ab, p_value.ba) < 0.1){
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

for (i in 1:length(table(dat.test$subgroup))){
  
  group.index.a <- as.numeric(names(table(dat.test$subgroup))[i])
  
  Y.a <- dat.test[dat.test$subgroup == group.index.a,]$Y
  W.a <- dat.test[dat.test$subgroup == group.index.a,]$W
  Z.a <- dat.test[dat.test$subgroup == group.index.a,]$Z
  X.a <- dat.test[dat.test$subgroup == group.index.a, -c("Y", "W", "Z", "subgroup")]
  
  ivreg.res.test.a <- instrumental_forest(X.a, Y.a, W.a, Z.a, 
                                          num.threads = 30, num.trees = 10000, 
                                          sample.fraction = 0.1, min.node.size = 10, ci.group.size = 100)
  late.a <- average_treatment_effect(ivreg.res.test.a)
  
  mean.cf <- c(mean.cf, late.a[1])
  std.error.cf <- c(std.error.cf, late.a[2])
  class.cf <- c(class.cf, names(table(dat.test$subgroup))[i])
}

sample.df <- data.frame(mean = mean.cf, std.error = std.error.cf, class = class.cf)

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
    y_position_vector_tmp[which(order(y_position_vector) == i)] <- max(y_position_vector_tmp) + 0.1
  } else {
    y_position_vector_tmp[which(order(y_position_vector) == i)] <- max(y_position_vector_tmp) + 0.15
  }
}

subgroup_matching_list <- list()
group_index = 1
for (i in as.character(sort(as.numeric(names(table(sample.df$class)))))){
  subgroup_matching_list[[i]] <- letters[group_index]
  group_index = group_index + 1
}

p.validation <- ggplot(sample.df, aes(factor(class), y=mean, ymin=mean - 1.96 * std.error, ymax = mean + 1.96 * std.error))+
  geom_pointrange() + geom_signif(comparisons = comparsion_list, annotations = annotations_vector, y_position = y_position_vector_tmp) +
  theme_classic() + xlab("Subgroup") + ylab("LATE") + ggtitle("Coronary Artery Disease") + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 20),
        axis.text = element_text(size = 25), axis.line=element_line(color='black')) + scale_x_discrete(labels = subgroup_matching_list)
p.validation
