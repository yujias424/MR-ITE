# Confidence interval

library(data.table)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(stringr)

# Set up a UKBB dictionary
UKBB_id <- c("whr",
             "31-0.0",
             "48-0.0",
             "1249-0.0",
             "1249-0.1",
             "1249-0.2",
             "1249-0.3",
             "1558-0.0",
             "1558-0.1",
             "1558-0.2",
             "1558-0.3",
             "1558-0.4",
             "1558-0.5",
             "1558-0.6",
             "2443-0.0",
             "4079-0.0",
             "4080-0.0",
             "20116-0.0",
             "20117-0.0",
             "20117-0.1",
             "20117-0.3",
             "21000-0.0",
             "21000-0.1",
             "21000-0.3",
             "21000-0.4",
             "21001-0.0",
             "21022-0.0",
             "22009-0.1",
             "22009-0.2",
             "22009-0.3",
             "22009-0.4",
             "22009-0.5",
             "22009-0.6",
             "22009-0.7",
             "22009-0.8",
             "22009-0.9",
             "22009-0.10",
             "22127-0.0",
             "23099-0.0",
             "30690-0.0",
             "30710-0.0",
             "30740-0.0",
             "30750-0.0",
             "30760-0.0",
             "30780-0.0",
             "30870-0.0",
             "22127-0.0",
             "22130-0.0")

UKBB_name <- c("Waist–hip ratio",
               "Sex", 
               "Waist circumference",
               "Past tobacco smoking level 0",
               "Past tobacco smoking level 1",
               "Past tobacco smoking level 2",
               "Past tobacco smoking level 3",
               "Alcohol intake frequency level 0",
               "Alcohol intake frequency level 1",
               "Alcohol intake frequency level 2",
               "Alcohol intake frequency level 3",
               "Alcohol intake frequency level 4",
               "Alcohol intake frequency level 5",
               "Alcohol intake frequency level 6",
               "Diabetes diagnosed by doctor",
               "Diastolic blood pressure",
               "Systolic blood pressure",
               "Smoking status",
               "Alcohol drinker status level 0",
               "Alcohol drinker status level 1",
               "Alcohol drinker status level 3",
               "Ethnic background type A",
               "Ethnic background type B",
               "Ethnic background type D",
               "Ethnic background type E",
               "Body mass index (BMI)",
               "Age at recruitment",
               "Genetic principal components (PC1)",
               "Genetic principal components (PC2)",
               "Genetic principal components (PC3)",
               "Genetic principal components (PC4)",
               "Genetic principal components (PC5)",
               "Genetic principal components (PC6)",
               "Genetic principal components (PC7)",
               "Genetic principal components (PC8)",
               "Genetic principal components (PC9)",
               "Genetic principal components (PC10)",
               "Doctor diagnosed asthma",
               "Body fat percentage",
               "Cholesterol",
               "C-reactive protein",
               "Glucose",
               "Glycated haemoglobin (HbA1c)",
               "HDL cholesterol",
               "LDL direct",
               "Triglycerides",
               "Doctor diagnosed asthma",
               "Doctor diagnosed COPD (chronic obstructive pulmonary disease)")

subfolder_name <- c("death_AUG07",
                   "death_AUG20",
                   "death_diabetes_related_AUG09",
                   "death_smoking",
                   "death_statin",
                   "severity_AUG20",
                   "severity_alcohol",
                   "severity_diabetes",
                   "severity_statin",
                   "severity_smoking",
                   "severity_diabetes_related_AUG08",
                   "severity_AUG05",
                   "severity_diabetes_related",
                   "death_diabetes")

subfolder_mean <- c(
  "Outcome: mortality, Covariates: cholesterol+, glucose+, HAb1c-",
  "Outcome: mortality, Covariates: cholesterol-, glucose-, HAb1c+",
  "Outcome: mortality, Covariates: cholesterol-, glucose-",
  "Outcome: mortality, Covariates: smoking+",
  "Outcome: mortality, Covariates: cholesterol+, glucose+, tx: statin",
  "Outcome: severity, Covariates: cholesterol-, glucose-",
  "Outcome: severity, Covariates: alcohol+",
  "Outcome: severity, Covariates: HAb1c+, tx: diabetes_diagnosed",
  "Outcome: severity, Covariates: smoking+",
  "Outcome: severity, Covariates: cholesterol+, glucose+, tx: statin",
  "Outcome: severity, Covariates: cholesterol-, glucose-",
  "Outcome: severity, Covariates: cholesterol+, glucose+, HAb1c-",
  "Outcome: severity, Covariates: cholesterol-, glucose-",
  "Outcome: mortality, Covariates: HAb1c+, tx: diabetes_diagnosed"
)

# Look Up Dictionary
UKBB_dict <- as.list(UKBB_name)
names(UKBB_dict) <- UKBB_id

subfolder_dict <- as.list(subfolder_mean)
names(subfolder_dict) <- subfolder_name

# File name manipulation
file_list <- list.files(path = "/home/kai/data/UKBB/covid19_clinical/",
                        pattern = "_tau_stats.csv$", recursive = TRUE)

# Exclude the backup folder
file_list <- file_list[!grepl("backups", file_list)]

# Need to replace / with \ first.
file_list_replace <- c()
for (i in file_list) {
  file_list_replace <- c(gsub("/", "\\", i, fixed = TRUE), file_list_replace)
}

# Build a large row data and make the plot
count_i <- 0
for (i in file_list_replace) {
  # Get different labels
  file_name <- str_match(i, "[^\\\\]+$")[1]
  factor <- str_match(file_name, ".*(?=_tau_stats.csv)")[1]
  organ <- str_match(i, ".{4,6}(?=\\\\)")[1]
  subfolder <- str_match(i, "(?<=\\\\).*(?=\\\\)")[1]

  # Get the file path.
  file_path <- paste0("/home/kai/data/UKBB/covid19_clinical/",
                      gsub("\\", "/", i, fixed = TRUE))

  factor_t <- factor
  if (factor %in% UKBB_id) {
    factor <- UKBB_dict[[factor]]
  }

  tau_df <- fread(file_path) %>% as.data.frame()

  if (length(unique(tau_df$treatment)) > 2) {
    next
  } else {
    # Data manipulation
    tau_df$treatment <- as.factor(tau_df$treatment) # factorized the treatment
    tau_df$group <- rep(paste0(subfolder, "-", factor),
                        length(tau_df$treatment))
    tau_df_t <- tau_df
    tau_df <- cbind.data.frame(group = tau_df$group,
                               treatment = tau_df$treatment,
                               tau = tau_df[, 4],
                               tau_se = tau_df[, 4]/tau_df$tau.zval)

    # Calculate the 95% confidence interval
    tau_df$min95 <- tau_df$tau - 1.96 * tau_df$tau_se
    tau_df$max95 <- tau_df$tau + 1.96 * tau_df$tau_se

    tau_df <- tau_df[order(tau_df$treatment, tau_df$min95), ]
    rownames(tau_df) <- NULL

    # Plot the errorbar plot
    q <- ggplot(tau_df, aes(x = as.numeric(rownames(tau_df)), y = tau)) +
      geom_pointrange(aes(ymin = min95,
                          ymax = max95,
                          colour = treatment), fatten = .5) +
      geom_point(aes(x = as.numeric(rownames(tau_df)),
                     y = tau,
                     shape = treatment)) +
      theme_classic() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(y = expression(paste(hat(tau)))) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(title = paste0("Organ: ", organ, "\n",
                          subfolder_dict[[subfolder]], "\n",
                          "Factor: ", factor)) +
      theme(plot.title = element_text(hjust = 0.5))

    save_path <- paste0("/home/yujia/Project/xgboost_HTE",
                        "/Tau_analysis/plots/patient_conf_interval/",
                        organ, "_", subfolder, "_", factor,
                        "_conf_interval.pdf")

    # Used for testing the plot
    # save_path <- paste0("/home/yujia/Project/xgboost_HTE",
    #                     "/Tau_analysis/",
    #                     organ, "_", subfolder, "_", factor,
    #                     "_conf_interval.pdf")

    ggsave(save_path, device = "pdf",
            plot = q, width = 20, height = 8.5, units = "in")
  }
  # break
}