#' This code is to generate the time variable for heart disease/stroke/t2dm
#' We set the heart disease happened after LDL measurement as 1 and the heart disease happened before LDL measurement and the rest as 0.
#' 
#' @author Yujia Shi
#' @date 2022.09.01

library(dplyr)
library(data.table)

# outcome <- fread("/home/yujia/Project/2023-07-20-individual_MR/dat/03_outcome_data/First-diagnoses-before-20211112-included-GP-and-hos-and-death-tidied-at-20230915.csv.gz")
outcome.james.latest <- read.csv("/home/yujia/Project/2023-07-20-individual_MR/dat/03_outcome_data/First-diagnoses-before-20221219-included-GP-till20210930-and-hos-and-death-till20221219-tidied-at-20230915.csv")
colnames(outcome.james.latest) <- c("IID", "af", "cad", "ischemic_stroke", "copd", "heart_failure", "hemorrhage_stroke", "htn", "stroke", "t2dm", "vte")
outcome.james.latest[outcome.james.latest == ""] <- NA

#' =======================
#' Get outcome for disease
#' =======================

outcome <- read.csv("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/outcome_with_date", sep = "\t")
outcome[outcome == ""] <- NA

# outcome <- outcome[, .SD, .SDcols = !c("cad", "stroke", "copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "ischemic_stroke")] # used for data.table
outcome <- outcome[, !(colnames(outcome) %in% c("cad", "stroke", "copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "ischemic_stroke"))] 
head(outcome.james.latest)
head(outcome)

# merge the outcome with james' latest outcome file.
outcome.latest <- merge(outcome, outcome.james.latest, all.x = T, by="IID")
outcome.latest$CAD_year <- as.integer(stringr::str_extract(outcome.latest$cad, "^\\d{4}"))
outcome.latest$attend_year <- as.integer(stringr::str_extract(outcome.latest$included_date, "^\\d{4}"))
outcome.latest$age_CAD <- outcome.latest$CAD_year - outcome.latest$birth_date
outcome.latest$age_censor <- 2022 - outcome.latest$birth_date

outcome <- outcome.latest %>% 
  mutate(across(age_CAD, ~ coalesce(.x, age_censor)))

# include patients with valid biological index assay date
# outcome <- outcome[complete.cases(outcome[, 13:dim(outcome)[2]]),] # we only keep patients without missing. 
# dim(outcome)

outcome <- outcome %>%
            mutate(max_date = do.call(pmax, c(select(., starts_with('X')), na.rm = TRUE)),
                    min_date = do.call(pmin, c(select(., starts_with('X')), 
            na.rm = TRUE)))

head(outcome)

#' ===============
#' CAD def based on birth) (date1)
#' ===============

outcome.ctrl <- outcome[(is.na(outcome$cad)), ]
outcome.tret <- outcome[(!is.na(outcome$cad)), ]
outcome.tret[is.na(outcome.tret)] <- "" # change NA to "", it would be easier for us to compare date.
head(outcome.tret)

outcome.cad <- rbind(outcome.ctrl, outcome.tret)
outcome.cad <- outcome.cad[order(outcome.cad$IID), ]
row.names(outcome.cad) <- NULL

outcome.cad.res <- rep(0, dim(outcome.cad)[1])
outcome.cad.res[which(outcome.cad$cad != "")] <- 1
outcome.cad$res <- outcome.cad.res
# for (i in included.covariates){
#   outcome.cad[is.na(outcome.cad[,i]), i] <- 0
# }

colnames(outcome.cad)[dim(outcome.cad)[2]] <- "CAD"
outcome.cad$cad <- NULL
print(table(outcome.cad$CAD))
outcome.cad <- outcome.cad[, c("IID", "age_CAD", "CAD")]

data.table::fwrite(outcome.cad, paste0("/home/yujia/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date1_james2022.gz"), sep = "\t", row.names = F)

#' ===============
#' CAD def based on date attending the assessment center (53-0.0) (date2)
#' ===============

# Option 1: Only patients with CAD after the assessed date will be selected as treatment group. (based on the discussion with Prof So)
outcome.ctrl <- outcome[(is.na(outcome$cad)), ]
outcome.tret <- outcome[(!is.na(outcome$cad)), ]
outcome.tret[is.na(outcome.tret)] <- "" # change NA to "", it would be easier for us to compare date.

# CAD def based on date attending the assessment center (53-0.0) (date1)
outcome.tret <- outcome.tret[(outcome.tret$included_date < outcome.tret$cad), ]
print(dim(outcome.tret))

# modified stroke column, now to the unknown stroke type.
outcome.tret[outcome.tret$stroke == outcome.tret$ischemic_stroke, "stroke"] <- ""
outcome.tret[outcome.tret$stroke == outcome.tret$hemorrhage_stroke, "stroke"] <- ""

# get disease covariates
outcome.tret[outcome.tret == ""] <- NA
included.covariates <- c("stroke", "copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "ischemic_stroke")
for (i in included.covariates){
  stroke.index.0 <- which(outcome.tret[, i] > outcome.tret$cad)
  stroke.index.1 <- which(outcome.tret[, i] <= outcome.tret$cad)
  outcome.tret[stroke.index.0, i] <- 0
  outcome.tret[stroke.index.1, i] <- 1

  outcome.ctrl[!is.na(outcome.ctrl[, i]), i]  <- 1
}
outcome.tret[is.na(outcome.tret)] <- 0

outcome.cad <- rbind(outcome.ctrl, outcome.tret)
outcome.cad <- outcome.cad[order(outcome.cad$IID), ]
row.names(outcome.cad) <- NULL

outcome.cad.res <- rep(0, dim(outcome.cad)[1])
outcome.cad.res[which(outcome.cad$included_date < outcome.cad$cad)] <- 1 # we need CAD occur after the program attended date (53-0.0) as case.
outcome.cad$res <- outcome.cad.res
for (i in included.covariates){
  outcome.cad[is.na(outcome.cad[,i]), i] <- 0
}

colnames(outcome.cad)[dim(outcome.cad)[2]] <- "CAD"
outcome.cad$cad <- NULL
print(table(outcome.cad$CAD))
outcome.cad <- outcome.cad[, c("IID", "age_CAD", "htn", "t2dm", "heart_failure", "hemorrhage_stroke", "ischemic_stroke", "CAD")]

data.table::fwrite(outcome.cad, paste0("/home/yujia/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date2_james2022.gz"), sep = "\t", row.names = F)

#' ===============
#' CAD def based on date assessing biological indexes (date3)
#' ===============

# Option 1: Only patients with CAD after the assessed date will be selected as treatment group. (based on the discussion with Prof So)
outcome.ctrl <- outcome[(is.na(outcome$cad)), ]
outcome.tret <- outcome[(!is.na(outcome$cad)), ]
outcome.tret[is.na(outcome.tret)] <- "" # change NA to "", it would be easier for us to compare date.

# CAD def based on date of biological assay (max_date) (date2)
outcome.tret <- outcome.tret[(outcome.tret$max_date < outcome.tret$cad), ]
print(dim(outcome.tret))

# modified stroke column, now to the unknown stroke type.
outcome.tret[outcome.tret$stroke == outcome.tret$ischemic_stroke, "stroke"] <- ""
outcome.tret[outcome.tret$stroke == outcome.tret$hemorrhage_stroke, "stroke"] <- ""

# get disease covariates
outcome.tret[outcome.tret == ""] <- NA
included.covariates <- c("stroke", "copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "ischemic_stroke")
for (i in included.covariates){
  stroke.index.0 <- which(outcome.tret[, i] > outcome.tret$cad)
  stroke.index.1 <- which(outcome.tret[, i] <= outcome.tret$cad)
  outcome.tret[stroke.index.0, i] <- 0
  outcome.tret[stroke.index.1, i] <- 1

  outcome.ctrl[!is.na(outcome.ctrl[, i]), i]  <- 1
}
outcome.tret[is.na(outcome.tret)] <- 0

outcome.cad <- rbind(outcome.ctrl, outcome.tret)
outcome.cad <- outcome.cad[order(outcome.cad$IID), ]
row.names(outcome.cad) <- NULL

outcome.cad.res <- rep(0, dim(outcome.cad)[1])
outcome.cad.res[which(outcome.cad$max_date < outcome.cad$cad)] <- 1 # we need CAD occur after the biological indexes assay date (max_date) as case.
outcome.cad$res <- outcome.cad.res
for (i in included.covariates){
  outcome.cad[is.na(outcome.cad[,i]), i] <- 0
}

colnames(outcome.cad)[dim(outcome.cad)[2]] <- "CAD"
outcome.cad$cad <- NULL
print(table(outcome.cad$CAD))
outcome.cad <- outcome.cad[, c("IID", "age_CAD", "htn", "t2dm", "heart_failure", "hemorrhage_stroke", "ischemic_stroke", "CAD")]

data.table::fwrite(outcome.cad, paste0("/home/yujia/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date3_james2022.gz"), sep = "\t", row.names = F)

# #' =====================
# #' Ischemic Stroke (IST)
# #' =====================

# #' Option 1: 
# #' Only patients with CAD after the assessed date will be selected as treatment group. (based on the discussion with Prof So)
# outcome.ctrl <- outcome[(is.na(outcome$ischemic_stroke)), ]
# outcome.tret <- outcome[(!is.na(outcome$ischemic_stroke)), ]
# outcome.tret[is.na(outcome.tret)] <- "" # change NA to "", it would be easier for us to compare date.
# outcome.tret <- outcome.tret[(outcome.tret$included_date < outcome.tret$ischemic_stroke), ]

# #' Option 2: 
# #' Patients with IS after the assessed date will be selected as treatment group and we also need the CAD occur before the IS and HS occured after or at the same day as the IS.
# # outcome.tret <- outcome[!is.na(outcome$ischemic_stroke), ]
# # outcome.tret[is.na(outcome.tret)] <- "" # change NA to "", it would be easier for us to compare date.
# # outcome.tret <- outcome.tret[(outcome.tret$included_date < outcome.tret$ischemic_stroke), ] 
# # outcome.tret <- outcome.tret[(outcome.tret$ischemic_stroke < outcome.tret$hemorrhage_stroke) | (outcome.tret$ischemic_stroke != "" & outcome.tret$hemorrhage_stroke == "") | (outcome.tret$ischemic_stroke != "" & outcome.tret$hemorrhage_stroke != "" & outcome.tret$hemorrhage_stroke < outcome.tret$included_date), ] 
# # outcome.tret <- outcome.tret[(outcome.tret$ischemic_stroke < outcome.tret$cad & outcome.tret$cad != "") | (outcome.tret$cad == ""), ] 
# # outcome.tret[outcome.tret == ""] <- NA

# # modified stroke column, now to the unknown stroke type.
# outcome.tret[outcome.tret$stroke == outcome.tret$ischemic_stroke, "stroke"] <- ""
# outcome.tret[outcome.tret$stroke == outcome.tret$hemorrhage_stroke, "stroke"] <- ""

# # get disease covariates
# outcome.tret[outcome.tret == ""] <- NA
# included.covariates <- c("cad", "copd", "htn", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "stroke")
# for (i in included.covariates){
#   stroke.index.0 <- which(outcome.tret[, i] >= outcome.tret$ischemic_stroke)
#   stroke.index.1 <- which(outcome.tret[, i] < outcome.tret$ischemic_stroke)
#   outcome.tret[stroke.index.0, i] <- 0
#   outcome.tret[stroke.index.1, i] <- 1

#   outcome.ctrl[!is.na(outcome.ctrl[, i]), i]  <- 1
# }
# outcome.tret[is.na(outcome.tret)] <- 0

# outcome.is <- rbind(outcome.ctrl, outcome.tret)
# outcome.is <- outcome.is[order(outcome.is$IID), ]
# row.names(outcome.is) <- NULL

# outcome.is.res <- rep(0, dim(outcome.is)[1])
# outcome.is.res[which(outcome.is$included_date < outcome.is$ischemic_stroke)] <- 1
# outcome.is$res <- outcome.is.res
# for (i in included.covariates){
#   outcome.is[is.na(outcome.is[,i]), i] <- 0
# }

# outcome.is$ischemic_stroke <- NULL
# colnames(outcome.is)[12] <- "IST"
# table(outcome.is$IST)

# data.table::fwrite(outcome.is, paste0("/home/yujia/Project/2023-07-20-individual_MR/dat/03_outcome_data/IST_outcome_include_comorbid"), sep = "\t", row.names = F)

# #' =====================
# #' Hypertension (HTN)
# #' =====================

# #' Option 1: 
# #' Only patients with CAD after the assessed date will be selected as treatment group. (based on the discussion with Prof So)
# outcome.ctrl <- outcome[(is.na(outcome$htn)), ]
# outcome.tret <- outcome[(!is.na(outcome$htn)), ]
# outcome.tret[is.na(outcome.tret)] <- "" # change NA to "", it would be easier for us to compare date.
# outcome.tret <- outcome.tret[(outcome.tret$included_date < outcome.tret$htn), ]

# # modified stroke column, now to the unknown stroke type.
# outcome.tret[outcome.tret$stroke == outcome.tret$ischemic_stroke, "stroke"] <- ""
# outcome.tret[outcome.tret$stroke == outcome.tret$hemorrhage_stroke, "stroke"] <- ""

# # get disease covariates
# outcome.tret[outcome.tret == ""] <- NA
# included.covariates <- c("cad", "copd", "t2dm", "vte", "af", "heart_failure", "hemorrhage_stroke", "ischemic_stroke", "stroke")
# for (i in included.covariates){
#   stroke.index.0 <- which(outcome.tret[, i] >= outcome.tret$htn)
#   stroke.index.1 <- which(outcome.tret[, i] < outcome.tret$htn)
#   outcome.tret[stroke.index.0, i] <- 0
#   outcome.tret[stroke.index.1, i] <- 1

#   outcome.ctrl[!is.na(outcome.ctrl[, i]), i]  <- 1
# }
# outcome.tret[is.na(outcome.tret)] <- 0

# outcome.is <- rbind(outcome.ctrl, outcome.tret)
# outcome.is <- outcome.is[order(outcome.is$IID), ]
# row.names(outcome.is) <- NULL

# outcome.is.res <- rep(0, dim(outcome.is)[1])
# outcome.is.res[which(outcome.is$included_date < outcome.is$htn)] <- 1
# outcome.is$res <- outcome.is.res
# for (i in included.covariates){
#   outcome.is[is.na(outcome.is[,i]), i] <- 0
# }

# outcome.is$htn <- NULL
# colnames(outcome.is)[12] <- "HTN"
# table(outcome.is$HTN)

# data.table::fwrite(outcome.is, paste0("/home/yujia/Project/2023-07-20-individual_MR/dat/03_outcome_data/HTN_outcome_include_comorbid"), sep = "\t", row.names = F)
