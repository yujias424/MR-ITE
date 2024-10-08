# import module
import os
import sys

import dowhy
import econml
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import shap
import statsmodels.nonparametric.smoothers_lowess as sl
from econml.grf import CausalForest, CausalIVForest, RegressionForest
from econml.iv.dml import DMLIV, NonParamDMLIV, OrthoIV
from econml.iv.dr import DRIV, ForestDRIV, LinearDRIV, SparseLinearDRIV
from econml.sklearn_extensions.linear_model import WeightedLassoCV
from matplotlib import cm
from matplotlib.colors import Normalize
from scipy import special
from scipy.interpolate import interp1d, interpn
from scipy.stats import pearsonr
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import (LinearRegression, LogisticRegression,
                                  LogisticRegressionCV)
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from statsmodels.stats.multitest import multipletests

import plot_top_clinics as ptc
import copy
import pickle

def warn(*args, **kwargs):
    pass
import warnings

warnings.warn = warn

def save_object(obj, filename):
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)

#########################
# Data preprocessing
#########################

# load dataset and sample
X_set1 = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.traits.set1.gz", sep = "\t")
X_set2 = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.traits.set2.gz", sep = "\t")
X_set3 = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TG/CAD/ukbb.covariate.TG.set3.gz", sep = "\t")

Z_set1 = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/TG/CAD/set1/1e-08/TG_prs.best", sep = " ") # score_constd
Z_set2 = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/TG/CAD/set2/1e-08/TG_prs.best", sep = " ") # score_constd
Z_set3 = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/TG/CAD/set3/1e-08/TG_prs.best", sep = " ") # score_constd

W = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TG/CAD/ukbb.phenotype.TG.mgdL", sep = "\t")

Y_date1 = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date1_james2022.gz", sep = "\t")
Y_date2 = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date2_james2022.gz", sep = "\t")
Y_date3 = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date3_james2022.gz", sep = "\t")

# harmonize the data
selected_id_set3 = set.intersection(set(X_set3['IID']), set(Z_set3['IID']), set(W['IID']), set(Y_date3['IID']) )
selected_id_set3 = list(selected_id_set3)
selected_id_set3_arr = np.array(selected_id_set3)
selected_id_set3_arr.sort()

selected_id_set1b = set.intersection(set(X_set1['IID']), set(Z_set1['IID']), set(W['IID']), set(Y_date3['IID']) )
selected_id_set1b = list(selected_id_set1b)
selected_id_set1b_arr = np.array(selected_id_set1b)
selected_id_set1b_arr.sort()

# model 1b
X_model1b = X_set1.loc[X_set1['IID'].isin(selected_id_set1b)].reset_index(drop = True)
Z_model1b = Z_set1.loc[Z_set1['IID'].isin(selected_id_set1b)].reset_index(drop = True)
W_model1b = W.loc[W['IID'].isin(selected_id_set1b)].reset_index(drop = True)
Y_model1b = Y_date3.loc[Y_date3['IID'].isin(selected_id_set1b)].reset_index(drop = True)

# model 3
X_model3 = X_set3.loc[X_set3['IID'].isin(selected_id_set3)].reset_index(drop = True)
Z_model3 = Z_set3.loc[Z_set3['IID'].isin(selected_id_set3)].reset_index(drop = True)
W_model3 = W.loc[W['IID'].isin(selected_id_set3)].reset_index(drop = True)
Y_model3 = Y_date3.loc[Y_date3['IID'].isin(selected_id_set3)].reset_index(drop = True)
X_model3 = pd.concat([X_model3, Y_model3.loc[:, ["htn", "t2dm", "heart_failure", "hemorrhage_stroke", "ischemic_stroke"]]], axis=1)
Y_model3 = Y_model3.loc[:, ["IID", "CAD"]]

# generate mat file for three models
# model 1b
W_model1b_mat = W_model1b.loc[:, ["30870-0.0"]]["30870-0.0"]
X_model1b_mat = X_model1b.iloc[:, 2:]
Y_model1b_mat = Y_model1b.loc[:, ["CAD"]]["CAD"]
Z_model1b_mat = Z_model1b.iloc[:, 3]
# model 3
W_model3_mat = W_model3.loc[:, ["30870-0.0"]]["30870-0.0"]
X_model3_mat = X_model3.iloc[:, 2:]
Y_model3_mat = Y_model3.loc[:, ["CAD"]]["CAD"]
Z_model3_mat = Z_model3.iloc[:, 3]

# correct the data type
X_model1b_mat = X_model1b_mat.astype({
    "22001-0.0": 'int64', "21022-0.0": 'float64', "22000-0.0": 'float64',
    "22009-0.1": 'float64', "22009-0.2": 'float64', "22009-0.3": 'float64', "22009-0.4": 'float64', "22009-0.5": 'float64',
    "22009-0.6": 'float64', "22009-0.7": 'float64', "22009-0.8": 'float64', "22009-0.9": 'float64', "22009-0.10": 'float64',
})

X_model3_mat = X_model3_mat.astype({
    "22001-0.0": 'int64', "21022-0.0": 'float64', 
    "4079-0.0": 'float64', "4080-0.0": 'float64', "189-0.0": 'float64', 
    "22009-0.1": 'float64', "22009-0.2": 'float64', "22009-0.3": 'float64', "22009-0.4": 'float64', "22009-0.5": 'float64',
    "22009-0.6": 'float64', "22009-0.7": 'float64', "22009-0.8": 'float64', "22009-0.9": 'float64', "22009-0.10": 'float64', "22000-0.0": 'float64',
    "21001-0.0": 'float64', "23099-0.0": 'float64', "21002-0.0": 'float64', # "whr": 'float64', 
    "Blood_pressure_medication": 'int64', "Cholesterol_lowering_medication": 'int64', "Insulin": 'int64', # "No_medication": 'int64', 
    "Non_alcohol_drinker": 'int64' , "Previous_alcohol_drinker": 'int64', "Current_alcohol_drinker": 'int64',
    "Non_smoker": 'int64' , "Previous_smoker": 'int64', "Current_smoker": 'int64',
    "30630-0.0": "float64", "30640-0.0": "float64", "30790-0.0": "float64", "30780-0.0": "float64", "30760-0.0": "float64", "30690-0.0": "float64", # lipid-related covariates
    "30700-0.0": 'float64', "30710-0.0": 'float64', "30720-0.0": 'float64', "30730-0.0": 'float64', "30680-0.0": 'float64', 
    "30740-0.0": 'float64', "30750-0.0": 'float64', "30650-0.0": 'float64', "30660-0.0": 'float64', 
    "30670-0.0": 'float64', "30770-0.0": 'float64', "30810-0.0": 'float64', "30830-0.0": 'float64', "30850-0.0": 'float64', 
    "30880-0.0": 'float64', "30890-0.0": 'float64', "30840-0.0": 'float64', "30860-0.0": 'float64', 
    "t2dm": 'int64', "htn": 'int64', "heart_failure": 'int64', "hemorrhage_stroke": 'int64', "ischemic_stroke": 'int64'
})

# correct the data type
X_model1b_mat.rename({
    "22001-0.0": 'Gender', "21022-0.0": 'Age',
    "22009-0.1": 'PC1', "22009-0.2": 'PC2', "22009-0.3": 'PC3', "22009-0.4": 'PC4', "22009-0.5": 'PC5',
    "22009-0.6": 'PC6', "22009-0.7": 'PC7', "22009-0.8": 'PC8', "22009-0.9": 'PC9', "22009-0.10": 'PC10', "22000-0.0": 'Genotype batch',
}, inplace=True)

X_model3_mat.rename({
    "22001-0.0": 'Gender', "21022-0.0": 'Age', 
    "4079-0.0": 'diastolic blood pressure', "4080-0.0": 'systolic blood pressure', "189-0.0": 'Townsend deprivation index', 
    "22009-0.1": 'PC1', "22009-0.2": 'PC2', "22009-0.3": 'PC3', "22009-0.4": 'PC4', "22009-0.5": 'PC5',
    "22009-0.6": 'PC6', "22009-0.7": 'PC7', "22009-0.8": 'PC8', "22009-0.9": 'PC9', "22009-0.10": 'PC10', "22000-0.0": 'Genotype batch',
    "21001-0.0": 'BMI', "23099-0.0": 'Body Fat Percentage', "21002-0.0": 'Weight', # "whr": 'Waist-hip-ratio', 
    "Blood_pressure_medication": 'Blood pressure medication', "Cholesterol_lowering_medication": 'Cholesterol lowering medication', "Insulin": 'Insulin', # "No_medication": 'No medication', 
    "Non_alcohol_drinker": 'Non-alcohol drinker' , "Previous_alcohol_drinker": 'Previous alcohol drinker', "Current_alcohol_drinker": 'Current alcohol drinker',
    "Non_smoker": 'Non-smoker' , "Previous_smoker": 'Previous smoker', "Current_smoker": 'Current smoker',
    "30630-0.0": "Apolipoprotein A", "30640-0.0": "Apolipoprotein B", "30790-0.0": "Lipoprotein A", "30780-0.0": "LDL-C", "30760-0.0": "HDL-C", "30690-0.0": "Total Cholesterol", # lipid-related covariates
    "30700-0.0": 'Creatinine', "30710-0.0": 'C-reactive protein', "30720-0.0": 'Cystatin C', "30730-0.0": 'Gamma glutamyltransferase', "30680-0.0": 'Calcium', 
    "30740-0.0": 'Glucose', "30750-0.0": 'HbA1c', "30650-0.0": 'Aspartate aminotransferase', "30660-0.0": 'Direct bilirubin', 
    "30670-0.0": 'Urea', "30770-0.0": '30770-0.0', "30810-0.0": '30810-0.0', "30830-0.0": 'SHBG', "30850-0.0": 'Testosterone', 
    "30880-0.0": 'Urate', "30890-0.0": 'Vitamin D', "30840-0.0": 'Total bilirubin', "30860-0.0": 'Total protein', 
    "t2dm": 'Type 2 diabetes history', "htn": 'Hypertension history', "heart_failure": 'Heart failure history', "hemorrhage_stroke": 'Hemorrhage Stroke history', "ischemic_stroke": 'Ischemic Stroke history'
}, inplace=True)

W_model1b_mat_binary = np.where(W_model1b_mat <= 150, 1, 0)
W_model1b_mat_binary = pd.DataFrame({"30870-0.0": W_model1b_mat_binary}).iloc[:, 0]
W_model3_mat_binary = np.where(W_model3_mat <= 150, 1, 0)
W_model3_mat_binary = pd.DataFrame({"30870-0.0": W_model3_mat_binary}).iloc[:, 0]

# Z_model1_mat_binary = np.where(Z_model1_mat <= 0, 1, 0)
# Z_model1_mat_binary = pd.DataFrame({"PRS": Z_model1_mat_binary}).iloc[:, 0]
# Z_model1b_mat_binary = np.where(Z_model1b_mat <= 0, 1, 0)
# Z_model1b_mat_binary = pd.DataFrame({"PRS": Z_model1b_mat_binary}).iloc[:, 0]
# Z_model2_mat_binary = np.where(Z_model2_mat <= 0, 1, 0)
# Z_model2_mat_binary = pd.DataFrame({"PRS": Z_model2_mat_binary}).iloc[:, 0]
# Z_model2b_mat_binary = np.where(Z_model2b_mat <= 0, 1, 0)
# Z_model2b_mat_binary = pd.DataFrame({"PRS": Z_model2b_mat_binary}).iloc[:, 0]
# Z_model3_mat_binary = np.where(Z_model3_mat <= 0, 1, 0)
# Z_model3_mat_binary = pd.DataFrame({"PRS": Z_model3_mat_binary}).iloc[:, 0]

#########################
# DRIV estimator
#########################

#' =======================
#' Continuous W (Full Set)
#' =======================
# model 1b
est_driv_continuousW_model1b = ForestDRIV(projection=False, discrete_treatment=False, discrete_instrument=False,
                                         n_estimators=5000, min_samples_leaf=100, max_samples=0.02,
                                         random_state=309, n_jobs=15) 
est_driv_continuousW_model1b.fit(Y_model1b_mat, W_model1b_mat, Z=Z_model1b_mat, X=X_model1b_mat, cache_values=True)
point_driv_continuousW_model1b = est_driv_continuousW_model1b.effect(X_model1b_mat)
print(pd.DataFrame({"dat": point_driv_continuousW_model1b}).describe())

# model 3
est_driv_continuousW_model3 = ForestDRIV(projection=False, discrete_treatment=False, discrete_instrument=False, 
                                         n_estimators=5000, min_samples_leaf=100, max_samples=0.02, 
                                         random_state=309, cov_clip = 5, n_jobs=15) 
est_driv_continuousW_model3.fit(Y_model3_mat, W_model3_mat, Z=Z_model3_mat, X=X_model3_mat, cache_values=True)
point_driv_continuousW_model3 = est_driv_continuousW_model3.effect(X_model3_mat)
print(pd.DataFrame({"dat": point_driv_continuousW_model3}).describe())

#' ================================
#' Binary W Continuous Z (Full Set)
#' ================================
# model 1b
est_driv_binaryW_continuousZ_model1b = ForestDRIV(projection=False, discrete_treatment=True, discrete_instrument=False, 
                                     n_estimators=5000, min_samples_leaf=100, max_samples=0.02, 
                                     random_state=309, cov_clip = 0.1, n_jobs=15) 
est_driv_binaryW_continuousZ_model1b.fit(Y_model1b_mat, W_model1b_mat_binary, Z=Z_model1b_mat, X=X_model1b_mat, cache_values=True)
point_driv_binaryW_continuousZ_model1b = est_driv_binaryW_continuousZ_model1b.effect(X_model1b_mat)
print(pd.DataFrame({"dat": point_driv_binaryW_continuousZ_model1b}).describe())

# model 3
est_driv_binaryW_continuousZ_model3 = ForestDRIV(projection=False, discrete_treatment=True, discrete_instrument=False, \
                                                n_estimators=5000, min_samples_leaf=100, max_samples=0.02, 
                                                random_state=309, cov_clip = 0.1, n_jobs=15) 
est_driv_binaryW_continuousZ_model3.fit(Y_model3_mat, W_model3_mat_binary, Z=Z_model3_mat, X=X_model3_mat, cache_values=True)
point_driv_binaryW_continuousZ_model3 = est_driv_binaryW_continuousZ_model3.effect(X_model3_mat)
print(pd.DataFrame({"dat": point_driv_binaryW_continuousZ_model3}).describe())

# calculate the upper bound and lower bound (Continuous W)
# full set
# model 1b
point_driv_lb_continuousW_model1b, point_driv_ub_continuousW_model1b = est_driv_continuousW_model1b.effect_interval(X_model1b_mat, alpha=0.1) # type: ignore
z_value_driv_continuousW_model1b = point_driv_continuousW_model1b/((point_driv_ub_continuousW_model1b-point_driv_lb_continuousW_model1b)/(2*1.645))
p_value_driv_continuousW_model1b = scipy.stats.norm.sf(abs(z_value_driv_continuousW_model1b)) # * 2
p_value_BH_driv_continuousW_model1b = multipletests(pvals = p_value_driv_continuousW_model1b, method = "fdr_bh", alpha=0.1)

full_results_continuousW_model1b = pd.DataFrame({"IID": selected_id_set1b_arr, "point": point_driv_continuousW_model1b, "upper_bound": point_driv_ub_continuousW_model1b, \
                                "lower_bound": point_driv_lb_continuousW_model1b, "z-value": z_value_driv_continuousW_model1b, \
                                "p_value": p_value_driv_continuousW_model1b, "p_value_corrected": p_value_BH_driv_continuousW_model1b[1]})

# model 3
point_driv_lb_continuousW_model3, point_driv_ub_continuousW_model3 = est_driv_continuousW_model3.effect_interval(X_model3_mat, alpha=0.1) # type: ignore
z_value_driv_continuousW_model3 = point_driv_continuousW_model3/((point_driv_ub_continuousW_model3-point_driv_lb_continuousW_model3)/(2*1.645))
p_value_driv_continuousW_model3 = scipy.stats.norm.sf(abs(z_value_driv_continuousW_model3)) # * 2
p_value_BH_driv_continuousW_model3 = multipletests(pvals = p_value_driv_continuousW_model3, method = "fdr_bh", alpha=0.1)

full_results_continuousW_model3 = pd.DataFrame({"IID": selected_id_set3_arr, "point": point_driv_continuousW_model3, "upper_bound": point_driv_ub_continuousW_model3, \
                                "lower_bound": point_driv_lb_continuousW_model3, "z-value": z_value_driv_continuousW_model3, \
                                "p_value": p_value_driv_continuousW_model3, "p_value_corrected": p_value_BH_driv_continuousW_model3[1]})

# calculate the upper bound and lower bound (Binary W, Continuous Z)
# full set
# model 1b
point_driv_lb_binaryW_continuousZ_model1b, point_driv_ub_binaryW_continuousZ_model1b = est_driv_binaryW_continuousZ_model1b.effect_interval(X_model1b_mat, alpha=0.1) # type: ignore
z_value_driv_binaryW_continuousZ_model1b = point_driv_binaryW_continuousZ_model1b/((point_driv_ub_binaryW_continuousZ_model1b-point_driv_lb_binaryW_continuousZ_model1b)/(2*1.645))
p_value_driv_binaryW_continuousZ_model1b = scipy.stats.norm.sf(abs(z_value_driv_binaryW_continuousZ_model1b)) # * 2
p_value_BH_driv_binaryW_continuousZ_model1b = multipletests(pvals = p_value_driv_binaryW_continuousZ_model1b, method = "fdr_bh", alpha=0.1)

full_results_binaryW_continuousZ_model1b = pd.DataFrame({"IID": selected_id_set1b_arr, "point": point_driv_binaryW_continuousZ_model1b, "upper_bound": point_driv_ub_binaryW_continuousZ_model1b, \
                                "lower_bound": point_driv_lb_binaryW_continuousZ_model1b, "z-value": z_value_driv_binaryW_continuousZ_model1b, \
                                "p_value": p_value_driv_binaryW_continuousZ_model1b, "p_value_corrected": p_value_BH_driv_binaryW_continuousZ_model1b[1]})

# model 3
point_driv_lb_binaryW_continuousZ_model3, point_driv_ub_binaryW_continuousZ_model3 = est_driv_binaryW_continuousZ_model3.effect_interval(X_model3_mat, alpha=0.1) # type: ignore
z_value_driv_binaryW_continuousZ_model3 = point_driv_binaryW_continuousZ_model3/((point_driv_ub_binaryW_continuousZ_model3-point_driv_lb_binaryW_continuousZ_model3)/(2*1.645))
p_value_driv_binaryW_continuousZ_model3 = scipy.stats.norm.sf(abs(z_value_driv_binaryW_continuousZ_model3)) # * 2
p_value_BH_driv_binaryW_continuousZ_model3 = multipletests(pvals = p_value_driv_binaryW_continuousZ_model3, method = "fdr_bh", alpha=0.1)

full_results_binaryW_continuousZ_model3 = pd.DataFrame({"IID": selected_id_set3_arr, "point": point_driv_binaryW_continuousZ_model3, "upper_bound": point_driv_ub_binaryW_continuousZ_model3, \
                                "lower_bound": point_driv_lb_binaryW_continuousZ_model3, "z-value": z_value_driv_binaryW_continuousZ_model3, \
                                "p_value": p_value_driv_binaryW_continuousZ_model3, "p_value_corrected": p_value_BH_driv_binaryW_continuousZ_model3[1]})

print("Summary of lb rb results: ")
print("Continuous W: ")
print("Model 1b:")
print(full_results_continuousW_model1b.describe())
print("Model 3:")
print(full_results_continuousW_model3.describe())
print("                              ")
print("Binary W Continuou Z: ")
print("Model 1b:")
print(full_results_binaryW_continuousZ_model1b.describe())
print("Model 3:")
print(full_results_binaryW_continuousZ_model3.describe())
print("                              ")

full_results_continuousW_model1b.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/01_individual_treatment_effect/driv_full_continuousW_te_ul_model1b.csv", sep=",", index=False)
full_results_binaryW_continuousZ_model1b.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_model1b.csv", sep=",", index=False)

full_results_continuousW_model3.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/01_individual_treatment_effect/driv_full_continuousW_te_ul_model3.csv", sep=",", index=False)
full_results_binaryW_continuousZ_model3.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_model3.csv", sep=",", index=False)

#########################
# Variabel Importance
#########################

# # save important variables
# # model 1
# print("Start calculating SHAP value.")
# shap_values_driv_binaryW_model1 = est_driv_binaryW_model1.shap_values(X_model1_mat)
# shap_pd_binaryW_model1 = pd.DataFrame(shap_values_driv_binaryW_model1['CAD']['30870-0.0'].values)
# shap_pd_binaryW_model1.columns = shap_values_driv_binaryW_model1['CAD']['30870-0.0'].feature_names
# value_pd_binaryW_model1 = pd.DataFrame(shap_values_driv_binaryW_model1['CAD']['30870-0.0'].data)
# value_pd_binaryW_model1.columns = shap_values_driv_binaryW_model1['CAD']['30870-0.0'].feature_names
# shap_pd_binaryW_model1.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_binaryW_shap_model1.csv.gz", index=False)
# value_pd_binaryW_model1.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_binaryW_value_model1.csv.gz", index=False)
# print("Finished calculating SHAP value on binary W model1.")
# shap_values_driv_continuousW_model1 = est_driv_continuousW_model1.shap_values(X_model1_mat)
# shap_pd_continuousW_model1 = pd.DataFrame(shap_values_driv_continuousW_model1['CAD']['30870-0.0'].values)
# shap_pd_continuousW_model1.columns = shap_values_driv_continuousW_model1['CAD']['30870-0.0'].feature_names
# value_pd_continuousW_model1 = pd.DataFrame(shap_values_driv_continuousW_model1['CAD']['30870-0.0'].data)
# value_pd_continuousW_model1.columns = shap_values_driv_continuousW_model1['CAD']['30870-0.0'].feature_names
# shap_pd_continuousW_model1.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_continuousW_shap_model1.csv.gz", index=False)
# value_pd_continuousW_model1.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_continuousW_value_model1.csv.gz", index=False)
# print("Finished calculating SHAP value on continuous W model1.")

# # model 1b
# print("Start calculating SHAP value.")
# shap_values_driv_binaryW_model1b = est_driv_binaryW_model1b.shap_values(X_model1b_mat)
# shap_pd_binaryW_model1b = pd.DataFrame(shap_values_driv_binaryW_model1b['CAD']['30870-0.0'].values)
# shap_pd_binaryW_model1b.columns = shap_values_driv_binaryW_model1b['CAD']['30870-0.0'].feature_names
# value_pd_binaryW_model1b = pd.DataFrame(shap_values_driv_binaryW_model1b['CAD']['30870-0.0'].data)
# value_pd_binaryW_model1b.columns = shap_values_driv_binaryW_model1b['CAD']['30870-0.0'].feature_names
# shap_pd_binaryW_model1b.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_binaryW_shap_model1b.csv.gz", index=False)
# value_pd_binaryW_model1b.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_binaryW_value_model1b.csv.gz", index=False)
# print("Finished calculating SHAP value on binary W model1b.")
# shap_values_driv_continuousW_model1b = est_driv_continuousW_model1b.shap_values(X_model1b_mat)
# shap_pd_continuousW_model1b = pd.DataFrame(shap_values_driv_continuousW_model1b['CAD']['30870-0.0'].values)
# shap_pd_continuousW_model1b.columns = shap_values_driv_continuousW_model1b['CAD']['30870-0.0'].feature_names
# value_pd_continuousW_model1b = pd.DataFrame(shap_values_driv_continuousW_model1b['CAD']['30870-0.0'].data)
# value_pd_continuousW_model1b.columns = shap_values_driv_continuousW_model1b['CAD']['30870-0.0'].feature_names
# shap_pd_continuousW_model1b.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_continuousW_shap_model1b.csv.gz", index=False)
# value_pd_continuousW_model1b.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_continuousW_value_model1b.csv.gz", index=False)
# print("Finished calculating SHAP value on continuous W model1b.")

# # model 2
# print("Start calculating SHAP value.")
# shap_values_driv_binaryW_model2 = est_driv_binaryW_model2.shap_values(X_model2_mat)
# shap_pd_binaryW_model2 = pd.DataFrame(shap_values_driv_binaryW_model2['CAD']['30870-0.0'].values)
# shap_pd_binaryW_model2.columns = shap_values_driv_binaryW_model2['CAD']['30870-0.0'].feature_names
# value_pd_binaryW_model2 = pd.DataFrame(shap_values_driv_binaryW_model2['CAD']['30870-0.0'].data)
# value_pd_binaryW_model2.columns = shap_values_driv_binaryW_model2['CAD']['30870-0.0'].feature_names
# shap_pd_binaryW_model2.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_binaryW_shap_model2.csv.gz", index=False)
# value_pd_binaryW_model2.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_binaryW_value_model2.csv.gz", index=False)
# print("Finished calculating SHAP value on binary W model2.")
# shap_values_driv_continuousW_model2 = est_driv_continuousW_model2.shap_values(X_model2_mat)
# shap_pd_continuousW_model2 = pd.DataFrame(shap_values_driv_continuousW_model2['CAD']['30870-0.0'].values)
# shap_pd_continuousW_model2.columns = shap_values_driv_continuousW_model2['CAD']['30870-0.0'].feature_names
# value_pd_continuousW_model2 = pd.DataFrame(shap_values_driv_continuousW_model2['CAD']['30870-0.0'].data)
# value_pd_continuousW_model2.columns = shap_values_driv_continuousW_model2['CAD']['30870-0.0'].feature_names
# shap_pd_continuousW_model2.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_continuousW_shap_model2.csv.gz", index=False)
# value_pd_continuousW_model2.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_continuousW_value_model2.csv.gz", index=False)
# print("Finished calculating SHAP value on continuous W model2.")

# # model 2b
# print("Start calculating SHAP value.")
# shap_values_driv_binaryW_model2b = est_driv_binaryW_model2b.shap_values(X_model2b_mat)
# shap_pd_binaryW_model2b = pd.DataFrame(shap_values_driv_binaryW_model2b['CAD']['30870-0.0'].values)
# shap_pd_binaryW_model2b.columns = shap_values_driv_binaryW_model2b['CAD']['30870-0.0'].feature_names
# value_pd_binaryW_model2b = pd.DataFrame(shap_values_driv_binaryW_model2b['CAD']['30870-0.0'].data)
# value_pd_binaryW_model2b.columns = shap_values_driv_binaryW_model2b['CAD']['30870-0.0'].feature_names
# shap_pd_binaryW_model2b.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_binaryW_shap_model2b.csv.gz", index=False)
# value_pd_binaryW_model2b.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_binaryW_value_model2b.csv.gz", index=False)
# print("Finished calculating SHAP value on binary W model2b.")
# shap_values_driv_continuousW_model2b = est_driv_continuousW_model2b.shap_values(X_model2b_mat)
# shap_pd_continuousW_model2b = pd.DataFrame(shap_values_driv_continuousW_model2b['CAD']['30870-0.0'].values)
# shap_pd_continuousW_model2b.columns = shap_values_driv_continuousW_model2b['CAD']['30870-0.0'].feature_names
# value_pd_continuousW_model2b = pd.DataFrame(shap_values_driv_continuousW_model2b['CAD']['30870-0.0'].data)
# value_pd_continuousW_model2b.columns = shap_values_driv_continuousW_model2b['CAD']['30870-0.0'].feature_names
# shap_pd_continuousW_model2b.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_continuousW_shap_model2b.csv.gz", index=False)
# value_pd_continuousW_model2b.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_continuousW_value_model2b.csv.gz", index=False)
# print("Finished calculating SHAP value on continuous W model2b.")

# model 3
# print("Start calculating SHAP value.")
# shap_values_driv_binaryW_continuousZ_model3 = est_driv_binaryW_continuousZ_model3.shap_values(X_model3_mat)
# save_object(shap_values_driv_binaryW_continuousZ_model3, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_binaryW_continuousZ_shap_model3.pkl")
# print("Finished calculating SHAP value on binary W continuous Z model3.")
# shap_values_driv_continuousW_model3 = est_driv_continuousW_model3.shap_values(X_model3_mat)
# save_object(shap_values_driv_continuousW_model3, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_continuousW_shap_model3.pkl")
# print("Finished calculating SHAP value on continuous W model3.")

# shap_values_driv_binaryW_binaryZ_model3 = est_driv_binaryW_binaryZ_model3.shap_values(X_model3_mat)
# save_object(shap_values_driv_binaryW_binaryZ_model3, "/mnt/md0/yujia/project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TG/CAD/03_variable_importance/driv_binaryW_binaryZ_shap_model3.pkl")
# print("Finished calculating SHAP value on binary W binary Z model3.")

print("Finished the pipeline.")