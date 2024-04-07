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
X_set1 = pd.read_csv("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.traits.set1.gz", sep = "\t")
X_set2 = pd.read_csv("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/traits/ukbb.covariate.traits.set2.gz", sep = "\t")
X_set3 = pd.read_csv("/home/yujia/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TC/CAD/ukbb.covariate.TC.set3.gz", sep = "\t")

Z_set1 = pd.read_csv("/home/yujia/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/TC/CAD/set1/1e-08/TC_prs.best", sep = " ") # score_constd
Z_set2 = pd.read_csv("/home/yujia/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/TC/CAD/set2/1e-08/TC_prs.best", sep = " ") # score_constd
Z_set3 = pd.read_csv("/home/yujia/Project/2023-07-20-individual_MR/dat/06_PRS_calculation/TC/CAD/set3/1e-08/TC_prs.best", sep = " ") # score_constd

W = pd.read_csv("~/Project/2023-07-20-individual_MR/dat/04_pheno_covar_data/TC/CAD/ukbb.phenotype.TC.mgdL", sep = "\t")

Y_date1 = pd.read_csv("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date1_james2022.gz", sep = "\t")
Y_date2 = pd.read_csv("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date2_james2022.gz", sep = "\t")
Y_date3 = pd.read_csv("~/Project/2023-07-20-individual_MR/dat/03_outcome_data/CAD_outcome_include_comorbid_date3_james2022.gz", sep = "\t")

# harmonize the data
selected_id_set1 = set.intersection(set(X_set1['IID']), set(Z_set1['IID']), set(W['IID']), set(Y_date1['IID']) )
selected_id_set2 = set.intersection(set(X_set2['IID']), set(Z_set2['IID']), set(W['IID']), set(Y_date2['IID']) )
selected_id_set3 = set.intersection(set(X_set3['IID']), set(Z_set3['IID']), set(W['IID']), set(Y_date3['IID']) )
selected_id_set1 = list(selected_id_set1)
selected_id_set2 = list(selected_id_set2)
selected_id_set3 = list(selected_id_set3)
selected_id_set1_arr = np.array(selected_id_set1)
selected_id_set2_arr = np.array(selected_id_set2)
selected_id_set3_arr = np.array(selected_id_set3)
selected_id_set1_arr.sort()
selected_id_set2_arr.sort()
selected_id_set3_arr.sort()

selected_id_set1b = set.intersection(set(X_set1['IID']), set(Z_set1['IID']), set(W['IID']), set(Y_date3['IID']) )
selected_id_set2b = set.intersection(set(X_set2['IID']), set(Z_set2['IID']), set(W['IID']), set(Y_date3['IID']) )
selected_id_set1b = list(selected_id_set1b)
selected_id_set2b = list(selected_id_set2b)
selected_id_set1b_arr = np.array(selected_id_set1b)
selected_id_set2b_arr = np.array(selected_id_set2b)
selected_id_set1b_arr.sort()
selected_id_set2b_arr.sort()

# model 1
X_model1 = X_set1.loc[X_set1['IID'].isin(selected_id_set1)].reset_index(drop = True)
Z_model1 = Z_set1.loc[Z_set1['IID'].isin(selected_id_set1)].reset_index(drop = True)
W_model1 = W.loc[W['IID'].isin(selected_id_set1)].reset_index(drop = True)
Y_model1 = Y_date1.loc[Y_date1['IID'].isin(selected_id_set1)].reset_index(drop = True)

# model 1b
X_model1b = X_set1.loc[X_set1['IID'].isin(selected_id_set1b)].reset_index(drop = True)
Z_model1b = Z_set1.loc[Z_set1['IID'].isin(selected_id_set1b)].reset_index(drop = True)
W_model1b = W.loc[W['IID'].isin(selected_id_set1b)].reset_index(drop = True)
Y_model1b = Y_date3.loc[Y_date3['IID'].isin(selected_id_set1b)].reset_index(drop = True)

# model 2
X_model2 = X_set2.loc[X_set2['IID'].isin(selected_id_set2)].reset_index(drop = True)
Z_model2 = Z_set2.loc[Z_set2['IID'].isin(selected_id_set2)].reset_index(drop = True)
W_model2 = W.loc[W['IID'].isin(selected_id_set2)].reset_index(drop = True)
Y_model2 = Y_date2.loc[Y_date2['IID'].isin(selected_id_set2)].reset_index(drop = True)
X_model2 = pd.concat([X_model2, Y_model2.loc[:, ["htn", "t2dm", "heart_failure", "hemorrhage_stroke", "ischemic_stroke"]]], axis=1)
Y_model2 = Y_model2.loc[:, ["IID", "CAD"]]

# model 2b
X_model2b = X_set2.loc[X_set2['IID'].isin(selected_id_set2b)].reset_index(drop = True)
Z_model2b = Z_set2.loc[Z_set2['IID'].isin(selected_id_set2b)].reset_index(drop = True)
W_model2b = W.loc[W['IID'].isin(selected_id_set2b)].reset_index(drop = True)
Y_model2b = Y_date3.loc[Y_date3['IID'].isin(selected_id_set2b)].reset_index(drop = True)
X_model2b = pd.concat([X_model2b, Y_model2b.loc[:, ["htn", "t2dm", "heart_failure", "hemorrhage_stroke", "ischemic_stroke"]]], axis=1)
Y_model2b = Y_model2b.loc[:, ["IID", "CAD"]]

# model 3
X_model3 = X_set3.loc[X_set3['IID'].isin(selected_id_set3)].reset_index(drop = True)
Z_model3 = Z_set3.loc[Z_set3['IID'].isin(selected_id_set3)].reset_index(drop = True)
W_model3 = W.loc[W['IID'].isin(selected_id_set3)].reset_index(drop = True)
Y_model3 = Y_date3.loc[Y_date3['IID'].isin(selected_id_set3)].reset_index(drop = True)
X_model3 = pd.concat([X_model3, Y_model3.loc[:, ["htn", "t2dm", "heart_failure", "hemorrhage_stroke", "ischemic_stroke"]]], axis=1)
Y_model3 = Y_model3.loc[:, ["IID", "CAD"]]

# generate mat file for three models
# model 1
W_model1_mat = W_model1.loc[:, ["30690-0.0"]]["30690-0.0"]
X_model1_mat = X_model1.iloc[:, 2:]
Y_model1_mat = Y_model1.loc[:, ["CAD"]]["CAD"]
Z_model1_mat = Z_model1.iloc[:, 3]

# model 1b
W_model1b_mat = W_model1b.loc[:, ["30690-0.0"]]["30690-0.0"]
X_model1b_mat = X_model1b.iloc[:, 2:]
Y_model1b_mat = Y_model1b.loc[:, ["CAD"]]["CAD"]
Z_model1b_mat = Z_model1b.iloc[:, 3]

# model 2
W_model2_mat = W_model2.loc[:, ["30690-0.0"]]["30690-0.0"]
X_model2_mat = X_model2.iloc[:, 2:]
Y_model2_mat = Y_model2.loc[:, ["CAD"]]["CAD"]
Z_model2_mat = Z_model2.iloc[:, 3]

# model 2b
W_model2b_mat = W_model2b.loc[:, ["30690-0.0"]]["30690-0.0"]
X_model2b_mat = X_model2b.iloc[:, 2:]
Y_model2b_mat = Y_model2b.loc[:, ["CAD"]]["CAD"]
Z_model2b_mat = Z_model2b.iloc[:, 3]

# model 3
W_model3_mat = W_model3.loc[:, ["30690-0.0"]]["30690-0.0"]
X_model3_mat = X_model3.iloc[:, 2:]
Y_model3_mat = Y_model3.loc[:, ["CAD"]]["CAD"]
Z_model3_mat = Z_model3.iloc[:, 3]

# correct the data type
X_model1_mat = X_model1_mat.astype({
    "22001-0.0": 'int64', "21022-0.0": 'float64', "22000-0.0": 'float64',
    "22009-0.1": 'float64', "22009-0.2": 'float64', "22009-0.3": 'float64', "22009-0.4": 'float64', "22009-0.5": 'float64',
    "22009-0.6": 'float64', "22009-0.7": 'float64', "22009-0.8": 'float64', "22009-0.9": 'float64', "22009-0.10": 'float64',
})

X_model1b_mat = X_model1b_mat.astype({
    "22001-0.0": 'int64', "21022-0.0": 'float64', "22000-0.0": 'float64',
    "22009-0.1": 'float64', "22009-0.2": 'float64', "22009-0.3": 'float64', "22009-0.4": 'float64', "22009-0.5": 'float64',
    "22009-0.6": 'float64', "22009-0.7": 'float64', "22009-0.8": 'float64', "22009-0.9": 'float64', "22009-0.10": 'float64',
})

X_model2_mat = X_model2_mat.astype({
    "22001-0.0": 'int64', "21022-0.0": 'float64', 
    "4079-0.0": 'float64', "4080-0.0": 'float64', "189-0.0": 'float64', 
    "22009-0.1": 'float64', "22009-0.2": 'float64', "22009-0.3": 'float64', "22009-0.4": 'float64', "22009-0.5": 'float64',
    "22009-0.6": 'float64', "22009-0.7": 'float64', "22009-0.8": 'float64', "22009-0.9": 'float64', "22009-0.10": 'float64', "22000-0.0": 'float64',
    "21001-0.0": 'float64', "23099-0.0": 'float64', "21002-0.0": 'float64', # "whr": 'float64', 
    "Blood_pressure_medication": 'int64', "Cholesterol_lowering_medication": 'int64', "Insulin": 'int64', # "No_medication": 'int64',
    "Non_alcohol_drinker": 'int64' , "Previous_alcohol_drinker": 'int64', "Current_alcohol_drinker": 'int64',
    "Non_smoker": 'int64' , "Previous_smoker": 'int64', "Current_smoker": 'int64'
})

X_model2b_mat = X_model2b_mat.astype({
    "22001-0.0": 'int64', "21022-0.0": 'float64', 
    "4079-0.0": 'float64', "4080-0.0": 'float64', "189-0.0": 'float64', 
    "22009-0.1": 'float64', "22009-0.2": 'float64', "22009-0.3": 'float64', "22009-0.4": 'float64', "22009-0.5": 'float64',
    "22009-0.6": 'float64', "22009-0.7": 'float64', "22009-0.8": 'float64', "22009-0.9": 'float64', "22009-0.10": 'float64', "22000-0.0": 'float64',
    "21001-0.0": 'float64', "23099-0.0": 'float64', "21002-0.0": 'float64', # "whr": 'float64', 
    "Blood_pressure_medication": 'int64', "Cholesterol_lowering_medication": 'int64', # "No_medication": 'int64', "Insulin": 'int64', 
    "Non_alcohol_drinker": 'int64' , "Previous_alcohol_drinker": 'int64', "Current_alcohol_drinker": 'int64',
    "Non_smoker": 'int64' , "Previous_smoker": 'int64', "Current_smoker": 'int64'
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
    "30870-0.0": 'float64', # lipid-related covariates
    "30700-0.0": 'float64', "30710-0.0": 'float64', "30720-0.0": 'float64', "30730-0.0": 'float64', "30680-0.0": 'float64', 
    "30740-0.0": 'float64', "30750-0.0": 'float64', "30650-0.0": 'float64', "30660-0.0": 'float64', 
    "30670-0.0": 'float64', "30770-0.0": 'float64', "30810-0.0": 'float64', "30830-0.0": 'float64', "30850-0.0": 'float64', 
    "30880-0.0": 'float64', "30890-0.0": 'float64', "30840-0.0": 'float64', "30860-0.0": 'float64', 
    "t2dm": 'int64', "htn": 'int64', "heart_failure": 'int64', "hemorrhage_stroke": 'int64', "ischemic_stroke": 'int64'
})

# correct the data type
X_model1_mat.rename({
    "22001-0.0": 'Gender', "21022-0.0": 'Age',
    "22009-0.1": 'PC1', "22009-0.2": 'PC2', "22009-0.3": 'PC3', "22009-0.4": 'PC4', "22009-0.5": 'PC5',
    "22009-0.6": 'PC6', "22009-0.7": 'PC7', "22009-0.8": 'PC8', "22009-0.9": 'PC9', "22009-0.10": 'PC10', "22000-0.0": 'Genotype batch',
}, inplace=True)

X_model1b_mat.rename({
    "22001-0.0": 'Gender', "21022-0.0": 'Age',
    "22009-0.1": 'PC1', "22009-0.2": 'PC2', "22009-0.3": 'PC3', "22009-0.4": 'PC4', "22009-0.5": 'PC5',
    "22009-0.6": 'PC6', "22009-0.7": 'PC7', "22009-0.8": 'PC8', "22009-0.9": 'PC9', "22009-0.10": 'PC10', "22000-0.0": 'Genotype batch',
}, inplace=True)

X_model2_mat.rename({
    "22001-0.0": 'Gender', "21022-0.0": 'Age', 
    "4079-0.0": 'diastolic blood pressure', "4080-0.0": 'systolic blood pressure', "189-0.0": 'Townsend deprivation index', 
    "22009-0.1": 'PC1', "22009-0.2": 'PC2', "22009-0.3": 'PC3', "22009-0.4": 'PC4', "22009-0.5": 'PC5',
    "22009-0.6": 'PC6', "22009-0.7": 'PC7', "22009-0.8": 'PC8', "22009-0.9": 'PC9', "22009-0.10": 'PC10', "22000-0.0": 'Genotype batch',
    "21001-0.0": 'BMI', # "whr": 'Waist-hip-ratio', "23099-0.0": 'Body Fat Percentage', "21002-0.0": 'Weight',
    "Blood_pressure_medication": 'Blood pressure medication', "Cholesterol_lowering_medication": 'Cholesterol lowering medication', # "No_medication": 'No medication', "Insulin": 'Insulin', 
    "Non_alcohol_drinker": 'Non-alcohol drinker' , "Previous_alcohol_drinker": 'Previous alcohol drinker', "Current_alcohol_drinker": 'Current alcohol drinker',
    "Non_smoker": 'Non-smoker' , "Previous_smoker": 'Previous smoker', "Current_smoker": 'Current smoker'
}, inplace=True)

X_model2b_mat.rename({
    "22001-0.0": 'Gender', "21022-0.0": 'Age', 
    "4079-0.0": 'diastolic blood pressure', "4080-0.0": 'systolic blood pressure', "189-0.0": 'Townsend deprivation index', 
    "22009-0.1": 'PC1', "22009-0.2": 'PC2', "22009-0.3": 'PC3', "22009-0.4": 'PC4', "22009-0.5": 'PC5',
    "22009-0.6": 'PC6', "22009-0.7": 'PC7', "22009-0.8": 'PC8', "22009-0.9": 'PC9', "22009-0.10": 'PC10', "22000-0.0": 'Genotype batch',
    "21001-0.0": 'BMI', "23099-0.0": 'Body Fat Percentage', "21002-0.0": 'Weight', # "whr": 'Waist-hip-ratio', 
    "Blood_pressure_medication": 'Blood pressure medication', "Cholesterol_lowering_medication": 'Cholesterol lowering medication', "Insulin": 'Insulin', # "No_medication": 'No medication', 
    "Non_alcohol_drinker": 'Non-alcohol drinker' , "Previous_alcohol_drinker": 'Previous alcohol drinker', "Current_alcohol_drinker": 'Current alcohol drinker',
    "Non_smoker": 'Non-smoker' , "Previous_smoker": 'Previous smoker', "Current_smoker": 'Current smoker'
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
    "30870-0.0": 'Triglycerides', # lipid-related covariates
    "30700-0.0": 'Creatinine', "30710-0.0": 'C-reactive protein', "30720-0.0": 'Cystatin C', "30730-0.0": 'Gamma glutamyltransferase', "30680-0.0": 'Calcium', 
    "30740-0.0": 'Glucose', "30750-0.0": 'HbA1c', "30650-0.0": 'Aspartate aminotransferase', "30660-0.0": 'Direct bilirubin', 
    "30670-0.0": 'Urea', "30770-0.0": '30770-0.0', "30810-0.0": '30810-0.0', "30830-0.0": 'SHBG', "30850-0.0": 'Testosterone', 
    "30880-0.0": 'Urate', "30890-0.0": 'Vitamin D', "30840-0.0": 'Total bilirubin', "30860-0.0": 'Total protein', 
    "t2dm": 'Type 2 diabetes history', "htn": 'Hypertension history', "heart_failure": 'Heart failure history', "hemorrhage_stroke": 'Hemorrhage Stroke history', "ischemic_stroke": 'Ischemic Stroke history'
}, inplace=True)

W_model1_mat_binary = np.where(W_model1_mat <= 220, 1, 0)
W_model1_mat_binary = pd.DataFrame({"30690-0.0": W_model1_mat_binary}).iloc[:, 0]
W_model1b_mat_binary = np.where(W_model1b_mat <= 220, 1, 0)
W_model1b_mat_binary = pd.DataFrame({"30690-0.0": W_model1b_mat_binary}).iloc[:, 0]
W_model2_mat_binary = np.where(W_model2_mat <= 220, 1, 0)
W_model2_mat_binary = pd.DataFrame({"30690-0.0": W_model2_mat_binary}).iloc[:, 0]
W_model2b_mat_binary = np.where(W_model2b_mat <= 220, 1, 0)
W_model2b_mat_binary = pd.DataFrame({"30690-0.0": W_model2b_mat_binary}).iloc[:, 0]
W_model3_mat_binary = np.where(W_model3_mat <= 220, 1, 0)
W_model3_mat_binary = pd.DataFrame({"30690-0.0": W_model3_mat_binary}).iloc[:, 0]

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
# # model 1
# est_driv_continuousW_model1 = ForestDRIV(projection=False, discrete_treatment=False, discrete_instrument=False,
#                                          n_estimators=5000, min_samples_leaf=100, max_samples=0.02,
#                                          random_state=309, n_jobs=15) 
# est_driv_continuousW_model1.fit(Y_model1_mat, W_model1_mat, Z=Z_model1_mat, X=X_model1_mat, cache_values=True)
# point_driv_continuousW_model1 = est_driv_continuousW_model1.effect(X_model1_mat)
# print(pd.DataFrame({"dat": point_driv_continuousW_model1}).describe())

# model 1b
est_driv_continuousW_model1b = ForestDRIV(projection=False, discrete_treatment=False, discrete_instrument=False,
                                         n_estimators=5000, min_samples_leaf=100, max_samples=0.02,
                                         random_state=309, n_jobs=15) 
est_driv_continuousW_model1b.fit(Y_model1b_mat, W_model1b_mat, Z=Z_model1b_mat, X=X_model1b_mat, cache_values=True)
point_driv_continuousW_model1b = est_driv_continuousW_model1b.effect(X_model1b_mat)
print(pd.DataFrame({"dat": point_driv_continuousW_model1b}).describe())

# # model 2
# est_driv_continuousW_model2 = ForestDRIV(projection=False, discrete_treatment=False, discrete_instrument=False,
#                                          n_estimators=5000, min_samples_leaf=100, max_samples=0.02,
#                                          random_state=309, cov_clip = 10, n_jobs=15) 
# est_driv_continuousW_model2.fit(Y_model2_mat, W_model2_mat, Z=Z_model2_mat, X=X_model2_mat, cache_values=True)
# point_driv_continuousW_model2 = est_driv_continuousW_model2.effect(X_model2_mat)
# print(pd.DataFrame({"dat": point_driv_continuousW_model2}).describe())

# model 2b
est_driv_continuousW_model2b = ForestDRIV(projection=False, discrete_treatment=False, discrete_instrument=False,
                                         n_estimators=5000, min_samples_leaf=100, max_samples=0.02,
                                         random_state=309, cov_clip = 10, n_jobs=15) 
est_driv_continuousW_model2b.fit(Y_model2b_mat, W_model2b_mat, Z=Z_model2b_mat, X=X_model2b_mat, cache_values=True)
point_driv_continuousW_model2b = est_driv_continuousW_model2b.effect(X_model2b_mat)
print(pd.DataFrame({"dat": point_driv_continuousW_model2b}).describe())

# model 3
est_driv_continuousW_model3 = ForestDRIV(projection=False, discrete_treatment=False, discrete_instrument=False, 
                                         n_estimators=5000, min_samples_leaf=150, max_samples=0.02, 
                                         random_state=309, cov_clip = 10, n_jobs=15) 
est_driv_continuousW_model3.fit(Y_model3_mat, W_model3_mat, Z=Z_model3_mat, X=X_model3_mat, cache_values=True)
point_driv_continuousW_model3 = est_driv_continuousW_model3.effect(X_model3_mat)
print(pd.DataFrame({"dat": point_driv_continuousW_model3}).describe())

#' ================================
#' Binary W Continuous Z (Full Set)
#' ================================
# # model 1
# est_driv_binaryW_continuousZ_model1 = ForestDRIV(projection=False, discrete_treatment=True, discrete_instrument=False, 
#                                      n_estimators=5000, min_samples_leaf=100, max_samples=0.02, 
#                                      random_state=309, cov_clip = 0.1, n_jobs=15) 
# est_driv_binaryW_continuousZ_model1.fit(Y_model1_mat, W_model1_mat_binary, Z=Z_model1_mat, X=X_model1_mat, cache_values=True)
# point_driv_binaryW_continuousZ_model1 = est_driv_binaryW_continuousZ_model1.effect(X_model1_mat)
# print(pd.DataFrame({"dat": point_driv_binaryW_continuousZ_model1}).describe())

# model 1b
est_driv_binaryW_continuousZ_model1b = ForestDRIV(projection=False, discrete_treatment=True, discrete_instrument=False, 
                                     n_estimators=5000, min_samples_leaf=100, max_samples=0.02, 
                                     random_state=309, cov_clip = 0.1, n_jobs=15) 
est_driv_binaryW_continuousZ_model1b.fit(Y_model1b_mat, W_model1b_mat_binary, Z=Z_model1b_mat, X=X_model1b_mat, cache_values=True)
point_driv_binaryW_continuousZ_model1b = est_driv_binaryW_continuousZ_model1b.effect(X_model1b_mat)
print(pd.DataFrame({"dat": point_driv_binaryW_continuousZ_model1b}).describe())

# # model 2
# est_driv_binaryW_continuousZ_model2 = ForestDRIV(projection=False, discrete_treatment=True, discrete_instrument=False,
#                                      n_estimators=5000, min_samples_leaf=100, max_samples=0.02, 
#                                      random_state=309, cov_clip = 0.1, n_jobs=15) 
# est_driv_binaryW_continuousZ_model2.fit(Y_model2_mat, W_model2_mat_binary, Z=Z_model2_mat, X=X_model2_mat, cache_values=True)
# point_driv_binaryW_continuousZ_model2 = est_driv_binaryW_continuousZ_model2.effect(X_model2_mat)
# print(pd.DataFrame({"dat": point_driv_binaryW_continuousZ_model2}).describe())

# model 2b
est_driv_binaryW_continuousZ_model2b = ForestDRIV(projection=False, discrete_treatment=True, discrete_instrument=False,
                                     n_estimators=5000, min_samples_leaf=100, max_samples=0.02, 
                                     random_state=309, cov_clip = 0.1, n_jobs=15) 
est_driv_binaryW_continuousZ_model2b.fit(Y_model2b_mat, W_model2b_mat_binary, Z=Z_model2b_mat, X=X_model2b_mat, cache_values=True)
point_driv_binaryW_continuousZ_model2b = est_driv_binaryW_continuousZ_model2b.effect(X_model2b_mat)
print(pd.DataFrame({"dat": point_driv_binaryW_continuousZ_model2b}).describe())

# model 3
est_driv_binaryW_continuousZ_model3 = ForestDRIV(projection=False, discrete_treatment=True, discrete_instrument=False, \
                                                n_estimators=5000, min_samples_leaf=200, max_samples=0.02, 
                                                random_state=309, cov_clip = 0.1, n_jobs=15) 
est_driv_binaryW_continuousZ_model3.fit(Y_model3_mat, W_model3_mat_binary, Z=Z_model3_mat, X=X_model3_mat, cache_values=True)
point_driv_binaryW_continuousZ_model3 = est_driv_binaryW_continuousZ_model3.effect(X_model3_mat)
print(pd.DataFrame({"dat": point_driv_binaryW_continuousZ_model3}).describe())

# #' ================================
# #' Binary W Binary Z (Full Set)
# #' ================================
# # model 1
# est_driv_binaryW_binaryZ_model1 = ForestDRIV(projection=False, discrete_treatment=True, discrete_instrument=True, 
#                                      n_estimators=5000, min_samples_leaf=100, max_samples=0.02, 
#                                      random_state=309, cov_clip = 0.1, n_jobs=15) 
# est_driv_binaryW_binaryZ_model1.fit(Y_model1_mat, W_model1_mat_binary, Z=Z_model1_mat_binary, X=X_model1_mat, cache_values=True)
# point_driv_binaryW_binaryZ_model1 = est_driv_binaryW_binaryZ_model1.effect(X_model1_mat)
# print(pd.DataFrame({"dat": point_driv_binaryW_binaryZ_model1}).describe())

# # model 1b
# est_driv_binaryW_binaryZ_model1b = ForestDRIV(projection=False, discrete_treatment=True, discrete_instrument=True, 
#                                      n_estimators=5000, min_samples_leaf=100, max_samples=0.02, 
#                                      random_state=309, cov_clip = 0.1, n_jobs=15) 
# est_driv_binaryW_binaryZ_model1b.fit(Y_model1b_mat, W_model1b_mat_binary, Z=Z_model1b_mat_binary, X=X_model1b_mat, cache_values=True)
# point_driv_binaryW_binaryZ_model1b = est_driv_binaryW_binaryZ_model1b.effect(X_model1b_mat)
# print(pd.DataFrame({"dat": point_driv_binaryW_binaryZ_model1b}).describe())

# # model 2
# est_driv_binaryW_binaryZ_model2 = ForestDRIV(projection=False, discrete_treatment=True, discrete_instrument=True,
#                                      n_estimators=5000, min_samples_leaf=100, max_samples=0.02, 
#                                      random_state=309, cov_clip = 0.1, n_jobs=15) 
# est_driv_binaryW_binaryZ_model2.fit(Y_model2_mat, W_model2_mat_binary, Z=Z_model2_mat_binary, X=X_model2_mat, cache_values=True)
# point_driv_binaryW_binaryZ_model2 = est_driv_binaryW_binaryZ_model2.effect(X_model2_mat)
# print(pd.DataFrame({"dat": point_driv_binaryW_binaryZ_model2}).describe())

# # model 2b
# est_driv_binaryW_binaryZ_model2b = ForestDRIV(projection=False, discrete_treatment=True, discrete_instrument=True,
#                                      n_estimators=5000, min_samples_leaf=100, max_samples=0.02, 
#                                      random_state=309, cov_clip = 0.1, n_jobs=15) 
# est_driv_binaryW_binaryZ_model2b.fit(Y_model2b_mat, W_model2b_mat_binary, Z=Z_model2b_mat_binary, X=X_model2b_mat, cache_values=True)
# point_driv_binaryW_binaryZ_model2b = est_driv_binaryW_binaryZ_model2b.effect(X_model2b_mat)
# print(pd.DataFrame({"dat": point_driv_binaryW_binaryZ_model2b}).describe())

# # model 3
# est_driv_binaryW_binaryZ_model3 = ForestDRIV(projection=False, discrete_treatment=True, discrete_instrument=True, \
#                                      n_estimators=5000, min_samples_leaf=100, max_samples=0.02, 
#                                      random_state=309, cov_clip = 0.1, n_jobs=15) 
# est_driv_binaryW_binaryZ_model3.fit(Y_model3_mat, W_model3_mat_binary, Z=Z_model3_mat_binary, X=X_model3_mat, cache_values=True)
# point_driv_binaryW_binaryZ_model3 = est_driv_binaryW_binaryZ_model3.effect(X_model3_mat)
# print(pd.DataFrame({"dat": point_driv_binaryW_binaryZ_model3}).describe())

# calculate the upper bound and lower bound (Continuous W)
# full set
# # model 1
# point_driv_lb_continuousW_model1, point_driv_ub_continuousW_model1 = est_driv_continuousW_model1.effect_interval(X_model1_mat, alpha=0.1) # type: ignore
# z_value_driv_continuousW_model1 = point_driv_continuousW_model1/((point_driv_ub_continuousW_model1-point_driv_lb_continuousW_model1)/(2*1.645))
# p_value_driv_continuousW_model1 = scipy.stats.norm.sf(abs(z_value_driv_continuousW_model1)) # * 2
# p_value_BH_driv_continuousW_model1 = multipletests(pvals = p_value_driv_continuousW_model1, method = "fdr_bh", alpha=0.1)

# full_results_continuousW_model1 = pd.DataFrame({"IID": selected_id_set1_arr, "point": point_driv_continuousW_model1, "upper_bound": point_driv_ub_continuousW_model1, \
#                                 "lower_bound": point_driv_lb_continuousW_model1, "z-value": z_value_driv_continuousW_model1, \
#                                 "p_value": p_value_driv_continuousW_model1, "p_value_corrected": p_value_BH_driv_continuousW_model1[1]})

# model 1b
point_driv_lb_continuousW_model1b, point_driv_ub_continuousW_model1b = est_driv_continuousW_model1b.effect_interval(X_model1b_mat, alpha=0.1) # type: ignore
z_value_driv_continuousW_model1b = point_driv_continuousW_model1b/((point_driv_ub_continuousW_model1b-point_driv_lb_continuousW_model1b)/(2*1.645))
p_value_driv_continuousW_model1b = scipy.stats.norm.sf(abs(z_value_driv_continuousW_model1b)) # * 2
p_value_BH_driv_continuousW_model1b = multipletests(pvals = p_value_driv_continuousW_model1b, method = "fdr_bh", alpha=0.1)

full_results_continuousW_model1b = pd.DataFrame({"IID": selected_id_set1b_arr, "point": point_driv_continuousW_model1b, "upper_bound": point_driv_ub_continuousW_model1b, \
                                "lower_bound": point_driv_lb_continuousW_model1b, "z-value": z_value_driv_continuousW_model1b, \
                                "p_value": p_value_driv_continuousW_model1b, "p_value_corrected": p_value_BH_driv_continuousW_model1b[1]})

# # model 2
# point_driv_lb_continuousW_model2, point_driv_ub_continuousW_model2 = est_driv_continuousW_model2.effect_interval(X_model2_mat, alpha=0.1) # type: ignore
# z_value_driv_continuousW_model2 = point_driv_continuousW_model2/((point_driv_ub_continuousW_model2-point_driv_lb_continuousW_model2)/(2*1.645))
# p_value_driv_continuousW_model2 = scipy.stats.norm.sf(abs(z_value_driv_continuousW_model2)) # * 2
# p_value_BH_driv_continuousW_model2 = multipletests(pvals = p_value_driv_continuousW_model2, method = "fdr_bh", alpha=0.1)

# full_results_continuousW_model2 = pd.DataFrame({"IID": selected_id_set2_arr, "point": point_driv_continuousW_model2, "upper_bound": point_driv_ub_continuousW_model2, \
#                                 "lower_bound": point_driv_lb_continuousW_model2, "z-value": z_value_driv_continuousW_model2, \
#                                 "p_value": p_value_driv_continuousW_model2, "p_value_corrected": p_value_BH_driv_continuousW_model2[1]})

# model 2b
point_driv_lb_continuousW_model2b, point_driv_ub_continuousW_model2b = est_driv_continuousW_model2b.effect_interval(X_model2b_mat, alpha=0.1) # type: ignore
z_value_driv_continuousW_model2b = point_driv_continuousW_model2b/((point_driv_ub_continuousW_model2b-point_driv_lb_continuousW_model2b)/(2*1.645))
p_value_driv_continuousW_model2b = scipy.stats.norm.sf(abs(z_value_driv_continuousW_model2b)) # * 2
p_value_BH_driv_continuousW_model2b = multipletests(pvals = p_value_driv_continuousW_model2b, method = "fdr_bh", alpha=0.1)

full_results_continuousW_model2b = pd.DataFrame({"IID": selected_id_set2b_arr, "point": point_driv_continuousW_model2b, "upper_bound": point_driv_ub_continuousW_model2b, \
                                "lower_bound": point_driv_lb_continuousW_model2b, "z-value": z_value_driv_continuousW_model2b, \
                                "p_value": p_value_driv_continuousW_model2b, "p_value_corrected": p_value_BH_driv_continuousW_model2b[1]})

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
# # model 1
# point_driv_lb_binaryW_continuousZ_model1, point_driv_ub_binaryW_continuousZ_model1 = est_driv_binaryW_continuousZ_model1.effect_interval(X_model1_mat, alpha=0.1) # type: ignore
# z_value_driv_binaryW_continuousZ_model1 = point_driv_binaryW_continuousZ_model1/((point_driv_ub_binaryW_continuousZ_model1-point_driv_lb_binaryW_continuousZ_model1)/(2*1.645))
# p_value_driv_binaryW_continuousZ_model1 = scipy.stats.norm.sf(abs(z_value_driv_binaryW_continuousZ_model1)) # * 2
# p_value_BH_driv_binaryW_continuousZ_model1 = multipletests(pvals = p_value_driv_binaryW_continuousZ_model1, method = "fdr_bh", alpha=0.1)

# full_results_binaryW_continuousZ_model1 = pd.DataFrame({"IID": selected_id_set1_arr, "point": point_driv_binaryW_continuousZ_model1, "upper_bound": point_driv_ub_binaryW_continuousZ_model1, \
#                                 "lower_bound": point_driv_lb_binaryW_continuousZ_model1, "z-value": z_value_driv_binaryW_continuousZ_model1, \
#                                 "p_value": p_value_driv_binaryW_continuousZ_model1, "p_value_corrected": p_value_BH_driv_binaryW_continuousZ_model1[1]})

# model 1b
point_driv_lb_binaryW_continuousZ_model1b, point_driv_ub_binaryW_continuousZ_model1b = est_driv_binaryW_continuousZ_model1b.effect_interval(X_model1b_mat, alpha=0.1) # type: ignore
z_value_driv_binaryW_continuousZ_model1b = point_driv_binaryW_continuousZ_model1b/((point_driv_ub_binaryW_continuousZ_model1b-point_driv_lb_binaryW_continuousZ_model1b)/(2*1.645))
p_value_driv_binaryW_continuousZ_model1b = scipy.stats.norm.sf(abs(z_value_driv_binaryW_continuousZ_model1b)) # * 2
p_value_BH_driv_binaryW_continuousZ_model1b = multipletests(pvals = p_value_driv_binaryW_continuousZ_model1b, method = "fdr_bh", alpha=0.1)

full_results_binaryW_continuousZ_model1b = pd.DataFrame({"IID": selected_id_set1b_arr, "point": point_driv_binaryW_continuousZ_model1b, "upper_bound": point_driv_ub_binaryW_continuousZ_model1b, \
                                "lower_bound": point_driv_lb_binaryW_continuousZ_model1b, "z-value": z_value_driv_binaryW_continuousZ_model1b, \
                                "p_value": p_value_driv_binaryW_continuousZ_model1b, "p_value_corrected": p_value_BH_driv_binaryW_continuousZ_model1b[1]})

# # model 2
# point_driv_lb_binaryW_continuousZ_model2, point_driv_ub_binaryW_continuousZ_model2 = est_driv_binaryW_continuousZ_model2.effect_interval(X_model2_mat, alpha=0.1) # type: ignore
# z_value_driv_binaryW_continuousZ_model2 = point_driv_binaryW_continuousZ_model2/((point_driv_ub_binaryW_continuousZ_model2-point_driv_lb_binaryW_continuousZ_model2)/(2*1.645))
# p_value_driv_binaryW_continuousZ_model2 = scipy.stats.norm.sf(abs(z_value_driv_binaryW_continuousZ_model2)) # * 2
# p_value_BH_driv_binaryW_continuousZ_model2 = multipletests(pvals = p_value_driv_binaryW_continuousZ_model2, method = "fdr_bh", alpha=0.1)

# full_results_binaryW_continuousZ_model2 = pd.DataFrame({"IID": selected_id_set2_arr, "point": point_driv_binaryW_continuousZ_model2, "upper_bound": point_driv_ub_binaryW_continuousZ_model2, \
#                                 "lower_bound": point_driv_lb_binaryW_continuousZ_model2, "z-value": z_value_driv_binaryW_continuousZ_model2, \
#                                 "p_value": p_value_driv_binaryW_continuousZ_model2, "p_value_corrected": p_value_BH_driv_binaryW_continuousZ_model2[1]})

# model 2b
point_driv_lb_binaryW_continuousZ_model2b, point_driv_ub_binaryW_continuousZ_model2b = est_driv_binaryW_continuousZ_model2b.effect_interval(X_model2b_mat, alpha=0.1) # type: ignore
z_value_driv_binaryW_continuousZ_model2b = point_driv_binaryW_continuousZ_model2b/((point_driv_ub_binaryW_continuousZ_model2b-point_driv_lb_binaryW_continuousZ_model2b)/(2*1.645))
p_value_driv_binaryW_continuousZ_model2b = scipy.stats.norm.sf(abs(z_value_driv_binaryW_continuousZ_model2b)) # * 2
p_value_BH_driv_binaryW_continuousZ_model2b = multipletests(pvals = p_value_driv_binaryW_continuousZ_model2b, method = "fdr_bh", alpha=0.1)

full_results_binaryW_continuousZ_model2b = pd.DataFrame({"IID": selected_id_set2b_arr, "point": point_driv_binaryW_continuousZ_model2b, "upper_bound": point_driv_ub_binaryW_continuousZ_model2b, \
                                "lower_bound": point_driv_lb_binaryW_continuousZ_model2b, "z-value": z_value_driv_binaryW_continuousZ_model2b, \
                                "p_value": p_value_driv_binaryW_continuousZ_model2b, "p_value_corrected": p_value_BH_driv_binaryW_continuousZ_model2b[1]})

# model 3
point_driv_lb_binaryW_continuousZ_model3, point_driv_ub_binaryW_continuousZ_model3 = est_driv_binaryW_continuousZ_model3.effect_interval(X_model3_mat, alpha=0.1) # type: ignore
z_value_driv_binaryW_continuousZ_model3 = point_driv_binaryW_continuousZ_model3/((point_driv_ub_binaryW_continuousZ_model3-point_driv_lb_binaryW_continuousZ_model3)/(2*1.645))
p_value_driv_binaryW_continuousZ_model3 = scipy.stats.norm.sf(abs(z_value_driv_binaryW_continuousZ_model3)) # * 2
p_value_BH_driv_binaryW_continuousZ_model3 = multipletests(pvals = p_value_driv_binaryW_continuousZ_model3, method = "fdr_bh", alpha=0.1)

full_results_binaryW_continuousZ_model3 = pd.DataFrame({"IID": selected_id_set3_arr, "point": point_driv_binaryW_continuousZ_model3, "upper_bound": point_driv_ub_binaryW_continuousZ_model3, \
                                "lower_bound": point_driv_lb_binaryW_continuousZ_model3, "z-value": z_value_driv_binaryW_continuousZ_model3, \
                                "p_value": p_value_driv_binaryW_continuousZ_model3, "p_value_corrected": p_value_BH_driv_binaryW_continuousZ_model3[1]})

# # calculate the upper bound and lower bound (Binary W, Binary Z)
# # full set
# # model 1
# point_driv_lb_binaryW_binaryZ_model1, point_driv_ub_binaryW_binaryZ_model1 = est_driv_binaryW_binaryZ_model1.effect_interval(X_model1_mat, alpha=0.1) # type: ignore
# z_value_driv_binaryW_binaryZ_model1 = point_driv_binaryW_binaryZ_model1/((point_driv_ub_binaryW_binaryZ_model1-point_driv_lb_binaryW_binaryZ_model1)/(2*1.645))
# p_value_driv_binaryW_binaryZ_model1 = scipy.stats.norm.sf(abs(z_value_driv_binaryW_binaryZ_model1)) # * 2
# p_value_BH_driv_binaryW_binaryZ_model1 = multipletests(pvals = p_value_driv_binaryW_binaryZ_model1, method = "fdr_bh", alpha=0.1)

# full_results_binaryW_binaryZ_model1 = pd.DataFrame({"IID": selected_id_set1_arr, "point": point_driv_binaryW_binaryZ_model1, "upper_bound": point_driv_ub_binaryW_binaryZ_model1, \
#                                 "lower_bound": point_driv_lb_binaryW_binaryZ_model1, "z-value": z_value_driv_binaryW_binaryZ_model1, \
#                                 "p_value": p_value_driv_binaryW_binaryZ_model1, "p_value_corrected": p_value_BH_driv_binaryW_binaryZ_model1[1]})

# # model 1b
# point_driv_lb_binaryW_binaryZ_model1b, point_driv_ub_binaryW_binaryZ_model1b = est_driv_binaryW_binaryZ_model1b.effect_interval(X_model1b_mat, alpha=0.1) # type: ignore
# z_value_driv_binaryW_binaryZ_model1b = point_driv_binaryW_binaryZ_model1b/((point_driv_ub_binaryW_binaryZ_model1b-point_driv_lb_binaryW_binaryZ_model1b)/(2*1.645))
# p_value_driv_binaryW_binaryZ_model1b = scipy.stats.norm.sf(abs(z_value_driv_binaryW_binaryZ_model1b)) # * 2
# p_value_BH_driv_binaryW_binaryZ_model1b = multipletests(pvals = p_value_driv_binaryW_binaryZ_model1b, method = "fdr_bh", alpha=0.1)

# full_results_binaryW_binaryZ_model1b = pd.DataFrame({"IID": selected_id_set1b_arr, "point": point_driv_binaryW_binaryZ_model1b, "upper_bound": point_driv_ub_binaryW_binaryZ_model1b, \
#                                 "lower_bound": point_driv_lb_binaryW_binaryZ_model1b, "z-value": z_value_driv_binaryW_binaryZ_model1b, \
#                                 "p_value": p_value_driv_binaryW_binaryZ_model1b, "p_value_corrected": p_value_BH_driv_binaryW_binaryZ_model1b[1]})

# # model 2
# point_driv_lb_binaryW_binaryZ_model2, point_driv_ub_binaryW_binaryZ_model2 = est_driv_binaryW_binaryZ_model2.effect_interval(X_model2_mat, alpha=0.1) # type: ignore
# z_value_driv_binaryW_binaryZ_model2 = point_driv_binaryW_binaryZ_model2/((point_driv_ub_binaryW_binaryZ_model2-point_driv_lb_binaryW_binaryZ_model2)/(2*1.645))
# p_value_driv_binaryW_binaryZ_model2 = scipy.stats.norm.sf(abs(z_value_driv_binaryW_binaryZ_model2)) # * 2
# p_value_BH_driv_binaryW_binaryZ_model2 = multipletests(pvals = p_value_driv_binaryW_binaryZ_model2, method = "fdr_bh", alpha=0.1)

# full_results_binaryW_binaryZ_model2 = pd.DataFrame({"IID": selected_id_set2_arr, "point": point_driv_binaryW_binaryZ_model2, "upper_bound": point_driv_ub_binaryW_binaryZ_model2, \
#                                 "lower_bound": point_driv_lb_binaryW_binaryZ_model2, "z-value": z_value_driv_binaryW_binaryZ_model2, \
#                                 "p_value": p_value_driv_binaryW_binaryZ_model2, "p_value_corrected": p_value_BH_driv_binaryW_binaryZ_model2[1]})

# # model 2b
# point_driv_lb_binaryW_binaryZ_model2b, point_driv_ub_binaryW_binaryZ_model2b = est_driv_binaryW_binaryZ_model2b.effect_interval(X_model2b_mat, alpha=0.1) # type: ignore
# z_value_driv_binaryW_binaryZ_model2b = point_driv_binaryW_binaryZ_model2b/((point_driv_ub_binaryW_binaryZ_model2b-point_driv_lb_binaryW_binaryZ_model2b)/(2*1.645))
# p_value_driv_binaryW_binaryZ_model2b = scipy.stats.norm.sf(abs(z_value_driv_binaryW_binaryZ_model2b)) # * 2
# p_value_BH_driv_binaryW_binaryZ_model2b = multipletests(pvals = p_value_driv_binaryW_binaryZ_model2b, method = "fdr_bh", alpha=0.1)

# full_results_binaryW_binaryZ_model2b = pd.DataFrame({"IID": selected_id_set2b_arr, "point": point_driv_binaryW_binaryZ_model2b, "upper_bound": point_driv_ub_binaryW_binaryZ_model2b, \
#                                 "lower_bound": point_driv_lb_binaryW_binaryZ_model2b, "z-value": z_value_driv_binaryW_binaryZ_model2b, \
#                                 "p_value": p_value_driv_binaryW_binaryZ_model2b, "p_value_corrected": p_value_BH_driv_binaryW_binaryZ_model2b[1]})

# # model 3
# point_driv_lb_binaryW_binaryZ_model3, point_driv_ub_binaryW_binaryZ_model3 = est_driv_binaryW_binaryZ_model3.effect_interval(X_model3_mat, alpha=0.1) # type: ignore
# z_value_driv_binaryW_binaryZ_model3 = point_driv_binaryW_binaryZ_model3/((point_driv_ub_binaryW_binaryZ_model3-point_driv_lb_binaryW_binaryZ_model3)/(2*1.645))
# p_value_driv_binaryW_binaryZ_model3 = scipy.stats.norm.sf(abs(z_value_driv_binaryW_binaryZ_model3)) # * 2
# p_value_BH_driv_binaryW_binaryZ_model3 = multipletests(pvals = p_value_driv_binaryW_binaryZ_model3, method = "fdr_bh", alpha=0.1)

# full_results_binaryW_binaryZ_model3 = pd.DataFrame({"IID": selected_id_set3_arr, "point": point_driv_binaryW_binaryZ_model3, "upper_bound": point_driv_ub_binaryW_binaryZ_model3, \
#                                 "lower_bound": point_driv_lb_binaryW_binaryZ_model3, "z-value": z_value_driv_binaryW_binaryZ_model3, \
#                                 "p_value": p_value_driv_binaryW_binaryZ_model3, "p_value_corrected": p_value_BH_driv_binaryW_binaryZ_model3[1]})

print("Summary of lb rb results: ")
print("Continuous W: ")
# print("Model 1:")
# print(full_results_continuousW_model1.describe())
print("Model 1b:")
print(full_results_continuousW_model1b.describe())
# print("Model 2:")
# print(full_results_continuousW_model2.describe())
print("Model 2b:")
print(full_results_continuousW_model2b.describe())
print("Model 3:")
print(full_results_continuousW_model3.describe())
print("                              ")
print("Binary W Continuou Z: ")
# print("Model 1:")
# print(full_results_binaryW_continuousZ_model1.describe())
print("Model 1b:")
print(full_results_binaryW_continuousZ_model1b.describe())
# print("Model 2:")
# print(full_results_binaryW_continuousZ_model2.describe())
print("Model 2b:")
print(full_results_binaryW_continuousZ_model2b.describe())
print("Model 3:")
print(full_results_binaryW_continuousZ_model3.describe())
print("                              ")
# print("Binary W Binary Z: ")
# print("Model 1:")
# print(full_results_binaryW_binaryZ_model1.describe())
# print("Model 1b:")
# print(full_results_binaryW_binaryZ_model1b.describe())
# print("Model 2:")
# print(full_results_binaryW_binaryZ_model2.describe())
# print("Model 2b:")
# print(full_results_binaryW_binaryZ_model2b.describe())
# print("Model 3:")
# print(full_results_binaryW_binaryZ_model3.describe())
# print("                              ")

# full_results_continuousW_model1.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_continuousW_te_ul_model1.csv", sep=",", index=False)
# full_results_binaryW_continuousZ_model1.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_model1.csv", sep=",", index=False)
# full_results_binaryW_binaryZ_model1.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_binaryW_binaryZ_te_ul_model1.csv", sep=",", index=False)

full_results_continuousW_model1b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_continuousW_te_ul_model1b.csv", sep=",", index=False)
full_results_binaryW_continuousZ_model1b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_model1b.csv", sep=",", index=False)
# full_results_binaryW_binaryZ_model1b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_binaryW_binaryZ_te_ul_model1b.csv", sep=",", index=False)

# full_results_continuousW_model2.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_continuousW_te_ul_model2.csv", sep=",", index=False)
# full_results_binaryW_continuousZ_model2.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_model2.csv", sep=",", index=False)
# full_results_binaryW_binaryZ_model2.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_binaryW_binaryZ_te_ul_model2.csv", sep=",", index=False)

full_results_continuousW_model2b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_continuousW_te_ul_model2b.csv", sep=",", index=False)
full_results_binaryW_continuousZ_model2b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_model2b.csv", sep=",", index=False)
# full_results_binaryW_binaryZ_model2b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_binaryW_binaryZ_te_ul_model2b.csv", sep=",", index=False)

full_results_continuousW_model3.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_continuousW_te_ul_model3.csv", sep=",", index=False)
full_results_binaryW_continuousZ_model3.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_binaryW_continuousZ_te_ul_model3.csv", sep=",", index=False)
# full_results_binaryW_binaryZ_model3.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/01_individual_treatment_effect/driv_full_binaryW_binaryZ_te_ul_model3.csv", sep=",", index=False)

#########################
# Variabel Importance
#########################

# # save important variables
# # model 1
# print("Start calculating SHAP value.")
# shap_values_driv_binaryW_model1 = est_driv_binaryW_model1.shap_values(X_model1_mat)
# shap_pd_binaryW_model1 = pd.DataFrame(shap_values_driv_binaryW_model1['CAD']['30690-0.0'].values)
# shap_pd_binaryW_model1.columns = shap_values_driv_binaryW_model1['CAD']['30690-0.0'].feature_names
# value_pd_binaryW_model1 = pd.DataFrame(shap_values_driv_binaryW_model1['CAD']['30690-0.0'].data)
# value_pd_binaryW_model1.columns = shap_values_driv_binaryW_model1['CAD']['30690-0.0'].feature_names
# shap_pd_binaryW_model1.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_binaryW_shap_model1.csv.gz", index=False)
# value_pd_binaryW_model1.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_binaryW_value_model1.csv.gz", index=False)
# print("Finished calculating SHAP value on binary W model1.")
# shap_values_driv_continuousW_model1 = est_driv_continuousW_model1.shap_values(X_model1_mat)
# shap_pd_continuousW_model1 = pd.DataFrame(shap_values_driv_continuousW_model1['CAD']['30690-0.0'].values)
# shap_pd_continuousW_model1.columns = shap_values_driv_continuousW_model1['CAD']['30690-0.0'].feature_names
# value_pd_continuousW_model1 = pd.DataFrame(shap_values_driv_continuousW_model1['CAD']['30690-0.0'].data)
# value_pd_continuousW_model1.columns = shap_values_driv_continuousW_model1['CAD']['30690-0.0'].feature_names
# shap_pd_continuousW_model1.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_continuousW_shap_model1.csv.gz", index=False)
# value_pd_continuousW_model1.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_continuousW_value_model1.csv.gz", index=False)
# print("Finished calculating SHAP value on continuous W model1.")

# # model 1b
# print("Start calculating SHAP value.")
# shap_values_driv_binaryW_model1b = est_driv_binaryW_model1b.shap_values(X_model1b_mat)
# shap_pd_binaryW_model1b = pd.DataFrame(shap_values_driv_binaryW_model1b['CAD']['30690-0.0'].values)
# shap_pd_binaryW_model1b.columns = shap_values_driv_binaryW_model1b['CAD']['30690-0.0'].feature_names
# value_pd_binaryW_model1b = pd.DataFrame(shap_values_driv_binaryW_model1b['CAD']['30690-0.0'].data)
# value_pd_binaryW_model1b.columns = shap_values_driv_binaryW_model1b['CAD']['30690-0.0'].feature_names
# shap_pd_binaryW_model1b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_binaryW_shap_model1b.csv.gz", index=False)
# value_pd_binaryW_model1b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_binaryW_value_model1b.csv.gz", index=False)
# print("Finished calculating SHAP value on binary W model1b.")
# shap_values_driv_continuousW_model1b = est_driv_continuousW_model1b.shap_values(X_model1b_mat)
# shap_pd_continuousW_model1b = pd.DataFrame(shap_values_driv_continuousW_model1b['CAD']['30690-0.0'].values)
# shap_pd_continuousW_model1b.columns = shap_values_driv_continuousW_model1b['CAD']['30690-0.0'].feature_names
# value_pd_continuousW_model1b = pd.DataFrame(shap_values_driv_continuousW_model1b['CAD']['30690-0.0'].data)
# value_pd_continuousW_model1b.columns = shap_values_driv_continuousW_model1b['CAD']['30690-0.0'].feature_names
# shap_pd_continuousW_model1b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_continuousW_shap_model1b.csv.gz", index=False)
# value_pd_continuousW_model1b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_continuousW_value_model1b.csv.gz", index=False)
# print("Finished calculating SHAP value on continuous W model1b.")

# # model 2
# print("Start calculating SHAP value.")
# shap_values_driv_binaryW_model2 = est_driv_binaryW_model2.shap_values(X_model2_mat)
# shap_pd_binaryW_model2 = pd.DataFrame(shap_values_driv_binaryW_model2['CAD']['30690-0.0'].values)
# shap_pd_binaryW_model2.columns = shap_values_driv_binaryW_model2['CAD']['30690-0.0'].feature_names
# value_pd_binaryW_model2 = pd.DataFrame(shap_values_driv_binaryW_model2['CAD']['30690-0.0'].data)
# value_pd_binaryW_model2.columns = shap_values_driv_binaryW_model2['CAD']['30690-0.0'].feature_names
# shap_pd_binaryW_model2.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_binaryW_shap_model2.csv.gz", index=False)
# value_pd_binaryW_model2.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_binaryW_value_model2.csv.gz", index=False)
# print("Finished calculating SHAP value on binary W model2.")
# shap_values_driv_continuousW_model2 = est_driv_continuousW_model2.shap_values(X_model2_mat)
# shap_pd_continuousW_model2 = pd.DataFrame(shap_values_driv_continuousW_model2['CAD']['30690-0.0'].values)
# shap_pd_continuousW_model2.columns = shap_values_driv_continuousW_model2['CAD']['30690-0.0'].feature_names
# value_pd_continuousW_model2 = pd.DataFrame(shap_values_driv_continuousW_model2['CAD']['30690-0.0'].data)
# value_pd_continuousW_model2.columns = shap_values_driv_continuousW_model2['CAD']['30690-0.0'].feature_names
# shap_pd_continuousW_model2.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_continuousW_shap_model2.csv.gz", index=False)
# value_pd_continuousW_model2.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_continuousW_value_model2.csv.gz", index=False)
# print("Finished calculating SHAP value on continuous W model2.")

# # model 2b
# print("Start calculating SHAP value.")
# shap_values_driv_binaryW_model2b = est_driv_binaryW_model2b.shap_values(X_model2b_mat)
# shap_pd_binaryW_model2b = pd.DataFrame(shap_values_driv_binaryW_model2b['CAD']['30690-0.0'].values)
# shap_pd_binaryW_model2b.columns = shap_values_driv_binaryW_model2b['CAD']['30690-0.0'].feature_names
# value_pd_binaryW_model2b = pd.DataFrame(shap_values_driv_binaryW_model2b['CAD']['30690-0.0'].data)
# value_pd_binaryW_model2b.columns = shap_values_driv_binaryW_model2b['CAD']['30690-0.0'].feature_names
# shap_pd_binaryW_model2b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_binaryW_shap_model2b.csv.gz", index=False)
# value_pd_binaryW_model2b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_binaryW_value_model2b.csv.gz", index=False)
# print("Finished calculating SHAP value on binary W model2b.")
# shap_values_driv_continuousW_model2b = est_driv_continuousW_model2b.shap_values(X_model2b_mat)
# shap_pd_continuousW_model2b = pd.DataFrame(shap_values_driv_continuousW_model2b['CAD']['30690-0.0'].values)
# shap_pd_continuousW_model2b.columns = shap_values_driv_continuousW_model2b['CAD']['30690-0.0'].feature_names
# value_pd_continuousW_model2b = pd.DataFrame(shap_values_driv_continuousW_model2b['CAD']['30690-0.0'].data)
# value_pd_continuousW_model2b.columns = shap_values_driv_continuousW_model2b['CAD']['30690-0.0'].feature_names
# shap_pd_continuousW_model2b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_continuousW_shap_model2b.csv.gz", index=False)
# value_pd_continuousW_model2b.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_continuousW_value_model2b.csv.gz", index=False)
# print("Finished calculating SHAP value on continuous W model2b.")

# model 3
print("Start calculating SHAP value.")
shap_values_driv_binaryW_continuousZ_model3 = est_driv_binaryW_continuousZ_model3.shap_values(X_model3_mat)
save_object(shap_values_driv_binaryW_continuousZ_model3, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_binaryW_continuousZ_shap_model3.pkl")
print("Finished calculating SHAP value on binary W continuous Z model3.")
shap_values_driv_continuousW_model3 = est_driv_continuousW_model3.shap_values(X_model3_mat)
save_object(shap_values_driv_continuousW_model3, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_continuousW_shap_model3.pkl")
print("Finished calculating SHAP value on continuous W model3.")

# shap_values_driv_binaryW_binaryZ_model3 = est_driv_binaryW_binaryZ_model3.shap_values(X_model3_mat)
# save_object(shap_values_driv_binaryW_binaryZ_model3, "/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/TC/CAD/03_variable_importance/driv_binaryW_binaryZ_shap_model3.pkl")
# print("Finished calculating SHAP value on binary W binary Z model3.")

print("Finished the pipeline.")