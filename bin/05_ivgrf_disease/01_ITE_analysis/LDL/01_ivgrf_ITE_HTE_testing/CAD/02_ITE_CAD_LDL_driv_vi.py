import pandas as pd
import numpy as np
import pickle

f = open("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/03_variable_importance/driv_binaryW_continuousZ_shap_model3.pkl", "rb")
bbb = pickle.load(f)
f.close()

shap_pd_binaryW_model3 = pd.DataFrame(bbb['CAD']['30780-0.0_1'].values)
shap_pd_binaryW_model3.columns = bbb['CAD']['30780-0.0_1'].feature_names
value_pd_binaryW_model3 = pd.DataFrame(bbb['CAD']['30780-0.0_1'].data)
value_pd_binaryW_model3.columns = bbb['CAD']['30780-0.0_1'].feature_names
shap_pd_binaryW_model3.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/03_variable_importance/driv_binaryW_shap_model3.csv.gz", index=False)
value_pd_binaryW_model3.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/03_variable_importance/driv_binaryW_value_model3.csv.gz", index=False)
print("Finished calculating SHAP value on binary W model3.")

f = open("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/03_variable_importance/driv_continuousW_shap_model3.pkl", "rb")
aaa = pickle.load(f)
f.close()

shap_pd_continuousW_model3 = pd.DataFrame(aaa['CAD']['30780-0.0'].values)
shap_pd_continuousW_model3.columns = aaa['CAD']['30780-0.0'].feature_names
value_pd_continuousW_model3 = pd.DataFrame(aaa['CAD']['30780-0.0'].data)
value_pd_continuousW_model3.columns = aaa['CAD']['30780-0.0'].feature_names
shap_pd_continuousW_model3.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/03_variable_importance/driv_continuousW_shap_model3.csv.gz", index=False)
value_pd_continuousW_model3.to_csv("/home/yujia/Project/2023-07-20-individual_MR/res/02_ITE_analysis/01_table/LDL/CAD/03_variable_importance/driv_continuousW_value_model3.csv.gz", index=False)
print("Finished calculating SHAP value on continuous W model3.")