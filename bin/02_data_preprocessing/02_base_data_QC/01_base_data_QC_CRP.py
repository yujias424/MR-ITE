# base data quality control
# following Choi's QC pipeline (we use UK summary statistics)
# 
# Author: Yujia Shi
# Date: 2023.02.08 

import os
import re

import numpy as np
import pandas as pd

os.chdir("/home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/")

# I have already use GENOMA to check the heritability 
# The pipeline starts from standard GWAS QC

# =======================
# Step 1: format the data
# =======================

crp = pd.read_csv("/home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CRP/download/03_UK/GCST90029070_buildGRCh37.tsv.gz", 
                  sep="\t", compression = 'gzip')

print("Total SNP in the GWAS data: ")
print(crp.shape)

# rename the data
crp.rename(columns={"chromosome":"CHR", "base_pair_location":"BP", "variant_id":"SNP", 
                    "beta":"BETA", "standard_error":"SE", "p_value":"P",
                    "effect_allele":"A1", "other_allele":"A2"}, inplace=True)

# reorder the column
cols = ["CHR", "BP", "SNP", "A1", "A2", "SE", "P", "BETA"]
crp = crp[cols]

# Remove NaN row
crp = crp.dropna()

# convert to correct column type
crp = crp.astype({'CHR':'int', 'BP':'int', "SNP":'str', 'A1':'str', 'A2':'str', 'SE':'float', 'P':'float', 'BETA':'float'})

# remove SNP without SNPID
crp = crp[crp['SNP'].str.startswith('rs')]

# save the formatted results
crp.to_csv("/home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CRP/QC/CRP_formated_nadrop.tsv.gz", sep="\t", index=False)

print("SNP left after basic QC: ")
print(crp.shape)

# SNP left after basic QC: 
# (10600971, 8)

# # =============================
# # Step 2: remove SNPs with MAF < 0.01 (No MAF info, no need to run.)
# # =============================

# cmd = "gunzip -c /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CRP/QC/CRP_formated_nadrop.tsv.gz |" + \
#       "awk 'NR==1 || ($10 > 0.01) {print}' |" + \
#       "gzip > /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CRP/QC/CRP.standardGWASQC.gz"
# os.system(cmd)

# cmd = "zcat /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CRP/QC/CRP.standardGWASQC.gz | wc -l"
# print("SNP left after removing SNPs with MAF < 0.01: ")
# print(os.popen(cmd).read())

# =============================
# Step 3: remove duplicate SNPS
# =============================

cmd = "gunzip -c /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CRP/QC/CRP_formated_nadrop.tsv.gz |" + \
      "awk '{seen[$3]++; if(seen[$3]==1){ print}}' |" + \
      "gzip - > /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CRP/QC/CRP.nodup.gz"
os.system(cmd)

cmd = "zcat /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CRP/QC/CRP.nodup.gz | wc -l"
print("SNP left after removing duplicate SNPs: ")
print(os.popen(cmd).read())

# SNP left after removing duplicate SNPs: 
# 10600971

# =============================
# Step 4: keep non-ambiguous SNPS
# =============================

cmd = "gunzip -c /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CRP/QC/CRP.nodup.gz | " + \
      "awk '!( ($4==\"A\" && $5==\"T\") || ($4==\"T\" && $5==\"A\") || ($4==\"G\" && $5==\"C\") || ($4==\"C\" && $5==\"G\")) {print}' | " + \
      "gzip > /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CRP/QC/CRP.QC.gz"
os.system(cmd)

cmd = "zcat /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CRP/QC/CRP.QC.gz | wc -l"
print("SNP left after removing ambiguous SNPs: ")
print(os.popen(cmd).read())

# SNP left after removing ambiguous SNPs: 
# 9053050