# base data quality control
# following Choi's QC pipeline
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

cystatinC = pd.read_csv("~/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/download/formatted_round3meta_eGFRcys_overall_IV_2GC_b36_MAFget005_Nget50_20120725_b37.csv.gz", sep=",", compression = 'gzip')

print("Total SNP in the GWAS data: ")
print(cystatinC.shape)

# rename the data
cystatinC.rename(columns={"rsID":"SNP", "allele1":"A1", "allele2":"A2", "freqA1":"EAF", "beta":"BETA", "se":"SE", "pval":"P", "N":"N"}, inplace=True)

# convert EAF to MAF
cystatinC['MAF'] = cystatinC['EAF'].apply(lambda x : x if x <= 0.5 else 1-x)

# Uppercase the A1 and A2
cystatinC['A1'] = cystatinC['A1'].str.upper()
cystatinC['A2'] = cystatinC['A2'].str.upper()

# reorder the column
cols = ["SNP", "A1", "A2", "SE", "P", "BETA", "MAF", "N", "EAF"]
cystatinC = cystatinC[cols]

# Remove NaN row
cystatinC = cystatinC.dropna()

# convert to correct column type
cystatinC = cystatinC.astype({"SNP":'str', 'A1':'str', 'A2':'str', 'SE':'float', 'P':'float', 'BETA':'float', 'MAF':'float', 'N':'int', 'EAF':'float'})

# remove SNP without SNPID
cystatinC = cystatinC[cystatinC['SNP'].str.startswith('rs')]

# save the formatted results
cystatinC.to_csv("/home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC_formated_nadrop.tsv.gz", sep="\t", index=False)

print("SNP left after basic QC: ")
print(cystatinC.shape)

# =============================
# Step 2: remove SNPs with MAF < 0.01
# =============================

cmd = "gunzip -c /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC_formated_nadrop.tsv.gz |" + \
      "awk 'NR==1 || ($7 > 0.01) {print}' |" + \
      "gzip > /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC.standardGWASQC.gz"
os.system(cmd)

cmd = "zcat /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC.standardGWASQC.gz | wc -l"
print("SNP left after removing SNPs with MAF < 0.01: ")
print(os.popen(cmd).read())

# SNP left after removing SNPs with MAF < 0.01: 
# 2194591

# =============================
# Step 3: remove duplicate SNPS
# =============================

cmd = "gunzip -c /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC.standardGWASQC.gz |" + \
      "awk '{seen[$1]++; if(seen[$1]==1){ print}}' |" + \
      "gzip - > /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC.nodup.gz"
os.system(cmd)

cmd = "zcat /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC.nodup.gz | wc -l"
print("SNP left after removing duplicate SNPs: ")
print(os.popen(cmd).read())

# SNP left after removing duplicate SNPs: 
# 2194591

# =============================
# Step 4: keep non-ambiguous SNPS
# =============================

cmd = "gunzip -c /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC.nodup.gz | " + \
      "awk '!( ($2==\"A\" && $3==\"T\") || ($2==\"T\" && $3==\"A\") || ($2==\"G\" && $3==\"C\") || ($2==\"C\" && $3==\"G\")) {print}' | " + \
      "gzip > /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC.QC.gz"
os.system(cmd)

cmd = "zcat /home/yujia/Project/2023-07-20-individual_MR/dat/02_base_data/CystatinC/QC/CystatinC.QC.gz | wc -l"
print("SNP left after removing ambiguous SNPs: ")
print(os.popen(cmd).read())

# SNP left after removing ambiguous SNPs: 
# 1855067