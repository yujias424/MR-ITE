# base data quality control
# following Choi's QC pipeline (we use UK summary statistics)
# 
# Author: Yujia Shi
# Date: 2023.02.08 

import os
import re

import numpy as np
import pandas as pd

# os.chdir("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/")

# I have already use GENOMA to check the heritability 
# The pipeline starts from standard GWAS QC

# =======================
# Step 1: format the data
# =======================

IGF1 = pd.read_csv("~/Project/2023-07-20-individual_MR/dat/02_base_data/IGF1/download/GCST90102623.tsv.gz", sep="\t", compression = 'gzip')

print("Total SNP in the GWAS data: ")
print(IGF1.shape)

# rename the data
IGF1.rename(columns={"chromosome":"CHR", "base_pair_location":"BP", "rs_id":"SNP", "effect_allele":"A1", "other_allele":"A2", "effect_allele_frequency":"EAF", "beta":"BETA", "standard_error":"SE", "p_value":"P", "n":"N"}, inplace=True)

# convert EAF to MAF
IGF1['MAF'] = IGF1['EAF'].apply(lambda x : x if x <= 0.5 else 1-x)

# Uppercase the A1 and A2
IGF1['A1'] = IGF1['A1'].str.upper()
IGF1['A2'] = IGF1['A2'].str.upper()

# reorder the column
cols = ["CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "BETA", "MAF", "EAF"]
IGF1 = IGF1[cols]

# Remove NaN row
IGF1 = IGF1.dropna()

# convert to correct column type
IGF1 = IGF1.astype({'CHR':'int', 'BP':'int', "SNP":'str', 'A1':'str', 'A2':'str', 'N':'int', 'SE':'float', 'P':'float', 'BETA':'float', 'MAF':'float', 'EAF':'float'})

# remove SNP without SNPID
IGF1 = IGF1[IGF1['SNP'].str.startswith('rs')]

# save the formatted results
IGF1.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1_formated_nadrop.tsv.gz", sep="\t", index=False)

print("SNP left after basic QC: ")
print(IGF1.shape)

# SNP left after basic QC: 
# (2459688, 11)

# =============================
# Step 2: remove SNPs with MAF < 0.01 (No MAF info, no need to run.)
# =============================

cmd = "gunzip -c /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1_formated_nadrop.tsv.gz |" + \
      "awk 'NR==1 || ($10 > 0.01) {print}' |" + \
      "gzip > /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1.standardGWASQC.gz"
os.system(cmd)

cmd = "zcat /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1.standardGWASQC.gz | wc -l"
print("SNP left after removing SNPs with MAF < 0.01: ")
print(os.popen(cmd).read())

# SNP left after removing duplicate SNPs: 
# 2459688

# =============================
# Step 3: remove duplicate SNPS
# =============================

cmd = "gunzip -c /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1.standardGWASQC.gz |" + \
      "awk '{seen[$3]++; if(seen[$3]==1){ print}}' |" + \
      "gzip - > /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1.nodup.gz"
os.system(cmd)

cmd = "zcat /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1.nodup.gz | wc -l"
print("SNP left after removing duplicate SNPs: ")
print(os.popen(cmd).read())

# SNP left after removing duplicate SNPs: 
# 2459688

# =============================
# Step 4: keep non-ambiguous SNPS
# =============================

cmd = "gunzip -c /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1.nodup.gz | " + \
      "awk '!( ($4==\"A\" && $5==\"T\") || ($4==\"T\" && $5==\"A\") || ($4==\"G\" && $5==\"C\") || ($4==\"C\" && $5==\"G\")) {print}' | " + \
      "gzip > /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1.QC.gz"
os.system(cmd)

cmd = "zcat /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1.QC.gz | wc -l"
print("SNP left after removing ambiguous SNPs: ")
print(os.popen(cmd).read())

# SNP left after removing ambiguous SNPs: 
# 2080031