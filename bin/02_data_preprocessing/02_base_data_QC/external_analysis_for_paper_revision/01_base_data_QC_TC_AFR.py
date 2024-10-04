# base data quality control
# following Choi's QC pipeline
# 
# Author: Yujia Shi
# Date: 2023.02.08 

import os
import re

import numpy as np
import pandas as pd

os.chdir("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/")

# I have already use GENOMA to check the heritability 
# The pipeline starts from standard GWAS QC

# =======================
# Step 1: format the data
# =======================

TC = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TC/download/without_UKB_TC_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz", sep="\t", compression = 'gzip')

print("Total SNP in the GWAS data: ")
print(TC.shape)

# Total SNP in the GWAS data: 
# (31916286, 14)

# rename the data
TC.rename(columns={"CHROM":"CHR", "POS_b37":"BP", "rsID":"SNP", "REF":"A2", "ALT":"A1", "POOLED_ALT_AF":"EAF", "EFFECT_SIZE":"BETA", "SE":"SE", "pvalue":"P", "N":"N"}, inplace=True)

# convert EAF to MAF
TC['MAF'] = TC['EAF'].apply(lambda x : x if x <= 0.5 else 1-x)

# Uppercase the A1 and A2
TC['A1'] = TC['A1'].str.upper()
TC['A2'] = TC['A2'].str.upper()

# reorder the column
cols = ["CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "BETA", "MAF", "EAF"]
TC = TC[cols]

# Remove NaN row
TC = TC.dropna()

# convert to correct column type
TC = TC.astype({'CHR':'int', 'BP':'int', "SNP":'str', 'A1':'str', 'A2':'str', 'N':'int', 'SE':'float', 'P':'float', 'BETA':'float', 'MAF':'float', 'EAF':'float'})

# remove SNP without SNPID
TC = TC[TC['SNP'].str.startswith('rs')]

# save the formatted results
TC.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TC/QC_AFR/TC_formated_nadrop.tsv.gz", sep="\t", index=False)

print("SNP left after basic QC: ")
print(TC.shape)

# SNP left after basic QC: 
# (31763752, 11)

# =============================
# Step 2: remove SNPs with MAF < 0.01
# =============================

cmd = "gunzip -c /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TC/QC_AFR/TC_formated_nadrop.tsv.gz |" + \
      "awk 'NR==1 || ($10 > 0.01) {print}' |" + \
      "gzip > /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TC/QC_AFR/TC.standardGWASQC.gz"
os.system(cmd)

cmd = "zcat /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TC/QC_AFR/TC.standardGWASQC.gz | wc -l"
print("SNP left after removing SNPs with MAF < 0.01: ")
print(os.popen(cmd).read())

# SNP left after removing SNPs with MAF < 0.01: 
# 16298613

# =============================
# Step 3: remove duplicate SNPS
# =============================

cmd = "gunzip -c /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TC/QC_AFR/TC.standardGWASQC.gz |" + \
      "awk '{seen[$3]++; if(seen[$3]==1){ print}}' |" + \
      "gzip - > /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TC/QC_AFR/TC.nodup.gz"
os.system(cmd)

cmd = "zcat /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TC/QC_AFR/TC.nodup.gz | wc -l"
print("SNP left after removing duplicate SNPs: ")
print(os.popen(cmd).read())

# SNP left after removing duplicate SNPs: 
# 16242897

# =============================
# Step 4: keep non-ambiguous SNPS
# =============================

cmd = "gunzip -c /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TC/QC_AFR/TC.nodup.gz | " + \
      "awk '!( ($4==\"A\" && $5==\"T\") || ($4==\"T\" && $5==\"A\") || ($4==\"G\" && $5==\"C\") || ($4==\"C\" && $5==\"G\")) {print}' | " + \
      "gzip > /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TC/QC_AFR/TC.QC.gz"
os.system(cmd)

cmd = "zcat /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TC/QC_AFR/TC.QC.gz | wc -l"
print("SNP left after removing ambiguous SNPs: ")
print(os.popen(cmd).read())

# SNP left after removing ambiguous SNPs: 
# 13943507