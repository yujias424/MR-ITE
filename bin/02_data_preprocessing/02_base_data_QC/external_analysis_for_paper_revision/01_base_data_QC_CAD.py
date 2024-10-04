# base data quality control
# following Choi's QC pipeline
# 
# Author: Yujia Shi
# Date: 2023.02.08 

import os
import re

import numpy as np
import pandas as pd

os.chdir("/home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/")

# I have already use GENOMA to check the heritability 
# The pipeline starts from standard GWAS QC

# =======================
# Step 1: format the data
# =======================

CAD = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/download/26343387-GCST003116-EFO_0000378-build37.f.tsv.gz", sep="\t", compression = 'gzip')

print("Total SNP in the GWAS data: ")
print(CAD.shape)

# Total SNP in the GWAS data: 
# (9455778, 13)

# rename the data
CAD.rename(columns={"chromosome":"CHR", "base_pair_location":"BP", "variant_id":"SNP", "effect_allele":"A1", "other_allele":"A2", "effect_allele_frequency":"EAF", "beta":"BETA", "standard_error":"SE", "p_value":"P"}, inplace=True)

# convert EAF to MAF
CAD['MAF'] = CAD['EAF'].apply(lambda x : x if x <= 0.5 else 1-x)

# Uppercase the A1 and A2
CAD['A1'] = CAD['A1'].str.upper()
CAD['A2'] = CAD['A2'].str.upper()

# reorder the column
cols = ["CHR", "BP", "SNP", "A1", "A2", "SE", "P", "BETA", "MAF", "EAF"]
CAD = CAD[cols]

# Remove NaN row
CAD = CAD.dropna()

# convert to correct column type
CAD = CAD.astype({'CHR':'int', 'BP':'int', "SNP":'str', 'A1':'str', 'A2':'str', 'SE':'float', 'P':'float', 'BETA':'float', 'MAF':'float', 'EAF':'float'})

# remove SNP without SNPID
CAD = CAD[CAD['SNP'].str.startswith('rs')]

# save the formatted results
CAD.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/QC/CAD_formated_nadrop.tsv.gz", sep="\t", index=False)

print("SNP left after basic QC: ")
print(CAD.shape)

# SNP left after basic QC: 
# (8622850, 10)

# =============================
# Step 2: remove SNPs with MAF < 0.01
# =============================

cmd = "gunzip -c /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/QC/CAD_formated_nadrop.tsv.gz |" + \
      "awk 'NR==1 || ($9 > 0.01) {print}' |" + \
      "gzip > /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/QC/CAD.standardGWASQC.gz"
os.system(cmd)

cmd = "zcat /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/QC/CAD.standardGWASQC.gz | wc -l"
print("SNP left after removing SNPs with MAF < 0.01: ")
print(os.popen(cmd).read())

# SNP left after removing SNPs with MAF < 0.01: 
# 8390372

# =============================
# Step 3: remove duplicate SNPS
# =============================

cmd = "gunzip -c /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/QC/CAD.standardGWASQC.gz |" + \
      "awk '{seen[$3]++; if(seen[$3]==1){ print}}' |" + \
      "gzip - > /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/QC/CAD.nodup.gz"
os.system(cmd)

cmd = "zcat /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/QC/CAD.nodup.gz | wc -l"
print("SNP left after removing duplicate SNPs: ")
print(os.popen(cmd).read())

# SNP left after removing duplicate SNPs: 
# 8390372

# =============================
# Step 4: keep non-ambiguous SNPS
# =============================

cmd = "gunzip -c /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/QC/CAD.nodup.gz | " + \
      "awk '!( ($4==\"A\" && $5==\"T\") || ($4==\"T\" && $5==\"A\") || ($4==\"G\" && $5==\"C\") || ($4==\"C\" && $5==\"G\")) {print}' | " + \
      "gzip > /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/QC/CAD.QC.gz"
os.system(cmd)

cmd = "zcat /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/CAD/QC/CAD.QC.gz | wc -l"
print("SNP left after removing ambiguous SNPs: ")
print(os.popen(cmd).read())

# SNP left after removing ambiguous SNPs: 
# 7092561