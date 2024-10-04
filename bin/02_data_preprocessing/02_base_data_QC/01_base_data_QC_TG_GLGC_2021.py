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

tg = pd.read_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG/download/without_UKB_logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz", sep="\t", compression = 'gzip')

print("Total SNP in the GWAS data: ")
print(tg.shape)

# drop first column since it is from hg18, 2nd column is hg19
tg = tg.iloc[:, 1:]

# split first column into two column CHR and BP
tg = pd.concat([tg['SNP_hg19'].str.split(':', expand = True), tg.iloc[:, 1:]], axis=1)
tg.rename(columns={0:"CHR", 1:"BP", "rsid":"SNP", "Freq.A1.1000G.EUR":"EAF", "beta":"BETA", "se":"SE", "P-value":"P"}, inplace=True)

# convert EAF to MAF
tg['MAF'] = tg['EAF'].apply(lambda x : x if x <= 0.5 else 1-x)

# Uppercase the A1 and A2
tg['A1'] = tg['A1'].str.upper()
tg['A2'] = tg['A2'].str.upper()

# N row to integer
tg.N = tg.N.astype(int)
tg.CHR = tg.CHR.str.extract(r'(?<=chr)(.*)')

# reorder the column
cols = ["CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "BETA", "MAF", "EAF"]
tg = tg[cols]

# Remove NaN row
tg = tg.dropna()

# convert to correct column type
tg = tg.astype({'CHR':'int', 'BP':'int', "SNP":'str', 'A1':'str', 'A2':'str', 'N':'int', 'SE':'float', 'P':'float', 'BETA':'float', 'MAF':'float', 'EAF':'float'})

# remove SNP without SNPID
tg = tg[tg['SNP'].str.startswith('rs')]

# save the formatted results
tg.to_csv("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/TG_formated_nadrop.tsv.gz", sep="\t", index=False)

print("SNP left after basic QC: ")
print(tg.shape) # (34823002, 11)

# =============================
# Step 2: remove SNPs with MAF < 0.01
# =============================

cmd = "gunzip -c /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/TG_formated_nadrop.tsv.gz |" + \
      "awk 'NR==1 || ($10 > 0.01) {print}' |" + \
      "gzip > /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/TG.standardGWASQC.gz"
os.system(cmd)

cmd = "zcat /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/TG.standardGWASQC.gz | wc -l"
print("SNP left after removing SNPs with MAF < 0.01: ") # 9865674
print(os.popen(cmd).read())

# =============================
# Step 3: remove duplicate SNPS
# =============================

cmd = "gunzip -c /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/TG.standardGWASQC.gz |" + \
      "awk '{seen[$3]++; if(seen[$3]==1){ print}}' |" + \
      "gzip - > /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/TG.nodup.gz"
os.system(cmd)

cmd = "zcat /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/TG.nodup.gz | wc -l"
print("SNP left after removing duplicate SNPs: ") # 9842222
print(os.popen(cmd).read())

# =============================
# Step 4: keep non-ambiguous SNPS
# =============================

cmd = "gunzip -c /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/TG.nodup.gz | " + \
      "awk '!( ($4==\"A\" && $5==\"T\") || ($4==\"T\" && $5==\"A\") || ($4==\"G\" && $5==\"C\") || ($4==\"C\" && $5==\"G\")) {print}' | " + \
      "gzip > /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/TG.QC.gz"
os.system(cmd)

cmd = "zcat /mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/TG_GLGC_2021/QC/TG.QC.gz | wc -l"
print("SNP left after removing ambiguous SNPs: ") # 8484277
print(os.popen(cmd).read())
