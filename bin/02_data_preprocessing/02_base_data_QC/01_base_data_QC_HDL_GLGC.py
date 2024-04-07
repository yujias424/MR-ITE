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

hdl = pd.read_csv("~/Project/2022-09-01-individual_MR/dat/02_base_data/HDL/download/jointGwasMc_HDL.txt.gz", sep="\t", compression = 'gzip')

print("Total SNP in the GWAS data: ")
print(hdl.shape)

# drop first column since it is from hg18, 2nd column is hg19
hdl = hdl.iloc[:, 1:]

# split first column into two column CHR and BP
hdl = pd.concat([hdl['SNP_hg19'].str.split(':', expand = True), hdl.iloc[:, 1:]], axis=1)
hdl.rename(columns={0:"CHR", 1:"BP", "rsid":"SNP", "Freq.A1.1000G.EUR":"EAF", "beta":"BETA", "se":"SE", "P-value":"P"}, inplace=True)

# convert EAF to MAF
hdl['MAF'] = hdl['EAF'].apply(lambda x : x if x <= 0.5 else 1-x)

# Uppercase the A1 and A2
hdl['A1'] = hdl['A1'].str.upper()
hdl['A2'] = hdl['A2'].str.upper()

# N row to integer
hdl.N = hdl.N.astype(int)
hdl.CHR = hdl.CHR.str.extract(r'(?<=chr)(.*)')

# reorder the column
cols = ["CHR", "BP", "SNP", "A1", "A2", "N", "SE", "P", "BETA", "MAF", "EAF"]
hdl = hdl[cols]

# Remove NaN row
hdl = hdl.dropna()

# convert to correct column type
hdl = hdl.astype({'CHR':'int', 'BP':'int', "SNP":'str', 'A1':'str', 'A2':'str', 'N':'int', 'SE':'float', 'P':'float', 'BETA':'float', 'MAF':'float', 'EAF':'float'})

# remove SNP without SNPID
hdl = hdl[hdl['SNP'].str.startswith('rs')]

# save the formatted results
hdl.to_csv("~/Project/2022-09-01-individual_MR/dat/02_base_data/HDL/QC/HDL_GLGC_formated_nadrop.tsv.gz", sep="\t", index=False)

print("SNP left after basic QC: ")
print(hdl.shape)

# =============================
# Step 2: remove SNPs with MAF < 0.01
# =============================

cmd = "gunzip -c /home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/HDL/QC/HDL_GLGC_formated_nadrop.tsv.gz |" + \
      "awk 'NR==1 || ($10 > 0.01) {print}' |" + \
      "gzip > /home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/HDL/QC/HDL.GLGC.standardGWASQC.gz"
os.system(cmd)

cmd = "zcat /home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/HDL/QC/HDL.GLGC.standardGWASQC.gz | wc -l"
print("SNP left after removing SNPs with MAF < 0.01: ")
print(os.popen(cmd).read())

# =============================
# Step 3: remove duplicate SNPS
# =============================

cmd = "gunzip -c /home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/HDL/QC/HDL.GLGC.standardGWASQC.gz |" + \
      "awk '{seen[$3]++; if(seen[$3]==1){ print}}' |" + \
      "gzip - > /home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/HDL/QC/HDL.GLGC.nodup.gz"
os.system(cmd)

cmd = "zcat /home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/HDL/QC/HDL.GLGC.nodup.gz | wc -l"
print("SNP left after removing duplicate SNPs: ")
print(os.popen(cmd).read())

# =============================
# Step 4: keep non-ambiguous SNPS
# =============================

cmd = "gunzip -c /home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/HDL/QC/HDL.GLGC.nodup.gz | " + \
      "awk '!( ($4==\"A\" && $5==\"T\") || ($4==\"T\" && $5==\"A\") || ($4==\"G\" && $5==\"C\") || ($4==\"C\" && $5==\"G\")) {print}' | " + \
      "gzip > /home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/HDL/QC/HDL.GLGC.QC.gz"
os.system(cmd)

cmd = "zcat /home/yujia/Project/2022-09-01-individual_MR/dat/02_base_data/HDL/QC/HDL.GLGC.QC.gz | wc -l"
print("SNP left after removing ambiguous SNPs: ")
print(os.popen(cmd).read())
