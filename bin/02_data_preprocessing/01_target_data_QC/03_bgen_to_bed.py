# Code to convert the bgen file to bed file
# 
# Author: Yujia Shi
# Date: 2021.11.21 

import os

DATAPATH = "/mnt/data/yujia/GWAS_ukbb_maf0.001_info0.8/"
os.chdir("/mnt/data/yujia/GWAS_ukbb_maf0.001_info0.8/") # set working directiory

for i in range(1, 23):
    
    chr = str(i)
    print("Converting chr " + chr + ".")
    
    cmd = "plink2 -bgen " + DATAPATH + "/chr" + chr + "_QC.bgen ref-first " + \
        "--sample " + DATAPATH + "/chr" + chr + "_QC.sample " + \
        "--out ./bed/chr" + chr + " --threads 30 " + \
        "--make-bed"
    os.system(cmd)
    
print("Finished the pipeline.")