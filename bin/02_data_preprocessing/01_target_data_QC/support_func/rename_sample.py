# This code is to rename the sample file using python since sample files are all the same for all included chromosome.

# Author: Yujia Shi
# Date: 2021.11.21 

import os

for i in range(0, 23):
    chr = i+1
    if i != 22:
        cmd = "cp ~/Project/2021-11-10-individual_MR/dat/sample/ukb_imp_chr22_v3.sample ~/Project/2021-11-10-individual_MR/dat/sample/ukb_imp_chr" + str(chr) + "_v3.sample"
        os.system(cmd)