# Code to merge all chromosome's bed file to a single bed file
# 
# Author: Yujia Shi
# Date: 2021.11.26

nohup plink --bed chr1.bed --bim chr1.bim --fam chr1.fam --merge-list list_beds.txt --make-bed --out ./allchrs/ukb_allchrs &

