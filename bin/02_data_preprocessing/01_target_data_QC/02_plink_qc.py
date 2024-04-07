# target data quality control following Choi's QC pipeline
# 
# Author: Yujia
# Date: 2021.11.21 

import os
import subprocess

# Initial the parameters
MIN_MAF = 0.001 # we pick 0.001 since we have a large UKBB dataset, as mentioned by Choi in his tutorial.
MIN_IS = 0.8
N_JOB_THREAD = 5
N_PLINK_THREAD = 30

DATAPATH = "/home/sohc/UKBB2"
kennDATAPATH = "/mnt/data/share/UKBB/genotype/imp_info0.3_maf0.001"
step1DATAPATH = "/home/yujia/Project/2022-09-01-individual_MR/dat/01_target_data/GWAS/GWAS_dat"
sampleDATAPATH = "/home/yujia/Project/2022-09-01-individual_MR/dat/01_target_data/sample"

# Find the SNP ID with MAF > MIN_MAF (0.001) and INFO_SCORE > MIN_IS (0.8)
for i in range(0, 23):
    chr = i + 1
    maf = "MAF=" + str(MIN_MAF)
    is_param = "IS=" + str(MIN_IS)
    # mfi_path = "/mnt/data/share/UKBB/genotype/imp_info0.3_maf0.001/mfi/" + "ukb_mfi_chr" + str(chr) + "_v3.txt"
    mfi_path = "/home/yujia/Project/2022-09-01-individual_MR/dat/01_target_data/mfi/" + "ukb_mfi_chr" + str(chr) + "_v3.txt"
    output_mfi_path = "/home/yujia/Project/2022-09-01-individual_MR/dat/01_target_data/GWAS/mfi/"  + "chr" + str(chr) + "_snpid_maf" + str(MIN_MAF) + "_is" + str(MIN_IS) + ".txt"

    cmd = "awk -v OFS='\t' -v " + is_param + " -v " + maf + \
            " '$6>MAF && $8>IS {print $2}' " + \
                mfi_path + " > " + output_mfi_path + ";"
    os.system(cmd) # cmd code is from Kenneth.

# QC using PLINK
os.chdir("/home/yujia/Project/2022-09-01-individual_MR/dat/01_target_data/GWAS/GWAS_dat/") # set working directiory
for i in range(0, 23):

    os.chdir("/home/yujia/Project/2022-09-01-individual_MR/dat/01_target_data/GWAS/GWAS_dat/")
    chr = str(i + 1)
    output_mfi_path = "/home/yujia/Project/2022-09-01-individual_MR/dat/01_target_data/GWAS/mfi/"  + "chr" + str(chr) + "_snpid_maf" + str(MIN_MAF) + "_is" + str(MIN_IS) + ".txt"

    # Step 1 
    # select score based on the MAF, Imputation 'Info score' and hardy-weinberg cutoff
    cmd = "plink2 -bgen " + DATAPATH + "/ukb_imp_chr" + chr + "_v3.bgen ref-first --threads 30 " + \
            "--sample " + sampleDATAPATH + "/ukb_imp_chr" + chr + "_v3.sample " + \
            "--extract " + output_mfi_path + " " + \
            "--hwe 1e-6 " + \
            "--geno 0.01 " + \
            "--rm-dup force-first --write-snplist --make-just-fam " + \
            "--out chr" + chr + "_QC " + \
            "--export bgen-1.2"

    os.system(cmd)
    print("Finished Step1.\n")

    # Step 2
    # Peform pruning to remove highly correlated SNPs
    cmd = "plink2 -bgen " + step1DATAPATH + "/chr" + chr + "_QC.bgen ref-first " + \
            "--sample " + step1DATAPATH + "/chr" + chr + "_QC.sample " + \
            "--keep chr" + chr + "_QC.fam " + \
            "--extract chr" + chr + "_QC.snplist " + \
            "--indep-pairwise 200 50 0.25 --threads 30 " + \
            "--out chr" + chr + "_QC"
    os.system(cmd)
    print("Finished Step2.\n")
    # QC Tutorial recommended param: 50 5 0.2 + PRSice recommended param: 200 50 0.25

    # Step 3
    # Calculate the Heterozygosity rates using Plink
    cmd = "plink2 -bgen " + step1DATAPATH +  "/chr" + chr + "_QC.bgen ref-first " + \
        "--sample " + step1DATAPATH + "/chr" + chr + "_QC.sample " + \
        "--keep chr" + chr + "_QC.fam " + \
        "--extract chr" + chr + "_QC.prune.in " + \
        "--het --threads 30 " + \
        "--out chr" + chr + "_QC"
    os.system(cmd)
    print("Finished Step3.\n")

    # =================================================== temporarily deprecated ==============================================================
    # Step 4 (deprecated as we have the family sample information from the UKBB dataset.)
    # Check Relatedness with cutoff 0.125 (--rel-cutoff), if plink2 is applied, use --king-cutoff instead with cutoff setting 0.0884
    # cmd =  "plink2 -bgen " + step1DATAPATH + "/chr" + chr + "_maf" + str(MIN_MAF) + "_is" + str(MIN_IS) + ".bgen ref-first --threads 30 " + \
    #     "--sample " + step1DATAPATH + "/chr" + chr + "_maf" + str(MIN_MAF) + "_is" + str(MIN_IS) + ".sample " + \
    #     "--extract chr" + chr + "_QC.prune.in " + \
    #     "--remove family.sample " + \
    #     "--out chr" + chr + "_QC_king"
    # os.system(cmd)
    # print("Finished Step4.")
    # =========================================================================================================================================

    # Step 4
    # Generate final QC'ed target data file (The third degree family were directly removed from the dataset.)
    cmd =  "plink2 -bgen " + step1DATAPATH + "/chr" + chr + "_QC.bgen ref-first " + \
            "--sample " + step1DATAPATH + "/chr" + chr + "_QC.sample " + \
            "--remove /home/yujia/Project/2022-09-01-individual_MR/dat/01_target_data/relatedness_kinship/relatedness_invalid.sample " + \
            "--extract " + step1DATAPATH  + "/chr" + chr + "_QC.prune.in " + \
            "--out chr" + chr + "_QC --threads 30 " + \
            "--export bgen-1.2"
    os.system(cmd)
    print("Finished Step4.\n")

    # Step 5
    # remove patients with high heterozygosity rate
    # Using Rscript to get the sample match the criterion.
    # remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean
    cmd = "Rscript /home/yujia/Project/2022-09-01-individual_MR/bin/02_data_preprocessing/01_target_data_QC/support_func/het_filter.R" # this file will generate the QC_valid.sample and QC_invalid.sample file
    os.system(cmd)
    
    os.chdir("/mnt/data/yujia/GWAS_ukbb_maf0.001_info0.8")
    cmd =  "plink2 -bgen " + step1DATAPATH + "/chr" + chr + "_QC.bgen ref-first " + \
            "--sample " + step1DATAPATH + "/chr" + chr + "_QC.sample " + \
            "--remove " + step1DATAPATH + "/chr" + chr + "_QC_invalid.sample " + \
            "--extract " + step1DATAPATH  + "/chr" + chr + "_QC.prune.in " + \
            "--out chr" + chr + "_QC --threads 30 " + \
            "--export bgen-1.2"
    os.system(cmd)
    print("Finished Step5.\n")

    print("Finished chromosome " + chr + ".\n")
    print("====================================================================================================================================")
