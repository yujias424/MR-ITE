# https://github.com/privefl/bigstatsr/issues/90

suppressPackageStartupMessages({
    library(data.table)
    library(bigsnpr)
    library(dplyr)
})

# =================================================
# Step 1: Read in the phenotype and covariate files
# =================================================
phenotype <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/external_analysis_for_paper_revision/IGF1/CAD/ukbb.phenotype.IGF1.nmolL")
covariate <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/04_pheno_covar_data/external_analysis_for_paper_revision/IGF1/CAD/ukbb.covariate.prs.IGF1.set1")
pheno <- merge(phenotype, covariate)

# ===========================
# Step 2: Obtain HapMap3 SNPs
# ===========================
info <- readRDS("/mnt/md0/public_data/Genodata/HapMap3/map_hm3_ldpred2.rds") # another option: map_hm3_plus.rds

# =====================================================
# Step 3: Load and transform the summary statistic file
# =====================================================

# get valid SNPs (We first try P-value 0.1)
valid.snps <- fread("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/external_analysis_for_paper_revision/LDPred2/EUR/IGF1/CAD/selected_snps_0.2.csv")

# Load and transform the summary statistic file
# Read in the summary statistic file
sumstats <- bigreadr::fread2("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/02_base_data/IGF1/QC/IGF1.QC.withPos.gz") 

340242/2077390

# LDpred 2 require the header to follow the exact naming
names(sumstats) <-
    c("chr",
    "pos",
    "rsid",
    "a1",
    "a0",
    "n_eff",
    "beta_se",
    "p",
    "beta",
    "MAF",
    "EAF")

print(dim(sumstats))

# Filter out hapmap SNPs
sumstats <- sumstats[sumstats$rsid%in% info$rsid,]
print(dim(sumstats))

# Select valid SNPs
sumstats <- sumstats[sumstats$rsid %in% valid.snps$SNP, ]
print(dim(sumstats))

# ===============================
# Step 4: Calculate the LD matrix
# ===============================

# Get maximum amount of cores
NCORES <- 40

# Open a temporary file
tmp <- tempfile(tmpdir = "/mnt/md1/")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)

# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file 
info_snp <- NULL
fam.order <- NULL

# now attach the genotype object
obj.bigSNP <- snp_attach(paste0("/mnt/md1/share/ukbb_geno_rds_allchr/ukb_allchrs.rds"))
valid.snps.ind <- rownames(obj.bigSNP$map[obj.bigSNP$map$marker.ID %in% valid.snps$SNP, ])
valid.snps.ind <- as.integer(valid.snps.ind)
valid.samps.ind <- rownames(obj.bigSNP$fam[obj.bigSNP$fam$sample.ID %in% pheno$IID, ])
valid.samps.ind <- as.integer(valid.samps.ind)
obj.bigSNP.select <- snp_subset(obj.bigSNP, ind.col = valid.snps.ind, ind.row = valid.samps.ind)

# now attach the genotype object
obj.bigSNP <- snp_attach(obj.bigSNP.select)
# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
# perform SNP matching
info_snp <- snp_match(sumstats, map, match.min.prop=0)

# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes
# Rename the data structures
CHR <- map$chr
POS <- map$pos
# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = "/mnt/md1/share/genetic_map/", ncores = 40)

# calculate LD
for (chr in 1:22) {
    print(chr)
    # Extract SNPs that are included in the chromosome
    ind.chr <- which(info_snp$chr == chr)
    ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
    # Calculate the LD
    corr0 <- snp_cor(
            genotype,
            ind.col = ind.chr2,
            ncores = NCORES,
            infos.pos = POS2[ind.chr2],
            size = 3 / 1000
        )
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}

# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,
        c("family.ID", "sample.ID"),
        c("FID", "IID"))

# ===================================
# Step 5: Perform LD score regression
# ===================================

df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(ld, 
                 length(ld), 
                 chi2 = (df_beta$beta / df_beta$beta_se)^2,
                 sample_size = df_beta$n_eff, 
                 blocks = NULL)
h2_est <- ldsc[["h2"]]

# ===================================
# Step 6: Calculate the null R2
# ===================================

# Reformat the phenotype file such that y is of the same order as the 
# sample ordering in the genotype file
pheno <- pheno[pheno$IID %in% fam.order$IID, ]
y <- pheno[fam.order, on = c("FID", "IID")]
colnames(y)[3] <- "IGF1"

# Calculate the null R2
# use glm for binary trait 
# (will also need the fmsb package to calculate the pseudo R2)
null.model <- paste("PC", 1:10, sep = "", collapse = "+") %>%
    paste0("IGF1~Gender+Genotype_Batch+Age+", .) %>%
    as.formula %>%
    lm(., data = y) %>%
    summary
null.r2 <- null.model$r.squared
null.r2

# ===================================
# Step 7: Obtain model PRS
# ===================================

# Get adjusted beta from the auto model
multi_auto <- snp_ldpred2_auto(
    corr,
    df_beta,
    h2_init = h2_est,
    vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
    ncores = NCORES
)
beta_auto <- sapply(multi_auto, function(auto)
    auto$beta_est)
print(multi_auto[[1]]$p_est)

# calculate PRS for all samples
ind.test <- 1:nrow(genotype)
genotype <- snp_fastImputeSimple(genotype, ncores = 35)
pred_auto <-
    big_prodMat(genotype,
                beta_auto,
                ind.col = info_snp$`_NUM_ID_`,
                ncores = 35)

# scale the PRS generated from AUTO
pred_scaled <- apply(pred_auto, 2, sd)
final_beta_auto <-
    rowMeans(beta_auto[,
                abs(pred_scaled -
                    median(pred_scaled)) <
                    3 * mad(pred_scaled)])
pred_auto <-
    big_prodVec(genotype,
        final_beta_auto,
        ind.col = info_snp$`_NUM_ID_`)

PRS_final <- data.frame("IID" = fam.order$IID, "PRS" = pred_auto)
head(pred_auto)

dir.create("/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/external_analysis_for_paper_revision/LDPred2/EUR/IGF1/CAD/set1", recursive = T)
fwrite(PRS_final, "/mnt/md0/yujia/project/2023-07-20-individual_MR/dat/06_PRS_calculation/external_analysis_for_paper_revision/LDPred2/EUR/IGF1/CAD/set1/IGF1_prs_0.2.csv")

# ======================================================
# Step 8: Get the final performance of the LDpred models
# ======================================================

reg.formula <- paste("PC", 1:10, sep = "", collapse = "+") %>%
    paste0("IGF1~PRS+Gender+Genotype_Batch+Age+", .) %>%
    as.formula
reg.dat <- y
reg.dat$PRS <- pred_auto
auto.model <- lm(reg.formula, dat=reg.dat) %>%
    summary
print(auto.model)
(result <- data.table(
    auto = auto.model$r.squared - null.r2,
    null = null.r2
))
print(result)
