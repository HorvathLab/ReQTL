# BUILD_COV_MATRIX.R
# LAST UPDATED BY LIAM FLINN SPURR ON SEPTEMBER 9, 2019
# THIS SCRIPT BUILDS THE COVARIATE MATRIX FOR REQTL ANALYSIS

# install missing required packages and load packages
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("tidyverse"); load_package("data.table"); load_package("factoextra")

# take arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# load in VAF matrix
genotypes <- fread(args[1])

# remove sample identifier for PCA transformation
genotypes_w <- genotypes %>% select(-1)

# load in covariates
cov <- fread(args[2]) %>% filter(Run %in% colnames(genotypes))

# get the first three principal components of the sample genotypes
pca <- prcomp(data = genotypes_w, ~., na.action = 'na.omit', scale = T)
pcs <- pca$x
pcs <- pca$rotation
pcs <- pcs[,1:3]
pcs <- pcs %>% as.data.frame()
pcs$Run <- row.names(pcs)

# join with other covariates
cov <- left_join(cov, pcs, by = "Run") %>% select(-2,-3)

# save the sample ids
n <- cov$Run

# transpose all but the Run column and convert to a matrix
cov_t <- as.data.frame(t(cov[,-1]))
colnames(cov_t) <- n
cov_t$id <- names(cov[-1])
cov_t <- cov_t %>% select(id, everything()) %>% as.matrix()

# write the output to a file
output_prefix <- args[3]
write.table(cov_t, paste0(output_prefix, '_cov_matrix.txt'), quote = F, row.names = F, sep = '\t')
