# HARMONIZE_MATRICES.R
# LAST UPDATED BY LIAM FLINN SPURR ON SEPTEMBER 9, 2019
# THIS SCRIPT HARMONIZES THE SAMPLES PRESENT IN THE INPUT MATRICES FOR REQTL ANALYSIS

# install missing required packages and load packages
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("tidyverse"); load_package("data.table")

# read arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# load in matrices
vaf <- data.frame(fread(args[1]))
exp <- data.frame(fread(args[2]))
cov <- data.frame(fread(args[3]))

# get the samples present in each
v.names <- colnames(vaf)[-1]
e.names <- colnames(exp)[-1]
c.names <- colnames(cov)[-1]

# get the samples in all 3 matrices
conc <- intersect(intersect(v.names, e.names), c.names)

# select only the concordant samples from all 3 matrices
vaf.clean <- vaf[,c("SNV", conc)]
exp.clean <- exp[,c("gene_id", conc)]
cov.clean <- cov[,c("id", conc)]

# write the fixed outputs to a file
write.table(vaf.clean, paste0(gsub(".txt", "", args[1]), "_harmonized.txt"), quote = F, row.names = F, sep = '\t')
write.table(exp.clean, paste0(gsub(".txt", "", args[2]), "_harmonized.txt"), quote = F, row.names = F, sep = '\t')
write.table(cov.clean, paste0(gsub(".txt", "", args[3]), "_harmonized.txt"), quote = F, row.names = F, sep = '\t')
