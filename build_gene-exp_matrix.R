# BUILD_GENE-EXP_MATRIX.R
# LAST UPDATED BY LIAM FLINN SPURR ON SEPTEMBER 9, 2019
# THIS SCRIPT BUILDS THE GENE EXPRESSION MATRIX FOR REQTL ANALYSIS

# install missing required packages and load packages
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("data.table"); load_package("tidyverse")

# take arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# identify the location of the readcount files
# you may wish to modify the pattern to match the naming conventions of your samples
# the pattern identified here is the suffix of the samples that will later be removed to get the sample name 
exp_path <- args[1]
pattern <- "_gene_abund.tab"
exp_files <- list.files(path = exp_path, pattern = pattern)

# load in the readcounts
df <- bind_rows(lapply(exp_files, function(f) {
  # select the necessary columns, including gene position and TPM value
  # fix sample name formatting, and remove genes without an annotation
  fread(paste0(exp_path, f)) %>% 
    select(`Gene ID`,`Gene Name`, Reference, Start, End, TPM) %>% 
    mutate(sample = gsub(x = f, pattern = pattern, replacement = '')) %>% # remove sample name suffix specified in "pattern"
    filter(`Gene Name` != '-')
}))

# create gene coordinate table
df_loc <- df %>% select(`Gene ID`, Reference, Start, End) %>% 
  distinct() %>% 
  mutate(Reference = paste0('chr', Reference)) %>% 
  rename(GeneID = `Gene ID`)

df <- df %>% select(`Gene ID`, `Gene Name`, TPM, sample)

# get number of samples
num_samples <- length(levels(factor(df$sample)))

# retain only genes that are expressed at TPM >=1 in at least 20% of the samples
df <- df %>% group_by(`Gene ID`) %>% mutate(count_zero = length(which(TPM < 1)),
                                             perc_zero = count_zero / num_samples)

# remove duplicate entries for each gene if necessary
df <- df %>% filter(perc_zero < 0.8) %>% 
  select(sample, TPM, `Gene ID`) %>% 
  group_by(sample, `Gene ID`) %>% 
  filter(row_number() == 1) %>%
  spread(sample, TPM) %>% 
  rename(GeneID = `Gene ID`)


# quantile normalize gene expression values
# (code adapted from MatrixEQTL sample code)
df_w <- df[-1]

for(gene in 1:nrow(df_w)) {
  mat = df_w[gene,]
  mat = apply(mat, 1, rank, ties.method = "average")
  mat = qnorm(mat / (ncol(df_w) + 1))
  df_w[gene,] = mat
}
df_w$gene_id <- df$GeneID
df_w <- df_w %>% select(gene_id, everything())

# write outputs
output_prefix <- args[2]
write.table(df_w, paste0(output_prefix, '_gene-exp_matrix.txt'), quote = F, row.names = F, sep = '\t')
write.table(df_loc, paste0(output_prefix, '_gene-exp-loc_matrix.txt'), quote = F, row.names = F, sep = '\t')
