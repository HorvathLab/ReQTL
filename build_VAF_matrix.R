# BUILD_VAF_MATRIX.R
# LAST UPDATED BY LIAM FLINN SPURR ON SEPTEMBER 9, 2019
# THIS SCRIPT BUILDS THE VAF MATRIX FOR REQTL ANALYSIS

# install missing required packages and load packages
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("dplyr"); load_package("tidyr"); load_package("data.table");

# take arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# identify the location of the readcount files
snv_path <- args[1]
snv_files <- list.files(path = snv_path, pattern = 'csv')

# load in the readcounts
df <- bind_rows(lapply(snv_files, function(f) {
  # create a temporary file for each individual readcounts csv
  temp <- fread(paste0(snv_path, f))
  
  # select the necessary columns and format the sample column
  temp <- temp %>% select(CHROM, POS, REF, ALT, AlignedReads, R) %>% 
    mutate(AlignedReads = gsub(x = AlignedReads, pattern = '\\.sorted', replacement = ''))
}))

# get number of samples and create SNP identifier
num_samples <- length(levels(factor(df$AlignedReads)))
df <- df %>% mutate(SNP = paste0(CHROM, ':', POS, '_', REF, '>', ALT))

# remove variants that are homozygous variant in more than 80% of the samples
df <- df %>% group_by(SNP) %>% mutate(count_homo_var = length(which(R == 1)),
                                     perc_homo_var = count_homo_var / num_samples) %>% ungroup()
df <- df %>% filter(perc_homo_var < 0.8)

# remove variants that are homozygous reference in more than 80% of the samples
df <- df %>% group_by(SNP) %>% mutate(count_homo_ref = length(which(R == 0)),
                                     perc_homo_ref = count_homo_ref / num_samples) %>% ungroup()
df <- df %>% filter(perc_homo_ref < 0.8)

# remove variants that are NA in more than 80% of the samples
df <- df %>% group_by(SNP) %>% mutate(count = n(),
                                      perc_non_na = count / num_samples) %>% ungroup()
df <- df %>% filter(perc_non_na > 0.2)

# select relevant columns and create SNP location matrix
df_loc <- df %>% select(SNP, CHROM, POS) %>% distinct() %>% mutate(CHROM = paste0('chr', CHROM))
df_loc <- as.matrix(df_loc)

# remove excess columns
df <- df %>% select(SNP, AlignedReads, R)

# build the SNP matrix
df <- df %>% spread(AlignedReads, R)
df <- as.matrix(df)

# write outputs
output_prefix <- args[2]
write.table(df, paste0(output_prefix, '_VAF_matrix.txt'), quote = F, row.names = F, sep = '\t')
write.table(df_loc, paste0(output_prefix, '_VAF-loc_matrix.txt'), quote = F, row.names = F, sep = '\t')
