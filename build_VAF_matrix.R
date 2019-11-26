# BUILD_VAF_MATRIX.R
# LAST UPDATED BY LIAM FLINN SPURR ON NOVEMBER 26, 2019

# load packages
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(data.table, quietly = TRUE))

handle_command_args <- function(args) {
  # make sure all flags are paired
  if(length(args) %% 2 != 0) stop("Command line arguments supplied incorrectly!")
  
  # load the flags into a "dictionary"
  arg_df <- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>%
    mutate_all(as.character)
  
  # identify the location of the readcount files
  snv_path <<- arg_df$value[arg_df$flag == "-r"]
  snv_files <<- list.files(path = snv_path, pattern = 'csv')
  
  # specify output prefix
  output_prefix <<- arg_df$value[arg_df$flag == "-o"]
  
}

# read arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
handle_command_args(args)

# load in the readcounts
df <- bind_rows(lapply(snv_files, function(f) {
  # create a temporary file for each individual readcounts csv
  temp <- fread(paste0(snv_path, f))
  
  # select the necessary columns and format the sample column
  temp <- temp %>% select(CHROM, POS, REF, ALT, AlignedReads, R) %>% 
    mutate(AlignedReads = gsub(x = AlignedReads, pattern = '\\.sorted', replacement = ''))
}))

# get number of samples and create SNV identifier
num_samples <- length(levels(factor(df$AlignedReads)))
df <- df %>% mutate(SNV = paste0(CHROM, ':', POS, '_', REF, '>', ALT))

# remove variants that are homozygous variant in more than 80% of the samples
df <- df %>% group_by(SNV) %>% mutate(count_homo_var = length(which(R == 1)),
                                     perc_homo_var = count_homo_var / num_samples) %>% ungroup()
df <- df %>% filter(perc_homo_var < 0.8)

# remove variants that are homozygous reference in more than 80% of the samples
df <- df %>% group_by(SNV) %>% mutate(count_homo_ref = length(which(R == 0)),
                                     perc_homo_ref = count_homo_ref / num_samples) %>% ungroup()
df <- df %>% filter(perc_homo_ref < 0.8)

# remove variants that are NA in more than 80% of the samples
df <- df %>% group_by(SNV) %>% mutate(count = n(),
                                      perc_non_na = count / num_samples) %>% ungroup()
df <- df %>% filter(perc_non_na > 0.2)

# select relevant columns and create SNV location matrix
df_loc <- df %>% select(SNV, CHROM, POS) %>% distinct() %>% mutate(CHROM = paste0('chr', CHROM))
df_loc <- as.matrix(df_loc)

# remove excess columns
df <- df %>% select(SNV, AlignedReads, R)

# build the SNV matrix
df <- df %>% spread(AlignedReads, R)
df <- as.matrix(df)

# write outputs
if (!dir.exists("output")) {
  cat('Creating output directory...\n')
  dir.create('output')
} 

write.table(df, paste0("output/", output_prefix, '_VAF_matrix.txt'), quote = F, row.names = F, sep = '\t')
write.table(df_loc, paste0("output/", output_prefix, '_VAF-loc_matrix.txt'), quote = F, row.names = F, sep = '\t')
cat(paste0("VAF matrix and locations for MatrixEQTL saved to output/", output_prefix, "_VAF_matrix.txt and output/", output_prefix, "_VAF-loc_matrix.txt.\n"))
