# LIAM SPURR
# HORVATH LAB
# 24MAY18
# BUILD_SNV_MATRIX.R

# install missing required packages and load packages
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("dplyr"); load_package("tidyr"); load_package("data.table");

# get the command line arguments
args = commandArgs(trailingOnly = T)

# identify the location of the readcount files
snv_path = args[1]
snv_files = list.files(path = snv_path, pattern = 'csv')

# create an empty data frame to hold the combined readcounts
df = data.frame()

# load in the readcounts
for (file in snv_files) {
  # create a temporary file for each individual readcounts csv
  temp = fread(paste0(snv_path, file))
  
  # select the necessary columns and format the sample column
  temp = temp %>% select(1:5, 14) %>% mutate(AlignedReads = gsub(x = AlignedReads, pattern = '\\.sorted', replacement = ''))
  
  # combine all of the files together
  df = rbind(df, temp)
}

# create SNP identifier
df = df %>% mutate(SNP = paste0(CHROM, ':', POS, '_', REF, '>', ALT))

# create the SNP location file
df_loc = df %>% select(SNP, CHROM, POS) %>% distinct() %>% mutate(CHROM = paste0('chr', CHROM))

# remove excess columns
df = df %>% select(SNP, AlignedReads, R)

# build the SNP matrix
df = df %>% spread(AlignedReads, R)

# write the matrices to a file
write.table(df, args[2], quote = F, row.names = F, sep = '\t')
write.table(df_loc, args[3], quote = F, row.names = F, sep = '\t')
