# LIAM SPURR
# HORVATH LAB
# 24MAY18
# BUILD_GENE-EXP_MATRIX.R

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

# identify the location of the stringtie gene expression files
gene_abund_path = args[1]
abund_files = list.files(path = gene_abund_path)

# create an empty data frame to hold the combined gene expresion files
df = data.frame()

# load in the readcounts
for (file in abund_files) {
  # create a temporary file for each individual gene expression file
  temp = fread(paste0(gene_abund_path, file))
  
  # select the necessary columns and create the sample column
  temp = temp %>% select(1,2,3,5,6,9) %>% 
    mutate(sample = !!file) %>%
    filter(`Gene Name` != '-')
  
  # combine all of the files together
  df = rbind(df, temp)
}

# remove spaces from column name
df = df %>% rename(GeneID = `Gene ID`)

# create gene location file
df_loc = df %>% select(1,3,4,5) %>% distinct() %>% mutate(Reference = paste0('chr', Reference))

# remove duplicate identifiers
df = df %>% select(sample, TPM, `Gene ID`) %>% group_by(sample, `Gene ID`) %>% filter(row_number() == 1)

# spread the file into a matrix
df = df %>% spread(sample, TPM)

# transform gene expression values into quantile normalized form
df_w = df[-1]
for(gene in 1:nrow(df_w)) {
  mat = df_w[gene,]
  mat = apply(mat, 1, rank, ties.method = "average")
  mat = qnorm(mat / (ncol(df_w) + 1))
  df_w[gene,] = mat
}
df_w$gene_id = df$GeneID
df_w = df_w %>% select(gene_id, everything())

# write outputs to files
write.table(df_w, args[2], quote = F, row.names = F, sep = '\t')
write.table(df_loc, args[3], quote = F, row.names = F, sep = '\t')
