# HARMONIZE_MATRICES.R
# LAST UPDATED BY LIAM FLINN SPURR ON NOVEMBER 26, 2019

# load packages
suppressMessages(library(tidyverse, quietly = TRUE))
suppressMessages(library(data.table, quietly = TRUE))

handle_command_args <- function(args) {
  # make sure all flags are paired
  if(length(args) %% 2 != 0) stop("Command line arguments supplied incorrectly!")
  
  # load the flags into a "dictionary"
  arg_df <<- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>%
    mutate_all(as.character)
  
  # load in matrices
  vaf <<- data.frame(fread(arg_df$value[arg_df$flag == "-r"]))
  express <<- data.frame(fread(arg_df$value[arg_df$flag == "-e"]))
  covar <<- data.frame(fread(arg_df$value[arg_df$flag == "-c"]))
  
}

# take arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
handle_command_args(args)

# get the samples present in each
v.names <- colnames(vaf)[-1]
s.names <- colnames(express)[-1]
c.names <- colnames(covar)[-1]

# get the samples in all 3 matrices
conc <- intersect(intersect(v.names, s.names), c.names)

# select only the concordant samples from all 3 matrices
vaf.clean <- vaf[,c("SNV", conc)]
exp.clean <- express[,c("gene_id", conc)]
cov.clean <- covar[,c("id", conc)]

# write the outputs
if (!dir.exists("output")) {
  cat('Creating output directory...\n')
  dir.create('output')
} 

write.table(vaf.clean, paste0("output/", gsub(".*/|.txt", "", arg_df$value[arg_df$flag == "-r"]), "_harmonized.txt"), quote = F, row.names = F, sep = '\t')
write.table(exp.clean, paste0("output/", gsub(".*/|.txt", "", arg_df$value[arg_df$flag == "-e"]), "_harmonized.txt"), quote = F, row.names = F, sep = '\t')
write.table(cov.clean, paste0("output/", gsub(".*/|.txt", "", arg_df$value[arg_df$flag == "-c"]), "_harmonized.txt"), quote = F, row.names = F, sep = '\t')
cat(paste0("Matrices containing concordant sample lists saved to output/", gsub(".*/|.txt", "", arg_df$value[arg_df$flag == "-r"]), "_harmonized.txt, output/", 
    gsub(".*/|.txt", "", arg_df$value[arg_df$flag == "-e"]), "_harmonized.txt, and output/", gsub(".*/|.txt", "", arg_df$value[arg_df$flag == "-c"]), "_harmonized.txt.\n"))
