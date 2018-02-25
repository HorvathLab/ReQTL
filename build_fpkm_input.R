# LIAM SPURR
# HORVATH LAB
# 19FEB18
# BUILD_FPKM_INPUT.R

# install missing required packages and load packages
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("MASS"); load_package("dplyr"); load_package("tidyr"); load_package("data.table");

# take path argument from the command line
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  print("ERROR: Please enter a valid file path and output prefix")
  quit()
}

fpkm_path <- args[1]


# function to load the fpkm files in the specified directory into a data frame
load_fpkms <- function(path) {
  # create an empty data frame
  cufflinkstable <- data.frame()
  
  # suppress messages that come from loading in files with fread
  oldw <- getOption("warn")
  options(warn = -1)
  
  # get the list of files and create an empty data frame
  files = list.files(path = path, pattern = "fpkm")
  
  for (f in files) {
    # skip empty expression files
    if(file.info(paste0(path, f))$size == 0) {
      cat('WARNING: EMPTY EXPRESSION FILE, this may cause problems with MatrixEQTL')
      next()
    }
      
    # load each file
    temp = fread(paste0(path, f), data.table = F, header = T)
    
    # get the sample ID -- NOTE: sample ID should be a 3-digit identifier as per README file
    sample_num = substr(f, 1, 3)
    
    # select/rename the relevant columns, removing duplicates (retaining the highest FPKM value for each gene)
    temp = temp %>% dplyr::select(gene_short_name, FPKM, locus) %>%
      rename(Gene = gene_short_name) %>% 
      group_by(Gene) %>% 
      arrange(desc(Gene, FPKM)) %>% 
      ungroup() %>% 
      distinct(Gene, .keep_all = TRUE) %>%
      mutate(Sample = !!sample_num)
    
    # add the new file to the data frame
    cufflinkstable <- rbind(cufflinkstable, as.data.frame(temp))
  }
  options(warn = oldw)
  return(cufflinkstable)
}

cuff <- load_fpkms(fpkm_path)

# get the gene locations
genes <- cuff %>% dplyr::select(Gene, locus) %>% distinct()

# split up the location into a start and end coordinate and arrange for efficient computations
genes <- genes %>% separate(locus, into = c("chr", "s1", "s2"), sep = "[-:]") %>% arrange(Gene)

# turn the gene location data frame into a matrix
gene_locs <- as.matrix(genes)

# remove location column
cuff <- cuff %>% dplyr::select(-locus)

# turn into horizontal data structure
cuff <- cuff %>% spread(Sample, FPKM)

# replace all values below 0.001 with 0.001
cuff[cuff < 0.001] <- 0.001

# turn the expression data frame into a matrix
cuff_matr <- as.matrix(cuff)

# write the output files
output_pref = args[2]

write.table(cuff_matr, paste0(output_pref, "_FPKMs.txt"), sep = "\t", quote = F, row.names = F)
write.table(gene_locs, paste0(output_pref, "_gene-locs.txt"), sep = "\t", quote = F, row.names = F)