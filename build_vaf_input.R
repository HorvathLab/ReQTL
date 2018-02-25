# LIAM SPURR
# HORVATH LAB
# 19FEB18
# BUILD_VAF_INPUT.R

# install missing required packages and load packages
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("MASS"); load_package("dplyr"); load_package("tidyr"); load_package("data.table");

# take arguments from the command line
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
  print("ERROR: Please enter a valid set of arguments")
  quit()
}

# load the vcf files
load_vcfs <- function(path) {
  vcfs = data.frame()
  vcf_path = path
  files = list.files(path = vcf_path, pattern = "filtered.vcf")
  
  # suppress messages that come from loading in files with fread
  oldw <- getOption("warn")
  options(warn = -1)
  
  for (f in files) {
    # skip any empty vcfs
    if(file.info(paste0(vcf_path, f))$size == 0) {
      cat('WARNING: EMPTY VCF FILE, this may cause problems with MatrixEQTL')
      next()
    }
    
    # load each file -- NOTE: commented lines must be removed as per README file
    temp = fread(paste0(vcf_path, f), header = T, data.table = F)
    
    # create sample number variable -- NOTE: sample ID should be a 3-digit number in the format 001 as per README file
    sample_num = substr(f, 1, 3)
    
    # select relevant columns and set proper column names
    temp = temp %>% dplyr::select(1, 2, 4, 5, 10)
    names(temp) = c("Chrom", "Pos", "Ref", "Alt", "Info")
    
    # remove indels and non-standard chromosomes
    temp = temp %>% filter(nchar(as.character(Ref)) == 1, 
                           nchar(as.character(Alt)) == 1,
                           nchar(as.character(Chrom)) < 6)
    
    # create unique SNP id numbers as well as sample id and genotype information columns
    temp = temp %>% mutate(SNP = paste0(Chrom, ":", Pos, "_", Ref, ">", Alt)) %>% 
      mutate(Sample = as.integer(!!sample_num),
             GT = substr(Info, 1, 3), 
             sumGT = as.numeric(ifelse(GT == "0/0", "0", 
                                       ifelse(GT == "0/1", "0.5",
                                              ifelse(GT == "1/1", "1", "NA")))))
    
    # add the new file's data to the existing data frame
    vcfs = rbind(vcfs, temp)
  }
  
  options(warn = oldw)
  # return the combined data frame
  return(vcfs)
}

path_to_vcfs <- args[1]
vcf_df <- load_vcfs(path_to_vcfs)

# load the readcount file
readcounts <- fread(args[3], header = T) %>% 
  mutate(SNP = paste0("chr", Chrom, ":", Pos, "_", Ref, ">", AltA)) %>%
  dplyr::select(-Chrom, -Pos, -Ref, -AltA)

# remove any variants from the vcfs not present in the readcounts file
vcf_df <- vcf_df %>% filter(SNP %in% readcounts$SNP)

# join the vcfs to the readcounts
vcf_df <- vcf_df %>% mutate(Sample = as.integer(Sample))
vcfs_comb <- left_join(vcf_df, readcounts, by = c("SNP", "Sample"))

# choose the proper VAF
if (args[2] == "T") {
  vcfs_comb <- vcfs_comb %>% rename(vaf = Ttr)
} else if (args[2] == "N") {
  vcfs_comb <- vcfs_comb %>% rename(vaf = Ntr)
} else {
  cat('Invalid tissue type specified')
}


# remove NA values, select relavant columns, and calculate vaf data for matrix
vcfs_comb <- vcfs_comb %>% filter(!is.na(vaf)) %>% 
  dplyr::select(SNP, sumGT, Sample, vaf) %>%
  mutate(GT = ifelse(sumGT == 2, 1, 
                     ifelse(sumGT == 1, vaf, NA)))

# remove unneeded columns and turn data frame into a horizontal data structure
vcfs_comb <- vcfs_comb %>% dplyr::select(-sumGT, -vaf)
vcfs_spread <- vcfs_comb %>% spread(Sample, GT, fill = 0)

# turn data frame into a matrix
matrix_vcfs <- as.matrix(vcfs_spread)

# write the vaf output file
output_pref = args[4]
write.table(matrix_vcfs, paste0(output_pref, "_SNPs.txt"), sep = '\t', quote = F, row.names = F)

# get the SNP positions
pos <- vcf_df %>% dplyr::select(SNP, Chrom, Pos) %>% filter(SNP %in% vcfs_comb$SNP)

# write the output SNP position file
write.table(pos, paste0(output_pref, "_SNP-locs.txt"), sep = "\t", quote = F, row.names = F)