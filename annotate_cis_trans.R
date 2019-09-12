# ANNOTATE_CIS_TRANS.R
# LAST UPDATED BY LIAM FLINN SPURR ON SEPTEMBER 9, 2019
# THIS SCRIPT ANNOTATES THE CIS/TRANS STATUS OF THE RESULTS OF A REQTL ANALYSIS

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

# load in a ReQTL result file
reqtls <- fread(args[1])

# load in the gene locations file
gene_locs <- fread(args[2])

# annotate which gene the SNP resides in
# classify ReQTLs in which the two members of the pair are in the same gene as cis
# classify all others as trans
res <- left_join(reqtls, gene_locs, by = c("gene" = "GeneID")) 
res <- res %>%
  mutate(Reference = gsub("chr", "", Reference),
         newSNV = SNV) %>%
  separate(newSNV, into = c("chrom", "pos", "ref", "alt")) %>% 
  mutate(class = ifelse(Reference == chrom & (pos >= Start & pos <= End), "cis", "trans")) %>%
  select(SNV:FDR, class)

# write the results to a file
output_prefix <- args[3]
write.table(res, paste0(output_prefix, "_cistrans_ann.txt"), quote = F, row.names = F, sep = '\t')