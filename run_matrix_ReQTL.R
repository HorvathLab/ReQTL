# RUN_MATRIX_REQTL.R
# LAST UPDATED BY LIAM FLINN SPURR ON NOVEMBER 26, 2019
# THIS SCRIPT IS BASED OFF OF THE MATRIX EQTL SAMPLE CODE
# ORIGINAL CODE AVAILABLE AT http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html#sample

# load packages
suppressMessages(library(MatrixEQTL))
suppressMessages(library(tidyverse))

# get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# manages command line arguments
handle_command_args <- function(args) {
  # make sure all flags are paired
  if(length(args) %% 2 != 0) stop("Command line arguments supplied incorrectly!")
  
  # load the flags into a "dictionary"
  arg_df <- data.frame(cbind(flag = args[seq(1, length(args), by = 2)], value = args[seq(2, length(args), by = 2)])) %>%
    mutate_all(as.character)
  
  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useModel <<- modelLINEAR
  
  # Genotype and gene location file names
  SNV_file_name <<- arg_df$value[arg_df$flag == "-s"]
  snvs_location_file_name <<- arg_df$value[arg_df$flag == "-sl"]
  
  # Gene expression file name
  expression_file_name <<- arg_df$value[arg_df$flag == "-g"]
  gene_location_file_name <<- arg_df$value[arg_df$flag == "-gl"]
  
  # Covariates file name
  covariates_file_name <<- ifelse(length(arg_df$value[arg_df$flag == "-c"]) > 0, arg_df$value[arg_df$flag == "-c"], "")
  
  # Whether to split into cis and trans
  split_cis_trans <<- arg_df$value[arg_df$flag == "-ct"]
  
  # Error covariance matrix
  # Set to numeric() for identity.
  errorCovariance <<- numeric()
  
  # Filename for qqplot
  qqplot_filename <<- paste0("output/", arg_df$value[arg_df$flag == "-o"], "_qqplot.tiff")
  
  if(split_cis_trans == "T") {
    # Output file name
    output_file_name_cis <<- paste0("output/", arg_df$value[arg_df$flag == "-o"], "_cis_ReQTLs.txt")
    output_file_name_tra <<- paste0("output/", arg_df$value[arg_df$flag == "-o"], "_trans_ReQTLs.txt")
    
    # Only associations significant at this level will be saved
    pvOutputThreshold_cis <<- as.numeric(arg_df$value[arg_df$flag == "-pcis"])
    pvOutputThreshold_tra <<- as.numeric(arg_df$value[arg_df$flag == "-ptr"])
    
    # Distance for local gene-SNV pairs
    cisDist <<- 1e6
    
  } else {
    # Output file name
    output_file_name <<- paste0("output/", arg_df$value[arg_df$flag == "-o"], "_all_ReQTLs.txt")
    
    # Only associations significant at this level will be saved
    pvOutputThreshold <<- as.numeric(arg_df$value[arg_df$flag == "-p"])
  }
}

handle_command_args(args)

## Load genotype data
snvs <- SlicedData$new()
snvs$fileDelimiter = "\t"      # the TAB character
snvs$fileOmitCharacters = "NA" # denote missing values
snvs$fileSkipRows = 1          # one row of column labels
snvs$fileSkipColumns = 1       # one column of row labels
snvs$fileSliceSize = 2000      # read file in slices of 2,000 rows
snvs$LoadFile(SNV_file_name)

## Load gene expression data

gene <- SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

## Load covariates
cvrt <- SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
if(length(covariates_file_name) > 1) {
  cvrt$LoadFile(covariates_file_name)
}

## Run the analysis
snvspos <- read.table(snvs_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos <- read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

if(split_cis_trans == "T") {
  me <- Matrix_eQTL_main(
    snps = snvs, 
    gene = gene, 
    cvrt = cvrt,
    output_file_name = output_file_name_tra,
    pvOutputThreshold  = pvOutputThreshold_tra,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE, 
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snvspos, 
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
} else {
  me <- Matrix_eQTL_main(
    snps = snvs, 
    gene = gene, 
    cvrt = cvrt,
    output_file_name = output_file_name,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE, 
    pvOutputThreshold= pvOutputThreshold,
    snpspos = snvspos, 
    genepos = genepos,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
}

tiff(filename = qqplot_filename)
plot(me)
dev.off()

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
