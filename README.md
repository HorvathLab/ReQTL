# ReQTL: Identifying correlations between expressed SNVs and gene expression using RNA-sequencing data
### Updated November 26, 2019

This toolkit contains the required scripts to transform sequencing files into ReQTL input files and run the MatrixEQTL R package to identify significant variation-expression relationships.


&nbsp;


## Getting Started

These instructions will get you a copy of the scripts up and running on your machine for development and testing purposes. See *Running the scripts* for notes on how to use the project on a live system. We have provided sample data that can be used to test the pipeline. It also serves as an example of the data format the pipeline expects.

This package was developed on R version 3.6.1 on macOS High Sierra.

### Prerequisites

* The following R packages installed on your machine:

  ```
  tidyverse
  MatrixEQTL
  data.table
  GenomicRanges
  ```
  
  You may install the packages using the following commands:
  ```
  install.packages(c("tidyverse", "MatrixEQTL", "data.table", "BiocManager"))
  BiocManager::install("GenomicRanges")
  ```

* Each of the following scripts downloaded on your machine

	```
	build_gene-exp_matrix.R
	build_VAF_matrix.R
	harmonize_matrices.R
	run_matrix_ReQTL.R
	annotate_cis_trans.R
	```
	You can obtain the full toolkit [here.](https://github.com/HorvathLab/ReQTL/archive/master.zip)
* Output *.csv* files from our ReadCounts tool (https://github.com/HorvathLab/NGS/tree/master/readCounts) containing the read counts extracted per SNV for each sample

* *gene_abund.tab* expression files from *Stringtie* (the scripts may be modified to take input from other expression quanitification software)

* (OPTIONAL) a file containing sample covariate information (see sample data)


&nbsp;


## Running the scripts

The scripts are designed to be run from the *Unix command line* (Terminal on macOS) from the root directory (by default ReQTL-master if the toolkit was downloaded from the link above). Make sure to *cd* to this directory before beginning.

***

&nbsp;

### build\_gene-exp_matrix.R

Transforms the raw expression files into a matrix with information from all provided samples

#### Input
* A directory containing the gene expression files (one per sample) from Stringtie
* The desired prefix of the output gene expression and gene location files


#### Output
* One file (in the output directory) with the gene expression values for MatrixEQTL
* One file (in the output directory) with the gene locations for MatrixEQTL 

#### Sample Command
```
Rscript build_gene-exp_matrix.R -e data/ -o ReQTL_test
```

&nbsp;

***

&nbsp;

### build\_VAF_matrix.R

Transforms the read counts into a variant fraction matrix with information from all provided samples

#### Input
* A directory containing the *.csv* files from the output of Readcounts	
* The desired prefix of the output SNV matrix and SNV location files


#### Output
* One file (in the output directory) with the SNV locations for MatrixEQTL 
* One file (in the output directory) with the SNV variant allele fraction matrix for MatrixEQTL


#### Sample command
```
Rscript build_VAF_matrix.R -r data/ -o ReQTL_test
```
&nbsp;

***

&nbsp;

### harmonize\_matrices.R

Harmonizes matrices so that all inputs for run_matrix_ReQTL.R contain the same samples (optional, but helps to avoid errors)

#### Input

* The path to the VAF matrix created by build_VAF_matrix.R
* The path to the gene expression matrix created by build_gene-exp_matrix.R
* The path to a covariate matrix (if no covariate information is used, you will need to modify this script to remove the references to the covariate matrix)

#### Output
* Three matrices (in the output directory) corresponding to the three input files, but including only samples that were present in all three input files


#### Sample command
```
Rscript harmonize_matrices.R -r output/ReQTL_test_VAF_matrix.txt -e output/ReQTL_test_gene-exp_matrix.txt -c data/covariates_matrix.txt
```
&nbsp;

***

&nbsp;

### run\_matrix_ReQTL.R

*This script is based off of the sample code from Shabalin, et al (2012)*

Runs the ReQTL analysis using MatrixEQTL

#### Input

* Names of the SNV (-s), SNV location (-sl), gene expression (-g), and gene location (-gl) files from build_gene-exp_matrix.R and build_VAF_matrix.R (or harmonize_matrices.R)
* *OPTIONAL:* Name of the covariates matrix (-c); we include an example in the "data" folder
* Prefix for the output files (-o)
* Logical (T or F, -ct) specifying whether to split the output into *cis* and *trans*
* P-value thresholds for the *cis* (-pcis) and *trans* (-ptr) output files or the unified output file (-p) depending on the logical specified above


#### Output
* One file (in the output directory) with the *cis* ReQTLs and one file (in the script’s directory) with the *trans* ReQTLs OR one file (in the script’s directory) with all of the unified ReQTLs depending on the logical specified above
* One QQ plot of p-values (in the output directory)


#### Sample commands

Splitting *cis* and *trans*
```
Rscript run_matrix_ReQTL.R -s output/ReQTL_test_VAF_matrix_harmonized.txt -sl output/ReQTL_test_VAF-loc_matrix.txt -g output/ReQTL_test_gene-exp_matrix_harmonized.txt -gl output/ReQTL_test_gene-exp-loc_matrix.txt -c output/covariates_matrix_harmonized.txt -ct T -o ReQTL_test -pcis 0.001 -ptr 0.00001
```

Unified *cis* and *trans*
```
Rscript run_matrix_ReQTL.R -s output/ReQTL_test_VAF_matrix_harmonized.txt -sl output/ReQTL_test_VAF-loc_matrix.txt -g output/ReQTL_test_gene-exp_matrix_harmonized.txt -gl output/ReQTL_test_gene-exp-loc_matrix.txt -c output/covariates_matrix_harmonized.txt -ct F -o ReQTL_test -p 0.001
```
&nbsp;

### annotate\_cis_trans.R

Annotates the output of ReQTL as cis/trans based on whether the SNV resides within its paired gene

#### Input

* The path to a ReQTL results file
* The path to a file containing gene location annotations
* The desired prefix of the output annotated results file

#### Output
* One file (in the output directory) with the ReQTLs annotated as *cis* or *trans*


#### Sample command
```
Rscript annotate_cis_trans.R -r output/ReQTL_test_all_ReQTLs.txt -g data/gene_locations_hg38.txt -o ReQTL_test
```


&nbsp;


## Authors and Acknowledgements

*Liam Spurr*, Nawaf Alomran, Pavlos Bousounis, Dacian Reece-Stremtan, Prashant N M, Hongyu Liu, Piotr Słowiński, Muzi Li, Qianqian Zhang, Justin Sein, Gabriel Asher, Keith A. Crandall, Krasimira Tsaneva-Atanasova, and Anelia Horvath 

We would like to thank the Matrix EQTL team (Shabalin, et al. 2012) for their sample code and R package upon which *run\_matrix_ReQTL.R* is based.
