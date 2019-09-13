# ReQTL: Identifying SNV – gene expression correlations using RNA-sequencing data

This toolkit contains the required scripts to transform sequencing files into ReQTL input files and run the MatrixEQTL R package to identify significant variation-expression relationships.

## Getting Started

These instructions will get you a copy of the scripts up and running on your machine for development and testing purposes. See *Running the scripts* for notes on how to use the project on a live system. We have provided sample data that can be used to test the pipeline. It also serves as an example of the data format the pipeline expects.

### Prerequisites

* Each of the following scripts copied to a working directory on your machine

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

## Running the scripts


### build\_gene-exp_matrix.R

Transforms the raw expression files into a matrix with information from all provided samples

#### Input
* A directory containing the gene expression files (one per sample) from Stringtie
* The desired prefix of the output gene expression and gene location files


#### Output
* One file (in the script’s directory) with the gene expression values for MatrixEQTL
* One file (in the script’s directory) with the gene locations for MatrixEQTL 

#### Sample Command
```
Rscript build_gene-exp_matrix.R /home/expression_files/ my_file_prefix
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
* One file (in the script’s directory) with the SNV locations for MatrixEQTL 
* One file (in the script’s directory) with the SNV variant allele fraction matrix for MatrixEQTL


#### Sample command
```
Rscript build_VAF_matrix.R /home/readcounts/ my_file_prefix
```
&nbsp;

***

&nbsp;

### harmonize\_matrices.R

Harmonizes matrices so that all inputs for run_matrix_ReQTL.R contain the same samples (optional, but helps to avoid errors)

#### Input

* The path to the VAF matrix created by build_VAF_matrix.R
* The path to the gene expression matrix created by build_gene-exp_matrix.R
* The path to the covariate matrix created by build_cov_matrix.R (if no covariate information is used, you will need to modify this script to remove the references to the covariate matrix)

#### Output
* The three input matrices with only the samples contained in all three


#### Sample command
```
Rscript harmonize_matrices.R VAF_matrix.txt gene-exp-matrix.txt cov_matrix.txt
```
&nbsp;

### run\_matrix_ReQTL.R

*This script is based off of the sample code from Shabalin, et al (2012)*

Runs the ReQTL analysis using MatrixEQTL

#### Input

* Names of the SNV (-s), SNV location (-sl), gene expression (-g), and gene location (-gl) files from build_gene-exp_matrix.R and build_VAF_matrix.R (or harmonize_matrices.R)
* *OPTIONAL:* Name of the covariates matrix (-c); we include an example in the "data" folder
* Logical (T or F, -ct) specifying whether to split the output into *cis* and *trans*
* Name of the output file for the qq plot (-qq)
* Names of the *cis* (-cis) and *trans* (-tr) output files or the unified output file (-o) depending on the logical specified above
* P-value thresholds for the *cis* (-pcis) and *trans* (-ptr) output files or the unified output file (-p) depending on the logical specified above


#### Output
* One file (in the script’s directory) with the *cis* ReQTLs
* One file (in the script’s directory) with the *trans* ReQTLs
OR
* One file (in the script’s directory) with all of the unified ReQTLs depending on the logical specified above


#### Sample commands

Splitting *cis* and *trans*
```
Rscript run_matrix_ReQTL.R -s VAF_matrix_harmonized.txt -sl VAF-loc_matrix.txt -g gene-exp_matrix_harmonized.txt -gl gene-exp-loc_matrix.txt -c cov_matrix_harmonized.txt -ct T -qq test_qqplot -cis output_cis -tr output_tr -pcis 0.001 -ptr 0.00001
```

Unified *cis* and *trans*
```
Rscript run_matrix_ReQTL.R -s VAF_matrix_harmonized.txt -sl VAF-loc_matrix.txt -g gene-exp_matrix_harmonized.txt -gl gene-exp-loc_matrix.txt -c cov_matrix_harmonized.txt -ct F -qq test_qqplot -o output -p 0.0001
```
&nbsp;

### annotate\_cis_trans.R

Annotates the output of ReQTL as cis/trans based on whether the SNV resides within its paired gene

#### Input

* The path to a ReQTL results file
* The path to the gene location matrix created by build_gene-exp_matrix.R
* The desired prefix of the output annotated results file

#### Output
* One file (in the script’s directory) with the ReQTLs annotated as *cis* or *trans*


#### Sample command
```
Rscript annotate_cis_trans.R output.txt gene-exp-loc_matrix.txt my_file_prefix
```
&nbsp;

## Authors

* **Liam Flinn Spurr**

## Acknowledgements

* Anelia Horvath, Muzi Li, and Nawaf Alomran for support and assistance in the development of this toolkit
* MatrixEQTL team for their sample code and R package upon which *run\_matrix_ReQTL.R* is based
