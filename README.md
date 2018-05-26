# ReQTL

This toolkit contains the required scripts to transform sequencing files into ReQTL input files then run the MatrixEQTL R package to determine significant variation-expression relationships.

## Getting Started

These instructions will get you a copy of the scripts up and running on your machine for development and testing purposes. See *Running the scripts* for notes on how to use the project on a live system.

### Prerequisites

* Each of the following scripts copied to a working directory on your machine

	```
	build_SNV_matrix.R
	build_gene-exp_matrix.R
	run_matrix_eqtl.R
	```
	You can obtain the full toolkit [here.](https://github.com/HorvathLab/ReQTL/archive/master.zip)
* Output *.csv* files from our ReadCounts tool (https://github.com/HorvathLab/NGS/tree/master/readCounts) containing the read counts extracted per SNV for each sample

* *gene_abund.tab* expression files from *Stringtie* (the scripts may be modified to take input from other expression quanitification software)

## Running the scripts


### build\_gene-exp_matrix.R

Transforms the raw expression files into a matrix with information from all provided samples

#### Input
* A directory containing the gene expression files (one per sample) from Stringtie
* The names of the output gene expression and gene location files


#### Output
* One file (in the script’s directory) with the gene locations for MatrixEQTL 
* One file (in the script’s directory) with the quantile-normalized expression values for MatrixEQTL

#### Sample Command
```
Rscript build_gene-exp_matrix.R /home/expression_files/ BRCA_gene-exp-matrix.txt BRCA_gene-loc-matrix.txt
```

&nbsp;

***

&nbsp;

### build\_SNV_matrix.R

Transforms the read counts into a variant fraction matrix with information from all provided samples

#### Input
* A directory containing the *.csv* files from the output of Readcounts	
* The names of the output SNV matrix and SNV location files


#### Output
* One file (in the script’s directory) with the SNV locations for MatrixEQTL 
* One file (in the script’s directory) with the SNV variant allele fraction matrix for MatrixEQTL


#### Sample command
```
Rscript build_SNV_matrix.R /home/readcounts/ BRCA_SNV_matrix.txt BRCA_SNV-loc-matrix.txt
```
&nbsp;

***

&nbsp;

### run\_matrix_eqtl.R

#### Input

* Names of the SNP, SNP location, expression, and gene location files from build_gene-exp_matrix.R and build_SNV_matrix.R
* Names of the cis and trans output files
* Name of the output file for the qq plot (with .tiff extension)
* *NOTE:* Covariates file can be used (requires modification of the script as specified in the in-code documentation) 

#### Output
* One file (in the script’s directory) with the cis eQTLs with a p value < 0.00001
* One file (in the script’s directory) with the trans eQTLs with a p value < 0.00001


#### Sample command
```
Rscript run_matrix_eqtl.R BRCA_SNV_matrix.txt BRCA_SNV-loc-matrix.txt BRCA_gene-exp-matrix.txt BRCA_gene-loc-matrix.txt BRCA_ReQTL_cis_output.txt BRCA_ReQTL_trans_output.txt BRCA_ReQTL_qqplot.tiff
```
&nbsp;

## Authors

* **Liam Flinn Spurr**

## Acknowledgements

* Muzi Li and Nawaf Alomran for support and assistance in the development of this toolkit
* MatrixEQTL team for their sample code and R package upon which *run\_matrix_eqtl.R* is based
