# ReQTL

This toolkit contains the required scripts to transform sequencing files into ReQTL input files then run MatrixEQTL to determine significant variation-expression relationships.

## Getting Started

These instructions will get you a copy of the scripts up and running on your machine for development and testing purposes. See *Running the scripts* for notes on how to use the project on a live system.

### Prerequisites

* Each of the following scripts copied to a working directory on your machine

	```
	organize_readcounts.py
	build_vaf_input.R
	build_fpkm_input.R
	run_matrix_eqtl.R
	```
	You can obtain the full toolkit [here.](https://github.com/HorvathLab/ReQTL/archive/master.zip)
* Output files from [RNA2DNAlign](https://github.com/HorvathLab/RNA2DNAlign) in a directory structure like:

	```	
	+-- RNA2DNAlign_out/
	|   +-- 001
	|   +-- 002
	|   +-- 003
	|   +-- 004
	|   +-- ...
	```
* *.fpkm* expression files from *cufflinks*
* *.vcf* variant call files from Samtools *mpileup* filtered through RNA2DNAlign with commented lines removed using the following (or similar) shell command

	```
	for i in *.filtered.vcf; do cat $i | egrep –v ‘^##’ > $i”_fix”; done

	```
	&nbsp;

## Running the scripts

### organize_readcounts.py

Transforms the per-sample readcount files from *RNA2DNAlign* into one combined readcounts file

#### Input
* Path to folder containing per-sample *RNA2DNAlign* outputs (see prerequisites for correct structure)
* The desired prefix of the output combined readcounts file

#### Output
* One file with the read counts from all samples

#### Sample command
```
NOTE: must be run using python 2

python2.6 organize_readcounts.py /home/samples/ BRCA
```

&nbsp;

***
&nbsp;


### build\_fpkm_input.R

Transforms the raw expression files into a matrix with information from all provided samples

#### Input
* Path to a folder with FPKM files (one per sample) from Cufflinks, where samples have a unique three-digit identifier at the beginning of the file name

	```Ex:
		+-- cufflinks_out/
		|   +-- 001_expression.fpkm_tracking
		|   +-- 002_expression.fpkm_tracking
		|   +-- 003_expression.fpkm_tracking
		|   +-- ... 
	```
	
* A prefix for the output files

	*Ex. BRCA\_19Feb will become BRCA\_19Feb\_FPKMs.txt and BRCA\_19Feb\_gene-locs.txt*

#### Output
* One file (in the script’s directory) with the gene locations for MatrixEQTL 
* One file (in the script’s directory) with the expression values for MatrixEQTL

#### Sample Command
```
Rscript build_fpkm_input.R /home/expression_files/ BRCA_19Feb
```

&nbsp;

***

&nbsp;

### build\_vaf_input.R

Transforms the readcounts and filtered *.vcf* files into a variant fraction matrix with information from all provided samples

#### Input
* A directory containing the filtered *.vcf* files (with three digit identifier as described for *build_fpkm_input.R*) from the output for RNA2DNAlign with commented lines removed (as described in *Prerequisites*)

	```NOTE: make sure the directory only contains the cleaned .vcfs```
	
* Indiciation of whether the tumor or normal VAF should be used (T => tumor, N => normal)
* The name of the readcounts_all.txt (the output from *organize_readcounts.py*)
* A prefix for the output files


#### Output
* One file (in the script’s directory) with the SNP locations for MatrixEQTL 
* One file (in the script’s directory) with the SNP VAF genotype information for MatrixEQTL


#### Sample command
```
Rscript build_vaf_input.R /home/clean_vcfs/ T BRCA_readcounts_all.txt BRCA_19Feb
```
&nbsp;

***

&nbsp;

### run\_matrix_eqtl.R

#### Input

* Names of the SNP, SNP location, expression, and gene location files from build_fpkm_input.R and build_vaf_input.R
* Names of the cis and trans output files

#### Output
* One file (in the script’s directory) with the cis eQTLs with a p value < 0.00001
* One file (in the script’s directory) with the trans eQTLs with a p value < 0.00001


#### Sample command
```
Rscript run_matrix_eqtl.R BRCA_19Feb_SNPs.txt BRCA_19Feb_SNP-locs.txt BRCA_19Feb_FPKMs.txt BRCA_19Feb_FPKM-locs.txt BRCA_19Feb_cis_eqtls.txt BRCA_19Feb_trans.txt
```
&nbsp;

## Authors

* **Liam Flinn Spurr** - *Initial work, R scripts*
* **Muzi Li** - *Python script*

## Acknowledgements

* MatrixEQTL team for their example script and R package upon which *run\_matrix_eqtl.R* is based
