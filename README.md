# ASA
Ancestral Spectrum Analyzer (ASA) is a program for ancestry inference as well as analysis of population structure with population-specific SNPs. It comes with a utility, Population-Specific SNP Screener (PSNPS), which can be used for identifying and selecting population-specific SNPs in reference populations.

## Dependencies
+ [EIGEN] (https://eigen.tuxfamily.org)
+ [Intel MKL] (https://software.intel.com/)
+ [Boost] (https://www.boost.org/)
+ [zlib] (https://www.zlib.net/)

## Installation
The archive includes executable binary "asa" and "psnps" pre-compiled under CentOS8.0 (x86-64). To compile from source code, edit Makefile to point EIGEN_PATH and MKL_PATH to your own locations of EIGEN3 and MKL, then type "make".

## Usage
Type "./asa --help" from the command line to display program options of ASA:

	--bfile		Input genotype file in plink binary format.
	--ifile		Input file that defines a panel of population-specific SNPs. The file is expected to have five columns without headers, which are SNP ID, population that the SNP is specific to, MAF in the population, minor and major alleles in reference populations.
	--out		Output file name [default: asa].
	--aiv		Calculate and output ancestral information vectors by the method of moment estimate.
	--mle		Calculate and output ancestral information vectors by the approximate maximum likelihood estimate.
	--psv		Calculate and output principal score vectors.
	--pca		Perform PCA and output top n eigenvectors, eigenvalues and PCs [default: n = 20].
	--grm		Calculate and output genetic relationship matrix.
	--batch-size	Number of SNPs to be processed in a batch [default: 10000].
	--thread-num	Number of threads on which the program will be running [default: thread number in your machine - 1].
	--match-alleles	Both minor and major alleles of the population-specific SNPs have to match the two alleles in the BIM file [default: as least one allele of the population-specific SNPs has to match one of the two alleles in the BIM file].
    
Type "./psnps --help" from the command line to display program options of PSNPS:

	--bfile		Input genotype file in plink binary format.
	--ref-pop	Input file that describes reference populations. The file is expected to have two columns without headers: the first is individual ID and the second is the population that the individual belongs to.
	--pop		Screen SNPs that are specific to the specified population. If multiple populations are specified, SNPs polymorphic in all specified populations and monmophic in unspecified populations will be found.
	--snp-num	Number of population-specific SNPs to be saved [default: all SNPs specific to the specified population(s)].
	--random-seed	Set a random seed for selecting the population-specific SNPs to be saved.
	--freq		Calculate and output allele frequencies in the reference populations.
	--out		Output file for saving population-specific SNPs or allele frequencies [default: psnps].
	--maf-min	Exclude SNPs with MAFs smaller than the specified value in the population(s) specified by --pop [default: 0]. Minor alleles are determined by the total samples of the reference populations.
	--maf-max	Exclude SNPs with MAFs larger than the specified value in the population(s) specified by --pop [default: 0.5].
	--miss-max	Exclude SNPs with missing rates larger than the specified value in the reference populations [default: 0].
	--dist-min	Exclude SNPs with distances from previous ones less than or equal to the specified value [default: 0].
	--batch-size	Number of SNPs to be processed in a batch [default: 10000].
	--thread-num	Number of threads on which the program will be running [default: thread number in your machine - 1].

## Credits
Qingmin Kuang and Gang Shi developed the original version (v1.0.0) of the software [1].

Gang Shi updated the software to v1.1.0 by rearranging options in ASA and PSNPS, adding allele columns in the file that defines population-specific SNPs, and providing the panels of population-specific SNPs used in [2].

## Citation
[1] Shi G, Kuang Q. Ancestral spectrum analysis with population-specific variants. Front Genet. 2021;12:724638.

[2] Shi G. Insights from analysis of ancient and modern DNAs with population-specific SNPs. In submission.
