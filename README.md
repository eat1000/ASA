# ASA
Ancestral Spectrum Analyzer (ASA) is a program for ancestry inference as well as analysis of population structure with population specific SNPs.

## Authors
Qingmin Kuang and Gang Shi (gshi@xidian.edu.cn)

## Dependencies
+ [EIGEN] (https://eigen.tuxfamily.org)
+ [Intel MKL] (https://software.intel.com/)
+ [Boost] (https://www.boost.org/)
+ [zlib] (https://www.zlib.net/)

## Installation
The archive includes an executable binary "asa" pre-compiled under CentOS8.0 (x86-64). To compile from source code, edit Makefile to point EIGEN_PATH and MKL_PATH to your own locations of EIGEN3 and MKL, then type "make".

## Usage
Type "./asa --help" from the command line to display program options:

    --bfile		Input genotype file in plink binary format.
    --ifile		Input information file of a panel of population specific SNPs. The file is expected to have three columns without headers: the first is SNP ID, the second is population that the SNP is specific to, and the third is MAF of the SNP in the population.
    --out		Output file name (default: asa).
    --aiv		Calculate and output ancestral information vectors by the method of moment estimate.
    --mle		Calculate and output ancestral information vectors by the approximated maximum likelihood estimate.
    --psv		Calculate and output principal score vectors.
    --pca		Perform PCA and output top n eigenvectors, eigenvalues and PCs (default: n = 20).
    --freq		Calculate and output allele frequencies.
    --grm		Calculate and output genetic relationship matrix.
    --maf-min	Exclude SNPs with MAFs smaller than the specified value (default: 0).
    --maf-max	Exclude SNPs with MAFs larger than the specified value (default: 0.5).
    --dist-min	Exclude SNPs with distances from previous ones less or equal to the specified value (default: 0).
    --miss-max	Exclude SNPs with missing rates larger than the specified value (default: 1).
    --batch-size	Number of SNPs to be processed in a batch (default: 10000).
    --thread-num	Number of threads on which the program will be running (default: thread number in your machine).

## Citation
Gang Shi, Qingmin Kuang. Ancestral spectrum analysis with population specific SNPs. In submission.
