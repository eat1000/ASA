# ASA
Ancestral Spectrum Analyzer (ASA) is a package for ancestry inference, co-ancestry analysis, or analysis of population structure with population-specific SNPs. It includes "asa" for global ancestral spectrum analysis, "lai" for local ancestry inference and "psnps" for screening population-specific SNPs in reference populations.

## Dependencies
+ [EIGEN] (https://eigen.tuxfamily.org)
+ [Intel MKL] (https://software.intel.com/)
+ [Boost] (https://www.boost.org/)
+ [zlib] (https://www.zlib.net/)
+ [htslib] (https://github.com/samtools/htslib)
  
## Installation
The archive includes executable binary "asa", "lai", "psnps" and "kartag" pre-compiled under CentOS8.0 (x86-64). To compile from source code, edit Makefile to point library paths to your own locations, then type "make".

## Usage
Type "./asa --help" from the command line to display program options of "asa":

	--bfile		Input genotype file in plink binary format.
	--ifile		Input file that defines a panel of population-specific SNPs. The file is expected to have at least five columns without headers, which are SNP ID, population that the SNP is specific to, MAF in the population, minor and major alleles in reference populations.
	--out		Output file name [default: asa].
	--aiv		Calculate and output ancestral information vectors by the method of moment estimate.
	--mle		Calculate and output ancestral information vectors by the approximate maximum likelihood estimate.
	--psv		Calculate and output principal score vectors.
	--pca		Perform PCA and output top n eigenvectors, eigenvalues and PCs [default: n = 20].
	--grm		Calculate and output genetic relationship matrix.
	--batch-size	Number of SNPs to be processed in a batch [default: 10000].
	--thread-num	Number of threads on which the program will be running [default: thread number in your machine - 1].
	--match-alleles	Both minor and major alleles of the population-specific SNPs have to match the two alleles in the BIM file [default: as least one allele of the population-specific SNPs has to match one of the two alleles in the BIM file].

Type "./lai --help" from the command line to display program options of "lai":

	--vfile		Prefix of haplotype file in bgzipped VCF format [.vcf.gz] with the associated index file in tbi or csi format [.vcf.gz.tbi] or [.vcf.gz.csi].
	--ifile		Input file that defines a panel of population-specific SNPs. The file is expected to have seven columns without headers, which are in the order of SNP ID, population that the SNP is specific to, MAF in the population, minor and major alleles in reference populations, chromsome and position.
 	--out		Prefix of output file [default: lai].
  	--chr		Specify the chromosome for the analysis [default: 1].
   	--pop-allele-id	Output populaiton-specific alleles carried by the speficied individual.
   	--lai		Calculate and output local ancestral inference.
  	--laiv		Output local ancestral information vectors [default: no output].
	--laiv2lai	Read the specifed local ancestral infomation vector file and output the local ancestral inference.
	--match-alleles	Both minor and major alleles of the population-specific SNPs have to match the two alleles in the haplotype file [default: as least one allele of the population-specific SNPs has to match one of the two alleles in the haplotype file].
	--window-size	Window size in the local ancestral information calculation [default: 2000000] (bp).
	--laiv-min	Minimum value of laiv for calling ancestries [default: 0.000001].
	--allele-min	Minimum number of alleles in the window for calling the ancestry [default: 2].
  
Type "./psnps --help" from the command line to display program options of "psnps":

	--bfile		Input genotype file in plink binary format.
	--ref-pop	Input file that describes reference populations. The file is expected to have two columns without headers: the first is individual ID and the second is the population that the individual belongs to.
 	--out		Output file for saving population-specific SNPs or allele frequencies [default: psnps].
	--snp-panel	Screen SNPs that are specific to one reference population at a time and compose a SNP panel.
	--pop		Screen SNPs that are specific to the specified population. If multiple populations are specified, SNPs polymorphic in all specified populations and monmophic in unspecified populations will be found.
	--snp-num	Number of population-specific SNPs to be saved [default: all SNPs specific to the specified population(s)].
	--random-seed	Set a random seed for selecting the population-specific SNPs to be saved.
	--freq		Calculate and output allele frequencies in the reference populations.
	--maf-min	Exclude SNPs with MAFs smaller than the specified value in the population(s) specified by --pop [default: 0]. Minor alleles are determined by the total samples of the reference populations.
	--maf-max	Exclude SNPs with MAFs larger than the specified value in the population(s) specified by --pop [default: 0.5].
	--miss-max	Exclude SNPs with missing rates larger than the specified value in the reference populations [default: 0].
	--dist-min	Exclude SNPs with distances from previous ones less than or equal to the specified value [default: 0].
	--tv-only	Keep SNPs that are transversions only.
	--batch-size	Number of SNPs to be processed in a batch [default: 10000].
	--thread-num	Number of threads on which the program will be running [default: thread number in your machine - 1].
	--sort		Sort the specified SNP panel file by chromosome and position.

## Credits
Qingmin Kuang and Gang Shi developed the original version (v1.0.0) of the software [1].

Gang Shi upgraded the software to v1.1.0 by updating options in "asa" and "psnps", adding allele columns in the file that defines population-specific SNPs, and enclosing the panels of population-specific SNPs used in [2].

Haoyue Fu and Gang Shi implemented the local ancetry inference method in v1.2.0 as "lai" and updated "psnps" [3].

## Citations
[1] Shi G, Kuang Q. Ancestral spectrum analysis with population-specific variants. Front Genet. 2021;12:724638.

[2] Shi G. Insights from the analysis of ancient and modern DNA with population-specific SNPs. https://doi.org/10.21203/rs.3.rs-3447042/v1

[3] Fu H, Shi G. Local ancestry inference with population-specific SNPs —— A study of admixed populations in the 1000 Genomes Project. In submission. 
