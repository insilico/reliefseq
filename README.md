ReliefSeq
========================

#### A feature selection tool for biological SEQuence data ####

### Description ###
ReliefSeq is a free, open-source command-line tool for analysis of GWAS (SNP) 
and other types of biological data. Several modes are available for various 
types of analysis.

ReliefSeq is developed by the In Silico Research Group at the Tandy School
of Computer Science of the [University of Tulsa](http://www.utulsa.edu).  Our
research is sponsored by the NIH and William K. Warren foundation.  For more
details, visit our research [website](http://insilico.utulsa.edu).

### Dependencies ###
* GNU Scientific library (libgsl)
* [Boost](http://www.boost.org) system, filesystem, and program-options libraries 
* OpenMP is required to take advantage of the parallelized distance matrix 
  calculations for ReliefF. This library is typically installed alongside the
  compiler toolchain.

### Compilation Environment and Instructions ###
To compile this code, a GNU toolchain and suitable environment are required.
GNU g++ has been used to successfully compile the code. We have successfully 
built and run ReliefSeq on:

 * Linux (64-bit Ubuntu) (gcc-4.6)
 * Linux (64-bit) gcc (Debian 4.4.5-8) 4.4.5

To build ReliefSeq, first run the bootstrap script:

    ./bootstrap.sh

Ignore any extraneous warnings. This calls autoreconf and generates the 
configure script.  From this point, the standard build procedure:

    ./configure && make && sudo make install

will generate the `Makefile`, compile and link the code, and copy the 
executable to the installation directory prefix (default of `/usr/local`). 

### Usage ###

reliefseq:

	--help                                produce help message
	--verbose                             verbose output
	--convert                             convert data set to data set - does not
																				run reliefseq
	--write-best-k                        optimize k, write best k's
	--write-each-k-scores                 optimize k, write best scores for each 
																				k
	-c [ --config-file ] arg              read configuration options from file - 
																				command line overrides these
	-s [ --snp-data ] arg                 read SNP attributes from genotype 
																				filename: txt, ARFF, plink (map/ped, 
																				binary, raw)
	--snp-file-type arg                   Ignore file extension and use type: 
																				textwhitesp, wekaarff, plinkped, 
																				plinkbed, plinkraw, dge, birdseed
	-n [ --numeric-data ] arg             read continuous attributes from 
																				PLINK-style covar file
	-X [ --numeric-transform ] arg        perform numeric transformation: 
																				normalize, standardize, zscore, log, 
																				sqrt, anscombe
	-a [ --alternate-pheno-file ] arg     specifies an alternative 
																				phenotype/class label file; one value 
																				per line
	-g [ --algorithm-mode ] arg (=relieff)
																				Relief algorithm mode 
																				(relieff|reliefseq)
	--seq-algorithm-mode arg (=snr)       Relief algorithm mode (snr|tstat)
	--seq-snr-mode arg (=snr)             Seq interaction algorithm SNR mode 
																				(snr|relieff)
	--seq-tstat-mode arg (=pval)          Seq interaction algorithm t-statistic 
																				mode (pval|abst|rawt)
	--seq-algorithm-s0 arg (=0.050000000000000003)
																				Seq interaction algorithm s0 (0.0 <= s0
																				<= 1.0)
	-t [ --num-target ] arg               target number of attributes to keep 
																				after backwards selection
	-r [ --iter-remove-n ] arg            number of attributes to remove per 
																				iteration of backwards selection
	-p [ --iter-remove-percent ] arg      percentage of attributes to remove per 
																				iteration of backwards selection
	--normalize-scores arg (=0)           normalize ReliefF scores? (0|1)
	-O [ --out-dataset-filename ] arg     write a new tab-delimited data set with
																				EC filtered attributes
	-o [ --out-files-prefix ] arg (=reliefseq_default)
																				use prefix for all output files
	--snp-metric arg (=gm)                metric for determining the difference 
																				between subjects (gm|am|nca|nca6)
	-B [ --snp-metric-nn ] arg (=gm)      metric for determining the difference 
																				between subjects (gm|am|nca|nca6|km)
	-W [ --snp-metric-weights ] arg (=gm) metric for determining the difference 
																				between SNPs (gm|am|nca|nca6)
	-N [ --numeric-metric ] arg (=manhattan)
																				metric for determining the difference 
																				between numeric attributes 
																				(manhattan=|euclidean)
	-x [ --snp-exclusion-file ] arg       file of SNP names to be excluded
	-k [ --k-nearest-neighbors ] arg (=10)
																				set k nearest neighbors (0=optimize k)
	--kopt-begin arg (=1)                 optimize k starting with kopt-begin
	--kopt-end arg (=1)                   optimize k ending with kopt-end
	--kopt-step arg (=1)                  optimize k incrementing with kopt-step
	-m [ --number-random-samples ] arg (=0)
																				number of random samples (0=all|1 <= n 
																				<= number of samples)
	-b [ --weight-by-distance-method ] arg (=equal)
																				weight-by-distance method 
																				(equal|one_over_k|exponential)
	--weight-by-distance-sigma arg (=2)   weight by distance sigma
	-d [ --diagnostic-tests ] arg         performs diagnostic tests and sends 
																				output to filename without running EC
	-D [ --diagnostic-levels-file ] arg   write diagnostic attribute level counts
																				to filename
	--dge-counts-data arg                 read digital gene expression counts 
																				from text file
	--dge-norm-factors arg                read digital gene expression 
																				normalization factors from text file
	--birdseed-snps-data arg              read SNP data from a birdseed formatted
																				file
	--birdseed-phenos-data arg            read birdseed subjects phenotypes from 
																				a text file
	--birdseed-subjects-labels arg        read subject labels from filename to 
																				override names from data file
	--birdseed-include-snps arg           include the SNP IDs listed in the text 
																				file
	--birdseed-exclude-snps arg           exclude the SNP IDs listed the text 
																				file
	--distance-matrix arg                 create a distance matrix for the loaded
																				samples and exit
	--gain-matrix arg                     create a GAIN matrix for the loaded 
																				samples and exit
	--dump-titv-file arg                  file for dumping SNP 
																				transition/transversion information

All commands will include an input file (`-s/--snp-data`), and, optionally, 
an output file prefix (`-o/--output-files-prefix`).

To perform a standard, all-default-parameters analysis,

    ./reliefseq -s snpdata.ped -o result

This will use genotype/phenotype information from `snpdata.ped`, a PLINK
plaintext GWAS file, in the feature selection.  All of the output files 
produced will be prepended with 'result'.

This produces a file called `result.reliefseq`, in which the SNPs are ranked 
in descending order.

For additional examples, see the 
[ReliefSeq](http://insilico.utulsa.edu/ReliefSeq.php)
page on our website.

### Contributors ###
See [AUTHORS](https://github.com/insilico/reliefseq/blob/master/AUTHORS) file.

### References ###
