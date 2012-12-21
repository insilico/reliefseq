ReliefSeq
========================

#### A feature selection tool for biological SEQuence data ####

### Description ###
ReliefSeq is a free, open-source command-line tool for analysis of GWAS (SNP) and
other types of biological data. Several modes are available for various types
of analysis.

ReliefSeq is developed by the In Silico Research Group at the Tandy School
of Computer Science of the [University of Tulsa](http://www.utulsa.edu).  Our
research is sponsored by the NIH and William K. Warren foundation.  For more
details, visit our research [website](http://insilico.utulsa.edu).

### Dependencies ###
* GNU Scientific library (libgsl)
* [Boost](http://www.boost.org) system, filesystem, and program-options libraries 
* OpenMP is required to take advantage of the parallelized tree growth in 
Random Jungle and distance matrix calculations for ReliefF.  This is another 
library typically installed alongside the compiler toolchain.

### Compilation Environment and Instructions ###
To compile this code, a GNU toolchain and suitable environment are required.
GNU g++ has been used to successfully compile the code.

We have successfully built and run ReliefSeq on:

 * Linux (64-bit Ubuntu) (gcc-4.6)

To build ReliefSeq, first run the bootstrap script

    ./bootstrap.sh

Ignore any extraneous warnings. This calls autoreconf and generates the 
configure script.  From this point, a standard

    ./configure && make && sudo make install

will generate the `Makefile`, compile and link the code, and copy the objects to
the installation directory (default of `/usr/local`). 

The resulting binary src/ec_static.exe will run as a command-line tool.

### Usage ###

[INSERT LATEST USAGE]

All commands will include an input file (`-s/--snp-data`), and, optionally, 
an output file prefix (`-o/--output-files-prefix`).

To perform a standard, all-default-parameters analysis,

    ./reliefseq -s snpdata.ped -o result

This will use genotype/phenotype information from `snpdata.ped`, a PLINK
plaintext GWAS file, in the feature selection.  All of the output files 
produced will be prepended with 'result'.

This produces a file called `result.reliefseq`, in which the SNPs are ranked 
in descending order.

For additional examples, see the [ReliefSeq](http://insilico.utulsa.edu/reliefseq)
page on our research website.

### Contributors ###
See [AUTHORS](https://github.com/insilico/reliefseq/blob/master/AUTHORS) file.

### References ###
