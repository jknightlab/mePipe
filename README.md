# The mePipe package

The `mePipe` R package provides a convenient way to carry out several steps that are often required as part 
of an eQTL analysis. [Matrix-eQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) is used 
to carry out the association testing. This very fast and supports the use of many covariates without a 
substantial increase in run time. In addition `mePipe` provides facilities to compute principle components 
of the gene expression data and include a subset of these as covariates to account for confounding variation 
and maximise the number significant associations. Some functionality to account for LD structure is also available. 
All computations can either be carried out on a single processor or run in parallel (as far as possible) on an 
SGE compute cluster.


## Installation
### Dependencies
The mePipe package requires R version 2.15 or newer to run and the following R packages have to be available. 
If they are not already present they need to be installed before mePipe:
* [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/)
* [trio](http://cran.r-project.org/web/packages/trio/index.html)
* [optparse](https://github.com/trevorld/optparse)
* [XML](http://cran.r-project.org/web/packages/XML/index.html)
* [Rsge](https://github.com/humburg/Rsge)

The first four of these are available from CRAN and Bioconductor and can be installed from within R with

```r
source("http://bioconductor.org/biocLite.R")
biocLite(c("MatrixEQTL", "trio", "optparse", "XML", "snow"))
```

Unfortunately `Rsge` has recently been removed from CRAN. 
The latest version is available from [here](https://github.com/humburg/Rsge). Download the source 
distribution to your home directory and (on the machine on which you intend to run Matrix-eQTL) run

```sh
R CMD INSTALL Rsge_0.6.4.tar.gz
```

### Installing mePipe from source

The latest version of the R source package is available [here](https://github.com/jknightlab/mePipe/releases). 
Download this into your home directory and (on the machine on which you would like to add the R package) run

```sh
R CMD INSTALL mePipe_latest.tar.gz
```

### Installing the command-line script
The R package includes the `runMatrixEQTL.R` script that allows use of the 
[MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) package from the command line. 
By default this will be installed into `<R_library>/mePipe/exec/`. When installing `mePipe` for the first 
time it is recommended to create a link to this script in a more convenient place, e.g.

```sh
mkdir ~/bin
cp -s <R_library>/mePipe/exec/runMatrixEQTL.R ~/bin/
```
