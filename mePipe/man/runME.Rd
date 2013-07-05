\name{runME}
\alias{runME}
\title{Run Matrix-eQTL analysis}
\usage{
  runME(expression, genotype, covariate, error, snpsPos,
    genePos, output = "",
    cisOutput = paste(output, "cis", sep = "."),
    exprOpt = getOptions(), genoOpt = getOptions(),
    covOpt = getOptions(), threshold = 1e-05,
    cisThreshold = 0.001,
    model = c("linear", "anova", "cross"), cis = 1e+06,
    bins = 0, verbose = FALSE, qqplot = FALSE,
    cluster = sge.getOption("sge.use.cluster"))
}
\arguments{
  \item{expression}{File name of gene expression data or a
  \code{SlicedData} object with expression data.}

  \item{genotype}{Name of file with genotyping data or a
  \code{SlicedData} object with genotypes.}

  \item{covariate}{Covariates. This may be a \code{matrix},
  \code{data.frame} or vector of file names.}

  \item{error}{Name of file with error covariance matrix.}

  \item{snpsPos}{Name of file with SNP positions.}

  \item{genePos}{Name of file with gene position.}

  \item{output}{Name of output file.}

  \item{cisOutput}{Name of output file for cis
  associations.}

  \item{exprOpt}{List of options for reading of gene
  expression data.}

  \item{genoOpt}{List of options for reading of genotype
  data.}

  \item{covOpt}{List of options for reading of covariate
  files.}

  \item{threshold}{P-value threshold up to which
  associations should be listed in the output file.}

  \item{cisThreshold}{P-value threshold up to which cis
  associations should be listed in the output file.}

  \item{model}{The model to use in the analysis.}

  \item{cis}{The distance (in bp) up to which a SNP should
  be concidered to be in cis of a gene.}

  \item{bins}{Number of p-value bins to use for plotting.}

  \item{verbose}{Flag indicating whether additional
  progress information should be printed.}

  \item{qqplot}{Flag indicating whether a QQ plot of
  p-values should be generated rather than a histogram.}
}
\value{
  The object returned by
  \code{\link[MatrixEQTL]{Matrix_eQTL_main}}
}
\description{
  Run Matrix-eQTL analysis
}
\note{
  Before calling this function a call to
  \code{sge.setDefaultOptions} has to be made.
}
\author{
  Peter Humburg
}

