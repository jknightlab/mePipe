\name{chooseCov}
\alias{chooseCov}
\title{Determine the optimal number of covariates.}
\usage{
  chooseCov(expression, genotype, covariate,
    candidates = seq(5, 50, by = 5), covThreshold = 0.01,
    covOpt = getOptions(), output = "covSelect",
    useCluster = FALSE, ...)
}
\arguments{
  \item{expression}{Name of file containing the gene
  expression data.}

  \item{genotype}{Name of file containing the genotyping
  data.}

  \item{covariate}{Character vector of file names for
  covariate files (see details).}

  \item{candidates}{Numeric vector listing the numbers of
  covariates that should be evaluated.}

  \item{covThreshold}{FDR threshold to use when determining
  the number of significant eQTLs.}

  \item{covOpt}{List of options to use for reading of
  covariate data.}

  \item{output}{Name of output directory}

  \item{useCluster}{Logical indicating whether jobs should
  be submited to SGE or run locally.}

  \item{...}{Addition parameters for \code{\link{runME}}}
}
\value{
  A list with elements \item{covariates}{Numeric vector
  with the number of covariates used in each iteration.}
  \item{eqtls}{A list with components \code{significant}
  and \code{all}. Each consisting of a numeric vector with
  the number of genes that have significant associations or
  any association with SNPs respectively.} \item{best}{Base
  of file name(s) that contain the results of the run that
  produced the most eQTLs.}
}
\description{
  Determine the optimal number of covariates.
}
\details{
  The argument \code{covariates} may be a vector of several
  file names. In this case the first one will be used to
  select varying numbers of covariates and the remaining
  ones will always be used.
}
\author{
  Peter Humburg
}

