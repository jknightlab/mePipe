\name{.submitChooseCov}
\alias{.submitChooseCov}
\title{Carry out computations to choose covariates.
This function is called by \code{\link{chooseCov}} to do the actual work
(possibly distributed on an SGE cluster)}
\usage{
  .submitChooseCov(candidates, covOpt, covariate,
    expression, genotype, output, ..., covThreshold)
}
\arguments{
  \item{candidates}{Numeric vector listing the numbers of
  covariates that should be evaluated.}

  \item{covOpt}{List of options to use for reading of
  covariate data.}

  \item{covariate}{Character vector of file names for
  covariate files (see details).}

  \item{expression}{Name of file containing the gene
  expression data.}

  \item{genotype}{Name of file containing the genotyping
  data.}

  \item{output}{Name of output directory}

  \item{...}{Addition parameters for \code{\link{runME}}}

  \item{covThreshold}{FDR threshold to use when determining
  the number of significant eQTLs.}
}
\description{
  Carry out computations to choose covariates. This
  function is called by \code{\link{chooseCov}} to do the
  actual work (possibly distributed on an SGE cluster)
}
\author{
  Peter Humburg
}

