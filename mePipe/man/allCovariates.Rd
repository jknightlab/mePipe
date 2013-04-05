\name{allCovariates}
\alias{allCovariates}
\title{Run principle component analysis of gene expression data to generate
table of possible covariates.}
\usage{
  allCovariates(expression, exprOpt = getOptions(), output,
    sep = "\t")
}
\arguments{
  \item{expression}{Name of file containing the gene
  expression data.}

  \item{exprOpt}{Options for reading of gene expression
  data.}

  \item{output}{Name of file to which the principle
  components should be saved (optional).}

  \item{sep}{Character used to separate entries in output
  file.}
}
\value{
  \code{data.frame} with all principle components
}
\description{
  Run principle component analysis of gene expression data
  to generate table of possible covariates.
}
\author{
  Peter Humburg
}

