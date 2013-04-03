\name{loadCovariates}
\alias{loadCovariates}
\title{Load covariates from multiple files}
\usage{
  loadCovariates(files, options = getOptions())
}
\arguments{
  \item{files}{Character vector with names of files
  containing covariates for the analysis.}

  \item{options}{List of options for reading of the
  covariate files.}
}
\value{
  A \code{\link[MatrixEQTL]{SlicedData}} object containing
  all covariates.
}
\description{
  Load covariates from multiple files
}
\author{
  Peter Humburg
}

