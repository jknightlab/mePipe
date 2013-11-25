\name{fitStats}
\alias{fitStats}
\title{Obtain model summaries for Matrix-eQTL model}
\usage{
  fitStats(expression, genotype, covariates)
}
\arguments{
  \item{expression}{vector of expression values}

  \item{genotype}{vector of genotypes}

  \item{\code{data.frame}}{of covariates}
}
\value{
  A list as produced by \code{summary.lm} with additional
  entry \code{pr2} giving the partial r^2 for all effects
  in the model.
}
\description{
  Obtain model summaries for Matrix-eQTL model
}
\author{
  Peter Humburg
}

