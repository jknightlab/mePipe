\name{plotCov}
\alias{plotCov}
\title{Plot change of identified eQTLs for different number of covariates}
\usage{
  plotCov(candidates, eqtls, selected, ...)
}
\arguments{
  \item{candidates}{Numerical vector with the number of
  covariates used in each iteration}

  \item{eqtls}{\code{data.frame} with the number of eQTLs
  identified in each iteration. Different columns give
  results for different settings.}

  \item{selected}{Numerical vector with the number of
  covariates chosen for each column of \code{eqtls}.}

  \item{...}{Further arguments to plot.}
}
\value{
  Called for its side effect
}
\description{
  Plot change of identified eQTLs for different number of
  covariates
}
\author{
  Peter Humburg
}

