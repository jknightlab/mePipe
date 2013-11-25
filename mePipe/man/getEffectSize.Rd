\name{getEffectSize}
\alias{getEffectSize}
\title{Compute effect size for all SNPs with significant associations.}
\usage{
  getEffectSize(hits, expression, genotype, covariate,
    minFDR)
}
\arguments{
  \item{hits}{A \code{data.frame} with eQTL results as
  produced by \code{Matrix-eQTL}}

  \item{expression}{Gene expression values.}

  \item{genotype}{Genotypes.}

  \item{covariate}{Additional covariates to include in the
  model.}

  \item{minFDR}{FDR threshold used to determine whether an
  association is significant.}
}
\value{
  A \code{data.frame} similar to \code{hits} with two
  additional columns. Column \code{var.explained} gives the
  variance explained by the linear model that was fitted
  for the SNP and column \code{effect.size} gives the
  variance explained by the genotype within this model (as
  measured by the partial r^2).
}
\description{
  Compute effect size for all SNPs with significant
  associations.
}
\author{
  Peter Humburg
}

