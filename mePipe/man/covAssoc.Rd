\name{covAssoc}
\alias{covAssoc}
\title{Identify covariates that are significantly associated with genotype}
\usage{
  covAssoc(genotype, covariate, otherCov, output, covOut,
    genoOpt = getOptions(), covOpt = getOptions(),
    model = c("linear", "anova", "cross"),
    threshold = 1e-05, exclude = 0.001)
}
\arguments{
  \item{genotype}{Name of file with genotyping data.}

  \item{covariate}{Name of file with covariates.}

  \item{otherCov}{Name of file(s) with manually selected
  covariates}

  \item{output}{Name to use for output file. Results of the
  association analysis will be written to this file.}

  \item{covOut}{Output file for filtered covariates.}

  \item{genoOpt}{Options for reading of genotype file.}

  \item{covOpt}{Options for reading of covariate file.}

  \item{model}{Type of model to use.}

  \item{threshold}{p-value threshold up to which variants
  should be reported.}

  \item{exclude}{FDR threshold below which associations are
  considered to be significant.}
}
\value{
  A character vector with the names of covariates that had
  significant associations.
}
\description{
  Identify covariates that are significantly associated
  with genotype
}
\details{
  All covariates in the input file that produce significant
  associations with SNPs are removed and the remaining ones
  are written to \code{covOut}.
}
\author{
  Peter Humburg
}

