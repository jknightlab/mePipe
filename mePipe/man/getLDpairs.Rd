\name{getLDpairs}
\alias{getLDpairs}
\title{Identify all SNPs in LD with peak SNP
The resulting list of eSNPs for the same gene may be considered independent signals.}
\usage{
  getLDpairs(eqtls, genotype, minFDR = 0.05, minR = 0.85,
    genoOpt = getOptions())
}
\arguments{
  \item{eqtls}{A \code{data.frame} with eQTLs as generated
  by package \code{MatrixEQTL}.}

  \item{genotype}{A \code{SlicedData} object with
  genotypes.}

  \item{minFDR}{FDR threshold to use when filtering
  associations. Only SNPs with an FDR below this threshold
  will be considered.}

  \item{minR}{Minimum R-squared required between two eSNPs
  for them to be considered in LD.}

  \item{genoOpt}{List of options for reading of genotype
  data.}
}
\value{
  A \code{data.frame} containing all significant peak SNPs
  that are not in (high) LD with more significant SNPs.
  Columns are the same as in Matrix-eQTL output plus two
  additional columns (\code{others} and \code{Rsquared})
  with comma separated list of SNPs that had significant
  p-values but are considered to be proxies of the peak SNP
  and their R^2.
}
\description{
  Identify all SNPs in LD with peak SNP The resulting list
  of eSNPs for the same gene may be considered independent
  signals.
}
\note{
  Unlike the raw Matrix-eQTL output all associations with
  an FDR below the threshold indicated by \code{minFDR}
  will be absent from the output produced by this function.
}
\author{
  Peter Humburg
}

