\name{getLDpairs}
\alias{getLDpairs}
\title{Identify all SNPs in LD with peak SNP
The resulting list of eSNPs for the same gene may be considered independent signals.}
\usage{
  getLDpairs(eqtls, genotype, minFDR = 0.05, maxP = NULL,
    minR = 0.8, genoOpt = getOptions())
}
\arguments{
  \item{eqtls}{A \code{data.frame} with eQTLs as generated
  by package \code{MatrixEQTL}.}

  \item{genotype}{A \code{SlicedData} object with
  genotypes.}

  \item{minFDR}{FDR threshold to use when filtering
  associations. Only SNPs with an FDR below this threshold
  will be considered.}

  \item{maxP}{Maximum p-value for correlation between to
  eSNPs that should still be considered to be in
  significant LD.}

  \item{minR}{Minimum R-squared at which two eSNPs will be
  assumed to have no independent effect.}

  \item{genoOpt}{List of options for reading of genotype
  data.}
}
\value{
  A list with components \describe{ \item{groups}{A
  \code{data.frame} containing all significant peak SNPs
  that are not in (high) LD with more significant SNPs.
  Columns are the same as in Matrix-eQTL output plus two
  additional columns (\code{others} and \code{Rsquared})
  with comma separated list of SNPs that had significant
  p-values but are considered to be proxies of the peak SNP
  and their R^2.} \item{proxies}{A \code{data.frame}
  containing the results of all proxy SNP computations} }
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

