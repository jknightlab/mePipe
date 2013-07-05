\name{getMultiPeak}
\alias{getMultiPeak}
\title{Identify secondary peaks and proxy SNPs
Test the set of eSNPs associated with each gene to identify secondary peaks
that are due to an independent effect rather than LD}
\usage{
  getMultiPeak(hits, pvalue = 1e-06, expression, genotype,
    covariate, minFDR, minR, exprOpt = getOptions(),
    genoOpt = getOptions(), covOpt = getOptions(), output,
    ...)
}
\arguments{
  \item{hits}{\code{data.frame} with SNP-gene pairs and
  associsted p-value, as produced by Matrix-eQTL.
  Additional columns may be included. In particular, output
  from \code{getLDpairs} is accepted.}

  \item{pvalue}{P-value threshold for association between
  secondary SNP and gene expression.}

  \item{expression}{Gene expression data. This can either
  be the file name or a \code{slicedData} object.}

  \item{genotype}{Genotype data. This can either be the
  file name or a \code{slicedData} object.}

  \item{covariat}{List of covariates to use. These can
  either be given as file names or a \code{slicedData}
  objects.}

  \item{minFDR}{Minimum FDR for associations to be
  considered.}

  \item{minR}{Minimum R-squared at which two eSNPs will be
  assumed to have no independent effect.}

  \item{exprOpt}{Options for reading of gene expression
  file.}

  \item{genoOpt}{Options for reading of gene expression
  file.}

  \item{covOpt}{Options for reading of gene expression
  file.}

  \item{...}{Further arguments to runME.}
}
\value{
  A \code{data.frame} with eSNPs grouped into peaks.
}
\description{
  Identify secondary peaks and proxy SNPs Test the set of
  eSNPs associated with each gene to identify secondary
  peaks that are due to an independent effect rather than
  LD
}
\author{
  Peter Humburg
}

