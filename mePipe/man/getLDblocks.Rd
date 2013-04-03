\name{getLDblocks}
\alias{getLDblocks}
\title{Identify LD blocks around eQTLs}
\usage{
  getLDblocks(eqtls, genotype, pos, dist = 500,
    window = 200, genoOpt = getOptions(), verbose = FALSE,
    minFDR = 0.05, useCluster = FALSE)
}
\arguments{
  \item{eqtls}{A \code{data.frame} with eQTLs as generated
  by package \code{MatrixEQTL}.}

  \item{genotype}{A \code{SlicedData} object with
  genotypes.}

  \item{pos}{A \code{data.frame} with columns \code{chrom}
  and \code{pos} giving SNP positions.}

  \item{dist}{Maximum allowed distance between adjacent
  SNPs within an LD block.}

  \item{verbose}{Logical indicating whether additional
  progress information should be printed.}

  \item{window}{Size of window around each eQTL (in number
  of SNPs) in which LD calculations are carried out.}

  \item{genoOpt}{List of options for reading of genotype
  data.}

  \item{minFDR}{Minimum FDR required for a SNP to be
  reported as part of a block of associations for a given
  gene.}

  \item{useCluster}{Logical indicating whether jobs should
  be submited to SGE or run locally.}
}
\value{
  A \code{data.frame} with single SNP eQTLs replaced by
  associations between LD blocks and gene expression.
}
\description{
  Identify LD blocks around eQTLs
}
\author{
  Peter Humburg
}
\seealso{
  \code{\link{ldBlock}} \code{{\link{getOptions}}}
}

