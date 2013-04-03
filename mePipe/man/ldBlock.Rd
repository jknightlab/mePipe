\name{ldBlock}
\alias{ldBlock}
\title{Identify region of high LD around a given SNP}
\usage{
  ldBlock(snp, geno, pos, dist = 500, window = 200, ...)
}
\arguments{
  \item{snp}{Name (typically the rs ID) of the SNP}

  \item{geno}{An object of class \code{"SnpMatrix"} with
  genotypes.}

  \item{pos}{A \code{"data.frame"} with columns
  \code{'chrom'} and \code{'pos'} giving the SNP locations.
  Row names correspond to SNP IDs. This is assumed to be
  sorted by position.}

  \item{dist}{Maximum allowed distance between adjacent
  SNPs within an LD block.}

  \item{window}{Size of window around \code{snp} (in number
  of SNPs) in which LD calculations are carried out.}

  \item{...}{Further arguments to
  \code{\link[trio]{findLDBlocks}}}
}
\value{
  A list with components \code{chrom}, \code{start},
  \code{end} and \code{snps}.
}
\description{
  Identify region of high LD around a given SNP
}
\author{
  Peter Humburg
}

