\name{getOptions}
\alias{getOptions}
\title{Generate list of default options}
\usage{
  getOptions(sep = "\t", missing = "NA", rowskip = 1L,
    colskip = 1L, slice = 2000L, ...)
}
\arguments{
  \item{sep}{Character used to seperate entries in input
  file.}

  \item{missing}{Character sequence used to signal missing
  data in input file.}

  \item{rowskip}{Number of rows to skip at beginning of
  input file.}

  \item{colskip}{Number of columns to skip at the start of
  each row.}

  \item{slice}{Size of slices to read.}

  \item{...}{Additional options.}
}
\value{
  A list with one entry for each option.
}
\description{
  Generate list of default options
}
\author{
  Peter Humburg
}

