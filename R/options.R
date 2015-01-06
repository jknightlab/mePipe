# TODO: Add comment
# 
# Author: humburg
###############################################################################


#' Generate list of default options
#' @param sep Character used to seperate entries in input file.
#' @param missing Character sequence used to signal missing data in input file.
#' @param rowskip Number of rows to skip at beginning of input file.
#' @param colskip Number of columns to skip at the start of each row.
#' @param slice Size of slices to read.
#' @param ... Additional options.
#' @return A list with one entry for each option.
#' 
#' @author Peter Humburg
#' @export
getOptions <- function(sep = '\t', missing = 'NA',	rowskip = 1L, colskip = 1L, 
		slice = 2000L, ...){
	others <- list(...)
	c(list(sep = sep, missing = missing, rowskip = rowskip, colskip = colskip,
					slice = slice, others))
}
