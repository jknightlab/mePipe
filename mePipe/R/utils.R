#' Get data for a named SNP from SlicedData object
getSNP <- function(data, snp){
	loc <- data$FindRow(snp)
	ans <- NULL
	if(!is.null(loc)){
		ans <- loc$row
	}
	ans
}

#' Extract named SNPs from SlicedData object and return them as aSnpMatrix
#' @importClassesFrom snpStats
toSnpMatrix <- function(data, snps){
	ans <- lapply(snps, getSNP, data=data)
	ans <- t(Reduce(rbind, ans))
	## only keep samples with complete information
	ans <- ans[apply(ans,1, function(x) !any(is.na(x))),]
	mode(ans) <- "raw"
	new("SnpMatrix", ans)	
}