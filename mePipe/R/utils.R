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

#' Extract named SNPs from SlicedData object and return them in the format
#' expected by CubeX
#' @author Peter Humburg
#' @keywords internal
toCubeX <- function(data, snps){
	mat <- lapply(snps, getSNP, data=data)
	mat <- Reduce(rbind, mat)
	t(apply(mat[-1, , drop=FALSE], 1, function(x, y){
				c(
						sum(x==0 & y==0, na.rm=TRUE), 
						sum(x==0 & y==1, na.rm=TRUE),
						sum(x==0 & y==2, na.rm=TRUE),
						sum(x==1 & y==0, na.rm=TRUE),
						sum(x==1 & y==1, na.rm=TRUE),
						sum(x==1 & y==2, na.rm=TRUE),
						sum(x==2 & y==0, na.rm=TRUE),
						sum(x==2 & y==1, na.rm=TRUE),
						sum(x==2 & y==2, na.rm=TRUE)
				)
			}, mat[1, , drop=FALSE]))
}