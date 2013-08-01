#' Get data for a named SNP from SlicedData object
#' @author Peter Humburg
#' @keywords internal
getSNP <- function(data, snp){
	loc <- data$FindRow(snp)
	ans <- NULL
	if(!is.null(loc)){
		ans <- loc$row
	} else{
		warning("Genotypes for ", snp, " not found")
	}
	ans
}

#' Extract rows from SlicedData object by name
#' @author Peter Humburg
#' @keywords internal
subsetRows <- function(data, names){
	if(!is(data, "SlicedData"))
		stop("Argument `data` has to be of class `SlicedData`")
	names <- unique(names)
	result <- lapply(names, function(x){
				ans <- SlicedData$new()
				ans$CreateFromMatrix(data$FindRow(x)$row)
				ans
			}
	)
	names(result) <- names
	result
}

#' Extract named SNPs from SlicedData object and return them in the format
#' expected by CubeX
#' @author Peter Humburg
#' @keywords internal
toCubeX <- function(data, snps){
	if(is(data, 'SlicedData')){
		mat <- lapply(snps, getSNP, data=data)
		missing <- sapply(mat, is.null) 
		if(any(missing)){
			warning("No data found for ", paste(snps[which(missing)], collapse=", "), 
					" in genotyping file")
		}
		mat <- do.call(rbind, mat)
	} else{
		if(nrow(data) != length(snps)){
			stop("Expected matrix with ", length(snps), " rows, got ", nrow(data))
		}
		mat <- data
	}
	
	ans <- matrix(nrow=0, ncol=9)
	if(length(mat) > 0 && nrow(mat) > 1){
		ans <- t(apply(mat[-1, , drop=FALSE], 1, function(x, y){
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
	ans
}

#' Load Dataset into SlicedData object
#' @author Peter Humnburg
#' @export
#' @keywords internal
loadData <- function(file, options){
	data <- SlicedData$new()
	data$fileDelimiter <- options$sep
	data$fileOmitCharacters <- options$missing
	data$fileSkipRows <- options$rowskip
	data$fileSkipColumns <- options$colskip
	data$fileSliceSize <- options$slice 
	data$LoadFile(file)
	
	data
}

#' @author Peter Humburg
#' @keywords internal
.getPackageNames <- function(){
	rev(c("base", "stats", "utils","methods","MatrixEQTL", "Rsge"))
}