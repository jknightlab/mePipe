#' Plot change of identified eQTLs for different number of covariates
#' @param candidates Numerical vector with the number of covariates used in each iteration
#' @param eqtls \code{data.frame} with the number of eQTLs identified in each iteration. Different
#' columns give results for different settings.
#' @param selected Numerical vector with the number of covariates chosen for each column of \code{eqtls}.
#' @param ... Further arguments to plot.
#' 
#' @return Called for its side effect
#' @export
#' @author Peter Humburg
plotCov <- function(candidates, eqtls, selected, ...){
	matplot(candidates, eqtls, type="b", pch=1, col=1:ncol(eqtls), 
			xlab="Number of covariates",
			ylab = "Number of genes with eQTLs", ...)
	if(!missing(selected)){
		selectedIdx <- sapply(selected, function(x) which(candidates == x))
		selectedObs <- mapply('[', eqtls, selectedIdx)
		points(selected, selectedObs, col=1:ncol(eqtls), pch=16)
	}
	legend("topleft", legend=names(eqtls), fill=1:ncol(eqtls))
}