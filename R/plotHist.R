#' Plot histogram of p-values
#' @param me Object returned by Matrix_eQTL_main.
#' @return Called for its side effect.
#' 
#' @author Peter Humburg
#' @export
plotHist <- function(me){
	cnts <- me$all$hist.counts
	bins <- me$all$hist.bins
	ntst <- me$all$ntests
	centers <- 0.5 * (bins[-1L] + bins[-length(bins)])
	density <- 0.5 / (bins[-1L] - centers) * cnts / ntst
	ntext <- paste('All', formatC(ntst, big.mark=',', format = 'f', digits = 0), 'gene-SNP pairs')
	r <- structure(list(breaks = bins, counts = cnts, density = density,
					mids = centers, equidist = FALSE), class = "histogram")
	plot(r, main = ntext, ylab = 'Density', xlab = 'P-values')
}
