#' Generate Q-Q plot of p-values
#' @param me Object returned by Matrix_eQTL_main
#' @return Called for its side effect.
#' 
#' @author Peter Humburg
#' @export
plotQQ <- function(me){
	
	cols <- 4:6
	types <- c('all', 'cis', 'trans')
	values <- lapply(me[types], 
			function(x){
				if(!is.null(x$eqtls) && nrow(x$eqtls) > 0){
					getPval(x)
				}
			}
	)
	idx <- !sapply(values, is.null)
	if(any(idx)){
		legendText <- types[idx]
		values <- values[idx]
		
		ntext <- 'Q-Q plot of all p-values'
		xlim <- c(0, max(sapply(values, function(x) max(x$xpvs))))
		ylim <- c(0, max(sapply(values, function(x) max(x$ypvs))))
		plot(0, xlab = '-log10(p-value), theoretical', 
				ylab = '-log10(p-value), actual', main = ntext, 
				type='n', xlim=xlim, ylim=ylim)
		
		for(i in 1:length(values)){
			points(c(values[[i]]$xpvs, 0), c(values[[i]]$ypvs, 0), pch=20, col=cols[i])
			lines(values[[i]]$xpos, values[[i]]$ypos, col = cols[i])
			lines(values[[i]]$ypos, values[[i]]$ypos, col = 'gray')
		}
		legend("topleft", legend=legendText, col=cols[1:length(values)], pch=20, lty=1)
	} else{
		warning("No data provided, skipping QQ-plot")
	}
	
}

getPval <- function(x){
	cnts <- x$hist.counts
	bins <- x$hist.bins
	ntst <- x$ntests
	
	cusu <- -log10( cumsum(cnts) / ntst)
	ypos <- -log10( bins[-1][is.finite(cusu)])
	xpos <- cusu[is.finite(cusu)]
	
	ypvs <- -log10( x$eqtls[ , 'pvalue'] )
	xpvs <- -log10( 1:length(ypvs) / ntst )
	
	list(ntst=ntst, xpos=xpos, ypos=ypos, xpvs=xpvs, ypvs=ypvs)
}