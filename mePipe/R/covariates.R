# Generating table of covariates from PCA components
# 
# Author: Peter Humburg
###############################################################################


#' Run principle component analysis of gene expression data to generate 
#' table of possible covariates.
#' @param expression Name of file containing the gene expression data.
#' @param exprOpt Options for reading of gene expression data.
#' @param output Name of file to which the principle components should be saved (optional).
#' @param sep Character used to separate entries in output file.
#' @return \code{data.frame} with all principle components
#' 
#' @author Peter Humburg
#' @export
allCovariates <- function(expression, exprOpt = getOptions(), output, sep = "\t"){
	
	## load gene expression data
	if(!is.null(sge.getOption("sge.use.cluster")) && sge.getOption("sge.use.cluster")){
		Rsge::sge.run(.submitAllCovariates, exprOpt=exprOpt, expression=expression, 
				output=output, sep=sep, packages=.getPackageNames())
	} else{
		.submitAllCovariates(exprOpt = exprOpt, expression = expression, output = output, sep = sep)
	}
}

#' @author Peter Humburg
.submitAllCovariates <- function(exprOpt, expression, output, sep) {
	gene <- loadData(expression, exprOpt)
	gene$SetNanRowMean()
	
	## combine into one big matrix for PCA
	gene$CombineInOneSlice()
	gene <- gene$getSlice(1)
	
	## run pca
	pca <- t(prcomp(t(gene), center=TRUE, scale = TRUE)$x)
	names <- unlist(strsplit(readLines(expression, 1), exprOpt$sep))
	if(exprOpt$colskip > 0) {
		dummy <- matrix(nrow=nrow(pca), ncol=exprOpt$colskip - 1)
		pca <- cbind(rownames(pca),dummy,pca)
		if(ncol(pca) > length(names)){
			names <- c(rep("", ncol(pca) - length(names)))
		}
	}
	if(exprOpt$rowskip > 0) {
		dummy <- matrix(ncol=ncol(pca), nrow=exprOpt$rowskip - 1)
		pca <- rbind(names, dummy, pca)
	}
	
	if(!missing(output)) write.table(pca, file=output, quote = FALSE, sep = sep,
				row.names=FALSE, col.names=FALSE)
	invisible(pca)
}

#' Identify covariates that are significantly associated with genotype
#' @param genotype Name of file with genotyping data.
#' @param covariate Name of file with covariates.
#' @param otherCov Name of file(s) with manually selected covariates
#' @param output Name to use for output file. 
#' Results of the association analysis will be written to this file.
#' @param covOut Output file for filtered covariates.
#' @param genoOpt Options for reading of genotype file.
#' @param covOpt Options for reading of covariate file.
#' @param model Type of model to use.
#' @param threshold p-value threshold up to which variants should be reported.
#' @param exclude FDR threshold below which associations are considered to be
#' significant. 
#' @return A character vector with the names of covariates that had significant associations.
#' @details All covariates in the input file that produce significant associations with
#' SNPs are removed and the remaining ones are written to \code{covOut}.
#' 
#' @author Peter Humburg
#' @export
covAssoc <- function(genotype, covariate, otherCov, output, covOut, genoOpt = getOptions(),
		covOpt = getOptions(), model = c('linear', 'anova', 'cross'), 
		threshold = 1e-5, exclude = 1e-3){
	model <- match.arg(model)
	## the interaction test uses a linear model, use this as base for covariate selection
	if(model == "cross") model <- "linear"

	## identify all covariates that are signifcantly associated with genotypes
	if(!is.null(sge.getOption("sge.use.cluster")) && sge.getOption("sge.use.cluster")){
		Rsge::sge.run(.submitCovAssoc, covariate = covariate, genotype = genotype, 
				output = output, covOpt = covOpt, genoOpt = genoOpt, threshold = threshold, 
				model = model, exclude = exclude, otherCov = otherCov, covOut = covOut,
				packages=.getPackageNames())
	} else{
		.submitCovAssoc(covariate = covariate, genotype = genotype, output = output, 
				covOpt = covOpt, genoOpt = genoOpt, threshold = threshold, model = model, 
				exclude = exclude, otherCov = otherCov, covOut = covOut)
	}
	
}

#' @author Peter Humburg
.submitCovAssoc <- function(covariate, genotype, output, covOpt, genoOpt, threshold, model, 
		exclude, otherCov, covOut) {
	me <- runME(covariate, genotype, output = output, exprOpt=covOpt,
			genoOpt = genoOpt, threshold = threshold, model = model, cis = -1)
	assoc <- unique(as.character(me$all$eqtls$gene[me$all$eqtls$FDR <= exclude]))
	
	## identify automatically generated (PCA) covariates that are correlated with
	## pre-determined covariates
	if(!missing(otherCov) && !is.null(otherCov)){
		if(model != "linear"){
			warning("Cannot test for correlation between covariates when using model ", model, 
					". Using ", sQuote("linear"), " instead.")
		}
		me <- runME(covariate, loadCovariates(otherCov), output = output, exprOpt=covOpt,
				genoOpt = genoOpt, threshold = threshold, model = "linear", cis = -1)
	}
	assoc <- c(assoc, unique(as.character(me$all$eqtls$gene[me$all$eqtls$FDR <= exclude])))
	
	if(!missing(covOut)){
		cov <- SlicedData$new()
		cov$fileDelimiter <- covOpt$sep
		cov$fileOmitCharacters = covOpt$missing
		cov$fileSkipRows = covOpt$rowskip
		cov$fileSkipColumns = covOpt$colskip
		cov$fileSliceSize = covOpt$slice
		cov$LoadFile(covariate)
		cov$CombineInOneSlice()
		
		covName <- cov$GetAllRowNames()
		sampleName <- unlist(strsplit(readLines(covariate, 1), covOpt$sep))
		cov <- cov$getSlice(1)
		
		if(covOpt$colskip > 0) {
			dummy <- matrix(nrow=nrow(cov), ncol=covOpt$colskip - 1)
			cov <- cbind(covName,dummy,cov)
			if(ncol(cov) > length(sampleName)){
				sampleName <- c(rep("", ncol(cov) - length(sampleName)))
			}
		}
		if(covOpt$rowskip > 0) {
			dummy <- matrix(ncol=ncol(cov), nrow=covOpt$rowskip - 1)
			cov <- rbind(sampleName, dummy, cov)
		}
		
		idx <- c(1:covOpt$rowskip, which(!covName %in% assoc) + covOpt$rowskip)
		cov <- cov[idx,]
		write.table(cov, file=covOut, quote=FALSE, sep=covOpt$sep, row.names=FALSE, col.names=FALSE)
	}
	assoc
}

#' Load covariates from multiple files
#' @param files Character vector with names of files containing covariates for the analysis. 
#' @param options List of options for reading of the covariate files.
#' @return A \code{\link[MatrixEQTL]{SlicedData}} object containing all covariates.
#' 
#' @author Peter Humburg
#' @export
loadCovariates <- function(files, options = getOptions()){
	files <- files[!sapply(files, is.null)]
	cvrt <- lapply(files, 
			function(f, opt) if(is(f, "SlicedData")) f else loadData(f, opt), 
			options)
	Reduce(combineSlicedData, cvrt)
}

#' Combine two \code{SlicedData} objects
#' @param x A \code{SlicedData} object
#' @param y A \code{SlicedData} object
#' @return A \code{SlicedData} object containg the data from \code{x} and \code{y}
#' 
#' @author Peter Humburg
combineSlicedData <- function(x, y){
	ans <- SlicedData$new()
	ans$CreateFromMatrix(rbind(as.matrix(x), as.matrix(y)))
}

#' Determine the optimal number of covariates.
#' @param expression Name of file containing the gene expression data.
#' @param genotype Name of file containing the genotyping data.
#' @param covariate Character vector of file names for covariate files (see details).
#' @param candidates Numeric vector listing the numbers of covariates that should be evaluated.
#' @param covThreshold FDR threshold to use when determining the number of significant eQTLs.
#' @param covOpt List of options to use for reading of covariate data. 
#' @param output Name of output directory
#' @param doCis Logical indicating whether covariates should be selected based on cis-associations.
#' @param doTrans Logical indicating whether covariates should be selected based on trans-associations.
#' @param ... Addition parameters for \code{\link{runME}}
#' @return A list with elements
#' \item{covariates}{Numeric vector with the number of covariates used in each iteration.}
#' \item{eqtls}{A list with components \code{significant} and \code{all}. Each consisting
#' of a numeric vector with the number of genes that have significant associations or any 
#' association with SNPs respectively.}
#' \item{best}{Base of file name(s) that contain the results of the run that produced the
#' most eQTLs.}
#' @details The argument \code{covariates} may be a vector of several file names. 
#' In this case the first one will be used to select varying numbers of covariates and
#' the remaining ones will always be used.
#' 
#' @author Peter Humburg
#' @import MatrixEQTL
#' @import Rsge
#' @export
chooseCov <- function(expression, genotype, covariate, candidates = seq(5, 50, by = 5), 
		covThreshold = 0.01, covOpt = getOptions(), output = "covSelect", doCis=FALSE, doTrans=FALSE,
		...){
	dir.create(output, showWarnings = FALSE)
	
	covCount <- Rsge::sge.parLapply(candidates, .submitChooseCov, covOpt = covOpt, covariate = covariate, 
			expression = expression, genotype = genotype, output = output, 
			covThreshold = covThreshold, ..., 
			packages=.getPackageNames(),
			njobs=length(candidates))
	covCount <- Reduce(function(x,y) 
				list(significant=rbind(x$significant, y$significant), total=rbind(x$total, y$total)), 
			covCount)
	
	bestCovIdx <- lapply(covCount, sapply, which.max)
	bestCov <- sapply(bestCovIdx, function(i) candidates[i])
	rownames(bestCov) <- c("all", "cis", "trans")
	bestCount <- mapply(function(x,i) {
				ans <- integer(length(i))
				names(ans) <- names(x)
				for(j in 1:length(i)){
					ans[j] <- x[i[j], j]
				}
				ans
			}, covCount, bestCovIdx)
	## flag indicating whether the total number of associations (rather than only the significant ones)
	## should be used
	useTotal <- apply(bestCount, 1, function(x) x[1] == 0)
	
	selection <- numeric()
	if(!doCis && !doTrans){
		selection <- c(all=bestCov["all", useTotal[1]+1])
	} 
	if(doCis){
		selection <- c(cis=bestCov["cis", useTotal[2] + 1])
	} 
	if(doTrans){
		selection <- c(selection, trans=bestCov["trans", useTotal[3] + 1])
	}
		
	best <- file.path(output, paste0("me_", selection, "_covariates.eQTL"))
	names(best) <- names(selection)
	if(doCis) best["cis"] <- paste(best["cis"], "cis", sep=".")
	list(covariates=candidates, eqtls=covCount, best=best, selected=selection)
}

#' Carry out computations to choose covariates.
#' This function is called by \code{\link{chooseCov}} to do the actual work
#' (possibly distributed on an SGE cluster)
#' @param n Numeric vector listing the numbers of covariates that should be evaluated.
#' @param covOpt List of options to use for reading of covariate data.
#' @param covariate Character vector of file names for covariate files (see details).
#' @param expression Name of file containing the gene expression data.
#' @param genotype Name of file containing the genotyping data.
#' @param output Name of output directory
#' @param ... Addition parameters for \code{\link{runME}}
#' @param covThreshold FDR threshold to use when determining the number of significant eQTLs.
#' @rdname submitChooseCov
#' @keywords internal
#' @author Peter Humburg
.submitChooseCov <- function(n, covOpt, covariate, expression, genotype, output, 
		..., covThreshold) {
	
	covCount <- list(significant=data.frame(all=integer(), cis=integer(), trans=integer()),
			total=data.frame(all=integer(), cis=integer(), trans=integer()))
	selectFrom <- as.matrix(loadCovariates(covariate[1], covOpt))
	if(n > nrow(selectFrom)){
		warning("Cannot use more than ", nrow(selectFrom), " covariates.")
		return(covCount)
	}
	if(length(n) == 0){
		stop("Invalid number of covariates selected. Please choose numbers between 0 and ", 
				nrow(selectFrom))
	}
	otherCov <- SlicedData$new()
	if(length(covariate) > 1){
		otherCov <- loadCovariates(covariate[-1], covOpt)
	}
	
	selectedCov <- SlicedData$new()
	if(n > 0){
		selectedCov <- selectedCov$CreateFromMatrix(selectFrom[1:n, , drop=FALSE])
	}
	if(length(covariate) > 1){
		if(n > 0){
			selectedCov <- combineSlicedData(selectedCov, otherCov)
		} else {
			selectedCov <- otherCov
		}
	}
	me <- runME(expression, genotype, selectedCov,
			output=file.path(output, paste("me_", n, "_covariates.eQTL", sep="")), ...)
	count <- c(all=0, cis=0, trans=0)
	fullCount <- c(all=0, cis=0, trans=0)
	if(!is.null(me$all$eqtls)){
		count["all"] <- length(unique((subset(me$all$eqtls, FDR < covThreshold)$gene)))
		fullCount["all"] <- length(unique((me$all$eqtls$gene)))
	}
	if(!is.null(me$trans$eqtls)){
		count["trans"] <- length(unique(subset(me$trans$eqtls, FDR < covThreshold)$gene))
		fullCount["trans"] <- length(unique((me$trans$eqtls$gene)))
	}
	if(!is.null(me$cis$eqtls)){
		count["cis"] <- length(unique((subset(me$cis$eqtls, FDR < covThreshold)$gene)))
		fullCount["cis"] <- length(unique((me$cis$eqtls$gene)))
	}
	if(count["all"] == 0) count["all"] <- sum(count)
	if(fullCount["all"] == 0) fullCount["all"] <- sum(fullCount)
	covCount$significant <- rbind(covCount$significant, as.integer(count))
	covCount$total <- rbind(covCount$total, as.integer(fullCount))
	colnames(covCount$significant) <- colnames(covCount$total) <- c("all", "cis", "trans")
	covCount
}