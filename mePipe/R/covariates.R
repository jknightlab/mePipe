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
		Rsge::sge.run(.submitAllCovariates, exprOpt = exprOpt, expression = expression, 
				output = output, sep = sep)
	} else{
		.submitAllCovariates(exprOpt = exprOpt, expression = expression, output = output, sep = sep)
	}
}

#' @author Peter Humburg
.submitAllCovariates <- function(exprOpt, expression, output, sep) {
	gene <- SlicedData$new()
	gene$fileDelimiter <- exprOpt$sep
	gene$fileOmitCharacters = exprOpt$missing
	gene$fileSkipRows = exprOpt$rowskip
	gene$fileSkipColumns = exprOpt$colskip
	gene$fileSliceSize = exprOpt$slice
	gene$LoadFile(expression)
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
				model = model, exclude = exclude, otherCov = otherCov, covOut = covOut)
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
	cvrt <- lapply(files, function(f){
				cvrt <- SlicedData$new();
				cvrt$fileDelimiter <- options$sep
				cvrt$fileOmitCharacters <- options$missing
				cvrt$fileSkipRows <- options$rowskip
				cvrt$fileSkipColumns <- options$colskip
				cvrt$fileSliceSize <- options$slice
				cvrt$LoadFile(f)
			})
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
		covThreshold = 0.01, covOpt = getOptions(), output = "covSelect", ...){
	dir.create(output, showWarnings = FALSE)
	
	covCount <- Rsge::sge.parLapply(candidates, .submitChooseCov, covOpt = covOpt, covariate = covariate, 
			expression = expression, genotype = genotype, output = output, 
			covThreshold = covThreshold, ..., 
			njobs=length(candidates))
	covCount <- Reduce(cbind, covCount)
	
	if(max(covCount[,1]) > 0){
		bestCov <- candidates[which.max(covCount[1,])]
	} else{
		bestCov <- candidates[which.max(covCount[2,])]
	}
	
	list(covariates=candidates, eqtls=list(significant=covCount[1,], all=covCount[2,]), 
			best=file.path(output, paste("me_", bestCov, "_covariates.eQTL", sep="")))
}

#' Carry out computations to choose covariates.
#' This function is called by \code{\link{chooseCov}} to do the actual work
#' (possibly distributed on an SGE cluster)
#' @param candidates Numeric vector listing the numbers of covariates that should be evaluated.
#' @param covOpt List of options to use for reading of covariate data.
#' @param covariate Character vector of file names for covariate files (see details).
#' @param expression Name of file containing the gene expression data.
#' @param genotype Name of file containing the genotyping data.
#' @param output Name of output directory
#' @param ... Addition parameters for \code{\link{runME}}
#' @param covThreshold FDR threshold to use when determining the number of significant eQTLs.
#' @rdname submitChooseCov
#' @author Peter Humburg
.submitChooseCov <- function(candidates, covOpt, covariate, expression, genotype, output, 
		..., covThreshold) {
	
	selectFrom <- as.matrix(loadCovariates(covariate[1], covOpt))
	if(any(candidates > nrow(selectFrom))){
		warning("Cannot use more than ", nrow(selectFrom), " covariates.")
		candidates <- candidates[candidates <= nrow(selectFrom)]
		if(length(candidates) == 0){
			stop("Invalid number of covariates selected. Please choose numbers between 0 and ", 
					nrow(selectFrom))
		}
	}
	otherCov <- SlicedData$new()
	if(length(covariate) > 1){
		otherCov <- loadCovariates(covariate[-1], covOpt)
	}
	maxCount <- 0
	covCount <- data.frame(all=integer(), total=integer())
	bestCov <- candidates[1]
	for(n in candidates){
		selectedCov <- SlicedData$new()
		if(n > 0){
			selectedCov <- selectedCov$CreateFromMatrix(selectFrom[1:n, , drop=FALSE])
		}
		if(length(covariate) > 1)
			selectedCov <- combineSlicedData(selectedCov, otherCov)
		me <- runME(expression, genotype, selectedCov, 
				output=file.path(output, paste("me_", n, "_covariates.eQTL", sep="")), ...)
		count <- 0
		fullCount <- 0
		if(!is.null(me$all$eqtls)){
			count <- length(unique((subset(me$all$eqtls, FDR < covThreshold)$gene)))
			fullCount <- length(unique((me$all$eqtls$gene)))
		}else{
			if(!is.null(me$trans$eqtls)){
				count <- length(unique(subset(me$trans$eqtls, FDR < covThreshold)$gene))
				fullCount <- length(unique((me$trans$eqtls$gene)))
			}
			if(!is.null(me$cis$eqtls)){
				count <- count + length(unique((subset(me$cis$eqtls, FDR < covThreshold)$gene)))
				fullCount <- fullCount + length(unique((me$cis$eqtls$gene)))
			}
		}
		covCount <- rbind(as.integer(count), as.integer(fullCount))
	}
	covCount
}