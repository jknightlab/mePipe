#' Run Matrix-eQTL analysis
#' @param expression File name of gene expression data or a \code{SlicedData} object with expression data.
#' @param genotype Name of file with genotyping data or a \code{SlicedData} object with genotypes.
#' @param covariate Covariates. This may be a \code{matrix}, 
#' \code{data.frame} or vector of file names.
#' @param error Name of file with error covariance matrix. 
#' @param snpsPos Name of file with SNP positions.
#' @param genePos Name of file with gene position.
#' @param output Name of output file.
#' @param cisOutput Name of output file for cis associations. 
#' @param exprOpt List of options for reading of gene expression data. 
#' @param genoOpt List of options for reading of genotype data.
#' @param covOpt List of options for reading of covariate files.
#' @param threshold P-value threshold up to which associations should be listed
#' in the output file.
#' @param cisThreshold P-value threshold up to which cis associations should be listed
#' in the output file.  
#' @param model The model to use in the analysis.
#' @param cis The distance (in bp) up to which a SNP should be concidered to be in cis 
#' of a gene.
#' @param bins Number of p-value bins to use for plotting. 
#' @param verbose Flag indicating whether additional progress information should be printed.
#' @param qqplot Flag indicating whether a QQ plot of p-values should be generated 
#' rather than a histogram.
#' @return The object returned by \code{\link[MatrixEQTL]{Matrix_eQTL_main}}
#' @note Before calling this function a call to \code{sge.setDefaultOptions} has to be made.
#' @author Peter Humburg
#' @import MatrixEQTL
#' @import Rsge
#' @export
runME <- function(expression, genotype, covariate, error, snpsPos, genePos,
		output="", cisOutput=paste(output, "cis", sep="."), 
		exprOpt=getOptions(), genoOpt=getOptions(),	covOpt=getOptions(), 
		threshold=1e-5, cisThreshold=1e-3, model=c('linear', 'anova', 'cross'), cis=1e6, bins=0, 
		verbose=FALSE, qqplot=FALSE){
	model <- match.arg(model)
	model <- switch(model,
			linear = modelLINEAR,
			anova = modelANOVA,
			cross = modelLINEAR_CROSS)
	if(missing(snpsPos) || missing(genePos)){
		if(cis > 0){
			warning("SNP and gene positions are required to distinguish between cis- and trans-effects.")
		}
		snpsPos <- NULL
		genePos <- NULL
		cis <- -1
	}
	if(cis < 0) cisThreshold <- 0
	
	if(sge.getOption("sge.use.cluster")){
		Rsge::sge.run(.submitRunME, error = error, expression = expression, exprOpt = exprOpt, 
				genotype = genotype, genoOpt = genoOpt, covariate = covariate, 
				combineSlicedData = combineSlicedData, cis = cis, snpsPos = snpsPos, 
				genePos = genePos, output = output, threshold = threshold, model = model, 
				verbose = verbose, cisOutput = cisOutput, cisThreshold = cisThreshold, 
				bins = bins, qqplot = qqplot)
	} else{
		.submitRunME(error = error, expression = expression, exprOpt = exprOpt, genotype = genotype, 
				genoOpt = genoOpt, covariate = covariate, combineSlicedData = combineSlicedData, 
				cis = cis, snpsPos = snpsPos, genePos = genePos, output = output, threshold = threshold, 
				model = model, verbose = verbose, cisOutput = cisOutput, cisThreshold = cisThreshold, 
				bins = bins, qqplot = qqplot)
	}
}

#' @author Peter Humburg
.submitRunME <- function(error, expression, exprOpt, genotype, genoOpt, covariate, 
		combineSlicedData, cis, snpsPos, genePos, output, threshold, model, verbose, 
		cisOutput, cisThreshold, bins, qqplot) {
	errorCov <- numeric()
	if(!missing(error) && error != ""){
		errorCov <- as.matrix(read.table(error))
	}
	
	## create data objects
	## gene expression
	gene <- NULL
	if(is(expression, "SlicedData")){
		gene <- expression
	} else{
		gene <- SlicedData$new();
		gene$fileDelimiter <- exprOpt$sep
		gene$fileOmitCharacters = exprOpt$missing
		gene$fileSkipRows = exprOpt$rowskip
		gene$fileSkipColumns = exprOpt$colskip
		gene$fileSliceSize = exprOpt$slice
		gene$LoadFile(expression)
	}
	## genotypes
	snps <- NULL
	if(is(genotype, "SlicedData")){
		snps <- genotype
	} else{
		snps <- SlicedData$new();
		snps$fileDelimiter <- genoOpt$sep
		snps$fileOmitCharacters <- genoOpt$missing
		snps$fileSkipRows <- genoOpt$rowskip
		snps$fileSkipColumns <- genoOpt$colskip
		snps$fileSliceSize <- genoOpt$slice
		snps$LoadFile(genotype)
	}
	## covariates (if any)
	if(!missing(covariate)){
		if(is(covariate, "SlicedData")){
			cvrt <- covariate
		} else{
			cvrt <- lapply(covariate,
					function(covariate){
						cvrt <- SlicedData$new()
						if(!is.null(covariate)) {
							if(is.data.frame(covariate)){
								covariate <- as.matrix(covariate)
							}
							if(is.matrix(covariate) && nrow(covariate) > 0){
								cvrt$CreateFromMatrix(covariate)
							} else if(is.character(covariate) && covariate != ""){
								cvrt <- loadCovariates(covariate, covOpt)
							}
							if(is(covariate, "SlicedData")){
								cvrt <- covariate
							}
						}
						cvrt
					}
			)
			cvrt <- Reduce(combineSlicedData, cvrt)
		}
	} else{
		cvrt <- SlicedData$new()
	}
	
	if(cis > 0){
		snpspos <- read.table(snpsPos, header=TRUE, stringsAsFactors=FALSE)
		genepos <- read.table(genePos, header=TRUE, stringsAsFactors=FALSE)
	}
	
	me <- Matrix_eQTL_main(
			snps = snps,
			gene = gene,
			cvrt = cvrt,
			output_file_name = output,
			pvOutputThreshold = threshold,
			useModel = model, 
			errorCovariance = errorCov, 
			verbose = verbose,
			output_file_name.cis = cisOutput,
			pvOutputThreshold.cis = cisThreshold,
			snpspos = snpspos,
			genepos = genepos,
			cisDist = cis,
			pvalue.hist = if(length(bins) > 1 || bins > 0) bins else FALSE
	)
	
	save(me, file=paste(output, 'rdata', sep='.'))
	if(length(bins) > 1 || bins > 0){
		pdf(paste(output, 'pdf', sep='.'))
		if(qqplot){
			plotQQ(me)
		}
		else{
			plotHist(me)
		}
		dev.off()
	}
	me
}