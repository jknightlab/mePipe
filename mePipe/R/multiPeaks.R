## Functions to identify secondary peaks
## Author: Peter Humburg, 2013
#########################################################

#' Identify secondary peaks and proxy SNPs
#' Test the set of eSNPs associated with each gene to identify secondary peaks
#' that are due to an independent effect rather than LD
#' @param hits \code{data.frame} with SNP-gene pairs and associsted p-value, 
#' as produced by Matrix-eQTL. Additional columns may be included. In particular,
#' output from \code{getLDpairs} is accepted.
#' @param pvalue P-value threshold for association between secondary SNP and 
#' gene expression.
#' @param expression Gene expression data. This can either be the file name 
#' or a \code{slicedData} object. 
#' @param genotype Genotype data. This can either be the file name 
#' or a \code{slicedData} object.
#' @param covariat List of covariates to use. These can either be given as file names 
#' or a \code{slicedData} objects.
#' @param minR Minimum R-squared at which two eSNPs will be assumed to have no independent effect.
#' @param exprOpt Options for reading of gene expression file. 
#' @param genoOpt Options for reading of gene expression file.
#' @param covOpt Options for reading of gene expression file.
#' @param ... Further arguments to runME.
#' @return A \code{data.frame} with eSNPs grouped into peaks.
#' 
#' @author Peter Humburg
#' @export
getMultiPeak <- function(hits, pvalue=1e-6, expression, genotype, covariate,
		minR, exprOpt=getOptions(), genoOpt=getOptions(), covOpt=getOptions(), ...){
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
	cvrt <- NULL
	if(!missing(covariate)){
		if(is(covariate, "SlicedData")){
			cvrt <- Reduce(combineSlicedData, covariate)
		} else{
			cvrt <- loadCovariates(covariate, covOpt)
		}
	}
	
	depth <- 1
	if(!"Rsquared" %in% names(hits)) hits <- getLDpairs(hits, snps, minR=minR)$group
	candidates <- hits
	candidates$finalPvalue <- as.numeric(NA)
	complete <- subset(candidates, FALSE)
	repeat{
		## only proceed with genes that have more than one SNP associated with them
		snpCount <- table(candidates$gene)
		complete <- rbind(complete, 
				subset(candidates, gene %in% names(snpCount)[snpCount == depth]))
		candidates <- subset(candidates, !gene %in% complete$gene)
		if(nrow(candidates) == 0) break
		ans <- Rsge::parLapply(unique(candidates$gene), .submitMultiPeak, candidates, 
				depth, pvalue, gene, snps, cvrt, ...)
		primary <- mapply(function(x,y) subset(x, secondary==y$snps[1])[,-ncol(x)], 
				primary, secondary, SIMPLIFY=FALSE)
		secondary <- Reduce(rbind, lapply(ans, '[[', "secondary"))
		secondary <- getLDpairs(secondary, snps, minR=minR)$group
		
		## update candidates
		for(p in primary){
			idx <- which(candidates$snps == p$snps[1] & candidates$gene == p$gene[1])
			candidates$finalPvalue[idx] <- p$pvalue[1]
		}
		for(g in unique(candidates$gene)){
			explained <- subset(secondary, gene == g)
			explained <- c(explained$snps, unlist(strsplit(explained$others, ",")))
			candidates <- subset(candidtates, gene == g & snps %in% explained)
		}
		candidates <- rbind(candidates, secondary)
		depth <- depth+1
	}
	
	## TODO: update FDR estimates
	
	complete
}

#' @author Peter Humburg
#' @keywords internal
.submitMultiPeak <- function(current, hits, depth, pvalue, expression, genotype, covariate, ...){
	hits <- subset(hits, gene == current)
	
	## restrict gene expression data to current gene
	expression$CreateFromMatrix(expression$FindRow(current)$row)
	
	ans <- data.frame(snps=character(), gene=character(), statistic=numeric(), 
			pvalue=numeric(), FDR=numeric(), stringsAsFactors=FALSE)
	
	## extract all candidate SNPs (including primary peak)
	snps <- lapply(unique(hits$snps), function(x){
				ans <- SlicedData$new()
				ans$CreateFromMatrix(genotype$FindRow(x)$row)
				ans
			}
	)
	names(snps) <- unique(hits$snps)
	secondary <- peakUpdate <- data.frame(snps=character(), gene=character(), 
			statistic=numeric(), pvalue=numeric(), FDR=numeric(), stringsAsFactors=FALSE)
	peakUpdate$secondary <- character()
	## Fit model including peak SNP and one other candidate
	tmp1 <- tempfile(pattern=paste(current, hits$snps[1], "secondaries", "", sep="_"), 
			tmpdir=".", fileext=".tmp")
	tmpCov <- Reduce(combineSlicedData, list(covariate, snps[1:depth]))
	me1 <- runME(expression, Reduce(combineSlicedData, snps[-(1:depth)]), 
			tmpCov, output=tmp1, threshold=pvalue, ...)
	unlink(paste0(tmp1, "*"))
	if(nrow(me1$all$eqtls)){
		snps <- c(snps[1:depth] , snps[as.character(me1$all$eqtls$snps)])
		for(i in (depth+1):length(snps)){
			for(j in 1:depth){
				tmp2 <- tempfile(pattern=paste(current, hits$snps[i], paste(hits$snps[j], 
										collapse="_"), "", sep="_"), 
						tmpdir=".", fileext=".tmp")
				tmpCov <- combineSlicedData(covariate, snps[c((1:depth)[-j], i)])
				me2 <- runME(expression, snps[[j]], tmpCov, output=tmp2, threshold=1, ...)
				if(me2$all$eqtls$pvalue > pvalue){
					warning("Peak SNP ", hits$snps[j], 
							" is not significant at requested significance threshold of ", pvalue)
				}
				eqtls <- me2$all$eqtls
				eqtls$secondary <- hits$snps[i]
				peakUpdate <- rbind(peakUpdate, eqtls)
				unlink(paste0(tmp2, "*"))
			}
		}
	}
	if(nrow(ans)){
		secondary <- me1$all$eqtls[order(me1$all$eqtls$pvalue),]
		primary <-subset(peakUpdate, secondary == secondary$snps[1])
	}
	list(primary=primary, secondary=secondary)
}
