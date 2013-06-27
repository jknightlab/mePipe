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
#' @param minFDR Minimum FDR for associations to be considered.
#' @param minR Minimum R-squared at which two eSNPs will be assumed to have no independent effect.
#' @param exprOpt Options for reading of gene expression file. 
#' @param genoOpt Options for reading of gene expression file.
#' @param covOpt Options for reading of gene expression file.
#' @param ... Further arguments to runME.
#' @return A \code{data.frame} with eSNPs grouped into peaks.
#' 
#' @author Peter Humburg
#' @export
getMultiPeak <- function(hits, pvalue=1e-6, expression, genotype, covariate, minFDR,
		minR, exprOpt=getOptions(), genoOpt=getOptions(), covOpt=getOptions(), output, ...){
	## create data objects
	## gene expression
	gene <- NULL
	if(is(expression, "SlicedData")){
		gene <- expression
	} else{
		gene <- loadData(expression, exprOpt)
	}
	## genotypes
	snps <- NULL
	if(is(genotype, "SlicedData")){
		snps <- genotype
	} else{
		snps <- loadData(genotype, genoOpt)
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
	hitsLD <- getLDpairs(hits, snps, minR=minR, minFDR=minFDR)
	hits <- hitsLD$groups
	ldTable <- hitsLD$proxies
	candidates <- hits
	candidates$finalPvalue <- as.numeric(NA)
	complete <- subset(candidates, FALSE)
	repeat{
		## only proceed with genes that have more than one SNP associated with them
		snpCount <- table(candidates$gene)
		complete <- rbind(complete, 
				subset(candidates, gene %in% names(snpCount)[snpCount <= depth]))
		candidates <- subset(candidates, !gene %in% complete$gene)
		if(nrow(candidates) == 0) break
		ans <- Rsge::sge.parLapply(unique(candidates$gene), .submitMultiPeak, candidates, 
				depth, pvalue, gene, snps, cvrt, ...)
		primary <- lapply(ans, '[[', "primary")
		secondary <- Reduce(rbind, lapply(ans, '[[', "secondary"))
		secondaryLD <- getLDpairs(secondary, snps, minR=minR, minFDR=1)
		secondary <- secondaryLD$groups
		ldTable <- rbind(ldTable, secondaryLD$proxies)
		
		## update candidates
		for(p in primary){
			idx <- which(candidates$snps == as.character(p$snps[1]) & 
							candidates$gene == as.character(p$gene[1]))
			candidates$finalPvalue[idx] <- p$pvalue[1]
		}
		if(!is.null(secondary) && nrow(secondary)){
			for(g in unique(candidates$gene)){
				secondPeak <- subset(secondary, gene == g)
				explained <- c(as.character(secondPeak$snps), 
						unlist(strsplit(as.character(secondPeak$others), ",")))
				candidates <- subset(candidates, !(gene == g & snps %in% explained))
				peak <- subset(hits, snps %in% secondPeak$snps & gene == g)
				peak$others <- secondPeak$others
				peak$Rsquared <- secondPeak$Rsquared
				peak$finalPvalue <- secondPeak$pvalue
				candidates <- rbind(candidates, peak)
			}
		} else {
			for(g in unique(candidates$gene)){
				remain <- subset(candidates, gene == g)
				peaks <- subset(remain, !is.na(finalPvalue))
				others <- subset(remain, is.na(finalPvalue))
				if(nrow(peaks) == 0){
					peaks <- others[1, ]
					others <- others[-1,]
				}
				if(nrow(others)){
					r2 <- sapply(as.character(others$snps), 
							function(s, p, x){
								tab <- subset(x, snp1 %in% p$snps & snp2==s)
								tab <- tab[match(p$snps, tab$snp1),]
								ans <- tab$Rsquared
								ans
							}, peaks, ldTable)
					if(length(dim(r2)) == 0){
						r2 <- matrix(r2, ncol=nrow(others))
					}
					idx <- apply(r2, 2, which.max)
					for(i in 1:length(idx)){
						candIdx <- which(candidates$gene == g & 
										candidates$snps == peaks$snps[idx[i]])
						base <- candidates[candIdx, ]
						if(is.na(base$others)){
							candidates$others[candIdx] <- as.character(others$snps[i])
							candidates$Rsquared[candIdx] <- r2[idx[i], i] 
						} else {
							candidates$others[candIdx] <- paste(base$others[!is.na(base$others)], 
									others$snps[i], sep=",")
							candidates$Rsquared[candIdx] <- paste(base$Rsquared[!is.na(base$Rsquared)], 
									r2[idx[i], i], sep=",")
						}
					}
					exclude <- others$snps
					candidates <- subset(candidates, !(gene == g & snps %in% exclude))
				}
			}
		}
		depth <- depth+1
	}
	if(!missing(output)){
		write.table(ldTable, file=paste0(output, "_LDtable"), row.names=FALSE,
				quote=FALSE, sep="\t")
	}
	## TODO: update FDR estimates
	
	complete
}

#' @author Peter Humburg
#' @keywords internal
.submitMultiPeak <- function(current, hits, depth, pvalue, expression, genotype, covariate, ...){
	hits <- subset(hits, gene == current)
	
	## restrict gene expression data to current gene
	tmpExpr <- SlicedData$new()
	tmpExpr$CreateFromMatrix(expression$FindRow(current)$row)
	expression <- tmpExpr
	
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
	secondaryPeak <- peakUpdate <- primary <- data.frame(snps=character(), gene=character(), 
			statistic=numeric(), pvalue=numeric(), FDR=numeric(), stringsAsFactors=FALSE)
	peakUpdate$secondary <- character()
	## Fit model including peak SNP and one other candidate
	tmp1 <- tempfile(pattern=paste(current, hits$snps[1], "secondaries", "", sep="_"), 
			tmpdir=".", fileext=".tmp")
	tmpCov <- Reduce(combineSlicedData, c(covariate, snps[1:depth]))
	me1 <- runME(expression, Reduce(combineSlicedData, snps[-(1:depth)]), 
			tmpCov, output=tmp1, threshold=pvalue, cisThreshold=0, cis=0, ...)
	unlink(paste0(tmp1, "*"))
	if(nrow(me1$all$eqtls)){
		snps <- c(snps[1:depth] , snps[as.character(me1$all$eqtls$snps)])
		for(i in (depth+1):length(snps)){
			for(j in 1:depth){
				tmp2 <- tempfile(pattern=paste(current, hits$snps[i], paste(hits$snps[j], 
										collapse="_"), "", sep="_"), 
						tmpdir=".", fileext=".tmp")
				tmpCov <- Reduce(combineSlicedData, c(covariate, snps[c((1:depth)[-j], i)]))
				me2 <- runME(expression, snps[[j]], tmpCov, output=tmp2, threshold=1, 
						cisThreshold=0, cis=0, ...)
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
	secondaryPeak <- me1$all$eqtls[order(me1$all$eqtls$pvalue),]
	primary <-subset(peakUpdate, secondary == as.character(secondaryPeak$snps[1]))
	
	list(primary=primary, secondary=secondaryPeak)
}
