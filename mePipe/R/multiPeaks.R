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
	
	complete <- Rsge::sge.parLapply(unique(as.character(hits$gene)), .submitMultiPeak, 
			hits, pvalue, gene, snps, cvrt, minR, minFDR, ..., packages=.getPackageNames())
	
	complete <- do.call(rbind, complete)
	complete <- complete[order(complete$pvalue), ]
	## TODO: update FDR estimates
	
	complete
}

#' @author Peter Humburg
#' @keywords internal
.submitMultiPeak <- function(current, hits, pvalue, expression, genotype, covariate, 
		minR, minFDR, ...){
	hits <- subset(hits, gene == current & FDR <= minFDR)
	hits$others <- NA
	hits$Rsquared <- NA
	
	depth <- 1
	hitsLD <- .computeLD(hits, genotype, current, maxP=NULL, minR=minR)
	hits[1,] <- hitsLD$groups
	ldTable <- hitsLD$proxies
	rm(hitsLD)
	hits$finalPvalue <- NA
	hits <- subset(hits, !snps %in% subset(ldTable, Rsquared >= minR)$snp2)
	
	## restrict gene expression data to current gene
	tmpExpr <- SlicedData$new()
	tmpExpr$CreateFromMatrix(expression$FindRow(current)$row)
	expression <- tmpExpr
	
	## extract all candidate SNPs (including primary peak)
	geno <- lapply(unique(hits$snps), function(x){
				ans <- SlicedData$new()
				ans$CreateFromMatrix(genotype$FindRow(x)$row)
				ans
			}
	)
	names(geno) <- unique(hits$snps)
	
	while(nrow(hits) > depth){
		genoIdx <- which(names(geno) %in% hits$snps[1:depth])
		
		## Fit model including peak SNP(s) and one other candidate
		tmp1 <- tempfile(pattern=paste(current, hits$snps[1], "secondaries", "", sep="_"), 
				tmpdir=".", fileext=".tmp")
		tmpCov <- do.call(combineSlicedData, c(covariate, geno[genoIdx]))
		tryCatch(
				## obtain list of SNPs that are still significant when controlling for
				## `depth` peak SNPs
				me1 <- runME(expression, do.call(combineSlicedData, 
								geno[as.character(hits$snps[-(1:depth)])]), 
						tmpCov, output=tmp1, threshold=pvalue, cisThreshold=0, cis=0, 
						cluster=FALSE, ...),
				error=function(e){
					if(grepl("Colinear", e$message)){
						me1 <- list()
						me1$all <- list()
						me1$all$eqtls <- data.frame(snps=character(), gene=character(),
								statistic=numeric(), pvalue=numeric(), FDR=numeric())
					} else{
						stop(e$message)
					}
				},
				finally=unlink(paste0(tmp1, "*"))
		)
		if(nrow(me1$all$eqtls)){
			me1$all$eqtls$others <- NA
			me1$all$eqtls$Rsquared <- NA
			me1$all$eqtls$finalPvalue <- NA
			## only keep SNPs for which p-value doesn't increase
			idx <- match(as.character(me1$all$eqtls$snps), as.character(hits$snps[-(1:depth)]))
			me1$all$eqtls <- subset(me1$all$eqtls, me1$all$eqtls$pvalue <= hits$pvalue[idx])
		}
		newPeak <- FALSE
		while(nrow(me1$all$eqtls) && !newPeak){
			## ensure that inclusion of the next most significant eSNP will
			## maintain significance of previous peaks
			peakUpdate <- data.frame(snps=character(), gene=character(), 
					statistic=numeric(), pvalue=numeric(), FDR=numeric(), stringsAsFactors=FALSE)
			peakUpdate$secondary <- character()
			snps <- geno[genoIdx]
			tmpCov <- do.call(combineSlicedData, c(covariate, 
							geno[as.character(me1$all$eqtls$snps[1])]))
			peakOK <- logical(depth)
			for(j in 1:depth){
				tmp2 <- tempfile(pattern=paste(current, me1$all$eqtls$snps[1], 
								paste(hits$snps[j], collapse="_"), "", sep="_"), 
						tmpdir=".", fileext=".tmp")
				
				tryCatch(
						me2 <- runME(expression, geno[genoIdx][[j]], tmpCov, 
								output=tmp2, threshold=1, cisThreshold=0, cis=0, 
								cluster=FALSE, ...),
						error=function(e){
							if(grepl("Colinear", e$message)){
								me2 <- list()
								me2$all <- list()
								me2$all$eqtls <- data.frame(snps=character(), gene=character(),
										statistic=numeric(), pvalue=numeric(), FDR=numeric())
							}else{
								stop(e$message)
							}
						},
						finally=unlink(paste0(tmp2, "*")))
				if(me2$all$eqtls$pvalue > pvalue || 
						me2$all$eqtls$pvalue > me1$all$eqtls$pvalue[1] ||
						me2$all$eqtls$pvalue > hits$pvalue[as.character(hits$snps) == as.character(me2$all$eqtls$snps)]){
					me1$all$eqtls <- me1$all$eqtls[-1,]
					break
				}
				peakOK[j] <- TRUE
				peakUpdate <- rbind(peakUpdate, me2$all$eqtls)
			}
			newPeak <- all(peakOK)
		}
		
		## update p-values for all SNPs that remain significant
		if(nrow(me1$all$eqtls)){
			idx <- match(as.character(me1$all$eqtls$snps), as.character(hits$snps))
			hits$finalPvalue[idx] <- me1$all$eqtls$pvalue
			idx <- match(as.character(peakUpdate$snps), as.character(hits$snps))
			hits$finalPvalue[idx] <- peakUpdate$pvalue
		}
		hits <- hits[order(hits$finalPvalue), ]
		
		## SNPs that are no longer significant
		remove <- subset(hits, !snps %in% c(as.character(me1$all$eqtls$snps),
						as.character(hits$snps[1:depth])))
		
		## compute LD between remaining SNPs and new peak
		if(nrow(me1$all$eqtls) + nrow(remove) > 0){
			secondaryLD <- .computeLD(rbind(me1$all$eqtls, remove), genotype, current,
					maxP=NULL, minR=minR)
			ldTable <- rbind(ldTable, secondaryLD$proxies)
		}
		## remove all SNPs in high LD with new peak
		hits <- subset(hits, !as.character(snps) %in% remove$snps & 
						!as.character(snps) %in% subset(secondaryLD$proxies, 
								snp1==me1$all$eqtls$snps[1] & Rsquared >= minR)$snp2)
		rm(secondaryLD)
		
		depth <- depth + 1
	}
	## for each peak, get list of proxy SNPs
	ldTable <- subset(ldTable, !snp2 %in% hits$snps & snp1 %in% hits$snps)
	assignedPeak <- by(ldTable, ldTable$snp2, function(x) {
				i <- which.max(x$Rsquared) 
				list(peak=x[i, "snp1"], Rsquared=x[i, "Rsquared"])
			})
	peaks <- sapply(assignedPeak, '[[', "peak")
	proxies <- tapply(names(assignedPeak), peaks, paste, collapse=",")
	r2 <- sapply(assignedPeak, '[[', "Rsquared")
	r2 <- tapply(r2, peaks, paste, collapse=",")
	idx <- match(names(proxies), as.character(hits$snps))
	hits$others[idx] <- paste(ifelse(is.na(hits$others[idx]), "", hits$others[idx]), 
			proxies, sep=",")
	hits$Rsquared[idx] <- paste(ifelse(is.na(hits$Rsquared[idx]), "", hits$Rsquared[idx]), 
			r2, sep=",")
	
	hits
}

#' @author Peter Humburg
#' @keywords internal
.submitMultiUpdate <- function(primary, secondary, candidates, snps, finalPvalue, ldTable) {
	idx <- which(candidates$snps == as.character(primary$snps[1]))
	candidates$finalPvalue[idx] <- primary$pvalue[1]
	
	if(!is.null(secondary) && nrow(secondary)){
		explained <- as.character(secondary$snps)
		if(!all(is.na(secondary$others))){
			explained <- c(explained, 
					unlist(strsplit(as.character(secondary$others), ",")[!is.na(secondary$others)]))
		}
		peak <- subset(candidates, snps %in% secondary$snps)
		if(nrow(peak)){
			peak$others <- secondary$others
			peak$Rsquared <- secondary$Rsquared
			peak$finalPvalue <- secondary$pvalue
			
			for(i in 1:nrow(peak)){
				candIdx <- match(as.character(peak$snps[i]), candidates$snps)
				if(!is.na(candidates$others[candIdx])){
					if(is.na(peak$others[i])){
						peak$others[i] <- candidates$others[candIdx]
						peak$Rsquared[i] <- candidates$Rsquared[candIdx]
					} else{
						peak$others[i] <- paste(candidates$others[candIdx], peak$others[i], sep=",")
						peak$Rsquared[i] <- paste(candidates$Rsquared[candIdx], peak$Rsquared[i], sep=",")
					}
				}
			}
			candidates <- subset(candidates, !snps %in% explained)
			candidates <- rbind(candidates, peak)
		}
	} else {
		peaks <- subset(candidates, !is.na(finalPvalue))
		others <- subset(candidates, is.na(finalPvalue))
		if(nrow(peaks) == 0){
			peaks <- others[1, ]
			others <- others[-1,]
		}
		if(nrow(others)){
			r2 <- sapply(as.character(others$snps), 
					function(s, p, x){
						tab <- subset(x, (snp1 %in% p$snps & snp2==s) | (snp2 %in% p$snps & snp1==s))
						i <- match(p$snps, tab$snp1)
						i[is.na(i)] <- match(p$snps[is.na(i)], tab$snp2)
						tab <- tab[i,]
						ans <- tab$Rsquared
						ans
					}, peaks, ldTable)
			if(length(dim(r2)) == 0){
				r2 <- matrix(r2, ncol=nrow(others))
			}
			idx <- apply(r2, 2, which.max)
			for(i in 1:length(idx)){
				candIdx <- which(candidates$snps == peaks$snps[idx[i]])
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
			candidates <- subset(candidates, !snps %in% exclude)
		}
	}
	candidates
}