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
getMultiPeak <- function(hits, p.value=1e-6, expression, genotype, covariate, minFDR,
		minR, snppos, window, exprOpt=getOptions(), genoOpt=getOptions(), 
		covOpt=getOptions(), output, ...){
	## only use peaks that reach significance
	hits <- subset(hits, FDR <= minFDR)
	complete <- cbind(subset(hits, FALSE), others=character(), Rsquared=character(),
			finalPvalue=numeric())
	
	if(nrow(hits)){
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
		snppos <- read.table(snppos, header=TRUE, stringsAsFactors=FALSE)
		names(snppos) <- c("snp", "chrom", "pos")
		complete <- Rsge::sge.parLapply(unique(as.character(hits$gene)), .submitMultiPeak, 
				hits, p.value, gene, snps, cvrt, minR, minFDR, snppos, window, ..., 
				packages=.getPackageNames())
		complete <- do.call(rbind, complete)
		complete <- complete[order(complete[["gene"]], complete[["var.explained"]], 
						complete[["minPvalue"]]), ]
		
	}
	complete
}

#' @author Peter Humburg
#' @keywords internal
.submitMultiPeak <- function(current, hits, pvalue, expression, genotype, covariate, 
		minR, minFDR, snppos, window, verbose=FALSE, ...){
	hits <- subset(hits, gene == current)
	
	if(verbose) message("Processing gene ", current, " (", nrow(hits), " eSNPs)")
	
	if(nrow(hits) && any(hits$FDR <= minFDR)){
		peakPos <- subset(snppos, snp == as.character(hits$snps[1]))
		candidates <- subset(snppos, chrom == peakPos$chrom & pos <= peakPos$pos + window & 
						pos >= peakPos$pos - window & snp %in% genotype$GetAllRowNames())$snp
		
		## extract all candidate SNPs (including primary peak)
		geno <- subsetRows(genotype, candidates)
		
		## extract current gene
		expression <- subsetRows(expression, current)[[1]]
		
		## Fit model for all SNPs that were previously filtered out
		tmp1 <- tempfile(pattern=paste(current, hits$snps[1], "others", "", sep="_"), 
				tmpdir=".", fileext=".tmp")
		tryCatch(
				me1 <- runME(expression, do.call(combineSlicedData, 
								geno[candidates[!candidates %in% hits$snps]]), 
						covariate, output=tmp1, threshold=1, cisThreshold=0, cis=0, 
						cluster=FALSE, ...),
				finally=unlink(paste0(tmp1, "*"))
		)
		hits <- rbind(hits, me1$all$eqtls)
		hits$var.explained <- NA
		hits$improvement <- NA
		hits$adj.r.squared <- NA
		hits$others <- NA
		hits$Rsquared <- NA
		hits$finalStatistic <- hits$statistic
		hits$finalPvalue <- hits$pvalue
		hits$finalCoef <- NA
		hits$minPvalue <- hits$pvalue
		hits$maxCoef <- NA
		hits$explained <- FALSE
		
		hitsLD <- .computeLD(hits, genotype, current, maxP=NULL, minR=minR)
		hits[1,] <- hitsLD$groups
		ldTable <- hitsLD$proxies
		rm(hitsLD)
		
		hits$explained[hits$snps %in% subset(ldTable, Rsquared >= minR)$snp2] <- TRUE
		rownames(hits) <- hits$snps
		
		## restrict gene expression data to current gene
		tmpExpr <- SlicedData$new()
		tmpExpr$CreateFromMatrix(expression$FindRow(current)$row)
		expression <- tmpExpr
		
		## model without genetic effdects
		cov.df <- as.data.frame(t(as.matrix(covariate)))
		fit <- lm(t(expression[[1]]) ~ ., data=cov.df)
		baseline <- summary(fit)
		
		cov.df <- cbind(t(geno[[as.character(hits$snps[1])]][[1]]), cov.df)
		fit <- lm(t(expression[[1]]) ~ ., data=cov.df)
		fitSummary <- summary(fit)
		
		hits$var.explained[1] <- fitSummary$r.squared
		hits$improvement[1] <- fitSummary$r.squared - baseline$r.squared
		hits$adj.r.squared[1] <- fitSummary$adj.r.squared
		hits$maxCoef[1] <- coef(fit)[2]
		hits$explained[1] <- TRUE
		
		peaks <- as.character(hits$snps[1])
		
		## iterate until no more significant SNPs are found
		while(any(!hits$explained)){
			genoIdx <- which(names(geno) %in% peaks)
			genoRemain <- which(names(geno) %in% as.character(hits$snps)[!hits$explained])
			## Fit model including peak SNP(s) and one other candidate
			tmp1 <- tempfile(pattern=paste(current, hits$snps[1], "secondaries", "", sep="_"), 
					tmpdir=".", fileext=".tmp")
			tmpCov <- do.call(combineSlicedData, c(covariate, geno[genoIdx]))
			tryCatch(
					## obtain list of SNPs that are still significant when controlling for
					## `depth` peak SNPs
					me1 <- runME(expression, do.call(combineSlicedData, geno[genoRemain]), 
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
				me1$all$eqtls$var.explained <- NA
				me1$all$eqtls$improvement <- NA
				me1$all$eqtls$adj.r.squared <- NA
				me1$all$eqtls$others <- NA
				me1$all$eqtls$Rsquared <- NA
				me1$all$eqtls$finalStatistic <- me1$all$eqtls$statistic
				me1$all$eqtls$finalPvalue <- me1$all$eqtls$pvalue
				me1$all$eqtls$finalCoef <- NA
				me1$all$eqtls$minPvalue <- me1$all$eqtls$pvalue
				me1$all$eqtls$minCoef <- NA
				me1$all$eqtls$explained <- FALSE
				
				## update minimum p-value where appropriate
				idx <- match(as.character(hits$snps[-(1:length(peaks))]), 
						as.character(me1$all$eqtls$snps))
				hits$minPvalue[-(1:length(peaks))] <- pmin(
						me1$all$eqtls$pvalue[idx], hits$minPvalue[-(1:length(peaks))], 
						hits$pvalue[-(1:length(peaks))], na.rm=TRUE)
			} else {
				hits[names(geno)[genoRemain], "explained"] <- TRUE
				break
			}
			newPeak <- FALSE
			while(nrow(me1$all$eqtls) && !newPeak){
				## ensure that inclusion of the next most significant eSNP will
				## improve model fit
				## fit model including covariates, all previous peaks and 
				## the most significant remaining SNP
				tmpData <- cbind(t(geno[[as.character(me1$all$eqtls$snps[1])]][[1]]), 
						cov.df)
				fit <- lm(t(expression[[1]]) ~ ., data=tmpData)
				fitSummary <- summary(fit)
				
				changed <- hits[peaks, "minPvalue"] < fitSummary$coefficients[3:(length(peaks)+2), 4]
				testCoef <- 1
				if(any(changed)){
					testCoef <- car::linearHypothesis(fit, peaks[changed], 
							hits[peaks, "maxCoef"])[["Pr(>F)"]][2]
				}
				if(fitSummary$adj.r.squared > hits[peaks[length(peaks)], "adj.r.squared"] && 
						(!any(changed) || testCoef > 0.05)){
					peaks <- c(peaks, as.character(me1$all$eqtls$snps[1]))
					idx <- match(rev(peaks), as.character(hits$snps))
					
					hits$finalStatistic[idx] <- fitSummary$coefficients[2:(length(peaks)+1), 3]
					hits$finalPvalue[idx] <- fitSummary$coefficients[2:(length(peaks)+1), 4]
					hits$finalCoef[idx] <- coef(fit)[2:(length(peaks)+1)]
					hits$minPvalue[idx] <- pmin(hits$finalPvalue[idx], hits$minPvalue[idx], na.rm=TRUE)
					hits$maxCoef[idx] <- sign(hits$finalCoef[idx])*pmax(abs(hits$finalCoef[idx]), hits$maxCoef[idx], na.rm=TRUE)
					hits$var.explained[idx[1]] <- fitSummary$r.squared
					hits$adj.r.squared[idx[1]] <- fitSummary$adj.r.squared
					hits$improvement[idx[1]] <- hits$var.explained[idx[1]] - hits$var.explained[idx[2]]
					hits$explained[idx[1]] <- TRUE
					
					cov.df <- tmpData
					newPeak <- TRUE
				} else {
					hits$explained[hits$snps == as.character(me1$all$eqtls$snps[1])] <- TRUE 
				}
				me1$all$eqtls <- me1$all$eqtls[-1, ]
			}
			
			## update p-values for all SNPs that remain significant
			if(nrow(me1$all$eqtls)){
				idx <- match(as.character(me1$all$eqtls$snps), as.character(hits$snps))
				hits$finalPvalue[idx] <- me1$all$eqtls$pvalue
				hits$finalStatistic[idx] <- me1$all$eqtls$statistic
			}
			hits <- hits[order(hits$var.explained, hits$finalPvalue), ]
			
			## compute LD between remaining SNPs and new peak
			if(nrow(me1$all$eqtls) > 0){
				secondaryLD <- .computeLD(hits[-(1:(length(peaks)-1)), ], genotype, current,
						maxP=NULL, minR=minR)
				ldTable <- rbind(ldTable, secondaryLD$proxies)
				## mark all SNPs in high LD with new peak
				hits$explained[as.character(hits$snps) %in% subset(secondaryLD$proxies, 
								snp1==peaks[length(peaks)] & Rsquared >= minR)$snp2] <- TRUE
				rm(secondaryLD)
			}
		}
		## remove SNPs that never reach significance
		hits <- hits[hits$minPvalue <= pvalue,]
		
		## for each peak, get list of proxy SNPs
		if(length(peaks) > 1){
			ldTable <- subset(ldTable, snp1 %in% peaks & snp2 %in% hits$snps[-(1:length(peaks))])
			if(nrow(ldTable)){
				assignedPeak <- by(ldTable, ldTable$snp2, function(x) {
							if(nrow(x)){
								i <- which.max(x$Rsquared) 
								list(peak=x[i, "snp1"], Rsquared=x[i, "Rsquared"])
							} else NA
						})
				peak <- sapply(assignedPeak, '[[', "peak")
				proxies <- tapply(names(assignedPeak), peak, paste, collapse=",")
				r2 <- sapply(assignedPeak, '[[', "Rsquared")
				r2 <- tapply(r2, peak, paste, collapse=",")
				idx <- match(names(proxies), as.character(hits$snps))
				hits$others[idx][is.na(hits$others[idx])] <- proxies[is.na(hits$others[idx])]
				hits$others[idx][!is.na(hits$others[idx])] <- 
						paste(hits$others[idx][!is.na(hits$others[idx])], 
								proxies[!is.na(hits$others[idx])], sep=",")
				hits$Rsquared[idx][is.na(hits$Rsquared[idx])] <- r2[is.na(hits$Rsquared[idx])]
				hits$Rsquared[idx][!is.na(hits$Rsquared[idx])] <- 
						paste(hits$Rsquared[idx][!is.na(hits$Rsquared[idx])], 
								r2[!is.na(hits$Rsquared[idx])], sep=",")
			}
		}
	}
	hits <- subset(hits, !is.na(var.explained))
	hits <- hits[, -ncol(hits)]

	if(verbose) message(nrow(hits), " peaks found")
	hits
}
