## Functions to add effect size estimates to output
## Author: Peter Humburg, 2013
#############################################################

#' Compute effect size for all SNPs with significant associations.
#' @param hits A \code{data.frame} with eQTL results as produced by \code{Matrix-eQTL}
#' @param expression Gene expression values.
#' @param genotype Genotypes.
#' @param covariate Additional covariates to include in the model.
#' @param minFDR FDR threshold used to determine whether an association is significant.
#' @param geneID Restrict computations to the indicated gene.
#' @param exprOpt Options for reading of gene expression file. 
#' @param genoOpt Options for reading of gene expression file.
#' @param covOpt Options for reading of gene expression file.
#' @return A \code{data.frame} similar to \code{hits} with two additional columns. Column
#' \code{var.explained} gives the variance explained by the linear model that was fitted for 
#' the SNP and column \code{effect.size} gives the variance explained by the genotype within
#' this model (as measured by the partial r^2). 
#' 
#' @author Peter Humburg
#' @export
getEffectSize <- function(hits, expression, genotype, covariate, minFDR, geneID,
		exprOpt, genoOpt, covOpt, ...){
	## only compute effect size for SNPs that are deemed significant
	hits <- subset(hits, FDR <= minFDR)
	if(!missing(geneID) && !is.null(geneID)) hits <- subset(hits, gene == geneID)
	complete <- cbind(subset(hits, FALSE), var.explained=numeric(), effect.size=numeric())
	
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
		
		complete <- Rsge::sge.parLapply(unique(as.character(hits$gene)), 
				.submitEffectSize, hits, gene, snps, cvrt, minFDR)
		complete <- do.call(rbind, complete)
		complete <- complete[order(complete[["FDR"]], complete[["effect.size"]]), ]
		
	}
	complete
}

.submitEffectSize <- function(current, hits, expr, genotype, covariate, minFDR, verbose=FALSE){
	hits$var.explained <- NA
	hits$effect.size <- NA
	hits <- subset(hits, gene == current & FDR <= minFDR)
	if(verbose) message("Processing gene ", current, " (", nrow(hits), " eSNPs)")
	
	if(nrow(hits)){
		## extract current gene (or probe)
		expr <- subsetRows(expr, current)[[1]]
		
		## extract all relevant SNPs
		geno <- subsetRows(genotype, hits$snps)
		
		## fit model and get summary
		for(i in 1:nrow(hits)){
			fitSummary <- fitStats(expr, t(geno[[as.character(hits$snps[i])]][[1]]), 
					as.data.frame(t(as.matrix(covariate))))
			hits[i, "var.explained"] <- fitSummary$r.squared
			hits[i, "effect.size"] <- fitSummary$pr2
		}
	}
	hits
}