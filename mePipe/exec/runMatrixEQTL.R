#!/usr/bin/Rscript

## run Matrix eQTL analysis
## Based on example script by Andrey A. Shabalin provided as part of Matrix eQTL documentation.
## Changes were made to allow running from the command line and to generate covariates from
## PCA components.
##
## Author: Peter Humburg, 2012

library(base)
library(utils)
library(stats)
library(grDevices)
library(graphics)
library(methods)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("mePipe"))

options(warn=1)
##TODO: add option to permute columns

## set up command line options
option_list <- list(
		make_option(c("-o", "--output"), default = '',
				help = "Name of output file."),
		make_option(c("-c", "--covariate"), default = NULL,
				help = "Name of file containing covariates."),
		make_option(c("-x", "--interaction", default=NULL), 
				help="Name of file containing a single covariate that should be tested for significant interaction with genotype. This is incomatible with option '--anova'."),
		make_option(c("-t", "--pthreshold"), default = 1e-5,
				help = "P-value threshold. Only associations significant at this level will be reported. [default: %default]"),
		make_option(c("-a", "--anova"), default = FALSE, action = 'store_true',
				help = "Flag indicating whether the ANOVA rather than the default linear model should be used. [default: %default]"),
		make_option(c('-i', '--cisdist'), default = 1e6, 
				help = "Distance between SNPs and genes below which effects are considerd to be in cis. [default: %default]"),
		make_option(c('-p', '--cisthreshold'), default = 0,
				help = "P-value threshold for reported cis-associations. If this is 0 no position specific annalysis will be carried out. [default: %default]"),
		make_option(c('-n', '--snpspos'), default = NULL,
				help = "Name of file containing SNP positions. This is required if the analysis should distinguish between cis and trans associations."),
		make_option(c('-g', '--genepos'), default = NULL,
				help = "Name of file containing gene positions. This is required if the analysis should distinguish between cis and trans associations."),
		make_option(c('-e', '--error'), default = '',
				help = "Name of file containing error covariance matrix"),
		make_option(c("--cluster"), action="store_true", default=FALSE,
				help="Submit jobs to SGE cluster instead of running locally. Note that this will only work when running on a login node of the cluster."),
		make_option(c('--selectcov'), default = FALSE, action = 'store_true',
				help = "Flag indicating whether covariates should be selected from the file given by option 'pca'. [default: %default]"),
		make_option(c('--filtercov'), default = FALSE, action = 'store_true',
				help = "Remove PCA covariates that are significantly associated with genotype or other covariates."),
		make_option(c('--excludecov'), default = 1e-3,
						help = "FDR threshold below which associations between principle components and genotypes are considered significant. [default: %default]"),
		make_option(c('--pcacov'), default = "[output]_pca_covariates.txt", 
				help = "Name of file with principle components of gene expression (will be created if it doesn't exist). [default: %default]"),
		make_option(c('--filterpca'), default = "[output]_pca_covariates.filtered.txt",
				help = "Name of file with filtered PCA covariates (will be created if it doesn't exist). [default: %default]"),
		make_option(c('--filterthreshold'), default = 1e-3,
				help = "P-value theshold up to which associations should be reported in the PCA association analysis. [default: %default]"),
		make_option(c('--filterout'), default = "[output]_pca_covariates.assoc",
				help = "Name of output file for PCA association analysis. [default: %default]"),
		make_option(c('--mincov'), default = 5L, 
				help = "Minimum number of PCA covariates to include. [default: %default]"),
		make_option(c('--maxcov'), default = 50L,
				help = "Maximum number of PCA covariates to include. [default: %default]"),
		make_option(c("--stepcov"), default = 5L,
				help = "Step size to use when evaluating different numbers of PCA covariates"),
		make_option(c("--ciscov"), default=0L,
				help="The number of PCA covariates to use for cis-analysis. Ignored if '--selectcov' is used. [default: %default]"),
		make_option(c("--transcov"), default=0L,
				help="The number of PCA covariates to use for trans-analysis. Ignored if '--selectcov' is used. [default: %default]"),
		make_option(c("--allcov"), default=0L,
				help="The number of PCA covariates to use for non-position-aware analysis. Ignored if '--selectcov' is used. [default: %default]"),
		make_option(c("--covthreshold"), default = 1e-3,
				help = "FDR threshold below which associations should be considered significant when choosing the optimal number of covariates. [default: %default]"),
		make_option(c("--covout"), default = "[output]_covSelect",
				help = "Name of output directory for covariate selection. This is interpreted relative to the output directory. [default: %default]"),
		make_option(c("--sortedSNPs"), action="store_true", default=FALSE,
				help="Flag indicationg whether the SNPs in 'snpspos' are sorted by genomic coordinate. [default: %default]"),
		make_option(c("--multiPeak"), action="store_true", default=FALSE,
				help="Flag indicating whether an attempt should be made to resolve multiple independent eSNPs for the same gene via multiple regression. [default: %default]"),
		make_option(c("--multiPvalue"), default=1e-6,
				help="P-value threshold to use for associations between secondary eSNPs and gene expression measurements. [default: %default]"),
		make_option(c("--ldBlocks"), action="store_true", default=FALSE,
				help="Flag indicating whether eQTLs should be summarised by LD block. [default: %default]"),
		make_option(c("--ldPairs"), action="store_true", default=FALSE,
				help="Flag indicating whether eQTLs should be grouped based on pairwise measures of LD. This differs from `ldBlocks` in that all SNPs that have R^2 > `ldR2` with the peak SNP will be grouped together and no attempt is made to infer local LD patterns. [default: %default]"),
		make_option(c("--maxDist"), default=1e5L,
				help="Maximum distance allowed between two adjacent SNPs within an LD block. [default: %default]"),
		make_option(c("--maxSNPs"), default=200L,
				help="Maximum number of SNPs to consider for each LD block. [default: %default]"),
		make_option(c("--ldFDR"), default=0.05,
				help="Maximum FDR of eQTLs to be included in list of SNPs for each block. Only blocks with at least one SNP significant at this level will be reported. [default: %default]"),
		make_option(c("--ldPvalR2"),
				help="Maximum p-value for correlation between two eSNPs for that should be considered part of the same signal by `ldPairs`. Use this option to test for significant correlation between SNPs instead of using the fixed R^2 threshold. Note that this may group SNPs with relatively low LD together, especially if the samplesize is large."),
		make_option(c("--ldR2"), default=0.8,
				help="Cut-off for R^2 between two eSNPs. All SNPs that have a higher R^2 with a peak SNP will be considered part of the same signal. [default: %default]"),
		make_option(c("--ldOnly"), action="store_true", default=FALSE,
				help="Compute LD blocks for existing eQTL results. This assumes that previous results can be loaded from the file implied by '--output'. Implies '--ldBlocks'"),
		make_option(c("-d", "--delim"), default = '\t',
				help = "Delimiter used in input files. [default: %default]"),
		make_option(c("-m", "--missing"), default = 'NA',
				help = "Missing value indicator. [default: %default]"),
		make_option(c("-r", "--rowskip"), default = 1L,
				help = "Number of rows to skip (column labels). [default: %default]"),
		make_option(c("-l", "--colskip"), default = 1L,
				help = "Number of columns (row labels). [default: %default]"),
		make_option(c("-s", "--slice"), default = 2000L,
				help = "Slice size. Files will be read in slices of this many rows. [default: %default]"),
		make_option(c("-b", "--bins"), default = 0,
				help = "Number of bins to use for p-value histogram. Set to 0 to disable histogram generation. [default: %default]"),
		make_option(c("-q", "--qqplot"), action = 'store_true', default = FALSE,
				help = "Flag indicating that a Q-Q plot rather than a histogram should be produced. Use the '--bins' option to control the resolution of the plot."),
		make_option(c("-v", "--verbose"), default = FALSE, action = 'store_true',
				help = "Flag indicating whether additional progress messages should be printed. [default: %default]")
)
parser <- OptionParser(usage = "%prog [options] expression_file genotype_file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

if(opt$ldOnly && !opt$ldBlocks) opt$ldPairs <- TRUE

## check that positional arguments are present
if(length(arguments$args) < 2){
	message("\nGene expression and genotype files are required\n")
	print_help(parser)
	stop()
}
if(length(arguments$args) > 2){
	message("\nFound extra arguments on command line:\n  ", paste(arguments$args, collapse=" "))
	print_help(parser)
	stop()
}

## check that mandatory options are present
if(opt$output == ''){
	message("\nOutput file name is required\n")
	print_help(parser)
	stop()
}

## set model type
if(opt$anova){
	opt$model <- 'anova'
} else{
	opt$model <- 'linear'
}

## check for incompatible options
if(!is.null(opt$interaction)){
	if(opt$anova){
		warning("Ignoring option '--anova' because interaction analysis was requested.")
	}
	opt$model <- "cross"
}


## if we are going to produce a Q-Q plot we need to adjust the bins
if(opt$qqplot && opt$bins > 0){
	opt$bins <- c(0.1^seq(0, -log10(opt$pthreshold)+1, length.out = opt$bins), 0)
}
## warn if qqplot was requested but bins is 0
if(opt$qqplot && opt$bins == 0){
	warning("Q-Q plot was requested but 'bins' is 0. No plots will be produced.")
}

opt$cisoutput <- paste(opt$output, 'cis', sep = '.')
## check whether we are running a position aware analysis and print summary of options
if(opt$cisdist == 0 && opt$cisthreshold > 0){
	opt$cisthreshold <- 0
}
doCis <- FALSE
doTrans <- FALSE
if(opt$cisthreshold > 0){
	doCis <- TRUE
}
if(opt$pthreshold > 0 && doCis){
	doTrans <- TRUE
}
doAll <- !doCis

if(!doAll && !opt$ldOnly){
	 if(opt$snpspos == ''){
		stop("SNP positions are required for position aware analysis.")
	}
	if(opt$genepos == ''){
		stop("Gene positions are required for position aware analysis.")
	} 
	if(doTrans && doCis){
		message("Running position aware analysis for cis- and trans-associations with\n cis distance: ", 
				opt$cisdist, "\n cis p-value: ", opt$cisthreshold, "\n trans p-value: ", 
				opt$pthreshold)
	} else if(doCis){
		message("Running position aware analysis for cis associations only with\n cis distance: ", 
				opt$cisdist, "\n cis p-value: ", opt$cisthreshold)
	} 
	message(" SNP positions: ", opt$snpspos, "\n Gene positions: ", opt$genepos)
}
if(doAll && !opt$ldOnly){
	if(opt$pthreshold == 0){
		message("\nAt least one of 'pthreshold' and 'cisthreshold' has to be greater than 0\n")
		print_help(parser)
		stop()
	}
	message("Running analysis without positional information with\n p-value: ", opt$pthreshold)
	opt$cisdist <- 0
}
if(!opt$ldOnly || opt$multiPeak){
	if(opt$anova){
		message(" Model: ANOVA")
	} else if(opt$model == "cross"){
		message(" Model: interaction")
	} else{
		message(" Model: linear")
	}
}
if(opt$ldBlocks){
	if(is.null(opt$snpspos)){
		stop("SNP positions are required for LD computations.")
	}
	message(" eQTLs will be summarised by LD block")
	message("    maximum number of SNPs per block: ", opt$maxSNPs)
	message("    maximum gap between SNPs within block: ", opt$maxDist)
}
if(opt$ldPairs){
	message(" eQTLs will be grouped by pairwise LD with peak SNP")
	message("    FDR threshold for reported associations: ", opt$ldFDR)
	message("    p-value threshold for correlations between SNPs: ", opt$ldPvalR2)
}
if(opt$multiPeak){
	message(" Multiple independent eQTLs will be resolved via multiple regression")
	message("    P-value threshold for secondary peaks: ", opt$multiPvalue)
	message("    R^2 threshold above which SNPs are associated with peak eSNP without further testing: ", opt$ldR2)
}

message("\nInput files:")
if(opt$ldOnly){
	message(" eQTL results: ", paste(opt$output, "rdata", sep = '.'))
}else{
	message(" Gene expression: ", arguments$args[1])
	message(" Genotypes: ", arguments$args[2])
	
	if(length(opt$covariate) > 0 && opt$covariate[1] != ''){
		message(" Covariates: ", opt$covariate)
	}
	if(opt$error != ''){
		message(" Covariance matrix: ", opt$error)
	}
}

message("\nOutput files:")
if(doCis){
	if(!opt$ldOnly) message(" Cis associations: ", opt$cisoutput)
	if(opt$ldBlocks){
		message(" LD blocks for cis associations: ", paste(opt$cisoutput, "LD", sep="_"))
	}
	if(opt$ldPairs){
		message(" SNP groups based on pairwise LD for cis associations: ", paste(opt$cisoutput, 
						"LDpair", sep="_"))
		message(" R^2 for all tested SNP pairs: ", paste(opt$cisoutput, "LDtable", sep="_"))
	}
	if(opt$multiPeak){
		message(" Peak eSNPs identified by multiple regression for cis associations: ",
				paste(opt$cisoutput, "peaks", sep="_"))
	}
	if(doTrans){
		if(!opt$ldOnly) message(" Trans associations: ", opt$output)
		if(opt$ldBlocks){
			message(" LD blocks for trans associations: ", paste(opt$output, "LD", sep="_"))
		}
		if(opt$ldBlocks){
			message(" SNP groups based on pairwise for trans associations: ", paste(opt$output, 
							"LDpair", sep="_"))
		}
		if(opt$ldPairs){
			message(" SNP groups based on pairwise LD for trans associations: ", paste(opt$output, 
							"LDpair", sep="_"))
			message(" R^2 for all tested SNP pairs: ", paste(opt$output, "LDtable", sep="_"))
		}
		if(opt$multiPeak){
			message(" Peak eSNPs identified by multiple regression for trans associations: ",
					paste(opt$output, "peaks", sep="_"))
		}
	}
} else{
	if(!opt$ldOnly) message(" All associations: ", opt$output)
	if(opt$ldBlocks){
		message(" LD blocks for all associations: ", paste(opt$output, "LD", sep="_"))
	}
	if(opt$ldPair){
		message(" SNP groups based on pairwise LD for all associations: ", 
				paste(opt$output, "LDpair", sep="_"))
		message(" R^2 for all tested SNP pairs: ", paste(opt$output, "LDtable", sep="_"))
	}
	if(opt$multiPeak){
		message(" Peak eSNPs identified by multiple regression for trans associations: ",
				paste(opt$output, "peaks", sep="_"))
	}
}
if(!opt$ldOnly) message(" R objects: ", paste(opt$output, "rdata", sep = '.'))
message("")

## set-up environment for Rsge
sge.options(sge.use.cluster=opt$cluster, sge.user.options="-S /bin/bash -V", sge.trace=opt$verbose)

## initialise number of principle components to use
selected <- list(selected=list(all=opt$allcov, cis=opt$ciscov, trans=opt$transcov))

if(!opt$ldOnly){
	if(opt$selectcov){
		## add output prefix to file names (unless custom names were provided)
		opt$pcacov <- gsub("[output]", opt$output, opt$pcacov, fixed=TRUE)
		opt$filterpca <- gsub("[output]", opt$output, opt$filterpca, fixed=TRUE)
		opt$covout <- gsub("[output]", opt$output, opt$covout, fixed=TRUE)
		opt$filterout <- gsub("[output]", opt$output, opt$filterout, fixed=TRUE)
		
		message("Selecting optimal number of PCA covariates ...")
		if(!file.exists(opt$pcacov) && (!opt$filtercov || !file.exists(opt$filterpca)) ){
			## create PCA covariates
			message("  Running PCA ...")
			allCovariates(arguments$args[1], 
					getOptions(sep = opt$delim, missing = opt$missing, rowskip = opt$rowskip, 
							colskip = opt$colskip, slice = opt$slice), 
					opt$pcacov, opt$delim)
		}
		else{
			message("  Using PCA covariates in ", opt$pcacov)
		}
		if(opt$filtercov){
			if(!file.exists(opt$filterpca)){
				## remove PCA covariates that are associated with genotype or other covariates (if testing for interaction)
				message("  Computing associations between covariates and genotypes ...")
				pcaAssoc <- covAssoc(arguments$args[2], opt$pcacov, c(opt$covariate, opt$interaction),
						opt$filterout, opt$filterpca,
						getOptions(sep = opt$delim, missing = opt$missing, rowskip = opt$rowskip, 
								colskip = opt$colskip, slice = opt$slice),
						getOptions(sep = opt$delim, missing = opt$missing, rowskip = opt$rowskip, 
								colskip = opt$colskip, slice = opt$slice),
						opt$model, opt$filterthreshold, opt$excludecov)
				message("  Excluding ", length(pcaAssoc), " PCs associated with genotype (FDR < ", 
						opt$excludecov, ")")
				message("  Excluded PCs: ", paste(pcaAssoc, collapse=", "))
			}
			else{
				message("  Using filtered PCA covariates in ", opt$filterpca)
			}
		}
		## determine optimal number of covariates
		covariate <- character()
		if(opt$filtercov){
			covariate <- opt$filterpca
		}
		else{
			covariate <- opt$pcacov
		}
		if(!is.null(opt$covariate)){
			covariate <- c(covariate, opt$covariate)
		}
		if(!is.null(opt$interaction)){
			covariate <- c(covariate, opt$interaction)
		}
		
		candidates <- seq(opt$mincov, opt$maxcov, by = opt$stepcov)
		fileOpt <- getOptions(sep = opt$delim, 
				missing = opt$missing, rowskip = opt$rowskip, 
				colskip = opt$colskip, slice = opt$slice)
		selected <- chooseCov(expression=arguments$args[1], genotype=arguments$args[2], 
				covariate=covariate, candidates=candidates,
				covThreshold=opt$covthreshold, covOpt=fileOpt,
				exprOpt = fileOpt, genoOpt = fileOpt,
				output=opt$covout, threshold = opt$pthreshold,
				model = opt$model, 
				error = opt$error, 
				verbose = opt$verbose,
				cisThreshold = opt$cisthreshold,
				snpsPos = opt$snpspos,
				genePos = opt$genepos,
				cis = opt$cisdist,
				bins = opt$bins, qqplot = opt$qqplot, 
				doCis = opt$cisthreshold > 0, 
				doTrans = opt$cisthreshold > 0 && opt$pthreshold > 0)
		
		for(i in 1:length(selected$selected)){
			if(names(selected$selected)[i] == "all"){
				message("Using ", selected$selected[i], " covariates")
			} else{
				message("Using ", selected$selected[i], " covariates for ", 
						names(selected$selected)[i], " analysis")
			}
		}
		
		## create plot of number of covariates vs number of eQTLs
		pdf(file=file.path(opt$covout, "CovVsEQTL.pdf")) 
		plotCov(selected$covariates, selected$eqtls$total[names(selected$selected)], 
				selected$selected, main="All associations returned by Matrix-eQTL")
		plotCov(selected$covariates, selected$eqtls$significant[names(selected$selected)], 
				selected$selected, main=paste0("Significant associations (fdr < ", 
						opt$covthreshold, ")"))
		dev.off()
		
		## copy selected results to output
		results <- lapply(selected$best, function(x) 
					dir(opt$covout, paste0("^", sub("\\.cis$", "", basename(x)), "(\\.rdata)?$")))
		if("cis" %in% names(results)){
			results[["cis"]][1] <- paste(results[["cis"]][1], "cis", sep=".")
		}
		target <- mapply(function(x, y) 
					gsub(basename(x), 
							paste0(basename(opt$output)), y), sub("\\.cis$", "", selected$best), 
				results)
		
		status <- file.copy(sapply(results, function(x) file.path(opt$covout, x)), 
				sapply(target, function(x) file.path(dirname(opt$output), x)),
				overwrite = TRUE)
		
		## combine MatrixEQTL result objects for cis- and trans-analysis
		if(doCis && doTrans){
			load(file.path(opt$covout, results$cis[2]))
			me.cis <- me
			load(file.path(opt$covout, results$trans[2]))
			me$param$output_file_name.cis <- me.cis$param$output_file_name.cis
			me$cis <- me.cis$cis
		} else{
			if(doAll) load(file.path(opt$covout, results$all[2]))
			if(doCis) load(file.path(opt$covout, results$cis[2]))
			if(doTrans) load(file.path(opt$covout, results$trans[2]))
		}
		me$all$covariates <- selected$selected["all"]
		me$cis$covariates <- selected$selected["cis"]
		me$trans$covariates <- selected$selected["trans"]
		save(me, file=paste(opt$output, "rdata", sep="."))
		
	} else{
		## or just run a single analysis with pre-specified covariates
		covariates <- c(opt$covariate, opt$interaction)
		me <- runME(arguments$args[1], arguments$args[2], covariates, 
				opt$error, opt$snpspos,	opt$genepos, opt$output, opt$cisoutput, 
				getOptions(sep = opt$delim, missing = opt$missing, 
						rowskip = opt$rowskip, colskip = opt$colskip, slice = opt$slice),
				getOptions(sep = opt$delim, missing = opt$missing, 
						rowskip = opt$rowskip, colskip = opt$colskip, slice = opt$slice),
				getOptions(sep = opt$delim, missing = opt$missing, 
						rowskip = opt$rowskip, colskip = opt$colskip, slice = opt$slice),
				opt$pthreshold, opt$cisthreshold, opt$model, opt$cisdist, opt$bins, opt$verbose,
				opt$qqplot)
		me$all$covariates <- selected$selected[["all"]]
		me$cis$covariates <- selected$selected[["cis"]]
		me$trans$covariates <- selected$selected[["trans"]]
		save(me, file=paste(opt$output, "rdata", sep="."))
	}
} else{
	load(paste(opt$output, 'rdata', sep='.'))
	if(doAll && (is.null(me$all) || is.null(me$all$eqtls) || nrow(me$all$eqtls) == 0)){
		message("No associations found. Did you mean to run a position aware analysis?")
	}
	if(doCis && (is.null(me$cis) || is.null(me$cis$eqtls) || nrow(me$cis$eqtls) == 0)){
		message("No cis-associations found. Maybe you didn't run a position aware analysis?")				
	}
	if(doTrans && (is.null(me$trans) || is.null(me$trans$eqtls) || nrow(me$trans$eqtls) == 0)){
		message("No trans associations found. Maybe you didn't run a position aware analysis?")
	}
}
if(opt$ldPairs){
	if(doAll){
		message("Computing pairwise LD ...")
		ans <- getLDpairs(me$all$eqtls, arguments$args[2], minFDR=opt$ldFDR, maxP=opt$ldPvalR2,
				minR=opt$ldR2, genoOpt=getOptions(sep = opt$delim, missing = opt$missing, 
						rowskip = opt$rowskip, colskip = opt$colskip, slice = opt$slice))
		write.table(ans$groups, file=paste(opt$output, "LDpair", sep="_"), 
				quote=FALSE, row.names=FALSE, sep="\t")
		write.table(ans$proxies, file=paste(opt$output, "LDtable", sep="_"), 
				quote=FALSE, row.names=FALSE, sep="\t")
	} else{
		if(doCis){
			message("Computing pairwise LD for cis associations...")
			ans <- getLDpairs(me$cis$eqtls, arguments$args[2], minFDR=opt$ldFDR, maxP=opt$ldPvalR2,
					minR=opt$ldR2, genoOpt=getOptions(sep = opt$delim, missing = opt$missing, 
							rowskip = opt$rowskip, colskip = opt$colskip, slice = opt$slice))
			write.table(ans$groups, file=paste(opt$cisoutput, "LDpair", sep="_"), 
					quote=FALSE, row.names=FALSE, sep="\t")
			write.table(ans$proxies, file=paste(opt$cisoutput, "LDtable", sep="_"), 
					quote=FALSE, row.names=FALSE, sep="\t")
		}
		if(doTrans){
			message("Computing pairwise LD for trans associations...")
			ans <- getLDpairs(me$trans$eqtls, arguments$args[2], minFDR=opt$ldFDR, maxP=opt$ldPvalR2,
					minR=opt$ldR2, genoOpt=getOptions(sep = opt$delim, missing = opt$missing, 
							rowskip = opt$rowskip, colskip = opt$colskip, slice = opt$slice))
			write.table(ans$groups, file=paste(opt$output, "LDpair", sep="_"), 
					quote=FALSE, row.names=FALSE, sep="\t")
			write.table(ans$proxies, file=paste(opt$output, "LDtable", sep="_"), 
					quote=FALSE, row.names=FALSE, sep="\t")
		}
	}
}
if(opt$ldBlocks){
	pos <- read.table(opt$snpspos, header=TRUE, stringsAsFactors=FALSE, row.names=1)
	names(pos) <- c("chrom", "pos")
	
	## SNPs in 'pos' have to be sorted by genomic position
	if(!opt$sortedSNPs) pos <- pos[order(pos$chrom, pos$pos),]
	
	if(doAll){
		message("Computing LD structure...")
		blocks <- getLDblocks(me$all$eqtls, arguments$args[2], pos, dist=opt$maxDist, 
				window=opt$maxSNPs, verbose=opt$verbose, minFDR=opt$ldFDR,
				genoOpt=getOptions(sep = opt$delim, missing = opt$missing, 
						rowskip = opt$rowskip, colskip = opt$colskip, slice = opt$slice))
		write.table(blocks, file=paste(opt$output, "LD", sep="_"), 
				quote=FALSE, row.names=FALSE, sep="\t")
	} else{
		if(doCis){
			message("Computing LD structure for cis associations...")
			blocks <- getLDblocks(me$cis$eqtls, arguments$args[2], pos, dist=opt$maxDist, 
					window=opt$maxSNPs, verbose=opt$verbose, minFDR=opt$ldFDR,
					genoOpt=getOptions(sep = opt$delim, missing = opt$missing, 
							rowskip = opt$rowskip, colskip = opt$colskip, slice = opt$slice))
			write.table(blocks, file=paste(opt$cisoutput, "LD", sep="_"), 
					quote=FALSE, row.names=FALSE, sep="\t")
		}
		if(doTrans){
			message("Computing LD structure for trans associations...")
			blocks <- getLDblocks(me$trans$eqtls, arguments$args[2], pos, dist=opt$maxDist,
					window=opt$maxSNPs, verbose=opt$verbose, minFDR=opt$ldFDR,
					genoOpt=getOptions(sep = opt$delim, missing = opt$missing, 
							rowskip = opt$rowskip, colskip = opt$colskip, slice = opt$slice))
			write.table(blocks, file=paste(opt$output, "LD", sep="_"), 
					quote=FALSE, row.names=FALSE, sep="\t")
		}
	}
	
}
if(opt$multiPeak){
	covariates <- list(opt$covariate, opt$interaction)
	fileOptions <- getOptions(sep = opt$delim, missing = opt$missing, 
			rowskip = opt$rowskip, colskip = opt$colskip, slice = opt$slice)
	selected <- lapply(me[c("all", "cis", "trans")], '[[', "covariate")
	
	if(doAll){
		message("Resolving multiple peaks...")
		allCovariates <- covariates
		if(me$all$covariates > 0){
			pc <- loadData(if(opt$filtercov) opt$filterout else opt$pcacov, 
					getOptions(sep = opt$delim, missing = opt$missing, 
							rowskip = opt$rowskip, colskip = opt$colskip, 
							slice = me$all$covariates))
			pc <- pc$CreateFromMatrix(pc[[1]])
			allCovariates <- c(pc, covariates)
		}
		multi <- getMultiPeak(me$all$eqtls, opt$multiPvalue, 
				arguments$args[1], arguments$arg[2], allCovariates,
				opt$ldR2, fileOptions, fileOptions, fileOptions)
		write.table(multi, file=paste(opt$output, "peaks", sep="_"), row.names=FALSE,
				quote=FALSE, sep="\t")
	} else{
		if(doCis){
			message("Resolving multiple peaks for cis associations...")
			cisCovariates <- covariates
			if(me$cis$covariates > 0){
				pc <- loadData(if(opt$filtercov) opt$filterout else opt$pcacov, 
						getOptions(sep = opt$delim, missing = opt$missing, 
								rowskip = opt$rowskip, colskip = opt$colskip, 
								slice = me$cis$covariates))
				pc <- pc$CreateFromMatrix(pc[[1]])
				cisCovariates <- c(pc, covariates)
			}
			multi <- getMultiPeak(me$cis$eqtls, opt$multiPvalue, 
					arguments$args[1], arguments$arg[2], cisCovariates,
					opt$ldR2, fileOptions, fileOptions, fileOptions)
			write.table(multi, file=paste(opt$cisoutput, "peaks", sep="_"), row.names=FALSE,
					quote=FALSE, sep="\t")
		}
		if(doTrans){
			message("Resolving multiple peaks for trans-associations...")
			transCovariates <- covariates
			if(me$trans$covariates > 0){
				pc <- loadData(if(opt$filtercov) opt$filterout else opt$pcacov, 
						getOptions(sep = opt$delim, missing = opt$missing, 
								rowskip = opt$rowskip, colskip = opt$colskip, 
								slice = me$trans$covariates))
				pc <- pc$CreateFromMatrix(pc[[1]])
				transCovariates <- c(pc, covariates)
			}
			multi$trans <- getMultiPeak(me$trans$eqtls, opt$multiPvalue, 
					arguments$args[1], arguments$arg[2], transCovariates,
					opt$ldR2, fileOptions, fileOptions, fileOptions)
			write.table(multi, file=paste(opt$output, "peaks", sep="_"), row.names=FALSE,
					quote=FALSE, sep="\t")
		}
	}
}
message("done.")