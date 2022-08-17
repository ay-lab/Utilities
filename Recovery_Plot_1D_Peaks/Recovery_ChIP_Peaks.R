#!/usr/bin/env Rscript

suppressMessages(library(GenomicRanges))
library(data.table)
library(optparse)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

Overlap1D <- function(Inpdata1, Inpdata2, boundary=1, offset=0, uniqov=TRUE) {

	ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))
	if (uniqov == TRUE) {
		ov_idx_file1 <- unique(ov1[,1])
		ov_idx_file2 <- unique(ov1[,2])		
	} else {
		ov_idx_file1 <- ov1[,1]
		ov_idx_file2 <- ov1[,2]
	}
	nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1)), ov_idx_file1)
	nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2)), ov_idx_file2)

	# return the overlapping and non-overlapping set of indices
	newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2)
	return(newList)

}


#================================================
option_list = list(
	make_option(c("--RefFile"), type="character", default=NULL, help="Reference 1D peak file. Mandatory parameter."),
	make_option(c("--InpFile"), type="character", default=NULL, help="Input ChIP or HiChIP peak file. Mandatory parameter."),

	make_option(c("--headerRef"), type="logical", action="store_true", default=FALSE, help="If TRUE, reference file has header information. Default FALSE."),
	make_option(c("--headerInp"), type="logical", action="store_true", default=FALSE, help="If TRUE, input file has header information. Default FALSE."),

	make_option(c("--ascend"), type="logical", action="store_true", default=FALSE, help="If TRUE, sorting of peaks should be in ascending order of the q-value related column. By default, sorting of peaks is done in descending order since MACS2 stores -log10(q-value) which should be sorted in descending order. Default FALSE."),
  	
  	make_option(c("--QcolInp"), type="integer", action="store", default=9, help="Column number storing the q-values (FDR) of the input interaction file. Default = 9, means that the 9'th column stores the q-value."),

	make_option(c("--offset"), type="integer", action="store", default=0, help="Offset / slack (in bp). Default 0 means a pair of peaks should be in exact overlapping. If 1000, 1 Kb slack is provided."),
	
	make_option(c("--OutDir"), type="character", default=NULL, help="Output directory."),

	make_option(c("--LabelRef"), type="character", default=NULL, help="Label / category / method name of the reference peak file."),
	make_option(c("--LabelInp"), type="character", default=NULL, help="Label / category / method name of the input peak file.")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$RefFile)) {
	print_help(opt_parser)
	stop("Reference peak file is not provided - check the option --RefFile \n", call.=FALSE)
}

if (is.null(opt$InpFile)) {
	print_help(opt_parser)
	stop("Input peak file is not provided - check the option --InpFile \n", call.=FALSE)
}

if (opt$QcolInp < 0) {
	print_help(opt_parser)
	stop("q-value (FDR) column for the input peak file should be > 0 - check the option --QcolInp \n", call.=FALSE)
}

if (is.null(opt$LabelRef)) {
	LabelRef <- 'Ref'
} else {
	LabelRef <- opt$LabelRef
}

if (is.null(opt$LabelInp)) {
	LabelInp <- 'Inp'
} else {
	LabelInp <- opt$LabelInp
}

##==================
system(paste("mkdir -p", opt$OutDir))

## read reference peak file
RefData <- data.table::fread(opt$RefFile, header=opt$headerRef)
RefData <- RefData[, 1:3]

## read input peak file
InpData <- data.table::fread(opt$InpFile, header=opt$headerInp)
InpData <- InpData[, c(1:3, opt$QcolInp)]

## sort the input peaks
if (opt$ascend == TRUE) {
	InpData <- InpData[order(InpData[, 4]), ]
} else {
	InpData <- InpData[order(-InpData[, 4]), ]
}

## now compute the step-wise overlap between the input and reference peaks
bool_DF <- FALSE
stepsize <- 1000
for (startidx in seq(1, nrow(InpData), by=stepsize)) {
	endidx <- min((startidx + stepsize - 1), nrow(InpData))
	peakdata <- InpData[1:endidx, ]
	ov <- Overlap1D(peakdata, RefData, boundary=1, offset=opt$offset)
	fracval <- (length(ov$B_AND_A) * 1.0) / nrow(RefData)
	currDF <- data.frame(peakcnt=endidx, frac=fracval)
	cat(sprintf("\n Computing peak overlap - input number of peaks : %s fraction of reference peaks covered : %s ", endidx, fracval))
	if (bool_DF == FALSE) {
		finalDF <- currDF
		bool_DF <- TRUE
	} else {
		finalDF <- rbind.data.frame(finalDF, currDF)
	}
}

write.table(finalDF, paste0(opt$OutDir, '/Recovery_', LabelInp, '_', LabelRef, '_offset_', opt$offset, '.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)


