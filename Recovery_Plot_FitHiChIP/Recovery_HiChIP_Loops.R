#!/usr/bin/env Rscript

suppressMessages(library(GenomicRanges))
library(data.table)
library(optparse)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

ReadContactData <- function(inpfile, headerInp=TRUE, chrNUM=FALSE, chrSET=c(1,4), mid=FALSE, midCOL=c(2), binSize=5000) {

	outDF <- data.table::fread(inpfile, header=headerInp)
    
    # handle chromosome intervals midpoints 
    if ((mid==TRUE) & (length(midCOL) > 0) & (binSize > 0)) {
        n <- ncol(outDF)
        for (i in (1:n)) {
            if (i %in% midCOL) {
                s1 <- outDF[,i] - (binSize/2)
                e1 <- outDF[,i] + (binSize/2)
                if (i == 1) {
                    ModDF <- cbind.data.frame(s1, e1)
                } else {
                    ModDF <- cbind.data.frame(ModDF, s1, e1)
                }
            } else {
                if (i == 1) {
                    ModDF <- outDF[,i]
                } else {
                    ModDF <- cbind.data.frame(ModDF, outDF[,i])
                }
            }
        }
        # reassign the modified data frame
        outDF <- ModDF
    }

    # handle chromosome numbers instead of name
    if (chrNUM==TRUE) {
        for (i in (1:length(chrSET))) {
            colno <- chrSET[i]      
            outDF[,colno] <- paste0('chr', as.character(outDF[,colno])) 
        }
    }
        
    return (outDF)
}


#=====================
OverlapLoop <- function(Inpdata1, Inpdata2, boundary=1, offset=0, uniqov=TRUE, IDX=FALSE) {

	# first compute the overlap between Inpdata1[,1:3],  Inpdata2[,1:3]
	# and between Inpdata1[,4:6],  Inpdata2[,4:6]
	ov1 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))
	ov2 <- as.data.frame(findOverlaps(GRanges(Inpdata1[,4], IRanges(Inpdata1[,5]+boundary-offset, Inpdata1[,6]-boundary+offset)),GRanges(Inpdata2[,4], IRanges(Inpdata2[,5]+boundary-offset, Inpdata2[,6]-boundary+offset))))
	overlap_uniq_mat <- ov1[unique(which(paste(ov1[,1], ov1[,2], sep=".") %in% paste(ov2[,1], ov2[,2], sep="."))),]

	# then compute the overlap between Inpdata1[,4:6],  Inpdata2[,1:3]
	# and between Inpdata1[,1:3],  Inpdata2[,4:6]
	# this is because input data can have first interval coordinate > second interval coordinate
	ov1A <- as.data.frame(findOverlaps(GRanges(Inpdata1[,4], IRanges(Inpdata1[,5]+boundary-offset, Inpdata1[,6]-boundary+offset)),GRanges(Inpdata2[,1], IRanges(Inpdata2[,2]+boundary-offset, Inpdata2[,3]-boundary+offset))))
	ov2A <- as.data.frame(findOverlaps(GRanges(Inpdata1[,1], IRanges(Inpdata1[,2]+boundary-offset, Inpdata1[,3]-boundary+offset)),GRanges(Inpdata2[,4], IRanges(Inpdata2[,5]+boundary-offset, Inpdata2[,6]-boundary+offset))))
	overlap_uniq_mat_1 <- ov1A[unique(which(paste(ov1A[,1], ov1A[,2], sep=".") %in% paste(ov2A[,1], ov2A[,2], sep="."))),]

	# now concatenate these two outcomes to get the final overlap statistics
	if (uniqov == TRUE) {
		ov_idx_file1 <- unique(c(overlap_uniq_mat[,1], overlap_uniq_mat_1[,1]))
		ov_idx_file2 <- unique(overlap_uniq_mat[,2], overlap_uniq_mat_1[,2])		
	} else {
		ov_idx_file1 <- c(overlap_uniq_mat[,1], overlap_uniq_mat_1[,1])
		ov_idx_file2 <- c(overlap_uniq_mat[,2], overlap_uniq_mat_1[,2])
	}
	nonov_idx_file1 <- setdiff(seq(1, nrow(Inpdata1)), ov_idx_file1)
	nonov_idx_file2 <- setdiff(seq(1, nrow(Inpdata2)), ov_idx_file2)

	# return the overlapping and non-overlapping set of indices
	if (IDX == FALSE) {
		newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2, A_AND_B.df = Inpdata1[ov_idx_file1, ], B_AND_A.df = Inpdata2[ov_idx_file2, ], A_MINUS_B.df = Inpdata1[nonov_idx_file1, ], B_MINUS_A.df = Inpdata2[nonov_idx_file2, ])
	} else {
		newList <- list(A_AND_B = ov_idx_file1, B_AND_A = ov_idx_file2, A_MINUS_B = nonov_idx_file1, B_MINUS_A = nonov_idx_file2)
	}
	return(newList)
}


#================================================
option_list = list(
	make_option(c("--RefFile"), type="character", default=NULL, help="Reference loop file (significant interactions). Mandatory parameter."),
	make_option(c("--InpFile"), type="character", default=NULL, help="Input (FitHiChIP or other methods) HiChIP significant loop file. Mandatory parameter."),

	make_option(c("--headerRef"), type="logical", action="store_true", default=FALSE, help="If TRUE, reference file has header information. Default FALSE."),
	make_option(c("--headerInp"), type="logical", action="store_true", default=FALSE, help="If TRUE, input file has header information. Default FALSE."),
  	
  	make_option(c("--QcolInp"), type="integer", action="store", default=0, help="Column number storing the q-values (FDR) of the input interaction file. Default = 0, means that the last column stores the q-value. Otherwise, a positive value for the target column number needs to be specified."),

  	make_option(c("--chrRef"), type="logical", action="store_true", default=FALSE, help="If TRUE, reference loop file has numbers as the chromosomes, for example 1, 2, instead of chr1, chr2. Default FALSE."),
  	make_option(c("--chrInp"), type="logical", action="store_true", default=FALSE, help="If TRUE, input loop file has numbers as the chromosomes, for example 1, 2, instead of chr1, chr2. Default FALSE."),

  	make_option(c("--midRef"), type="logical", action="store_true", default=FALSE, help="If TRUE, reference interaction file has only the midpoints of an  interval (chr1, mid1, chr2, mid2). Otherwise (default = FALSE), the reference file columns are (chr1, start1, end1, chr2, start2, end2)."),
  	make_option(c("--midInp"), type="logical", action="store_true", default=FALSE, help="If TRUE, input interaction file has only the midpoints of an  interval (chr1, mid1, chr2, mid2). Otherwise (default = FALSE), the input file columns are (chr1, start1, end1, chr2, start2, end2)."),

  	make_option(c("--binsizeRef"), type="integer", action="store", default=5000, help="Bin size (in bp) for the reference interaction file. Required if --midRef is 1. Default = 5000 (5 Kb)"), 
  	make_option(c("--binsizeInp"), type="integer", action="store", default=5000, help="Bin size (in bp) for the input interaction file. Required if --midInp is 1. Default = 5000 (5 Kb)"), 

	make_option(c("--offset"), type="integer", action="store", default=5000, help="Offset / slack (in bp). Default 5000 (5 Kb) means a pair of loops with interacting bins within 5 Kb for both sides, will be considered as overlapping. If 0, exact overlap will be computed."),
	
	make_option(c("--OutDir"), type="character", default=NULL, help="Output directory."),

	make_option(c("--LabelRef"), type="character", default=NULL, help="Label / category / method name of the reference loop file."),
	make_option(c("--LabelInp"), type="character", default=NULL, help="Label / category / method name of the input loop file.")

); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$RefFile)) {
	print_help(opt_parser)
	stop("Reference interaction file is not provided - check the option --RefFile \n", call.=FALSE)
}

if (is.null(opt$InpFile)) {
	print_help(opt_parser)
	stop("Input interaction file is not provided - check the option --InpFile \n", call.=FALSE)
}

if (opt$QcolInp < 0) {
	print_help(opt_parser)
	stop("q-value (FDR) column for the input interaction file should be either 0 (last column) or positive - check the option --QcolInp \n", call.=FALSE)
}

if (opt$midRef && opt$binsizeRef <= 0) {
	print_help(opt_parser)
	stop("Reference interaction file has only midpoints of the chromosome intervals, but the bin size is <= 0 - check the options --midRef and --binsizeRef \n", call.=FALSE)
}

if (opt$midInp && opt$binsizeInp <= 0) {
	print_help(opt_parser)
	stop("Input interaction file has only midpoints of the chromosome intervals, but the bin size is <= 0 - check the options --midInp and --binsizeInp \n", call.=FALSE)
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

## read reference contact file
RefData <- ReadContactData(opt$RefFile, headerInp=opt$headerRef, chrNUM=opt$chrRef, chrSET=c(1,4), mid=opt$midRef, midCOL=c(2,4), binSize=opt$binsizeRef)

## read input contact file
InpData <- ReadContactData(opt$InpFile, headerInp=opt$headerInp, chrNUM=opt$chrInp, chrSET=c(1,4), mid=opt$midInp, midCOL=c(2,4), binSize=opt$binsizeInp)

## get the q-value column information for the input file
if (opt$QcolInp == 0) {
	# the last column is the q-value column
	QcolInp <- ncol(InpData)
} else {
	# q-value column is a positive and custom value
	if (opt$midInp) {
		# here two columns are already added in the input interaction data
		# so we also shift the q-value column
		QcolInp <- opt$QcolInp + 2
	} else {
		QcolInp <- opt$QcolInp
	}
}

## now filter the reference data - only 1st 6 columns
RefData <- RefData[, 1:6]

## now filter the InpData data - only 1st 6 columns, and the q-value column
InpData <- InpData[, c(1:6, QcolInp)]

## sort the input interactions by the q-value (7th column)
InpData <- InpData[order(InpData[, 7]), ]

## now compute the step-wise overlap between the input interactions and the 
## reference loops
bool_DF <- FALSE
stepsize <- 1000
for (startidx in seq(1, nrow(InpData), by=stepsize)) {
	endidx <- min((startidx + stepsize - 1), nrow(InpData))
	loopdata <- InpData[1:endidx, ]
	ov <- OverlapLoop(loopdata, RefData, boundary=1, offset=opt$offset)
	fracval <- (length(ov$B_AND_A) * 1.0) / nrow(RefData)
	currDF <- data.frame(loopcnt=endidx, frac=fracval)
	cat(sprintf("\n Computing loop overlap - input number of loops : %s fraction of reference loops covered : %s ", endidx, fracval))
	if (bool_DF == FALSE) {
		finalDF <- currDF
		bool_DF <- TRUE
	} else {
		finalDF <- rbind.data.frame(finalDF, currDF)
	}
}

write.table(finalDF, paste0(opt$OutDir, '/Recovery_', LabelInp, '_', LabelRef, '_offset_', opt$offset, '.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)


