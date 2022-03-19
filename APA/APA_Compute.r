#!/usr/bin/env Rscript

#===========================================================
# R script for computing the Aggregate Peak Analysis (APA)
# given a set of (significant) interactions

#Author: Sourya Bhattacharyya
#Vijay-Ay lab, LJI
#===========================================================

suppressMessages(library(GenomicRanges))
library(optparse)
library(tools)
library(ggplot2)
library(reshape2)	# used for melt function
library(RColorBrewer)

# variables depicting the dimension of plots
PlotWidth <- 6	#8
PlotHeight <- 4	#6

# font size used in texts
FontSize <- 20

# Color vector used for plotting the matrices
colorvec <- colorRampPalette(brewer.pal(9,"YlGnBu"))(10)

#================================================
ExtractChrData <- function(InpFile, chrName, OutFile=NULL, header=TRUE, dist=c(-1,-1), mid=FALSE) {
	if (is.null(OutFile)) {
		OutFile <- paste0(dirname(InpFile), "/temp_Chr_data.bed")
	}

	if ((dist[1] > 0) & (dist[2] > 0) & (dist[2] > dist[1])) {
		distthrlow <- dist[1]
		distthrhigh <- dist[2]		
	} else {
		distthrlow <- -1
		distthrhigh <- -1
	}

	if (file_ext(InpFile) == "gz") {
		if (header == TRUE) {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		} else {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))						
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("zcat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("zcat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		}
	} else {
		if (header == TRUE) {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if (NR>1 && $1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		} else {
			if (mid == TRUE) {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $3==\"", chrName, "\" && (abs($2-$4)>=dt) && (abs($2-$4)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $3==\"", chrName, "\") {print $0}}\' -  > ", OutFile))						
				}
			} else {
				if ((distthrlow > 0) & (distthrhigh > 0)) {
					system(paste0("cat ", InpFile, " | awk -v dt=", distthrlow, " -v dth=", distthrhigh, " \' function abs(v) {return v < 0 ? -v : v} {if ($1==\"", chrName, "\" && $4==\"", chrName, "\" && (abs($2-$5)>=dt) && (abs($2-$5)<=dth)) {print $0}}\' -  > ", OutFile))
				} else {
					system(paste0("cat ", InpFile, " | awk \' {if ($1==\"", chrName, "\" && $4==\"", chrName, "\") {print $0}}\' -  > ", OutFile))
				}
			}
		}
	}

}	# end function

ReadContactData <- function(inpfile, headerInp=TRUE, chrNUM=FALSE, chrSET=c(1,4), mid=FALSE, midCOL=c(2), binSize=5000) {
    
    # handle gzipped file
	if (file_ext(inpfile) == "gz") {
    	outDF <- read.table(gzfile(inpfile), header=headerInp, sep="\t", stringsAsFactors=F)
    } else {
        outDF <- read.table(inpfile, header=headerInp, sep="\t", stringsAsFactors=F)
    }

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
  	make_option(c("--InpFile"), type="character", default=NULL, help="Input interaction / Loop file"),
  	make_option(c("--headerInp"), type="logical", action="store_true", default=FALSE, help="If TRUE, input interaction file has header. Default FALSE."),
  	make_option(c("--midInp"), type="logical", action="store_true", default=FALSE, help="If TRUE, input interaction file has the bins in (chr1, mid1, chr2, mid2) instead of the conventional (chr1, start1, end1, chr2, start2, end2) format (similar to FitHiC). Default FALSE."),  
  	make_option(c("--chrInp"), type="logical", action="store_true", default=FALSE, help="If TRUE, input interaction file has numbers as the chromosomes, for example 1, 2, instead of chr1, chr2...Default FALSE."),

  	make_option(c("--RefFile"), type="character", default=NULL, help="Reference interaction / Loop file, preferably from Hi-C data"),
  	make_option(c("--headerRef"), type="logical", action="store_true", default=FALSE, help="Similar to the option --headerInp for the reference loop file. Default FALSE."),
  	make_option(c("--midRef"), type="logical", action="store_true", default=FALSE, help="Similar to the option --midInp for the reference loop file. Default FALSE."),
  	make_option(c("--chrRef"), type="logical", action="store_true", default=FALSE, help="Similar to the option --chrInp for the reference loop file. Default FALSE."),  	

  	make_option(c("--cccol"), type="integer", action="store", default=7, help="Column number storing the contact counts in the REFERENCE interaction file. Default = 7."),  	

  	make_option(c("--binsize"), type="integer", action="store", default=5000, help="Bin size. Default 5000 (5 Kb). Both the input and reference loop files should have identical bin sizes."),
  	
  	make_option(c("--window"), type="integer", action="store", default=50000, help="Window size (+/- in both directions). Default 50000: 50 Kb in either directions). Should be a multiple of the parameter --binsize. Set as recommended in the MANGO paper."),

  	make_option(c("--distthrlow"), type="integer", action="store", default=150000, help="Interactions having distance lower than this threshold are discarded from analysis. Default 150000 (150 Kb). Set as recommended in the MANGO paper."),  
  	make_option(c("--distthrhigh"), type="integer", action="store", default=1000000, help="Interactions having distance higher than this threshold are discarded from analysis. Default = 1000000 (1 Mb). Set as recommended in the MANGO paper."),

  	make_option(c("--low"), type="integer", action="store", default=15000, help="Lower distance threshold from the anchor segment, which is used to compute the APA score and the plot range. Default 15000, means that pixels with distance >= 15 Kb will be considered for APA mean score and pixel range computation. Set as recommended in the MANGO paper."),
  	
  	make_option(c("--high"), type="integer", action="store", default=30000, help="Upper distance threshold from the anchor segment, which is used to compute the APA score and the plot range. Default 30000, means that pixels with distance <= 30 Kb will be considered for APA mean score and pixel range computation. Set as recommended in the MANGO paper."),

  	make_option(c("--cclim"), type="integer", action="store", default=50, help="Limit of the contact count for color scale display. Default = 50. Can be adjusted according to the visual APA color scale."),
  	  	
  	make_option(c("--OutPrefix"), type="character", default=NULL, help="Output prefix (including the output directory) used to name the output files. Mandatory parameter."),

  	make_option(c("--overwrite"), type="logical", action="store_true", default=FALSE, help="If TRUE, overwrites existing results. Default FALSE.")  	
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);  		

if (is.null(opt$InpFile)) {
	print_help(opt_parser)
	stop("Input interaction file is not provided. Check the option --InpFile \n", call.=FALSE)
}
if (is.null(opt$RefFile)) {
	print_help(opt_parser)
	stop("Input interaction file is not provided. Check the option --RefFile \n", call.=FALSE)
}
if (is.null(opt$OutPrefix)) {
	print_help(opt_parser)
	stop("The output prefix (including the output directory). Check the option --OutPrefix \n", call.=FALSE)
}

#===============
# create the output directory
#===============
outbasedir <- dirname(opt$OutPrefix)
onlyprefix <- basename(opt$OutPrefix)
outdir <- paste0(outbasedir, '/b', opt$binsize, '_w', opt$window, '_L', opt$distthrlow, '_H', opt$distthrhigh, '_l', opt$low, '_h', opt$high)
system(paste('mkdir -p', outdir))

#===============
## get the chromosome names and numbers
#===============
temploopfile <- paste0(outdir, '/temp_loop_chr.txt')
if (opt$headerInp == FALSE) {
	system(paste0("awk \'{print $1}\' ", opt$InpFile, " | sort -k1,1 | uniq > ", temploopfile))
} else {
	system(paste0("awk \'{if (NR>1) {print $1}}\' ", opt$InpFile, " | sort -k1,1 | uniq > ", temploopfile))
}

tempLoopData <- read.table(temploopfile, header=F, sep="\t", stringsAsFactors=F)
if (opt$chrInp) {
	## chromosome numbers are provided
	ChrList_Num <- sort(as.vector(tempLoopData[,1]))
	ChrList_NameNum <- paste0("chr", ChrList_Num)
} else {
	## "chr" + numbers
	ChrList_NameNum <- sort(as.vector(tempLoopData[,1]))
	ChrList_Num <- gsub("chr", "", ChrList_NameNum)
}
system(paste("rm", temploopfile))

#===============
# output files storing the matrix outputs
#===============
MatrixFinalOutFile <- paste0(outdir, '/', onlyprefix, '_DumpAPAMatrix.bed')
Matrix_Count_FinalOutFile <- paste0(outdir, '/', onlyprefix, '_Dump_APA_Count_Matrix.bed')
NormMatrixFinalOutFile <- paste0(outdir, '/', onlyprefix, '_Dump_Norm_APAMatrix.bed')

if ((file.exists(MatrixFinalOutFile) == FALSE) | (file.exists(Matrix_Count_FinalOutFile) == FALSE) | (file.exists(NormMatrixFinalOutFile) == FALSE) | (opt$overwrite == TRUE)) {

	# text file storing the results
	textfile <- paste0(outdir, '/', onlyprefix, '_Result_Summary.txt')
	cat(sprintf("\n **** textfile : %s ", textfile))
	fp_out <- file(textfile, "w")

	# these are temporary files, which store chromosome specific (and maintaining distthrlow and distthrhigh)
	# interaction from input and reference files
	InpTempFile = paste0(outdir, '/', onlyprefix, '_Temp_Inp.bed')
	RefTempFile = paste0(outdir, '/', onlyprefix, '_Temp_Ref.bed')

	#==========================
	# first check if there is partial results already computed (upto certain chromosome)
	#==========================
	PartialRes <- FALSE
	for (chr_idx in (length(ChrList_NameNum):1)) {
		chrName <- ChrList_NameNum[chr_idx]
		MatrixOutFile <- paste0(outdir, '/', onlyprefix, '_DumpAPAMatrix_after_', chrName ,'.bed')
		Matrix_Count_OutFile <- paste0(outdir, '/', onlyprefix, '_Dump_APA_Count_Matrix_after_', chrName ,'.bed')
		if ((file.exists(MatrixOutFile) == TRUE) & (file.exists(Matrix_Count_OutFile) == TRUE) & (opt$overwrite == FALSE)) {
			# partial result upto current chromosome is available
			PartialRes <- TRUE
			break
		}
	}
	if (PartialRes == TRUE) {
		# store the chromosome index upto which partial results are computed
		Partial_Res_Chr_Idx <- chr_idx
		# read the results and store them in matrix structures
		APA_Matrix <- as.matrix(read.table(MatrixOutFile, header=F))
		Count_APA_Mat <- as.matrix(read.table(Matrix_Count_OutFile, header=F))
		# dimension of the matrix
		nelem <- ncol(APA_Matrix)
		outtext <- paste0("\n\n **** Partial result is already computed upto the chromosome index : ", Partial_Res_Chr_Idx, " Chromosome name ", ChrList_NameNum[Partial_Res_Chr_Idx], "****\n")
		writeLines(outtext, con=fp_out, sep="\n")

	} else {
		# indicates that we have to start from the beginning
		Partial_Res_Chr_Idx <- 0
		# dimension of output matrix
		nelem <- (as.integer(opt$window / opt$binsize) * 2) + 1
		# the matrix storing the absolute contact count information (from the given Hi-C matrix)
		APA_Matrix <- matrix(0, nelem, nelem)
		# this matrix stores the count of interactions (overlap) for individual matrix elements
		Count_APA_Mat <- matrix(0, nelem, nelem)
		outtext <- paste0("\n\n **** There is no prior computation - start from scratch ****\n")
		writeLines(outtext, con=fp_out, sep="\n")		
	}

	# the middle position of the matrix, corresponding to the interaction region
	Central_row <- as.integer((nelem - 1) / 2) + 1
	Central_col <- Central_row

	# This matrix stores the normalized contact count of binned interactions
	# computed by dividing APA_Matrix with Count_APA_Mat
	Norm_APA_Matrix <- matrix(0, nelem, nelem)

	outtext <- paste0("\n\n **** The matrix dimension : ", nelem, " X ", nelem, "****\n")
	writeLines(outtext, con=fp_out, sep="\n")

	#======================================
	# Loop through individual chromosomes,
	# extract data (of that particular chromosome) from the input and reference file
	# append the overlap information in the final matrix
	#======================================

	# checkpoint condition
	if (Partial_Res_Chr_Idx < length(ChrList_NameNum)) {

		# loop for individual chromosomes	
		for (chr_idx in ((Partial_Res_Chr_Idx+1):length(ChrList_NameNum))) {
			# output files storing the matrix after processing individual chromosomes
			chrName <- ChrList_NameNum[chr_idx]
			MatrixOutFile <- paste0(outdir, '/', onlyprefix, '_DumpAPAMatrix_after_', chrName ,'.bed')
			Matrix_Count_OutFile <- paste0(outdir, '/', onlyprefix, '_Dump_APA_Count_Matrix_after_', chrName ,'.bed')

			#===================================================
			# first extract input data of this chromosome
			#===================================================

			# depending on the chromosome name style provided in the input interaction file
			# either the chromosome name is a string, or a number
			if (opt$chrInp == TRUE) {
				curr.chr <- ChrList_Num[chr_idx]
			} else {
				curr.chr <- ChrList_NameNum[chr_idx]
			}

			# dump the current chromosome specific CIS interactions from the input file
			# maintaining the lower distance threshold between the interactions 
			# (specified in the option --distthrlow and --distthrhigh)
			ExtractChrData(InpFile=opt$InpFile, chrName=curr.chr, OutFile=InpTempFile, header=opt$headerInp, mid=opt$midInp, dist=c(opt$distthrlow, opt$distthrhigh))			
			cat(sprintf("\n Extracted input interaction file - current chromosome : %s", curr.chr))

			# check the extracted file size
			# continue if file size = 0
			if (file.info(InpTempFile)$size == 0) {
				next
			}

			#===================================================
			# now extract the reference data for this chromosome
			#===================================================

			# depending on the chromosome name style provided in the reference interaction file
			if (opt$chrRef == TRUE) {
				curr.chr <- ChrList_Num[chr_idx]
			} else {
				curr.chr <- ChrList_NameNum[chr_idx]
			}

			# dump the current chromosome specific CIS interactions from the reference file
			# no distance threshold is applied - important
			# since we want to get the overlapping interactions 
			# (at any distance) from the reference file
			ExtractChrData(InpFile=opt$RefFile, chrName=curr.chr, OutFile=RefTempFile, header=opt$headerRef, mid=opt$midRef)
			cat(sprintf("\n Extracted reference interaction file  - current chromosome : %s", curr.chr))

			# check the extracted file size
			# continue if file size = 0
			if (file.info(RefTempFile)$size == 0) {
				next
			}	

			#==================================
			# read the input interaction file (with respect to the current chromosome)
			# by construction, the extracted file has no header information
			#==================================
			InpData <- ReadContactData(InpTempFile, headerInp=FALSE, chrNUM=opt$chrInp, mid=opt$midInp, midCOL=c(2,4), binSize=opt$binsize)
			if (nrow(InpData) == 0) {
				next
			}
			cat(sprintf("\n Read input interaction file - current chromosome"))

			#==================================
			# read the reference interaction file
			# by construction, the extracted file has no header information
			#==================================
			RefData <- ReadContactData(RefTempFile, headerInp=FALSE, chrNUM=opt$chrRef, mid=opt$midRef, midCOL=c(2,4), binSize=opt$binsize)
			if (nrow(RefData) == 0) {
				next
			}
			cat(sprintf("\n Read reference interaction file  - current chromosome"))

			#==================================
			# if input data does not follow the boundary of bin size 
			# (i.e. have variable bin size like CHICAGO or hichipper loops)
			# adjust the input data to have bins aligned with the specified bin size
			# and also extend the boundaries to cover the allowed windows
			#==================================
			# creating vector for all the loop entries

			# interacting bin 1 (for all the loops) - start - aligned w.r.t the bin size			
			bin1_start_vec <- (round((InpData[, 2] * 1.0) / opt$binsize) * opt$binsize)
			# interacting bin 1 (for all the loops) - end - aligned w.r.t the bin size
			bin1_end_vec <- pmax((bin1_start_vec + opt$binsize), (round((InpData[, 3] * 1.0) / opt$binsize) * opt$binsize))
			# interacting bin 2 (for all the loops) - start - aligned w.r.t the bin size
			bin2_start_vec <- (round((InpData[, 5] * 1.0) / opt$binsize) * opt$binsize)
			# interacting bin 2 (for all the loops) - end - aligned w.r.t the bin size
			bin2_end_vec <- pmax((bin2_start_vec + opt$binsize), (round((InpData[, 6] * 1.0) / opt$binsize) * opt$binsize))

			# computing offset (window) for all the loop entries (as a vector)
			# applicable for each side of the bin
			bin1_offset_each_side_vec <- ((opt$binsize + (2 * opt$window) - (bin1_end_vec - bin1_start_vec)) / 2)
			bin2_offset_each_side_vec <- ((opt$binsize + (2 * opt$window) - (bin2_end_vec - bin2_start_vec)) / 2)

			# modify the input data to include the window (offset as well)
			InpData[, 2] <- pmax(0, (bin1_start_vec - bin1_offset_each_side_vec))
			InpData[, 3] <- bin1_end_vec + bin1_offset_each_side_vec
			InpData[, 5] <- pmax(0, (bin2_start_vec - bin2_offset_each_side_vec))
			InpData[, 6] <- bin2_end_vec + bin2_offset_each_side_vec

			#======================================
			# compute the overlap of the interaction file with the reference data
			# subject to the extension of window size in either direction
			#======================================
			LoopOv <- OverlapLoop(InpData, RefData, uniqov=FALSE, IDX=TRUE)

			# overlapping indices of input loops 
			idx_InpData_List <- LoopOv$A_AND_B
			
			# overlapping indices of reference loops 
			idx_RefData_List <- LoopOv$B_AND_A

			# compare the overlapping reference interaction and 
			# find the X and Y indices correponding to it
			# here we process the whole list
			x_idx_List <- nelem - as.integer(pmax(0,(InpData[idx_InpData_List, 3] - RefData[idx_RefData_List, 3])) / opt$binsize)
			y_idx_List <- nelem - as.integer(pmax(0,(InpData[idx_InpData_List, 6] - RefData[idx_RefData_List, 6])) / opt$binsize)
			if (opt$midRef == TRUE) {
				currcc_List <- RefData[idx_RefData_List, (opt$cccol + 2)]
			} else {
				currcc_List <- RefData[idx_RefData_List, opt$cccol]
			}

			# fill the APA matrix
			for (x_idx in (1:nelem)) {
				for (y_idx in (1:nelem)) {
					target_idx_list <- which((x_idx_List == x_idx) & (y_idx_List == y_idx))
					if (length(target_idx_list) > 0) {
						APA_Matrix[x_idx, y_idx] <- APA_Matrix[x_idx, y_idx] + sum(currcc_List[target_idx_list])
						Count_APA_Mat[x_idx, y_idx] <- Count_APA_Mat[x_idx, y_idx] + length(target_idx_list)
					}
				}
			}

			outtext <- paste0("\n **** Processed all interactions of input data for the chromosome: ", ChrList_NameNum[chr_idx])
			writeLines(outtext, con=fp_out, sep="\n")

			# clear the memories:
			InpData <- InpData[c(), ]
			RefData <- RefData[c(), ]

			# dump intermediate matrix
			zz <- file(description=MatrixOutFile, "w") 
			write.table(APA_Matrix, file=zz, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)
			close(zz)	

			# dump intermediate overlap count matrix
			zz <- file(description=Matrix_Count_OutFile,"w") 
			write.table(Count_APA_Mat, file=zz, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)
			close(zz)	

			# as the partial result for the latest chromosome is dumped, no need to 
			# dump the partial result involving the last chromosome
			oldMatrixOutFile <- paste0(outdir, '/', onlyprefix, '_DumpAPAMatrix_after_', ChrList_NameNum[chr_idx-1],'.bed')
			oldMatrix_Count_OutFile <- paste0(outdir, '/', onlyprefix, '_Dump_APA_Count_Matrix_after_', ChrList_NameNum[chr_idx-1],'.bed')
			if (file.exists(oldMatrixOutFile) == TRUE) {
				system(paste0('rm ', oldMatrixOutFile))	
				outtext <- paste0("\n **** Removed dumped interaction matrix for the chromosome: ", ChrList_NameNum[chr_idx-1])
				writeLines(outtext, con=fp_out, sep="\n")		
			}
			if (file.exists(oldMatrix_Count_OutFile) == TRUE) {
				system(paste0('rm ', oldMatrix_Count_OutFile))
				outtext <- paste0("\n **** Removed dumped overlap count matrix for the chromosome: ", ChrList_NameNum[chr_idx-1])
				writeLines(outtext, con=fp_out, sep="\n")			
			}

		}	# end chromosome main loop
	}	# end check point

	# close the text file connection
	close(fp_out)	

	# delete the temporary files
	if (file.exists(InpTempFile) == TRUE) {
		system(paste0('rm ', InpTempFile))	
	}
	if (file.exists(RefTempFile) == TRUE) {
		system(paste0('rm ', RefTempFile))
	}

	#========================================================
	# dump the generated matrix (raw contact count information)
	zz <- file(description=MatrixFinalOutFile,"w") 
	write.table(APA_Matrix, file=zz, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)
	close(zz)

	# dump the generated overlap count matrix 
	zz <- file(description=Matrix_Count_FinalOutFile,"w") 
	write.table(Count_APA_Mat, file=zz, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)
	close(zz)

	# now delete the partial result for the last chromosome
	oldMatrixOutFile <- paste0(outdir, '/', onlyprefix, '_DumpAPAMatrix_after_', ChrList_NameNum[length(ChrList_NameNum)],'.bed')
	oldMatrix_Count_OutFile <- paste0(outdir, '/', onlyprefix, '_Dump_APA_Count_Matrix_after_', ChrList_NameNum[length(ChrList_NameNum)],'.bed')
	if (file.exists(oldMatrixOutFile) == TRUE) {
		system(paste0('rm ', oldMatrixOutFile))	
		outtext <- paste0("\n **** Removed dumped interaction matrix for the chromosome: ", ChrList_NameNum[length(ChrList_NameNum)])
		writeLines(outtext, con=fp_out, sep="\n")		
	}
	if (file.exists(oldMatrix_Count_OutFile) == TRUE) {
		system(paste0('rm ', oldMatrix_Count_OutFile))
		outtext <- paste0("\n **** Removed dumped overlap count matrix for the chromosome: ", ChrList_NameNum[length(ChrList_NameNum)])
		writeLines(outtext, con=fp_out, sep="\n")			
	}	

	# now derive the normalized (mean) contact count matrix
	# by dividing APA_Matrix with the Count_APA_Mat contents
	for (i in (1:nelem)) {
		for (j in (1:nelem)) {
			if (Count_APA_Mat[i,j] > 0) {
				Norm_APA_Matrix[i,j] <- (APA_Matrix[i,j] * 1.0) / Count_APA_Mat[i,j]
			} else {
				Norm_APA_Matrix[i,j] <- (APA_Matrix[i,j] * 1.0 + 1) / (Count_APA_Mat[i,j] + 1)
			}
		}
	}

	# dump the normalized matrix contents
	zz <- file(description=NormMatrixFinalOutFile,"w") 
	write.table(Norm_APA_Matrix, file=zz, row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE, append=FALSE)	
	close(zz)

} else {

	cat(sprintf("\n All the matrices are pre-computed !!! just plotting !!!"))

	APA_Matrix <- as.matrix(read.table(MatrixFinalOutFile, header=F))
	Count_APA_Mat <- as.matrix(read.table(Matrix_Count_FinalOutFile, header=F))
	Norm_APA_Matrix <- as.matrix(read.table(NormMatrixFinalOutFile, header=F))

	# dimension of output matrix
	nelem <- ncol(APA_Matrix)

	# the middle position of the matrix, corresponding to the interaction region
	Central_row <- as.integer((nelem - 1) / 2) + 1
	Central_col <- Central_row

}
#========================================================
# find APA related score values from both raw contact count and normalized contact count matrices

# assuming coordinates of bins are positive in X axis right direction
# and positive Y axis in top direction

# top right location (both X and Y are positive)
pos_x_pos_y_low_bin_x = Central_col + as.integer(nelem / 4)
pos_x_pos_y_high_bin_x = nelem
pos_x_pos_y_low_bin_y = Central_row + as.integer(nelem / 4)
pos_x_pos_y_high_bin_y = nelem
APA_score_pos_x_pos_y = (APA_Matrix[Central_col, Central_row] * 1.0) / mean(APA_Matrix[pos_x_pos_y_low_bin_x:pos_x_pos_y_high_bin_x, pos_x_pos_y_low_bin_y:pos_x_pos_y_high_bin_y])
APA_score_norm_pos_x_pos_y = (Norm_APA_Matrix[Central_col, Central_row] * 1.0) / mean(Norm_APA_Matrix[pos_x_pos_y_low_bin_x:pos_x_pos_y_high_bin_x, pos_x_pos_y_low_bin_y:pos_x_pos_y_high_bin_y])

# bottom right location (X positive, Y negative)
pos_x_neg_y_low_bin_x = Central_col + as.integer(nelem / 4)
pos_x_neg_y_high_bin_x = nelem
pos_x_neg_y_low_bin_y = 1
pos_x_neg_y_high_bin_y = as.integer(nelem / 4)
APA_score_pos_x_neg_y = (APA_Matrix[Central_col, Central_row] * 1.0) / mean(APA_Matrix[pos_x_neg_y_low_bin_x:pos_x_neg_y_high_bin_x, pos_x_neg_y_low_bin_y:pos_x_neg_y_high_bin_y])
APA_score_norm_pos_x_neg_y = (Norm_APA_Matrix[Central_col, Central_row] * 1.0) / mean(Norm_APA_Matrix[pos_x_neg_y_low_bin_x:pos_x_neg_y_high_bin_x, pos_x_neg_y_low_bin_y:pos_x_neg_y_high_bin_y])

# top left location (X negative, Y positive)
neg_x_pos_y_low_bin_x = 1
neg_x_pos_y_high_bin_x = as.integer(nelem / 4)
neg_x_pos_y_low_bin_y = Central_row + as.integer(nelem / 4)
neg_x_pos_y_high_bin_y = nelem
APA_score_neg_x_pos_y = (APA_Matrix[Central_col, Central_row] * 1.0) / mean(APA_Matrix[neg_x_pos_y_low_bin_x:neg_x_pos_y_high_bin_x, neg_x_pos_y_low_bin_y:neg_x_pos_y_high_bin_y])
APA_score_norm_neg_x_pos_y = (Norm_APA_Matrix[Central_col, Central_row] * 1.0) / mean(Norm_APA_Matrix[neg_x_pos_y_low_bin_x:neg_x_pos_y_high_bin_x, neg_x_pos_y_low_bin_y:neg_x_pos_y_high_bin_y])

# bottom left (both X and Y negative)
neg_x_neg_y_low_bin_x = 1
neg_x_neg_y_high_bin_x = as.integer(nelem / 4)
neg_x_neg_y_low_bin_y = 1
neg_x_neg_y_high_bin_y = as.integer(nelem / 4)
APA_score_neg_x_neg_y = (APA_Matrix[Central_col, Central_row] * 1.0) / mean(APA_Matrix[neg_x_neg_y_low_bin_x:neg_x_neg_y_high_bin_x, neg_x_neg_y_low_bin_y:neg_x_neg_y_high_bin_y])
APA_score_norm_neg_x_neg_y = (Norm_APA_Matrix[Central_col, Central_row] * 1.0) / mean(Norm_APA_Matrix[neg_x_neg_y_low_bin_x:neg_x_neg_y_high_bin_x, neg_x_neg_y_low_bin_y:neg_x_neg_y_high_bin_y])

# in the MANGO paper, 15-30 Kb (or in general the low and high thresholds provided in this command line option)
# in positive X and Y directions are used for computing the APA score
low_pixel_offset <- as.integer(opt$low / opt$binsize)
high_pixel_offset <- as.integer(opt$high / opt$binsize)

# implementation of MANGO paper based scores
# we check 15-30 Kb (or the specified distance thresholds) on both side (positive and negative X and Y directions)
# The MANGO paper uses positive x 15-30 Kb offset and negative y 15-30 Kb offset
# that is, APA_score_MANGO_paper_pos_x_neg_y and Norm_APA_score_MANGO_paper_pos_x_neg_y
APA_score_MANGO_paper_pos_x_pos_y <- (APA_Matrix[Central_col, Central_row] * 1.0) / mean(APA_Matrix[(Central_col + low_pixel_offset):(Central_col + high_pixel_offset), (Central_row + low_pixel_offset):(Central_row + high_pixel_offset)])
Norm_APA_score_MANGO_paper_pos_x_pos_y <- (Norm_APA_Matrix[Central_col, Central_row] * 1.0) / mean(Norm_APA_Matrix[(Central_col + low_pixel_offset):(Central_col + high_pixel_offset), (Central_row + low_pixel_offset):(Central_row + high_pixel_offset)])

APA_score_MANGO_paper_pos_x_neg_y <- (APA_Matrix[Central_col, Central_row] * 1.0) / mean(APA_Matrix[(Central_col + low_pixel_offset):(Central_col + high_pixel_offset), (Central_row - low_pixel_offset):(Central_row - high_pixel_offset)])
Norm_APA_score_MANGO_paper_pos_x_neg_y <- (Norm_APA_Matrix[Central_col, Central_row] * 1.0) / mean(Norm_APA_Matrix[(Central_col + low_pixel_offset):(Central_col + high_pixel_offset), (Central_row - low_pixel_offset):(Central_row - high_pixel_offset)])

APA_score_MANGO_paper_neg_x_pos_y <- (APA_Matrix[Central_col, Central_row] * 1.0) / mean(APA_Matrix[(Central_col - low_pixel_offset):(Central_col - high_pixel_offset), (Central_row + low_pixel_offset):(Central_row + high_pixel_offset)])
Norm_APA_score_MANGO_paper_neg_x_pos_y <- (Norm_APA_Matrix[Central_col, Central_row] * 1.0) / mean(Norm_APA_Matrix[(Central_col - low_pixel_offset):(Central_col - high_pixel_offset), (Central_row + low_pixel_offset):(Central_row + high_pixel_offset)])

APA_score_MANGO_paper_neg_x_neg_y <- (APA_Matrix[Central_col, Central_row] * 1.0) / mean(APA_Matrix[(Central_col - low_pixel_offset):(Central_col - high_pixel_offset), (Central_row - low_pixel_offset):(Central_row - high_pixel_offset)])
Norm_APA_score_MANGO_paper_neg_x_neg_y <- (Norm_APA_Matrix[Central_col, Central_row] * 1.0) / mean(Norm_APA_Matrix[(Central_col - low_pixel_offset):(Central_col - high_pixel_offset), (Central_row - low_pixel_offset):(Central_row - high_pixel_offset)])

# compute the mean APA score (w.r.t MANGO paper definition) for all four corners
APA_score_MANGO_paper_mean <- mean(c(APA_score_MANGO_paper_pos_x_pos_y, APA_score_MANGO_paper_pos_x_neg_y, APA_score_MANGO_paper_neg_x_pos_y, APA_score_MANGO_paper_neg_x_neg_y))
Norm_APA_score_MANGO_paper_mean <- mean(c(Norm_APA_score_MANGO_paper_pos_x_pos_y, Norm_APA_score_MANGO_paper_pos_x_neg_y, Norm_APA_score_MANGO_paper_neg_x_pos_y, Norm_APA_score_MANGO_paper_neg_x_neg_y))

# compute the ratio of central bin to the remaining matrix

APA_Ratio_Central_Rest <- (APA_Matrix[Central_col, Central_row] * 1.0) / mean(c(mean(APA_Matrix[1:(Central_col - 1), 1:(Central_row - 1)]), mean(APA_Matrix[(Central_col + 1):nelem, 1:(Central_row - 1)]), mean(APA_Matrix[1:(Central_col - 1), (Central_row + 1):nelem]), mean(APA_Matrix[(Central_col + 1):nelem, (Central_row + 1):nelem])))

Norm_APA_Ratio_Central_Rest <- (Norm_APA_Matrix[Central_col, Central_row] * 1.0) / mean(c(mean(Norm_APA_Matrix[1:(Central_col - 1), 1:(Central_row - 1)]), mean(Norm_APA_Matrix[(Central_col + 1):nelem, 1:(Central_row - 1)]), mean(Norm_APA_Matrix[1:(Central_col - 1), (Central_row + 1):nelem]), mean(Norm_APA_Matrix[(Central_col + 1):nelem, (Central_row + 1):nelem])))

# compute the Z score of central bin, using individual portions of the matrix

# with respect to top right segment
APA_zscore_pos_x_pos_y <- ((APA_Matrix[Central_col, Central_row] - mean(APA_Matrix[pos_x_pos_y_low_bin_x:pos_x_pos_y_high_bin_x, pos_x_pos_y_low_bin_y:pos_x_pos_y_high_bin_y])) * 1.0) / sd(APA_Matrix[pos_x_pos_y_low_bin_x:pos_x_pos_y_high_bin_x, pos_x_pos_y_low_bin_y:pos_x_pos_y_high_bin_y])

Norm_APA_zscore_pos_x_pos_y <- ((Norm_APA_Matrix[Central_col, Central_row] - mean(Norm_APA_Matrix[pos_x_pos_y_low_bin_x:pos_x_pos_y_high_bin_x, pos_x_pos_y_low_bin_y:pos_x_pos_y_high_bin_y])) * 1.0) / sd(Norm_APA_Matrix[pos_x_pos_y_low_bin_x:pos_x_pos_y_high_bin_x, pos_x_pos_y_low_bin_y:pos_x_pos_y_high_bin_y])

# with respect to top left segment
APA_zscore_neg_x_pos_y <- ((APA_Matrix[Central_col, Central_row] - mean(APA_Matrix[neg_x_pos_y_low_bin_x:neg_x_pos_y_high_bin_x, neg_x_pos_y_low_bin_y:neg_x_pos_y_high_bin_y])) * 1.0) / sd(APA_Matrix[neg_x_pos_y_low_bin_x:neg_x_pos_y_high_bin_x, neg_x_pos_y_low_bin_y:neg_x_pos_y_high_bin_y])

Norm_APA_zscore_neg_x_pos_y <- ((Norm_APA_Matrix[Central_col, Central_row] - mean(Norm_APA_Matrix[neg_x_pos_y_low_bin_x:neg_x_pos_y_high_bin_x, neg_x_pos_y_low_bin_y:neg_x_pos_y_high_bin_y])) * 1.0) / sd(Norm_APA_Matrix[neg_x_pos_y_low_bin_x:neg_x_pos_y_high_bin_x, neg_x_pos_y_low_bin_y:neg_x_pos_y_high_bin_y])

# with respect to bottom right segment
APA_zscore_pos_x_neg_y <- ((APA_Matrix[Central_col, Central_row] - mean(APA_Matrix[pos_x_neg_y_low_bin_x:pos_x_neg_y_high_bin_x, pos_x_neg_y_low_bin_y:pos_x_neg_y_high_bin_y])) * 1.0) / sd(APA_Matrix[pos_x_neg_y_low_bin_x:pos_x_neg_y_high_bin_x, pos_x_neg_y_low_bin_y:pos_x_neg_y_high_bin_y])

Norm_APA_zscore_pos_x_neg_y <- ((Norm_APA_Matrix[Central_col, Central_row] - mean(Norm_APA_Matrix[pos_x_neg_y_low_bin_x:pos_x_neg_y_high_bin_x, pos_x_neg_y_low_bin_y:pos_x_neg_y_high_bin_y])) * 1.0) / sd(Norm_APA_Matrix[pos_x_neg_y_low_bin_x:pos_x_neg_y_high_bin_x, pos_x_neg_y_low_bin_y:pos_x_neg_y_high_bin_y])

# with respect to bottom left segment
APA_zscore_neg_x_neg_y <- ((APA_Matrix[Central_col, Central_row] - mean(APA_Matrix[neg_x_neg_y_low_bin_x:neg_x_neg_y_high_bin_x, neg_x_neg_y_low_bin_y:neg_x_neg_y_high_bin_y])) * 1.0) / sd(APA_Matrix[neg_x_neg_y_low_bin_x:neg_x_neg_y_high_bin_x, neg_x_neg_y_low_bin_y:neg_x_neg_y_high_bin_y])

Norm_APA_zscore_neg_x_neg_y <- ((Norm_APA_Matrix[Central_col, Central_row] - mean(Norm_APA_Matrix[neg_x_neg_y_low_bin_x:neg_x_neg_y_high_bin_x, neg_x_neg_y_low_bin_y:neg_x_neg_y_high_bin_y])) * 1.0) / sd(Norm_APA_Matrix[neg_x_neg_y_low_bin_x:neg_x_neg_y_high_bin_x, neg_x_neg_y_low_bin_y:neg_x_neg_y_high_bin_y])

#================================
# write down the raw APA count and associated metrics

# first the raw contact count based APA scores
textfile <- paste0(outdir, '/', onlyprefix, '_APA_Result_Summary.txt')
fp_out <- file(textfile, "w")

outtext <- paste0("\n\n **** Raw APA scores **** \n\n top right (+x, +y): ", APA_score_pos_x_pos_y, "\n Bottom right (+x, -y): ", APA_score_pos_x_neg_y, "\n Top left (-x, +y): ", APA_score_neg_x_pos_y, "\n Bottom left (-x, -y): ", APA_score_neg_x_neg_y)
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n\n **** Raw APA scores (with respect to MANGO paper) ***** \n\n top right (+x, +y): ", APA_score_MANGO_paper_pos_x_pos_y, "\n bottom right (original) (+x, -y): ", APA_score_MANGO_paper_pos_x_neg_y, "\n top left (-x, +y): ", APA_score_MANGO_paper_neg_x_pos_y, "\n bottom left (-x, -y): ", APA_score_MANGO_paper_neg_x_neg_y, "\n Mean with respect to all four corners: ", APA_score_MANGO_paper_mean)
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n\n **** Raw APA scores **** \n\n Ratio of the central bin to the rest of the matrix: ", APA_Ratio_Central_Rest, " \n\n zscore top right (+x, +y): ", APA_zscore_pos_x_pos_y, "\n bottom right (+x, -y): ", APA_zscore_pos_x_neg_y, "\n top left (-x, +y): ", APA_zscore_neg_x_pos_y, "\n bottom left (-x, -y): ", APA_zscore_neg_x_neg_y)
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n\n **** Normalized APA (raw APA / occurrence) scores **** \n\n top right (+x, +y): ", APA_score_norm_pos_x_pos_y, "\n Bottom right (+x, -y): ", APA_score_norm_pos_x_neg_y, "\n Top left (-x, +y): ", APA_score_norm_neg_x_pos_y, "\n Bottom left (-x, -y): ", APA_score_norm_neg_x_neg_y)
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n\n **** Normalized APA (raw APA / occurrence) scores (with respect to MANGO paper) ***** \n\n top right (+x, +y): ", Norm_APA_score_MANGO_paper_pos_x_pos_y, "\n bottom right (original) (+x, -y): ", Norm_APA_score_MANGO_paper_pos_x_neg_y, "\n top left (-x, +y): ", Norm_APA_score_MANGO_paper_neg_x_pos_y, "\n bottom left (-x, -y): ", Norm_APA_score_MANGO_paper_neg_x_neg_y, "\n Mean with respect to all four corners: ", Norm_APA_score_MANGO_paper_mean)
writeLines(outtext, con=fp_out, sep="\n")

outtext <- paste0("\n\n **** Normalized APA (raw APA / occurrence) scores **** \n\n Ratio of the central bin to the rest of the matrix: ", Norm_APA_Ratio_Central_Rest, " \n\n zscore top right (+x, +y): ", Norm_APA_zscore_pos_x_pos_y, "\n bottom right (+x, -y): ", Norm_APA_zscore_pos_x_neg_y, "\n top left (-x, +y): ", Norm_APA_zscore_neg_x_pos_y, "\n bottom left (-x, -y): ", Norm_APA_zscore_neg_x_neg_y)
writeLines(outtext, con=fp_out, sep="\n")

close(fp_out)

#========================================================
# in all of these plots, value of Y axis increases to the bottom (according to the graphics convention)
# the legends are adjusted accordingly
#========================================================
# print the raw APA matrix as a colored plot
plotfile <- paste0(outdir, '/', onlyprefix, '_Raw_APA_Plot.pdf')

# comment - sourya
melt_matrix <- as.data.frame(melt(APA_Matrix))
colnames(melt_matrix) <- c("X", "Y", "CC")

titlestr <- paste0("APA= ", round(APA_score_MANGO_paper_pos_x_neg_y, 2))

pos_x_neg_y_label <- paste0('APA: ', round(APA_score_pos_x_neg_y, 2))
pos_x_pos_y_label <- paste0('APA: ', round(APA_score_pos_x_pos_y, 2))
neg_x_neg_y_label <- paste0('APA: ', round(APA_score_neg_x_neg_y, 2))
neg_x_pos_y_label <- paste0('APA: ', round(APA_score_neg_x_pos_y, 2))

centerlabel <- paste0('R: ', round(APA_Ratio_Central_Rest, 2))

curr_plot1 <- ggplot(melt_matrix, aes(x = X, y = Y)) + geom_raster(aes(fill=CC)) + scale_fill_gradientn(colours=colorvec) + annotate("text", x=0, y=1, label= ((-1) * (opt$window / 1000)), size=8) + annotate("text", x=(nelem+1), y=1, label= (opt$window / 1000), size=8) + annotate("text", x=0, y=nelem, label= (opt$window / 1000), size=8) + xlab("Kb from U Anchor") + ylab("Kb from D Anchor") + theme(text = element_text(size=FontSize)) + ggtitle(titlestr) + annotate("text", x=(nelem-3), y=3, label= pos_x_neg_y_label, size=8) + annotate("text", x=(nelem-3), y=(nelem-3), label = pos_x_pos_y_label, size=8) + annotate("text", x=3, y=3, label = neg_x_neg_y_label, size=8) + annotate("text", x=3, y=(nelem-3), label = neg_x_pos_y_label, size=8) + annotate("text", x=(Central_col-1), y=Central_row, label = centerlabel, size=8) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

ggsave(plotfile, plot = curr_plot1, width=PlotWidth, height=PlotHeight)

# also add one plot which do not have any labels
plotfile_without_label <- paste0(outdir, '/', onlyprefix, '_Raw_APA_Plot_without_label.pdf')
curr_plot1A <- ggplot(melt_matrix, aes(x = X, y = Y)) + geom_raster(aes(fill=CC)) + scale_fill_gradientn(colours=colorvec) + theme(text = element_text(size=FontSize)) + ggtitle(titlestr) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
ggsave(plotfile_without_label, plot = curr_plot1A, width=PlotWidth, height=PlotHeight)

# #========================================================
# print the normalized APA matrix as a colored plot
plotfile <- paste0(outdir, '/', onlyprefix, '_Normalized_APA_Plot.pdf')

# comment - sourya
melt_matrix <- as.data.frame(melt(Norm_APA_Matrix))
colnames(melt_matrix) <- c("X", "Y", "CC")

titlestr <- paste0("APA= ", round(Norm_APA_score_MANGO_paper_pos_x_neg_y, 2))

pos_x_neg_y_label <- paste0('APA: ', round(APA_score_norm_pos_x_neg_y, 2))
pos_x_pos_y_label <- paste0('APA: ', round(APA_score_norm_pos_x_pos_y, 2))
neg_x_neg_y_label <- paste0('APA: ', round(APA_score_norm_neg_x_neg_y, 2))
neg_x_pos_y_label <- paste0('APA: ', round(APA_score_norm_neg_x_pos_y, 2))

centerlabel <- paste0('R: ', round(Norm_APA_Ratio_Central_Rest, 2))

# limit the normalized contact (color scale) from 1 to 20
curr_plot2 <- ggplot(melt_matrix, aes(x = X, y = Y)) + geom_raster(aes(fill=CC)) + scale_fill_gradientn(colours=colorvec, limits=c(1, opt$cclim)) + annotate("text", x=0, y=1, label= ((-1) * (opt$window / 1000)), size=8) + annotate("text", x=(nelem+1), y=1, label= (opt$window / 1000), size=8) + annotate("text", x=0, y=nelem, label= (opt$window / 1000), size=8) + xlab("Kb from U Anchor") + ylab("Kb from D Anchor") + theme(text = element_text(size=FontSize)) + ggtitle(titlestr) + annotate("text", x=(nelem-3), y=3, label= pos_x_neg_y_label, size=8) + annotate("text", x=(nelem-3), y=(nelem-3), label = pos_x_pos_y_label, size=8) + annotate("text", x=3, y=3, label = neg_x_neg_y_label, size=8) + annotate("text", x=3, y=(nelem-3), label = neg_x_pos_y_label, size=8) + annotate("text", x=(Central_col-1), y=Central_row, label = centerlabel, size=8) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

ggsave(plotfile, plot = curr_plot2, width=PlotWidth, height=PlotHeight)

# also add one plot which do not have any labels
plotfile_without_label <- paste0(outdir, '/', onlyprefix, '_Normalized_APA_Plot_without_label.pdf')

curr_plot2A <- ggplot(melt_matrix, aes(x = X, y = Y)) + geom_raster(aes(fill=CC)) + scale_fill_gradientn(colours=colorvec, limits=c(1, opt$cclim)) + theme(text = element_text(size=FontSize)) + ggtitle(titlestr) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

ggsave(plotfile_without_label, plot = curr_plot2A, width=PlotWidth, height=PlotHeight)
