#!/usr/bin/env Rscript

#===========================================================
# R script which takes input of one or more peak files
# computes their relative overlap
# and finally plots a venn diagram showing the relative overlap

# Author: Sourya Bhattacharyya
# Vijay-Ay lab, LJI
#===========================================================

library(UtilRPckg)
library(optparse)
library(data.table)

# venn diagram 
library(venn)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

#================================================
# get the union of indices of the merged peak set
Get_Union_Set <- function(Ov_Idx_MergedRef, set_idx) {
	for (i in 1:length(set_idx)) {
		idx <- set_idx[i]
		if (i == 1){
			UnionSet <- Ov_Idx_MergedRef[[idx]]
		} else {
			UnionSet <- union(UnionSet, Ov_Idx_MergedRef[[idx]])
		}
	}
	return(UnionSet)
}

#================================================
option_list = list(
	make_option(c("--FileList"), type="character", default=NULL, help="Comma or colon separated list of peak files (default %default)"),
	make_option(c("--Labels"), type="character", default=NULL, help="Comma or colon separated list of strings, associated with individual peak files. Used as the labels of the generated venn diagram as well. Default: category1, category2, etc.."),
  	make_option(c("--OutPrefix"), type="character", default=NULL, help="Prefix (including the output directory) of the output files."),
  	make_option(c("--Summit"), type="logical", action="store_true", default=FALSE, help="If TRUE, peaks relative to the summit position (and of a specified offset) are only considered."),
  	make_option(c("--Offset"), type="integer", action="store", default=500, help="If Summit option is TRUE, peaks relative to the summit position and spanning this Offset number of base pairs are onky considered. Default = 500, means the peaks spanning 500 bp around the summit would be considered."),
  	make_option(c("--Dump"), type="logical", action="store_true", default=FALSE, help="If TRUE, peaks in different partitions are dumped.")
); 

parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args

InpFileList <- as.character(unlist(strsplit(opt$FileList,"[,:]")))

# categorical information of input interaction / loop files
if (is.null(opt$Labels)) {
	InpLabels <- paste0('category', seq(1:length(InpFileList)))
} else {
	InpLabels <- as.character(unlist(strsplit(opt$Labels, "[,:]")))
	if (length(InpLabels) < length(InpFileList)) {
		stop("Number of labels provided is not same as the number of input files - check the option --Labels \n", call.=FALSE)
	}
}

if (is.null(opt$OutPrefix)) {
	print_help(opt_parser)
	stop("Output prefix is not given - check the option --OutPrefix \n", call.=FALSE)
}

# create the output directory
OutDir <- dirname(opt$OutPrefix)
system(paste('mkdir -p', OutDir))

# only the prefix string of the output file names
OnlyPrefix <- basename(opt$OutPrefix)

#===================================
# first create a master file containing peaks
MasterFile <- paste0(OutDir, '/', OnlyPrefix, 'MasterFile_Peak.bed')
MasterTempFile <- paste0(OutDir, '/', OnlyPrefix, 'MasterFile_Peak_Temp.bed')
MasterTempFile1 <- paste0(OutDir, '/', OnlyPrefix, 'MasterFile_Peak_Temp1.bed')

tempfile <- paste0(OutDir, '/temp_extracted_peak.bed')
if (file.exists(MasterFile) == FALSE) {
	# first merge individual peaks
	for (i in (1:length(InpFileList))) {
		currfile <- InpFileList[i]
		# extract current peaks
		UtilRPckg::ExtractPeaksOrSummits(currfile, chrFilt=TRUE, OutFile=tempfile, Summit=opt$Summit, Offset=opt$Offset)
		if (i == 1) {
			system(paste("cat", tempfile, ">", MasterTempFile))
		} else {
			system(paste("cat", tempfile, ">>", MasterTempFile))
		}
	}
	# sort and remove duplicates
	system(paste("cut -f1-3 ", MasterTempFile, " | sort -k1,1 -k2,2n -k3,3n | uniq >", MasterTempFile1))
	system(paste("bedtools merge -i", MasterTempFile1, ">", MasterFile))
	cat(sprintf("\n\n\n *** extracted peaks for all input files and merged ****** \n\n\n"))
}

# read the merged peak file
Merged_PeakData <- data.table::fread(MasterFile, header=F)

#===================================
# then compute overlap of peaks for this master file
# and for individual peak files
vennPlotFile <- paste0(OutDir, '/', OnlyPrefix, 'Venn_Peak_Overlap_NEW.pdf')
vennPlotFile_NoCount <- paste0(OutDir, '/', OnlyPrefix, 'Venn_Peak_Overlap_NEW_NoCount.pdf')
# if (file.exists(vennPlotFile) == FALSE) {

	# list of indices - 
	# overlapping indices with respect to merged peak files
	Ov_Idx_MergedRef <- list()
	# overlapping indices with respect to individual input peak files
	Ov_Idx_Inp <- list()

	for (i in (1:length(InpFileList))) {
		InpData <- data.table::fread(InpFileList[i], header=F)
		# some peak intervals may be duplicated - so get the unique peaks
		InpData <- unique(InpData[, 1:3])
		cat(sprintf("\n *** computing peak overlap - read input file idx : %s \n", i))

		# compute peak overlap using the custom function
		CurrOv <- UtilRPckg::Overlap1D(Merged_PeakData, InpData)
		
		Ov_Idx_MergedRef[[i]] <- CurrOv$A_AND_B
		Ov_Idx_Inp[[i]] <- CurrOv$B_AND_A

		cat(sprintf("\n **** overlap between merged peak to category: %s length overlap w.r.t merged peak: %s  length overlap w.r.t input: %s  ", InpLabels[i], length(Ov_Idx_MergedRef[[i]]), length(Ov_Idx_Inp[[i]])))
	}

	# plot the venn diagram with these overlapping indices
	# and place a legend at the bottom
	names(Ov_Idx_MergedRef) <- seq(1, length(InpLabels))
	legendstr <- c()
	for (i in (1:length(InpLabels))) {
		legendstr <- c(legendstr, paste(i, ":", InpLabels[i]))
	}

	#===================================
	# plot venn diagram displaying the counts
	# output file storing the plot
	pdf(vennPlotFile, width=8, height=5)

	if (length(InpLabels)<=2) {
    	cexil_val <- 1.5
    } else if (length(InpLabels)<=4) {
    	cexil_val <- 1
    } else if (length(InpLabels)<=5) {
		cexil_val <- 0.6
    } else if (length(InpLabels)<=6) {
		cexil_val <- 0.5
    } else {
		cexil_val <- 0.45
    }
    # if (length(InpLabels)<=3) {
        venn(Ov_Idx_MergedRef, counts=TRUE, zcolor="style", cexil=cexil_val)
    # } else {
    #     venn(Ov_Idx_MergedRef, counts=TRUE, ellipse=TRUE, zcolor="style", cexil=cexil_val)
    # }

	legend("bottomright",legend=legendstr, lty=1, lwd=1, cex=1.2)
	dev.off()

	#===================================
	# plot venn diagram without displaying the counts
	pdf(vennPlotFile_NoCount, width=8, height=5)
 	# if (length(InpLabels)<=3) {
        venn(Ov_Idx_MergedRef, counts=FALSE, zcolor="style")
    # } else {
    #     venn(Ov_Idx_MergedRef, ellipse=TRUE, zcolor="style")
    # }
	legend("bottomright",legend=legendstr, lty=1, lwd=1, cex=1.2)
	dev.off()
	#===================================

# }

#***********************************
# create a pie chart representing the 
# peaks present in different cell types
# peaks present in one cell, two cells, etc.
#***********************************

# matrix of dimension : row = nrow(Merged_PeakData) and col = length(InpFileList)
# stores binary values regarding the overlap of peaks for individual input data
Peak_Mat <- matrix(0, nrow=nrow(Merged_PeakData), ncol=length(InpFileList))

# place those matrix elements with 1
# which overlap with input peak file
for (j in (1:length(InpFileList))) {
	Peak_Mat[Ov_Idx_MergedRef[[j]], j] <- 1
}

# derive the row wise sum of this matrix
RowSumPeakMat <- apply(Peak_Mat, 1, sum)
NumPeak_All <- c()
for (cell_count in c(1:length(InpFileList))) {
	cat(sprintf("\n Analyzing peaks in %s no of cells -- ", cell_count))
	num_Peaks <- length(which(RowSumPeakMat==cell_count))
	cat(sprintf("\n Number of peaks belonging to %s cell types is : %s ", cell_count, num_Peaks))
	NumPeak_All <- c(NumPeak_All, num_Peaks)
}

# plot the pie chart
lbls <- paste0(c(1:length(InpFileList)), "_Cell")
pct <- round(NumPeak_All/sum(NumPeak_All)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 

PlotFile <- paste0(OutDir, '/', OnlyPrefix, 'Pie_Chart_Cell_Specific_Peaks.pdf')
pdf(PlotFile, width=6, height=4)
pie(NumPeak_All, labels = lbls, col=rainbow(length(lbls)), main="Peaks in different cell types")
dev.off()

#==========================
# dump some of the peak sets
#==========================
if (opt$Dump == TRUE) {

	# generate exclusive peaks for each categories and dump into separate output file
	for (i in (1:length(InpFileList))) {
		# generate exclusive peaks for the category i, i.e. the set InpLabels[i]
		cat(sprintf("\n\n *** Generating exclusive peaks for the category -- : %s no of peaks of this category (w.r.t merged data) : %s ", i, length(Ov_Idx_MergedRef[[i]])))

		# generate the union of peaks covered by all other input categories (w.r.t the merged peak data)
		other_cat_set <- setdiff(seq(1, length(InpFileList)), i)
		UnionPeakSet_Other_Cat <- Get_Union_Set(Ov_Idx_MergedRef, other_cat_set)
		cat(sprintf("\n *** no of peaks covered by all other categories : %s ", length(UnionPeakSet_Other_Cat)))

		# exclusive peaks only to the current category i
		CurrCat_ExclPeak_Set <- setdiff(Ov_Idx_MergedRef[[i]], UnionPeakSet_Other_Cat)
		cat(sprintf("\n *** no of peaks exclusive to the category : %s label : %s is : %s ", i, InpLabels[i], length(CurrCat_ExclPeak_Set)))

		# dump the exclusive peaks
		if (length(CurrCat_ExclPeak_Set) > 0) {
			CurrOutFile <- paste0(dirname(opt$OutPrefix), '/Peaks_Exclusive_', InpLabels[i], '.bed')
			write.table(Merged_PeakData[CurrCat_ExclPeak_Set, ], CurrOutFile, row.names=F, col.names=F, sep="\t", quote=F, append=F)
		}
	}

	# now generate the peaks common to two input sets
	for (n1 in (1:(length(InpFileList) - 1))) {
		for (n2 in ((n1+1):length(InpFileList))) {
			cat(sprintf("\n\n *** Dump peaks common to two categories - processing indices %s and %s ", n1, n2))
			Curr_Common_Peak_Set <- intersect(Ov_Idx_MergedRef[[n1]], Ov_Idx_MergedRef[[n2]])
			# dump the peaks common to these two categories
			if (length(Curr_Common_Peak_Set) > 0) {
				CurrOutFile <- paste0(dirname(opt$OutPrefix), '/Peaks_Common_', InpLabels[n1], '_and_', InpLabels[n2], '.bed')
				write.table(Merged_PeakData[Curr_Common_Peak_Set, ], CurrOutFile, row.names=F, col.names=F, sep="\t", quote=F, append=F)
			}
		}
	}

	# now generate the peaks common to three input sets
	if (length(InpFileList) >= 3) {
		for (n1 in (1:(length(InpFileList) - 2))) {
			for (n2 in ((n1+1):(length(InpFileList) - 1))) {
				for (n3 in ((n2+1):length(InpFileList))) {
					cat(sprintf("\n\n *** Dump peaks common to three categories - processing indices %s and %s and %s ", n1, n2, n3))
					Curr_Common_Peak_Set <- intersect(intersect(Ov_Idx_MergedRef[[n1]], Ov_Idx_MergedRef[[n2]]), Ov_Idx_MergedRef[[n3]])
					# dump the peaks common to these three categories
					if (length(Curr_Common_Peak_Set) > 0) {
						CurrOutFile <- paste0(dirname(opt$OutPrefix), '/Peaks_Common_', InpLabels[n1], '_and_', InpLabels[n2], '_and_', InpLabels[n3], '.bed')
						write.table(Merged_PeakData[Curr_Common_Peak_Set, ], CurrOutFile, row.names=F, col.names=F, sep="\t", quote=F, append=F)
					}
				}
			}
		}
	}

}	# end dump condition




