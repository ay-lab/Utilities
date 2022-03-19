#!/usr/bin/env Rscript

# ===========================================================
# script to analyze two or more interaction (loop) files
# and generate a venn diagram among those interactions
# ===>> upto 5 interaction files are supported

# Author: Sourya Bhattacharyya
# Vijay-Ay lab, LJI
# ===================================================

suppressMessages(library(GenomicRanges))
library(data.table)
library(optparse)
library(venn)
library(VennDiagram)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

fillcolors <- c("lightsteelblue", "saddlebrown", "lightpink", "burlywood", "darkgoldenrod", "seagreen")

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
	make_option(c("--FileList"), type="character", default=NULL, help="Comma or colon separated list of interaction / Loop files which are used for mutual comparison. Mandatory parameter."),
	make_option(c("--HeaderList"), type="character", default=NULL, help="A list of 1's or 0's separated by comma or colon, indicating whether the corresponding input file has a header (1) or not (0). Default: a list of all 1's (header = TRUE)"),	
	make_option(c("--Labels"), type="character", default=NULL, help="Comma or colon separated list of strings, indicating the label (class / cell / tissue) of individual input files"),

	make_option(c("--offset"), type="integer", action="store", default=0, help="Offset / slack (in bp). Default 0. If 0, means exact overlap is computed. If 5000, a pair of loops with interacting bins within 5 Kb for both sides will be considered as overlapping."),
  	make_option(c("--OutDir"), type="character", default=NULL, help="Output directory."),

	make_option(c("--RefSample"), type="integer", action="store", default=0, help="Integer which can be 0 (default) or between 1 to the number of samples (input files) provided. If > 0, corresponding input file will be treated as reference; overlap with respect to that file will be computed."),  	
  	make_option(c("--Dump"), type="logical", action="store_true", default=FALSE, help="If TRUE, loops belonging to different sets are dumped. Works when the number of input files <= 3.")
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
}

# boolean information of header for individual inputs
if (is.null(opt$HeaderList)) {
	HeaderList <- rep(1, length(InpFileList))
} else {
	HeaderList <- as.integer(unlist(strsplit(opt$HeaderList, "[,:]")))
	if (length(HeaderList) < length(InpFileList)) {
		HeaderList <- c(HeaderList, rep(0, (length(InpFileList) - length(HeaderList))))
	}
}

# create the output directory
system(paste("mkdir -p", opt$OutDir))
cat(sprintf("\n\n\n ==>>> InpFileList : %s ", paste(InpFileList, collapse=" ")))

# dump the input file list and the categories in a 
# separate text file
if (opt$Dump == TRUE) {
	textfile <- paste0(opt$OutDir, '/FileList.txt')
	fp_out <- file(textfile, "w")
	for (i in (1:length(InpFileList))) {
		outtext <- paste0("\n Input file number : ", i, "\n Name: ", InpFileList[i], "\n Label: ", InpLabels[i])
		writeLines(outtext, con=fp_out, sep="\n")
	}	
	close(fp_out)
}


#============================
# check if the value of RefSample is between 1 to length(InpFileList)
# in such case, overlap is computed only with respect to that sample
if ((opt$RefSample > 0) & (opt$RefSample <= length(InpFileList))) {

	# extract only the loops of the referred sample
	MasterFile <- paste0(opt$OutDir, '/Master_Interactions.bed')
	i <- opt$RefSample
	if (HeaderList[i] == 1) {
		system(paste("awk \'{if (NR>1) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}\'", InpFileList[i], ">", MasterFile))	
	} else {
		system(paste("awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}\'", InpFileList[i], ">", MasterFile))	
	}	

	# now read the reference loops
	Merged_IntData <- data.table::fread(MasterFile, header=F, sep="\t", stringsAsFactors=F)
	cat(sprintf("\n **** Number of rows in Merged_IntData : %s ", nrow(Merged_IntData)))

	# then compute overlap of loops for this master file
	# and for individual input files

	# list of indices - 
	# overlapping indices with respect to reference loops
	Ov_Idx_MergedRef <- list()

	# compute the overlap of merged loops
	# vs individual input files
	for (i in (1:length(InpFileList))) {
		InpData <- data.table::fread(InpFileList[i], header=as.logical(HeaderList[i]), sep="\t", stringsAsFactors=F)
		cat(sprintf("\n *** computing loop overlap - read input file idx : %s \n", i))

		# compute loop overlap using the custom function
		CurrOv <- OverlapLoop(Merged_IntData, InpData, boundary=1, offset=opt$offset)
		Ov_Idx_MergedRef[[i]] <- CurrOv$A_AND_B
		frac_ov_merged <- round((length(CurrOv$A_AND_B) * 1.0) / nrow(Merged_IntData), 2)
		cat(sprintf("\n **** overlap between merged loop to category: %s length overlap w.r.t merged loop: %s  fraction overlap merged ref : %s ", InpLabels[i], length(Ov_Idx_MergedRef[[i]]), frac_ov_merged))
	}

	# add - sourya
	if (length(InpFileList) <= 5) {
		# for categories upto 5, use the VennDiagram package
		vennPlotFile <- paste0(opt$OutDir, '/Venn_Loop_Overlap_Slack', opt$offset, '_', opt$RefSample, '.png')
		width_val <- 500 * length(InpFileList)
		height_val <- 500 * length(InpFileList)
		VennDiagram::venn.diagram(Ov_Idx_MergedRef, vennPlotFile, imagetype = "png", width=width_val, height=height_val, category.names = InpLabels, print.mode = c("raw", "percent"), fill=fillcolors[1:length(InpFileList)], sigdigs=3, euler.d=FALSE, scaled=FALSE, cat.cex=rep(1, length(InpFileList)))
	
	} else {
		# otherwise, use the venn package
		vennPlotFile <- paste0(opt$OutDir, '/Venn_Loop_Overlap_Slack', opt$offset, '_', opt$RefSample, '.pdf')
		vennPlotFile_NoCount <- paste0(opt$OutDir, '/Venn_Loop_Overlap_Slack', opt$offset, '_', opt$RefSample, '_NoCount.pdf')
		names(Ov_Idx_MergedRef) <- seq(1, length(InpLabels))
		legendstr <- c()
		for (i in (1:length(InpLabels))) {
			legendstr <- c(legendstr, paste(i, ":", InpLabels[i]))
		}		
		if (length(InpLabels)<=5) {
			cexil_val <- 1.5
		} else if (length(InpLabels)<=6) {
			cexil_val <- 0.5
		} else {
			cexil_val <- 0.45
		}
		pdf(vennPlotFile, width=8, height=6)
		venn::venn(Ov_Idx_MergedRef, zcolor="style", cexil=cexil_val, borders=FALSE)		
		legend("bottomright",legend=legendstr, lty=1, lwd=1, cex=1.2)
		dev.off()

		# plot pie chart as well
		pieChartPlotFile <- paste0(opt$OutDir, '/Pie_Chart_Loops.pdf')
		OvMat <- matrix(0, nrow=nrow(Merged_IntData), ncol=(length(InpFileList) - 1))
		c <- 0
		for (i in (1:length(InpFileList))) {
			if (i != opt$RefSample) {
				c <- c + 1
				OvMat[Ov_Idx_MergedRef[[i]], c] <- 1
			}
		}
		RowSums <- apply(OvMat, 1, sum)
		NumLoops_All <- c()
		for (cell_count in c(1:(length(InpFileList) - 1))) {
			cat(sprintf("\n cell_count : %s ", cell_count))
			num_Loops <- length(which(RowSums==cell_count))
			cat(sprintf("\n Number of loops belonging to %s cell types is : %s ", cell_count, num_Loops))
			NumLoops_All <- c(NumLoops_All, num_Loops)
		}

		# plot the pie chart
		lbls <- paste0(seq(1,(length(InpFileList) - 1)), "_Cell")
		pct <- round(NumLoops_All / sum(NumLoops_All)*100)
		lbls <- paste(lbls, NumLoops_All, pct) # add percents to labels 
		lbls <- paste(lbls,"%",sep="") # add % to labels 
		pdf(pieChartPlotFile, width=6, height=4)
		pie(NumLoops_All, labels = lbls, col=rainbow(length(lbls)), main="Loops in different cell types")
		dev.off()

	}
	# end add - sourya

	# dump the loops exclusive to the reference samples
	# i.e. which do not overlap with any other sample
	if (opt$Dump == TRUE)  {
		# prepare the set of loop indices exclusive to the reference loops
		excl_Ref <- seq(1, nrow(Merged_IntData))
		for (i in (1:length(InpFileList))) {
			if (i == opt$RefSample) {
				next
			}
			excl_Ref <- setdiff(excl_Ref, Ov_Idx_MergedRef[[i]])
		}
		# dump those loops
		write.table(Merged_IntData[excl_Ref, ], paste0(opt$OutDir, '/excl_Ref.bed'), row.names = FALSE, col.names = TRUE, sep = "\t", quote=FALSE, append=FALSE)
	}

} else {

	# first create a temporary interaction file
	# storing interactions from all the input files
	# (first 6 columns)
	MasterFile <- paste0(opt$OutDir, '/Master_Interactions.bed')
	if (file.exists(MasterFile) == FALSE) {
		for (i in (1:length(InpFileList))) {
			if (i == 1) {
				if (HeaderList[i] == 1) {
					system(paste("awk \'{if (NR>1) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}\'", InpFileList[i], ">", MasterFile))	
				} else {
					system(paste("awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}\'", InpFileList[i], ">", MasterFile))	
				}
			} else {
				if (HeaderList[i] == 1) {
					system(paste("awk \'{if (NR>1) {print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}\'", InpFileList[i], ">>", MasterFile))
				} else {
					system(paste("awk \'{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}\'", InpFileList[i], ">>", MasterFile))
				}
			}
		}	
	}


	# now sort the master interaction file 
	# according to the columns 1, 2 and 5
	MasterUniqFile <- paste0(opt$OutDir, '/Master_Interactions_Unique.bed')
	if (file.exists(MasterUniqFile) == FALSE) {
		system(paste("sort -k1,1 -k2,2n -k5,5n", MasterFile, "| awk -F\"\t\" \'!seen[$1, $2, $3, $4, $5, $6]++\' - >", MasterUniqFile))
	}

	# now read the complete set of (merged interactions)
	Merged_IntData <- data.table::fread(MasterUniqFile, header=F, sep="\t", stringsAsFactors=F)
	cat(sprintf("\n **** Number of rows in Merged_IntData : %s ", nrow(Merged_IntData)))

	#===================================
	# then compute overlap of loops for this master file
	# and for individual input files

	# list of indices - 
	# overlapping indices with respect to merged loops
	Ov_Idx_MergedRef <- list()

	# overlapping indices with respect to individual input loops
	Ov_Idx_Inp <- list()

	# count vector
	CountVec <- c()

	# compute the overlap of merged loops
	# vs individual input files
	for (i in (1:length(InpFileList))) {
		InpData <- data.table::fread(InpFileList[i], header=as.logical(HeaderList[i]), sep="\t", stringsAsFactors=F)
		cat(sprintf("\n *** computing loop overlap - read input file idx : %s \n", i))

		# compute peak overlap using the custom function
		CurrOv <- OverlapLoop(Merged_IntData, InpData, boundary=1, offset=opt$offset)
		
		Ov_Idx_MergedRef[[i]] <- CurrOv$A_AND_B
		Ov_Idx_Inp[[i]] <- CurrOv$B_AND_A

		frac_ov_merged <- round((length(CurrOv$A_AND_B) * 1.0) / nrow(Merged_IntData), 2)
		CountVec <- c(CountVec, frac_ov_merged)

		cat(sprintf("\n **** overlap between merged loop to category: %s length overlap w.r.t merged loop: %s  length overlap w.r.t input: %s  fraction overlap merged ref : %s ", InpLabels[i], length(Ov_Idx_MergedRef[[i]]), length(Ov_Idx_Inp[[i]]), frac_ov_merged))
	}

	# venn diagram
	if (length(InpFileList) <= 5) {
		# for categories upto 5, use the VennDiagram package		
		vennPlotFile <- paste0(opt$OutDir, '/Venn_Loop_Overlap_Slack', opt$offset, '.png')
		width_val <- 500 * length(InpFileList)
		height_val <- 500 * length(InpFileList)	
		VennDiagram::venn.diagram(Ov_Idx_MergedRef, vennPlotFile, imagetype = "png", width=width_val, height=height_val, category.names = InpLabels, print.mode = c("raw", "percent"), fill=fillcolors[1:length(InpFileList)], sigdigs=3, euler.d=FALSE, scaled=FALSE, cat.cex=rep(0.8, length(InpFileList)))
	} else {

		# otherwise, use the venn package
		vennPlotFile <- paste0(opt$OutDir, '/Venn_Loop_Overlap_Slack', opt$offset, '.pdf')
		vennPlotFile_NoCount <- paste0(opt$OutDir, '/Venn_Loop_Overlap_Slack', opt$offset, '_NoCount.pdf')
		names(Ov_Idx_MergedRef) <- seq(1, length(InpLabels))
		legendstr <- c()
		for (i in (1:length(InpLabels))) {
			legendstr <- c(legendstr, paste(i, ":", InpLabels[i]))
		}
		if (length(InpLabels)<=5) {
			cexil_val <- 1.5
		} else if (length(InpLabels)<=6) {
			cexil_val <- 0.5
		} else {
			cexil_val <- 0.45
		}

		pdf(vennPlotFile, width=8, height=6)
		venn::venn(Ov_Idx_MergedRef, zcolor="style", cexil=cexil_val, borders=FALSE)
		# venn::venn(Ov_Idx_MergedRef, zcolor=fillcolors, cexil=cexil_val, borders=FALSE)		
		legend("bottomright",legend=legendstr, lty=1, lwd=1, cex=0.5)
		dev.off()

		# plot pie chart as well
		pieChartPlotFile <- paste0(opt$OutDir, '/Pie_Chart_Loops.pdf')
		OvMat <- matrix(0, nrow=nrow(Merged_IntData), ncol=length(InpFileList))
		c <- 0
		for (i in (1:length(InpFileList))) {
			OvMat[Ov_Idx_MergedRef[[i]], i] <- 1
		}
		RowSums <- apply(OvMat, 1, sum)
		NumLoops_All <- c()
		for (cell_count in c(1:length(InpFileList))) {
			cat(sprintf("\n cell_count : %s ", cell_count))
			num_Loops <- length(which(RowSums==cell_count))
			cat(sprintf("\n Number of loops belonging to %s cell types is : %s ", cell_count, num_Loops))
			NumLoops_All <- c(NumLoops_All, num_Loops)
		}

		# plot the pie chart
		lbls <- paste0(seq(1,length(InpFileList)), "_Cell")
		pct <- round(NumLoops_All / sum(NumLoops_All)*100)
		lbls <- paste(lbls, NumLoops_All, pct) # add percents to labels 
		lbls <- paste(lbls,"%",sep="") # add % to labels 
		pdf(pieChartPlotFile, width=6, height=4)
		pie(NumLoops_All, labels = lbls, col=rainbow(length(lbls)), main="Loops in different cell types")
		dev.off()

	}

}	# end RefSample condition

