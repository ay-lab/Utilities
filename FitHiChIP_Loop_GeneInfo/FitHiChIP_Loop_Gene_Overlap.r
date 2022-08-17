#!/usr/bin/env Rscript

suppressMessages(library(GenomicRanges))
library(data.table)
options(scipen = 10)
options(datatable.fread.datatable=FALSE)

##====================
## parameters - users need to edit them
##====================

## FitHiChIP loop file
## first 6 fields should be chr1, start1, end1, chr2, start2, end2
## that is, two interacting intervals
LoopFile <- 'FitHiChIP.interactions_FitHiC_Q0.01.bed'

## GTF file (tab delimited) containing at least the following fields:
## gene names, chromosome and TSS information
GTFFile <- 'gEncode_input.gtf'

## column containing the chromosome information in the GTF file
chrcol <- 1 ## user can edit depending on the input GTF file

## column containing the TSS information in the GTF file
TSScol <- 2 ## user can edit depending on the input GTF file

## column containing the gene name information in the GTF file
geneNamecol <- 3 ## user can edit depending on the input GTF file

## offset (in bp) between TSS and the FitHiChIP loop anchors 
## considered for overlap
## 5000 means if the TSS is within 5 Kb of an interacting bin, it will be considered as overlapping
OffsetValue <- 5000

## output file
OutFile <- 'FitHiChIP_loops_overlap_gene.bed'

##====================
## important functions
##====================

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

##====================
## main code
##====================

## read FitHiChIP loops - assume it has header
LoopData <- data.table::fread(LoopFile, header=T)
CN <- colnames(LoopData)

## read GTF file - assume it has header
GTFData <- data.table::fread(GTFFile, header=T)

ov1 <- Overlap1D(GTFData[,c(chrcol, TSScol, TSScol)], LoopData[,1:3], boundary=0, offset=OffsetValue, uniqov=FALSE)
ov2 <- Overlap1D(GTFData[,c(chrcol, TSScol, TSScol)], LoopData[,4:6], boundary=0, offset=OffsetValue, uniqov=FALSE)

if (length(ov1$A_AND_B) > 0) {
	df1 <- cbind.data.frame(GTFData[ov1$A_AND_B, c(chrcol, TSScol, geneNamecol)], LoopData[ov1$B_AND_A, 1:3], LoopData[ov1$B_AND_A, 4:6], LoopData[ov1$B_AND_A, 7:ncol(LoopData)])
	colnames(df1) <- c('chr', 'TSS', 'geneName', 'Bin1Chr', 'Bin1Start', 'Bin1End', 'Bin2Chr', 'Bin2Start', 'Bin2End', CN[7:length(CN)])
}

if (length(ov2$A_AND_B) > 0) {
	df2 <- cbind.data.frame(GTFData[ov2$A_AND_B, c(chrcol, TSScol, geneNamecol)], LoopData[ov2$B_AND_A, 4:6], LoopData[ov2$B_AND_A, 1:3], LoopData[ov2$B_AND_A, 7:ncol(LoopData)])
	colnames(df2) <- c('chr', 'TSS', 'geneName', 'Bin1Chr', 'Bin1Start', 'Bin1End', 'Bin2Chr', 'Bin2Start', 'Bin2End', CN[7:length(CN)])
}

if (length(ov1$A_AND_B) > 0) & (length(ov2$A_AND_B) > 0) {
	finalDF <- rbind.data.frame(df1, df2)
} else if (length(ov1$A_AND_B) > 0) {
	finalDF <- df1
} else {
	finalDF <- df2
}
write.table(finalDF, OutFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

