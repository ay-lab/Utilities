
The code is for identifying the genes overlapping FitHiChIP loop endpoints.

In the parameter section, specify:

LoopFile: containing the FitHiChIP loops
	User needs to put the input loop file.

GTFFile:
	Reference GTF files containing gene annotations.
	The file should contain information of chromosome, gene TSS and gene names

chrcol:
	Column containing the chromosome in the GTFFile

TSScol:
	Column containing the gene TSS in the GTFFile

geneNamecol:
	Column containing the gene name in the GTFFile	

OffsetValue:
	The offset (slack) allowed to compute the overlap between gene TSS and the interacting bins. For example, offset = 5000 means 5 Kb slack is allowed.


OutFile:
	User can put a custom output file name.

	Format of this output file:
		First 3 columns: 'chr', 'TSS', 'geneName'
			Overlapping gene information.
		Next 3 columns: 'Bin1Chr', 'Bin1Start', 'Bin1End'
			FitHiChIP loop anchor overlapping the gene.
		Next 3 columns: 'Bin2Chr', 'Bin2Start', 'Bin2End'
			Other interacting bin of this FitHiChIP loop.
		Subsequent columns: FitHiChIP information related to this loop.




