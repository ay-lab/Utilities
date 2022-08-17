Overlap of various loop files (upto 5 loop files)
and their venn diagram
============================

Usage: Venn_Interactions.r [options]


Options:
	
	--FileList=FILELIST
		Comma or colon separated list of interaction / Loop files which are used for mutual comparison. First 6 fields of the interaction files should be chr1, start1, end1, chr2, start2, end2. That is, interacting chromosomes and the corresponding bin intervals should be provided in the first 6 fields. Mandatory parameter.

	--HeaderList=HEADERLIST
		A list of 1's or 0's separated by comma or colon, indicating whether the corresponding input file has a header (1) or not (0). Default: a list of all 1's (header = TRUE for all input interaction files)

	--Labels=LABELS
		Comma or colon separated list of strings, indicating the label (class / cell / tissue) of individual input interaction files.

	--offset=OFFSET
		Offset / slack (in bp). Default 0. If 0, means exact overlap is computed between the loops. If 5000, a pair of loops with interacting bins within 5 Kb for both sides will be considered as overlapping. In FitHiChIP, we used 5000 (i.e. 5 Kb slack).

	--OutDir=OUTDIR
		Output directory to contain the overlapping interactions and the final venn diagram.

	--RefSample=REFSAMPLE
		Integer which can be 0 (default) or between 1 to the number of samples (input files) provided. If 0 (default), overlap between all the set of loops are computed, like a classic venn diagram. Hoowever, if > 0, corresponding input file will be treated as reference. Here, all overlap will be computed with respect to that input file. See Fig. 2g in the FitHiChIP manuscript, where we have used GM12878 Hi-C HiCCUPS loops as the reference, and computed overlap with respect to it.
		
	--Dump
		If TRUE, loops belonging to different sets (i.e. section in Venn diagram) are dumped. Applicable when the number of input files <= 3.

	-h, --help
		Show this help message and exit


