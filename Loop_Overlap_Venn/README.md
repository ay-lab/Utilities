Overlap of various loop files (upto 5 loop files)
and their venn diagram
============================

Usage: Venn_Interactions.r [options]


Options:
	--FileList=FILELIST
		Comma or colon separated list of interaction / Loop files which are used for mutual comparison. Mandatory parameter.

	--HeaderList=HEADERLIST
		A list of 1's or 0's separated by comma or colon, indicating whether the corresponding input file has a header (1) or not (0). Default: a list of all 1's (header = TRUE)

	--Labels=LABELS
		Comma or colon separated list of strings, indicating the label (class / cell / tissue) of individual input files

	--offset=OFFSET
		Offset / slack (in bp). Default 0. If 0, means exact overlap is computed. If 5000, a pair of loops with interacting bins within 5 Kb for both sides will be considered as overlapping.

	--OutDir=OUTDIR
		Output directory.

	--RefSample=REFSAMPLE
		Integer which can be 0 (default) or between 1 to the number of samples (input files) provided. If > 0, corresponding input file will be treated as reference; overlap with respect to that file will be computed.

	--Dump
		If TRUE, loops belonging to different sets are dumped. Works when the number of input files <= 3.

	-h, --help
		Show this help message and exit


