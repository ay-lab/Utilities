Overlap of peak (or 1D bed) files
==================================


Usage: Peak_Intersect.r [options]


Options:
	--FileList=FILELIST
		Comma or colon separated list of peak files (default NULL)

	--Labels=LABELS
		Comma or colon separated list of strings, associated with individual peak files. Used as the labels of the generated venn diagram as well. Default: category1, category2, etc..

	--OutPrefix=OUTPREFIX
		Prefix (including the output directory) of the output files.

	--Summit
		If TRUE, peaks relative to the summit position (and of a specified offset) are only considered.

	--Offset=OFFSET
		If Summit option is TRUE, peaks relative to the summit position and spanning this Offset number of base pairs are onky considered. Default = 500, means the peaks spanning 500 bp around the summit would be considered.

	--Dump
		If TRUE, peaks in different partitions are dumped.

	-h, --help
		Show this help message and exit


