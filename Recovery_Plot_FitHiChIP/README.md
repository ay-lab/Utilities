
Recovery of reference loops by FitHiChIP (or HiChIP loops from other methdos)
=================

Computes the percentage of reference loops (like HiCCUPS loops) which are captured by increasing number of FitHiChIP (or other HiChIP loop caller) loops, when these loops are sorted by statistical significance.

Rscript Recovery_HiChIP_Loops.R [options]

Options:

	--RefFile=REFFILE
		Reference loop file (significant interactions). Mandatory parameter.

	--InpFile=INPFILE
		Input (FitHiChIP or other methods) HiChIP significant loop file. Mandatory parameter.

	--headerRef
		If TRUE, reference file has header information. Default FALSE.

	--headerInp
		If TRUE, input file has header information. Default FALSE.

	--QcolInp=QCOLINP
		Column number storing the q-values (FDR) of the input interaction file. Default = 0, means that the last column stores the q-value. Otherwise, a positive value for the target column number needs to be specified.

	--chrRef
		If TRUE, reference loop file has numbers as the chromosomes, for example 1, 2, instead of chr1, chr2. Default FALSE.

	--chrInp
		If TRUE, input loop file has numbers as the chromosomes, for example 1, 2, instead of chr1, chr2. Default FALSE.

	--midRef
		If TRUE, reference interaction file has only the midpoints of an  interval (chr1, mid1, chr2, mid2). Otherwise (default = FALSE), the reference file columns are (chr1, start1, end1, chr2, start2, end2).

	--midInp
		If TRUE, input interaction file has only the midpoints of an  interval (chr1, mid1, chr2, mid2). Otherwise (default = FALSE), the input file columns are (chr1, start1, end1, chr2, start2, end2).

	--binsizeRef=BINSIZEREF
		Bin size (in bp) for the reference interaction file. Required if --midRef is 1. Default = 5000 (5 Kb)

	--binsizeInp=BINSIZEINP
		Bin size (in bp) for the input interaction file. Required if --midInp is 1. Default = 5000 (5 Kb)

	--offset=OFFSET
		Offset / slack (in bp). Default 5000 (5 Kb) means a pair of loops with interacting bins within 5 Kb for both sides, will be considered as overlapping. If 0, exact overlap will be computed.

	--OutDir=OUTDIR
		Output directory.

	--LabelRef=LABELREF
		Label / category / method name of the reference loop file.

	--LabelInp=LABELINP
		Label / category / method name of the input loop file.

	-h, --help
		Show this help message and exit





A sample script is provided to show how this routine is executed.


Output
==========

An output file of the following name will be created under the specified output directory.

Recovery_[LabelInp]_[LabelRef]_offset_[offset].txt

where,

LabelInp: Parameter --LabelInp

LabelRef: Parameter --LabelRef

offset: Parameter --offset


This output file will store two columns:

loopcnt: number of loops from FitHiChIP (or any HiChIP loop caller) sorted by q-value, used to compute the recovery of reference loops.

frac: fraction of reference loops overlapped.


Plotting for a number of different methods
==============

Once this script is run keeping the reference loops fixed and varying the input loops for FitHiChIP and other competing approaches, each script will generate corresponding output files.

Use excel / other plotting tools to plot the trend of overlap according to increasing number of HiChIP loops used.








