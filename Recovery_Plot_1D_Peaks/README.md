Recovery plot of reference ChIP-seq peaks
====================

See FitHiChIP manuscript, supplementary Fig. 11

Computes the percentage of reference peaks (like ChIP-seq peaks) which are captured by increasing number of HiChIP peaks or any other peaks, sorted by statistical significance.


Usage: Recovery_ChIP_Peaks.R [options]


Options:

	--RefFile=REFFILE
		Reference 1D peak file. Mandatory parameter.

	--InpFile=INPFILE
		Input ChIP or HiChIP peak file. Mandatory parameter.

	--headerRef
		If TRUE, reference file has header information. Default FALSE.

	--headerInp
		If TRUE, input file has header information. Default FALSE.

	--ascend
		If TRUE, sorting of peaks should be in ascending order of the q-value related column. By default, sorting of peaks is done in descending order since MACS2 stores -log10(q-value) which should be sorted in descending order. Default FALSE.

	--QcolInp=QCOLINP
		Column number storing the q-values (FDR) of the input interaction file. Default = 9, means that the 9'th column stores the q-value.

	--offset=OFFSET
		Offset / slack (in bp). Default 0 means a pair of peaks should be in exact overlapping. If 1000, 1 Kb slack is provided.

	--OutDir=OUTDIR
		Output directory.

	--LabelRef=LABELREF
		Label / category / method name of the reference peak file.

	--LabelInp=LABELINP
		Label / category / method name of the input peak file.

	-h, --help
		Show this help message and exit


A sample script is provided to show how this routine is executed.

Output
======

An output file of the following name will be created under the specified output directory.

Recovery_[LabelInp]_[LabelRef]offset[offset].txt

where,


LabelInp: Parameter --LabelInp

LabelRef: Parameter --LabelRef

offset: Parameter --offset

This output file will store two columns:

peakcnt: number of peaks from the input file sorted by statistical significance, used to compute the recovery of reference peaks.

frac: fraction of reference peaks overlapped.

Plotting recovery for a number of different methods
=========================

Once this script is run keeping the reference peaks fixed and varying the input peaks for various ChIP-seq or HiChIP specific peak files, each script will generate corresponding output files.

Use excel / other plotting tools to plot the trend of overlap according to increasing number of HiChIP peaks used.

