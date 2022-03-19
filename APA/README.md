Aggregate Peak Analysis (APA)
=========================

Enrichment of HiC or HiChIP (or any other 3C assays like ChIA-PET) signficant contacts 
with respect to reference HiC background contacts.

Defined in the paper:
Rao et. al., A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping, Cell 2014.


Usage: APA_Compute.r [options]


Options:
	--InpFile=INPFILE
		Input interaction / Loop file

	--headerInp
		If TRUE, input interaction file has header. Default FALSE.

	--midInp
		If TRUE, input interaction file has the bins in (chr1, mid1, chr2, mid2) instead of the conventional (chr1, start1, end1, chr2, start2, end2) format (similar to FitHiC). Default FALSE.

	--chrInp
		If TRUE, input interaction file has numbers as the chromosomes, for example 1, 2, instead of chr1, chr2...Default FALSE.

	--RefFile=REFFILE
		Reference interaction / Loop file, preferably from Hi-C data

	--headerRef
		Similar to the option --headerInp for the reference loop file. Default FALSE.

	--midRef
		Similar to the option --midInp for the reference loop file. Default FALSE.

	--chrRef
		Similar to the option --chrInp for the reference loop file. Default FALSE.

	--cccol=CCCOL
		Column number storing the contact counts in the REFERENCE interaction file. Default = 7.

	--binsize=BINSIZE
		Bin size. Default 5000 (5 Kb). Both the input and reference loop files should have identical bin sizes.

	--window=WINDOW
		Window size (+/- in both directions). Default 50000: 50 Kb in either directions). Should be a multiple of the parameter --binsize. Set as recommended in the MANGO paper.

	--distthrlow=DISTTHRLOW
		Interactions having distance lower than this threshold are discarded from analysis. Default 150000 (150 Kb). Set as recommended in the MANGO paper.

	--distthrhigh=DISTTHRHIGH
		Interactions having distance higher than this threshold are discarded from analysis. Default = 1000000 (1 Mb). Set as recommended in the MANGO paper.

	--low=LOW
		Lower distance threshold from the anchor segment, which is used to compute the APA score and the plot range. Default 15000, means that pixels with distance >= 15 Kb will be considered for APA mean score and pixel range computation. Set as recommended in the MANGO paper.

	--high=HIGH
		Upper distance threshold from the anchor segment, which is used to compute the APA score and the plot range. Default 30000, means that pixels with distance <= 30 Kb will be considered for APA mean score and pixel range computation. Set as recommended in the MANGO paper.

	--cclim=CCLIM
		Limit of the contact count for color scale display. Default = 50. Can be adjusted according to the visual APA color scale.

	--OutPrefix=OUTPREFIX
		Output prefix (including the output directory) used to name the output files. Mandatory parameter.

	--overwrite
		If TRUE, overwrites existing results. Default FALSE.

	-h, --help
		Show this help message and exit


