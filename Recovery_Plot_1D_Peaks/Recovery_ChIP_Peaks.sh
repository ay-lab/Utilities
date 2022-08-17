#!/bin/bash

CurrDir='Utility_Scripts/Recovery_Plot_1D_Peaks'
cd $CurrDir

CodeExec=$CurrDir'/Recovery_ChIP_Peaks.R'

InpPeakFile='Inp_ChIP_Peaks.bed'

RefPeakFile='Ref_ChIP_Peaks.bed'

Rscript $CodeExec --RefFile $RefLoopFile --InpFile $InpLoopFile --headerInp --headerRef --QcolInp 0 --binsizeInp 5000 --binsizeRef 10000 --offset 5000 --OutDir $CurrDir --LabelRef 'HiCCUPS_Rao2014_GM12878_HiC_H3K27ac' --LabelInp 'FitHiChIP_GM12878_H3K27ac'

