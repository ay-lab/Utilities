#!/bin/bash

CurrDir='Utility_Scripts/Recovery_Plot_FitHiChIP'
cd $CurrDir

CodeExec=$CurrDir'/Recovery_HiChIP_Loops.R'

InpLoopFile='FitHiChIP.interactions_FitHiC_Q0.01.bed'

RefLoopFile='HiC/RAO_2014/GSE63525_GM12878_primary_replicate_HiCCUPS_looplist_20Kb_2Mb.txt'

Rscript $CodeExec --RefFile $RefLoopFile --InpFile $InpLoopFile --headerInp --headerRef --QcolInp 0 --binsizeInp 5000 --binsizeRef 10000 --offset 5000 --OutDir $CurrDir --LabelRef 'HiCCUPS_Rao2014_GM12878_HiC_H3K27ac' --LabelInp 'FitHiChIP_GM12878_H3K27ac'
