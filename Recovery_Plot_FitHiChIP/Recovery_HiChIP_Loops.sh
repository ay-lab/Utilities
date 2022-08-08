#!/bin/bash -ex
#PBS -l nodes=1:ppn=1
#PBS -l mem=20GB
#PBS -l walltime=04:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
source ~/.bashrc
#source ~/.bash_profile
hostname
TMPDIR=/scratch
cd $PBS_O_WORKDIR

CurrDir='/home/sourya/proj/DEVELOPED_PACKAGES/Utility_Scripts/Recovery_Plot_FitHiChIP'
cd $CurrDir

CodeExec=$CurrDir'/Recovery_HiChIP_Loops.R'

InpLoopFile='/mnt/BioAdHoc/Groups/vd-vijay/sourya/Projects/2017_FitHiChIP/Results/HiChIP_Output_Loops/Mumbach_2017_paper/GM_HiChIP_H3K27ac_MERGED/FitHiChIP_HiCPro_Feb2018/FitHiChIP_Peak2ALL_b5000_L20000_U2000000/P2PBckgr_1/Coverage_Bias/FitHiC_BiasCorr_Resid_0_EqOcc_1/BinomDistr/FitHiChIP.interactions_FitHiC_Q0.01.bed'

RefLoopFile='/mnt/BioAdHoc/Groups/vd-vijay/sourya/DATA/HiC/RAO_2014/GSE63525_GM12878_primary_replicate_HiCCUPS_looplist_20Kb_2Mb.txt'

Rscript $CodeExec --RefFile $RefLoopFile --InpFile $InpLoopFile --headerInp --headerRef --QcolInp 0 --binsizeInp 5000 --binsizeRef 10000 --offset 5000 --OutDir $CurrDir --LabelRef 'HiCCUPS_Rao2014_GM12878_HiC_H3K27ac' --LabelInp 'FitHiChIP_GM12878_H3K27ac'

