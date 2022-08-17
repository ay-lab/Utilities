#!/bin/bash

# jar file of Juicer tool
# downloaded from https://github.com/aidenlab/juicer/wiki/Download
JuicerJarExec='juicer_tools.jar'

## FitHiChIP interaction file
fithichipfilename='FitHiChIP.interactions_FitHiC_Q0.01.bed'

## reference genome
RefGenome='hg19'

## output directory
OutDir='Motif_Out'

## file to be used for motif calling
## uses first 6 fields of FitHiChIP interactions
MotifInpFile='CTCF_Motif_Input.bed'
cut -f1-6 $fithichipfilename > $MotifInpFile

# call the CTCF motif orientation command
java -jar ${JuicerJarExec} motifs ${RefGenome} ${OutDir} ${MotifInpFile}


