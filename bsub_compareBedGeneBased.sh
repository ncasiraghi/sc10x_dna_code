#!/bin/sh 
#BSUB -J compbedgene
#BSUB -n 50
#BSUB -q night 
#BSUB -e compbedgene.log 
#BSUB -o compbedgene.txt 

module load R/3.5.1

Rscript compareBedGeneBased.R
