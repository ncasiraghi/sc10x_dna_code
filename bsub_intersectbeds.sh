#!/bin/sh 
#BSUB -J intbed
#BSUB -n 50
#BSUB -q long 
#BSUB -e intbed.log 
#BSUB -o intbed.txt 

module load R/3.5.1 bedtools

Rscript intersectsegsdists.GeneBased.parallel.R
