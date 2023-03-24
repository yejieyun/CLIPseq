#!/bin/sh

#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#This script will load all the modules needed for the pipeline.
#Run before running the actual pipeline.
#''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

#############
#Load package
#############
module load cutadapt/2.10
module load STAR/2.7.6a
module load samtools/1.13.0
module load anaconda3/2020.11
module load bedtools/2.29.2
module load MEME-suite/5.4.1
module load FastQC
module load homer
module load bedtools 

######################
#Activate anaconda env
######################
conda activate py2