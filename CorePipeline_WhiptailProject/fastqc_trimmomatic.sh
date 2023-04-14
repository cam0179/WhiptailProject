#!/usr/bin/env bash

#SBATCH --job-name=FastQC-trimm
#SBATCH --output=logs/FastQC-trimm_%A.log
#SBATCH --ntasks=20
#SBATCH --mem=15gb
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --qos=class
#SBATCH --constraint=dmc

module load fastqc/0.11.9

workDir=/scratch/WhiptailProject

sampleDir=${workDir}/Cleaned/Trimmomatic
outDir=${workDir}/QC/Trimmomatic

mkdir -p $outDir

fastqc ${sampleDir}/* --outdir=$outDir -t 20
