#!/usr/bin/env bash

#SBATCH --job-name=FastQC
#SBATCH --output=logs/FastQC_%A.log
#SBATCH --ntasks=20
#SBATCH --mem=15gb
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --qos=class
#SBATCH --constraint=dmc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gas0042@auburn.edu

module load fastqc/0.11.9

workDir="/scratch/WhiptailProject"

sampleDir="${workDir}/raw"
outDir="${workDir}/qc"

mkdir -p $outDir

fastqc ${sampleDir}/* --outdir=$outDir -t 20
