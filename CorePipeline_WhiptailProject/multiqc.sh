#!/usr/bin/env bash

#SBATCH --job-name=MultiQC
#SBATCH --output=logs/MultiQC_%A.log
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --qos=class
#SBATCH --constraint=dmc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gas0042@auburn.edu

source activate /home/aubclsb0305/.conda/envs/Multiqc

workDir=/scratch/WhiptailProject
fqc_rawDir=${workDir}/QC/Raw
fqc_trimmDir=${workDir}/QC/Trimmomatic
trimmomaticDir=${workDir}/scripts/logs/Trimmomatic_1021277.log
hisat2Log=${workDir}/Mapped/Gabriel/hisat-log/logs

outDir=${workDir}/QC/Multiqc

mkdir -p $outDir

multiqc $fqc_rawDir $fqc_trimmDir $trimmomaticDir $hisat2Log --outdir $outDir