#!/usr/bin/env bash

#SBATCH --job-name=StringTie
#SBATCH --output=logs/Stringtie_%A.log
#SBATCH --cpus-per-task=20
#SBATCH --mem=10gb
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --qos=class
#SBATCH --constraint=dmc

source activate /home/aubclsb0305/.conda/envs/Stringtie

workDir=/scratch/WhiptailProject

bamDir=${workDir}/Mapped/P_raffonei
outDir=${workDir}/Transcripts
samplesFile=${workDir}/samplesID.txt
refFile=${workDir}/Reference/P_raffonei/rPodRaf1_trasformed.gtf

mkdir -p $outDir

while read id; do

stringtie -p 20 -e -B -G  ${refFile} -l ${id} -o ${outDir}/${id}/${id}.gtf ${bamDir}/${id}/${id}_sorted.bam

done < ${samplesFile}

prepDE.py -i ${outDir} -g ${outDir}/gene_count_matrix.csv -t ${outDir}/transcript_count_matrix.csv
