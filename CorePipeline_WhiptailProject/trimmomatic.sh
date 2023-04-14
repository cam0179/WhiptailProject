#!/usr/bin/env bash

#SBATCH --job-name=Trimmomatic
#SBATCH --output=logs/Trimmomatic_%A.log
#SBATCH --ntasks=20
#SBATCH --mem=25gb
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --qos=class
#SBATCH --constraint=dmc

module load trimmomatic/0.39

workDir=/scratch/WhiptailProject

fastqDir=${workDir}/raw
outDir=${workDir}/Cleaned/Trimmomatic
samplesFile=${workDir}/samplesID.txt

mkdir -p $outDir/logs

while read id; do

    R1_in=${fastqDir}/${id}_1.fq.gz
    R2_in=${fastqDir}/${id}_2.fq.gz

    R1_out=${outDir}/${id}_1.fq.gz
    R2_out=${outDir}/${id}_2.fq.gz
    UP1_out=${outDir}/${id}_U1.fq.gz
    UP2_out=${outDir}/${id}_U2.fq.gz

    log=${outDir}/logs/${id}_trim_out.log
    adapters=${workDir}/AdaptersToTrim_All.fa

    trimmomatic PE -threads 15 -phred33 -trimlog $log \
        $R1_in $R2_in \
        $R1_out $UP1_out $R2_out $UP2_out \
        ILLUMINACLIP:${adapters}:2:35:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36

done < $samplesFile
