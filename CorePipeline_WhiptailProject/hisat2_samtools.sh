#!/usr/bin/env bash

#SBATCH --job-name=HISAT2
#SBATCH --output=logs/HISAT2_%A.log
#SBATCH --cpus-per-task=20
#SBATCH --mem=40gb
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --qos=class
#SBATCH --constraint=dmc
#SBATCH --mail-user=gas0042@auburn.edu
#SBATCH --mail-type=ALL

source activate /home/aubclsb0305/.conda/envs/Hisat2

workDir=/scratch/WhiptailProject

fastqDir=${workDir}/Cleaned/Trimmomatic
outDir=${workDir}/Mapped/P_raffonei
samplesFile=${workDir}/samplesID.txt
indexFiles=${workDir}/Reference/P_raffonei/Hisat2_idx/rPodRaf1

mkdir -p $outDir

while read id; do

    mkdir -p ${outDir}/${id}

    R1=${fastqDir}/${id}_1.fq.gz
    R2=${fastqDir}/${id}_2.fq.gz
    
    samOut=${outDir}/${id}/${id}.sam
    bamSorted=${outDir}/${id}/${id}_sorted.bam
    bamStats=${outDir}/${id}/${id}_stats.txt
    

    hisat2 -p 20 --dta --phred33 \
        -x $indexFiles \
        -1 $R1 -2 $R2 \
        -S $samOut
    
    samtools sort -@ 20 -o $bamSorted $samOut
    samtools flagstat $bamSorted > $bamStats

done < $samplesFile