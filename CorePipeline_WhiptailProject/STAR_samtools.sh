#!/usr/bin/env bash

#SBATCH --job-name=STAR_sam
#SBATCH --output=logs/STARS-log_%A.log
#SBATCH --cpus-per-task=20
#SBATCH --mem=40gb
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --qos=class
#SBATCH --constraint=dmc
#SBATCH --mail-user=mag0099@auburn.edu
#SBATCH --mail-type=ALL

module load star

workDir=/scratch/WhiptailProject

fastqDir=${workDir}/Cleaned
outDir=${workDir}/Mapped/Melissa/Star-log
samplesFile=${workDir}/samplesID.txt
indexFiles=${workDir}/Reference/P_raffonei/Star_idx/rPodRaf1

mkdir -p $outDir


STAR --genomeDir $indexFiles \
--runThreadN 20 \
--readFilesIn $fastqDir \
--outFileNamePrefix $workDir/Mapped/Melissa/STAR_ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard
--readFilesCommand zcat
