#!/usr/bin/bash

#SBATCH --job-name=Star_idx
#SBATCH --output=logs/Star-idx_%A.log
#SBATCH --ntasks=10
#SBATCH --mem=100gb
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --qos=large
#SBATCH --constraint=dmc
#SBATCH --mail-user=mag0099@auburn.edu
#SBATCH --mail-type=ALL


module load star/2.7.0e
module load samtools/1.13

workDir=/scratch/WhiptailProject

baseName=${workDir}/Reference/P_raffonei/rPodRaf1
outIndex=${workDir}/Reference/P_raffonei/Star_idx

mkdir -p ${workDir}/Reference/P_raffonei/Star_idx

STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir $outIndex \
--genomeFastaFiles ${baseName}.fna \
--sjdbGTFfile ${baseName}.gtf \
--sjdbOverhang -1
